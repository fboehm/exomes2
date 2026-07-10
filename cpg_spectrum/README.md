# CpG mutation spectrum (gnomAD v4.1 & RGC Million Exomes)

Counts mutated C and G SNVs, split by CpG context, substitution type, gene, and
(for coding CpGs) whether the CpG sits within a single codon or is split across
two codons.

One pipeline handles two data sources, selected by `source:` in `config.yaml`
(or `--config source=…`):

- **`gnomadv4-1`** — v4.1 release sites VCFs, already VEP-annotated and chr-prefixed.
- **`rgc`** — RGC Million Exomes frequency VCFs, which carry **no functional
  annotation** and use Ensembl-style contigs (`1..22,X,Y`). For this source the
  pipeline renames contigs to `chr*` and runs Ensembl VEP itself.

Everything downstream of the `bcftools +split-vep` TSV (`cpg_spectrum.py`,
`refine_exon_boundary.py`, `summarize.R`, the exon map) is source-agnostic.

## Inputs

- **Source VCFs** — path template set per source in `config.yaml`
  (`sources.<source>.vcf_template`, with a `{chrom}` = `chr1..chrY` wildcard):
  - `gnomadv4-1`: v4.1 release **sites** VCFs, which carry the `vep` (VEP CSQ) INFO
    field used for gene/codon annotation.
  - `rgc`: `rgc_me_variant_frequencies_{chrom}_20231004.vcf.gz`. Whole-cohort
    allele count is `INFO/ALL_AC` (per-population `*_AC` are probabilistic
    floats); contigs inside are `1..22,X,Y`.
- **GRCh38 reference FASTA**, **uncompressed** and **chr-prefixed** (the `.fai`
  reader seeks by byte offset, so bgzipped won't work; the Snakefile rejects it).
  Shared by both sources — for `rgc` the contig rename (`1→chr1`, via
  `data/ens2ucsc_chrs.txt`) makes RGC records match this same reference and the
  exon map. The `faidx` rule builds the index (needs `samtools`); by hand,
  `samtools faidx GRCh38.fa` first.
- **`rgc` only:** Ensembl VEP + a GRCh38 offline cache (see the VEP section below).

Confirm the CSQ subfield names for a `gnomadv4-1` source once (the `rgc` source
generates them itself via VEP `--fields`, so no check needed):

```bash
bcftools +split-vep -l gnomad.exomes.v4.1.sites.chr1.vcf.bgz
```

You should see `Gene`, `SYMBOL`, `Feature`, `BIOTYPE`, `Consequence`,
`CANONICAL`, `MANE_SELECT`, `Codons`. If any differ, adjust `SPLIT_VEP_FMT` /
`sources.gnomadv4-1.csq_tag` in the Snakefile/config.

## Run

### Option A: Snakemake (recommended)

Set `source`, `ref`, `exon_map`, and the per-source `vcf_template` in
`config.yaml`, then:

```bash
# gnomAD (already annotated — no VEP needed).
# NOTE: the target (`all`) must come BEFORE --config. --config consumes every
# following token as key=value, so `--config source=… all` makes Snakemake choke
# on `all`.
snakemake -j 8 all --config source=gnomadv4-1

# RGC Million Exomes (runs VEP; --use-envmodules picks up the cluster's `ensembl` module)
snakemake -j 8 --use-envmodules all --config source=rgc
```

This runs each chromosome in parallel, merges, and applies the exon-boundary
refinement, producing `variant_level.tsv` and `gene_level.refined.tsv`. Outputs
are written under a per-source tree, `results/<source>/` (e.g.
`results/rgc/variant_level.tsv`, `results/gnomadv4-1/variant_level.tsv`), so the
two sources never overwrite each other and can be built side by side. For
`rgc` it first renames contigs, normalizes/filters, and annotates with VEP
(`prep_chrom` → `vep_chrom` → `classify_chrom`); for `gnomadv4-1` the VEP rule is
absent and `classify_chrom` reads the already-annotated prepped VCF.

### Option B: by hand (gnomAD)

The recipe below is the `gnomadv4-1` case (already annotated). For `rgc` you must
first rename contigs and run VEP (see the VEP section); prefer Snakemake there.

```bash
REF=GRCh38.fa
for vcf in gnomad.*.v4.1.sites.chr*.vcf.bgz; do
  bcftools norm -m- -f "$REF" "$vcf" \
  | bcftools view -f PASS -v snps -i 'REF="C" || REF="G"' \
  | bcftools +split-vep -a vep -d \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%Gene\t%SYMBOL\t%Feature\t%BIOTYPE\t%Consequence\t%CANONICAL\t%MANE_SELECT\t%Codons\n'
done | python3 scripts/cpg_spectrum.py \
        --fasta "$REF" \
        --variant-out variant_level.tsv \
        --gene-out gene_level.tsv

python3 scripts/refine_exon_boundary.py \
  --gene-tsv gene_level.tsv \
  --exon-map ../uniprot_exon_map/mane_exons_grch38.tsv \
  --out gene_level.refined.tsv
```

- `norm -m-` splits multiallelics so each ALT gets a scalar `AC`.
- `view -f PASS -v snps -i ...` keeps PASS biallelic C/G SNVs.
- `+split-vep -d` emits one row per transcript consequence; the Python script
  re-groups them per variant and picks one anchor transcript per gene
  (MANE Select, else Ensembl canonical).
- `refine_exon_boundary.py` re-checks coding CpGs against the exon map (below).

## VEP annotation & cache (rgc source)

RGC-ME frequency VCFs carry no gene/consequence annotation, so the `rgc` source
runs Ensembl VEP itself (rule `vep_chrom`) to reproduce gnomAD's CSQ subfields
(`Gene, SYMBOL, Feature, BIOTYPE, Consequence, CANONICAL, MANE_SELECT, Codons`).
Because those subfield names and the `Codons` field match gnomAD's, the
downstream scripts are unchanged.

VEP is expected on the cluster as a module; run Snakemake with
`--use-envmodules` (the module name is `vep.module` in `config.yaml`). A conda
fallback (`envs/vep.yaml`) is used with `--use-conda` instead.

**One-time cache download.** VEP needs a GRCh38 offline cache; `vep.cache_dir`
must contain it and `vep.release` must match its version. On a login node (which
usually has network), with the module loaded:

```bash
module load ensembl                      # provides `vep` and `vep_install`
vep --help | grep -i 'versions\|ensembl-vep'   # note the release, set vep.release to it
mkdir -p /scratch/jacks.local/frederick.boehm/exomes/data/vep_cache

# cache only (-a c); we pass our own --fasta, so skip the cache FASTA
vep_install -a c -s homo_sapiens -y GRCh38 \
  --CACHEDIR /scratch/jacks.local/frederick.boehm/exomes/data/vep_cache \
  --CACHE_VERSION 112 --NO_HTSLIB
```

This creates `…/vep_cache/homo_sapiens/112_GRCh38/`, which includes the
`chr_synonyms.txt` the rule passes via `--synonyms` so the Ensembl cache resolves
the `chr`-prefixed contigs. Match `--CACHE_VERSION` to the VEP binary's release
and to `vep.release`; align both with MANE v1.5 / GENCODE v49 (Ensembl ~112) so
transcript IDs line up with `mane_exons_grch38.tsv`.

To restrict to on-target exome calls, set `sources.rgc.include_extra:
"ON_TARGET=1"` (confirm the `INFO/ON_TARGET` Type in the VCF header first).

## Outputs

All output paths below are relative to the per-source tree `results/<source>/`
(e.g. `results/rgc/variant_level.tsv`). The `summary/` and `awk` examples in the
following sections show bare names for brevity — prefix them the same way.

`variant_level.tsv` (one row per variant allele, includes intergenic):

| chrom | pos | ref | alt | ac | is_cpg | sub_type |

`gene_level.refined.tsv` (one row per variant x anchor gene):

| chrom | pos | ref | alt | ac | is_cpg | sub_type | gene | symbol | transcript | biotype | consequence | coding_base | codon_pos | codon_class | partner_coding |

`sub_type` is forward-strand genomic (`C>T`, `G>A`, ...). `coding_base`/`codon_pos`
are coding-strand (from VEP `Codons`). `codon_class` is one of
`single_codon`, `split_codon`, `exon_boundary`, `not_cpg`, `noncoding`, `NA`.
`partner_coding` (`yes`/`no`/`unknown`/empty) records the exon-map check from the
refinement step (`unknown` = anchor transcript absent from the exon map, so the
codon label was left untouched).

(`gene_level.tsv` is the pre-refinement intermediate, same columns minus
`partner_coding`.)

## Summaries

Genome-wide spectrum (sites and AC-weighted), e.g. with awk:

```bash
# sites: C/G x CpG x substitution
tail -n +2 variant_level.tsv \
  | awk -F'\t' '{s[$6"\t"$7]++; a[$6"\t"$7]+=$5}
                END{for(k in s) print k"\t"s[k]"\t"a[k]}' \
  | sort   # columns: is_cpg  sub_type  n_sites  sum_AC
```

Per gene:

```bash
tail -n +2 gene_level.tsv \
  | awk -F'\t' '{n[$9]++; a[$9]+=$5} END{for(g in n) print g"\t"n[g]"\t"a[g]}' \
  | sort -k2,2nr   # symbol  n_sites  sum_AC
```

CpG codon class, per gene (coding CpGs only):

```bash
tail -n +2 gene_level.refined.tsv \
  | awk -F'\t' '$6==1 && $15 ~ /codon|exon_boundary/ {print $9"\t"$15}' \
  | sort | uniq -c   # counts of single_codon / split_codon / exon_boundary per gene
```

Or load both TSVs into R/pandas and pivot — they're tidy/long on purpose.

## Exon-boundary CpGs (handled)

`codon_class` initially assumes the CpG's partner base is the adjacent CDS base.
When a CpG straddles an exon/intron junction (partner base intronic), that
assumption fails. The refinement step (`refine_exon_boundary.py`, run
automatically by the Snakefile) checks each coding CpG's partner position
against the per-exon CDS spans in `../uniprot_exon_map/mane_exons_grch38.tsv`
and relabels such sites `exon_boundary`, recording the check in
`partner_coding`. For a clean intra-CDS analysis, filter to
`codon_class in {single_codon, split_codon}`.
