# CpG mutation spectrum (gnomAD v4.1)

Counts mutated C and G SNVs, split by CpG context, substitution type, gene, and
(for coding CpGs) whether the CpG sits within a single codon or is split across
two codons.

## Inputs

- gnomAD v4.1 release **sites** VCFs (`gnomad.*.v4.1.sites.chr*.vcf.bgz`), which
  carry the `vep` (VEP CSQ) INFO field used for gene/codon annotation.
- GRCh38 reference FASTA, **uncompressed** (the `.fai` reader seeks by byte
  offset, so bgzipped won't work; the Snakefile rejects it). Contig names must be
  `chr`-prefixed to match the VCF. The Snakemake `faidx` rule builds the `.fai`
  index automatically (needs `samtools`); for a by-hand run do `samtools faidx GRCh38.fa` first.

Confirm the CSQ subfield names in your files once:

```bash
bcftools +split-vep -l gnomad.exomes.v4.1.sites.chr1.vcf.bgz
```

You should see `Gene`, `SYMBOL`, `Feature`, `BIOTYPE`, `Consequence`,
`CANONICAL`, `MANE_SELECT`, `Codons`. If any differ, adjust the `-f` string below.

## Run

### Option A: Snakemake (recommended)

Edit `REF`, `VCF_TEMPLATE`, and `EXON_MAP` at the top of the `Snakefile`, then:

```bash
snakemake -j 8 all
```

This runs each chromosome in parallel, merges, and applies the exon-boundary
refinement, producing `variant_level.tsv` and `gene_level.refined.tsv`.

### Option B: by hand

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

## Outputs

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
