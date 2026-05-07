# UniProt per-exon mapping tables

Two TSVs, one row per exon of one canonical transcript per protein-coding gene:

- `mane_exons_grch38.tsv` — 171,315 rows, 19,769 genes (19,228 MANE Select + 541 Ensembl canonical fallback)
- `mane_exons_grch37.tsv` — 170,634 rows, 19,648 genes (19,133 MANE Select + 515 Ensembl canonical fallback)

Anchor transcript per gene: MANE Select where available, else Ensembl_canonical.
Both tables use the **same Ensembl transcript IDs** (no UCSC liftOver — GENCODE
publishes the lift37 GTF natively).

## Columns

| # | Column | Notes |
|---|---|---|
| 1 | `gene_id` | Ensembl gene ID, unversioned |
| 2 | `gene_name` | HGNC symbol |
| 3 | `transcript_id` | Ensembl transcript ID, unversioned |
| 4 | `transcript_id_versioned` | with `.N` version |
| 5 | `uniprot_acc` | UniProt canonical accession (Swiss-Prot) |
| 6 | `chrom` | e.g. `chr19` |
| 7 | `exon_start` | 1-based, inclusive (GTF convention) |
| 8 | `exon_end` | 1-based, inclusive |
| 9 | `strand` | `+` or `-` |
| 10 | `exon_number` | 1-based, in 5'→3' transcript order |
| 11 | `exon_id` | Ensembl exon ID, unversioned |
| 12 | `cds_start` | coding portion only; empty for UTR-only exons |
| 13 | `cds_end` | empty for UTR-only exons |
| 14 | `aa_start` | first UniProt AA position whose codon overlaps this exon's CDS; 1-based |
| 15 | `aa_end` | last AA position; AAs at exon boundaries appear on both flanking rows because their codon spans the splice site |
| 16 | `anchor_source` | `MANE_Select` or `Ensembl_canonical` |
| 17 | `aa_range_consistent` | `T` if total CDS length / 3 equals UniProt protein length; `F` for selenoproteins, readthrough cases, and lift37 partial mappings |
| 18 | `uniprot_seq` | full canonical protein sequence (repeated per row of the same transcript) |

## Caveats

- **Boundary codons**: When a codon spans an exon-exon junction, the AA position appears on both exons (so a query "which exon contains residue X?" never misses). If you want non-overlapping ranges, take `aa_end - 1` for all but the last coding exon when the last codon is split.
- **`aa_range_consistent=F` rows** (~0.3% GRCh38, ~1.3% GRCh37): the AA columns may not align with the UniProt sequence. Causes:
  - Selenoproteins (UGA stop codon read through as Sec)
  - GENCODE/UniProt isoform-call disagreements
  - Lift37 truncation at scaffold breaks (mostly the 193 GRCh37-only mismatches)
  Filter with `awk -F'\t' '$17=="T"'` to keep only fully-consistent transcripts.
- **UTR-only exons** have empty `cds_start`/`cds_end`/`aa_start`/`aa_end` (~9,200 such rows in GRCh38).
- **Sequences are repeated per row** (~88 MB of sequence × 9-exon avg). If that's a problem, deduplicate with `sort -u -k3,3` then join back, or extract a separate FASTA.

## Reproduction

The full pipeline (downloads + filtering + both builds) is wrapped in a Snakefile:

```
snakemake -j 4 all
```

That fetches MANE, GENCODE GRCh38 + GRCh37-lift, and UniProt into `data/`, filters
the UniProt idmapping down to Ensembl_TRS rows, and emits both TSVs. To pin
different upstream releases, edit `MANE_VERSION` / `GENCODE_VERSION` at the top
of the `Snakefile` (verify against the FTP listings first — Snakemake won't
auto-discover newer releases).

To run a single build without Snakemake:

```
python3 scripts/build_exon_map.py \
  --gtf data/gencode.v49.primary_assembly.annotation.gtf.gz \
  --mane data/MANE.GRCh38.v1.5.summary.txt.gz \
  --fasta data/UP000005640_9606.fasta.gz \
  --xref data/uniprot_ensembl_trs.tsv.gz \
  --out mane_exons_grch38.tsv
```

Swap the GTF for `gencode.v49lift37.annotation.gtf.gz` for the GRCh37 build.

## Data sources

- MANE Select v1.5 — `ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/`
- GENCODE v49 (GRCh38) and v49lift37 (GRCh37) — `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/`
- UniProt human reference proteome `UP000005640_9606` and `HUMAN_9606_idmapping.dat.gz` — `ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/`
