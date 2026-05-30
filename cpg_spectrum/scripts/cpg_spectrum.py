#!/usr/bin/env python3
"""Classify gnomAD C/G SNVs by CpG context, mutation spectrum, gene, and codon.

Reads a TSV stream produced by `bcftools +split-vep` (one row per
transcript consequence; see cpg_spectrum/README.md for the exact command) and
emits two tidy tables:

  --variant-out : one row per variant allele (deduped across transcripts).
                  Covers the genome-wide question "mutated Cs / Gs, CpG vs not,
                  broken down by substitution type". Includes intergenic sites.

  --gene-out    : one row per (variant allele x gene), using a single anchor
                  transcript per gene (MANE Select, else Ensembl canonical).
                  Adds the codon classification for coding CpGs.

CpG status is taken from the reference: a C is CpG if the next base (POS+1) is
G; a G is CpG if the previous base (POS-1) is C. Substitution type is reported
on the forward/reference strand (so C>T and G>A stay distinct).

Codon classification (coding CpGs only) is read from VEP's strand-corrected
`Codons` field. 5'-CG-3' is its own reverse complement, so the changed base's
codon position and coding-strand identity fully determine it:

    split across two codons  <=>  (coding base C at codon pos 3)
                              or   (coding base G at codon pos 1)
    otherwise                ->   within a single codon

Caveat: if a CpG straddles an exon/intron boundary (its partner base is
intronic) the "adjacent CDS base" assumption fails; such sites are the
apparent splits at exon edges and should be confirmed against an exon map.
"""
import argparse
import sys


# --- minimal random-access FASTA reader (no pysam) -------------------------

class Fasta:
    """Single-base lookups against an uncompressed, samtools-faidx'd FASTA."""

    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as f:
            for line in f:
                name, length, offset, linebases, linewidth = line.split("\t")[:5]
                self.idx[name] = (
                    int(length), int(offset), int(linebases), int(linewidth)
                )

    def base(self, chrom, pos):
        """Return the 1-based reference base at chrom:pos, or 'N' if out of range."""
        rec = self.idx.get(chrom)
        if rec is None:
            return "N"
        length, offset, linebases, linewidth = rec
        if pos < 1 or pos > length:
            return "N"
        byte = offset + (pos - 1) // linebases * linewidth + (pos - 1) % linebases
        self.fh.seek(byte)
        return self.fh.read(1).decode().upper()


# --- codon classification --------------------------------------------------

def classify_codon(codons, is_cpg):
    """Return (coding_base, codon_pos, codon_class) from a VEP Codons string.

    codons looks like 'Cgt/Tgt' (changed base uppercase). Only meaningful for
    coding SNVs; codon_class is one of single_codon / split_codon / not_cpg /
    noncoding / NA.
    """
    if not codons or codons == "." or "/" not in codons:
        return "", "", "noncoding"
    ref_codon = codons.split("/")[0]
    ups = [i for i, c in enumerate(ref_codon) if c.isupper()]
    if len(ups) != 1:
        return "", "", "NA"          # MNV/indel-style codon; not handled here
    i = ups[0]
    base = ref_codon[i].upper()
    pos = i + 1
    if not is_cpg:
        return base, str(pos), "not_cpg"
    split = (base == "C" and pos == 3) or (base == "G" and pos == 1)
    return base, str(pos), "split_codon" if split else "single_codon"


# --- per-gene anchor-transcript selection ----------------------------------

def pick_anchor_rows(rows):
    """One row per gene: MANE Select > Ensembl canonical > first seen."""
    best = {}
    for r in rows:
        gene = r["gene"]
        if not gene or gene == ".":
            continue
        rank = 2 if r["mane"] not in ("", ".") else (1 if r["canonical"] == "YES" else 0)
        cur = best.get(gene)
        if cur is None or rank > cur[0]:
            best[gene] = (rank, r)
    return [v[1] for v in best.values()]


COLS = ["chrom", "pos", "ref", "alt", "ac",
        "gene", "symbol", "feature", "biotype", "consequence",
        "canonical", "mane", "codons"]


def parse_row(line):
    f = line.rstrip("\n").split("\t")
    return dict(zip(COLS, f))


def is_cpg(fasta, chrom, pos, ref):
    if ref == "C":
        return fasta.base(chrom, pos + 1) == "G"
    if ref == "G":
        return fasta.base(chrom, pos - 1) == "C"
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True,
                    help="Reference FASTA (uncompressed, with .fai); contigs must match the VCF")
    ap.add_argument("--variant-out", required=True)
    ap.add_argument("--gene-out", required=True)
    ap.add_argument("--vcf-tsv", default="-",
                    help="bcftools +split-vep TSV stream (default: stdin)")
    args = ap.parse_args()

    fasta = Fasta(args.fasta)
    src = sys.stdin if args.vcf_tsv == "-" else open(args.vcf_tsv)
    vout = open(args.variant_out, "w")
    gout = open(args.gene_out, "w")

    vout.write("\t".join(["chrom", "pos", "ref", "alt", "ac",
                          "is_cpg", "sub_type"]) + "\n")
    gout.write("\t".join(["chrom", "pos", "ref", "alt", "ac", "is_cpg", "sub_type",
                          "gene", "symbol", "transcript", "biotype", "consequence",
                          "coding_base", "codon_pos", "codon_class"]) + "\n")

    def flush(rows):
        if not rows:
            return
        r0 = rows[0]
        chrom, pos, ref, alt = r0["chrom"], int(r0["pos"]), r0["ref"], r0["alt"]
        if len(ref) != 1 or len(alt) != 1 or ref not in ("C", "G"):
            return                                   # defensive; bcftools already filters
        cpg = is_cpg(fasta, chrom, pos, ref)
        sub = f"{ref}>{alt}"
        vout.write("\t".join([chrom, str(pos), ref, alt, r0["ac"],
                              "1" if cpg else "0", sub]) + "\n")
        for r in pick_anchor_rows(rows):
            cbase, cpos, cclass = classify_codon(r["codons"], cpg)
            gout.write("\t".join([chrom, str(pos), ref, alt, r0["ac"],
                                  "1" if cpg else "0", sub,
                                  r["gene"], r["symbol"], r["feature"],
                                  r["biotype"], r["consequence"],
                                  cbase, cpos, cclass]) + "\n")

    key = None
    rows = []
    for line in src:
        if not line.strip():
            continue
        r = parse_row(line)
        k = (r["chrom"], r["pos"], r["ref"], r["alt"])
        if k != key:
            flush(rows)
            rows = []
            key = k
        rows.append(r)
    flush(rows)

    for fh in (vout, gout):
        fh.close()


if __name__ == "__main__":
    main()
