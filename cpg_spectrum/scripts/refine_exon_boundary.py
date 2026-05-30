#!/usr/bin/env python3
"""Refine codon_class for coding CpGs against a per-exon CDS map.

cpg_spectrum.py labels a coding CpG single_codon / split_codon assuming the
CpG's partner base is the adjacent CDS base. That assumption fails when the CpG
straddles an exon/intron junction: the partner base is genomically adjacent but
intronic (or UTR), so it is not the next codon base at all.

This step re-checks each coding CpG. The CpG occupies genomic positions p and
its neighbour (p+1 for a ref C, p-1 for a ref G). If that partner position is
not inside any CDS exon of the row's anchor transcript, the codon-split call is
unreliable and codon_class is rewritten to `exon_boundary`. A `partner_coding`
column (yes/no/unknown) records the check; `unknown` means the transcript was
not found in the exon map, so the original label is kept.

Strand is irrelevant here: CpG bases are genomically adjacent and CDS spans are
genomic, so the containment test is the same for + and - strand genes.
"""
import argparse
import sys
from collections import defaultdict


def strip_version(tid):
    return tid.split(".")[0] if tid else tid


def load_cds_intervals(path):
    """transcript (unversioned) -> list of (cds_start, cds_end) genomic spans."""
    cds = defaultdict(list)
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        tcol, scol, ecol = col["transcript_id"], col["cds_start"], col["cds_end"]
        for line in f:
            r = line.rstrip("\n").split("\t")
            s, e = r[scol], r[ecol]
            if not s or not e or s in (".", "NA") or e in (".", "NA"):
                continue                       # UTR-only exon: no CDS span
            cds[strip_version(r[tcol])].append((int(s), int(e)))
    return cds


def partner_in_cds(intervals, partner_pos):
    return any(s <= partner_pos <= e for s, e in intervals)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-tsv", required=True, help="gene_level.tsv from cpg_spectrum.py")
    ap.add_argument("--exon-map", required=True, help="mane_exons_grch3*.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cds = load_cds_intervals(args.exon_map)

    with open(args.gene_tsv) as f, open(args.out, "w") as out:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        out.write("\t".join(header + ["partner_coding"]) + "\n")

        for line in f:
            r = line.rstrip("\n").split("\t")
            codon_class = r[col["codon_class"]]
            partner_coding = ""
            if r[col["is_cpg"]] == "1" and codon_class in ("single_codon", "split_codon"):
                ref = r[col["ref"]]
                pos = int(r[col["pos"]])
                partner = pos + 1 if ref == "C" else pos - 1
                tx = strip_version(r[col["transcript"]])
                if tx not in cds:
                    partner_coding = "unknown"          # transcript not in map; keep label
                elif partner_in_cds(cds[tx], partner):
                    partner_coding = "yes"
                else:
                    partner_coding = "no"
                    r[col["codon_class"]] = "exon_boundary"
            out.write("\t".join(r + [partner_coding]) + "\n")


if __name__ == "__main__":
    main()
