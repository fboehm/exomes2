#!/usr/bin/env python3
"""Build per-exon TSV mapping genomic coordinates to UniProt canonical sequences.

Anchor transcript per gene: MANE Select if available, else Ensembl_canonical.
One row per exon. Columns include exon and CDS coords plus the AA range that
exon's CDS encodes, and the full UniProt canonical protein sequence.
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def strip_version(tid):
    return tid.split(".")[0] if tid else tid


_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def parse_gtf_attrs(s):
    out = {}
    tags = []
    for m in _ATTR_RE.finditer(s):
        k, v = m.group(1), m.group(2)
        if k == "tag":
            tags.append(v)
        else:
            out[k] = v
    out["_tags"] = tags
    return out


def load_fasta(path):
    """Parse UniProt FASTA. Returns dict[accession] = sequence."""
    seq = {}
    cur = None
    buf = []
    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seq[cur] = "".join(buf)
                parts = line[1:].split("|")
                cur = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seq[cur] = "".join(buf)
    return seq


def load_mane(path):
    """Return dict[ensembl_transcript_unversioned] = ensembl_gene_unversioned for MANE Select transcripts.

    Note: MANE summary v1.5+ no longer carries a UniProt column. We resolve
    UniProt accessions for these transcripts via the Ensembl_TRS xref file.
    """
    out = {}
    with open_maybe_gz(path) as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        i_gene = header.index("Ensembl_Gene")
        i_trs = header.index("Ensembl_nuc")
        i_status = header.index("MANE_status")
        for line in f:
            row = line.rstrip("\n").split("\t")
            if "MANE Select" not in row[i_status]:
                continue
            out[strip_version(row[i_trs])] = strip_version(row[i_gene])
    return out


def load_uniprot_xref(path, fasta_keys):
    """Return dict[ensembl_transcript_unversioned] = uniprot_acc (canonical only).

    Accessions in the xref may carry an isoform suffix like "P31946-2".
    "P31946" or "P31946-1" both refer to the canonical isoform; anything else is
    an alternative isoform and is dropped (we want canonical-only).
    The base accession (suffix stripped) must exist in the reference-proteome FASTA.
    """
    out = {}
    with open_maybe_gz(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc, key, val = parts[0], parts[1], parts[2]
            if key != "Ensembl_TRS":
                continue
            if "-" in acc:
                base, suffix = acc.split("-", 1)
                if suffix != "1":
                    continue
                acc = base
            if acc not in fasta_keys:
                continue
            tid = strip_version(val)
            if tid not in out or acc < out[tid]:
                out[tid] = acc
    return out


def load_gtf(path):
    """Return dict[transcript_id_unversioned] = {meta, exons[], cds[]}."""
    transcripts = {}
    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature not in ("transcript", "exon", "CDS"):
                continue
            attrs = parse_gtf_attrs(cols[8])
            tid_raw = attrs.get("transcript_id", "")
            tid = strip_version(tid_raw)
            if not tid:
                continue
            t = transcripts.setdefault(tid, {
                "gene_id": "",
                "gene_name": "",
                "chrom": "",
                "strand": "",
                "tags": set(),
                "biotype": "",
                "transcript_id_versioned": tid_raw,
                "exons": [],
                "cds": [],
            })
            if feature == "transcript":
                t["gene_id"] = strip_version(attrs.get("gene_id", ""))
                t["gene_name"] = attrs.get("gene_name", "")
                t["chrom"] = cols[0]
                t["strand"] = cols[6]
                t["tags"] = set(attrs.get("_tags", []))
                t["biotype"] = attrs.get("transcript_type", "")
                t["transcript_id_versioned"] = tid_raw
            elif feature == "exon":
                t["exons"].append({
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "exon_number": int(attrs.get("exon_number", "0")),
                    "exon_id": strip_version(attrs.get("exon_id", "")),
                })
            elif feature == "CDS":
                t["cds"].append({
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "phase": int(cols[7]) if cols[7] != "." else 0,
                })
    return transcripts


def pick_anchor_transcripts(transcripts, mane_transcripts):
    """Return dict[gene_id] = (transcript_id, anchor_source).

    Preference order: MANE_Select tag in GTF, MANE membership from summary file
    (covers GRCh37 lift where the GTF tag may be absent), then Ensembl_canonical.
    """
    by_gene = defaultdict(list)
    for tid, t in transcripts.items():
        if t.get("biotype") != "protein_coding":
            continue
        if not t["cds"]:
            continue
        by_gene[t["gene_id"]].append((tid, t))
    chosen = {}
    for gene_id, trs in by_gene.items():
        mane_tagged = [tid for tid, t in trs if "MANE_Select" in t["tags"]]
        if mane_tagged:
            chosen[gene_id] = (mane_tagged[0], "MANE_Select")
            continue
        mane_summary = [tid for tid, _ in trs if tid in mane_transcripts]
        if mane_summary:
            chosen[gene_id] = (mane_summary[0], "MANE_Select")
            continue
        canon = [tid for tid, t in trs if "Ensembl_canonical" in t["tags"]]
        if canon:
            chosen[gene_id] = (canon[0], "Ensembl_canonical")
    return chosen


def transcript_order(records, strand):
    return sorted(records, key=lambda r: r["start"], reverse=(strand == "-"))


def build_rows(transcripts, anchors, xref, fasta):
    rows = []
    skipped = defaultdict(int)
    aa_mismatches = []
    for gene_id, (tid, source) in anchors.items():
        t = transcripts[tid]
        uniprot = xref.get(tid)
        if not uniprot:
            skipped["no_uniprot_id"] += 1
            continue
        if uniprot not in fasta:
            skipped["uniprot_id_not_in_fasta"] += 1
            continue
        seq = fasta[uniprot]

        cds_sorted = transcript_order(t["cds"], t["strand"])
        cds_aa = {}
        cumul = 0
        for c in cds_sorted:
            length = c["end"] - c["start"] + 1
            aa_start = cumul // 3 + 1
            aa_end = (cumul + length - 1) // 3 + 1
            cds_aa[(c["start"], c["end"])] = (aa_start, aa_end)
            cumul += length
        # Sanity check: total CDS length / 3 should equal protein length (start codon included; stop excluded).
        expected = len(seq)
        observed = cumul // 3
        aa_consistent = observed == expected
        if not aa_consistent:
            aa_mismatches.append((tid, uniprot, observed, expected))

        exons_sorted = transcript_order(t["exons"], t["strand"])
        for idx, e in enumerate(exons_sorted, 1):
            cds_start = ""
            cds_end = ""
            aa_start = ""
            aa_end = ""
            for c in cds_sorted:
                if c["start"] >= e["start"] and c["end"] <= e["end"]:
                    cds_start = c["start"]
                    cds_end = c["end"]
                    aa_start, aa_end = cds_aa[(c["start"], c["end"])]
                    break
            rows.append({
                "gene_id": gene_id,
                "gene_name": t["gene_name"],
                "transcript_id": tid,
                "transcript_id_versioned": t["transcript_id_versioned"],
                "uniprot_acc": uniprot,
                "chrom": t["chrom"],
                "exon_start": e["start"],
                "exon_end": e["end"],
                "strand": t["strand"],
                "exon_number": e["exon_number"] or idx,
                "exon_id": e["exon_id"],
                "cds_start": cds_start,
                "cds_end": cds_end,
                "aa_start": aa_start,
                "aa_end": aa_end,
                "anchor_source": source,
                "aa_range_consistent": "T" if aa_consistent else "F",
                "uniprot_seq": seq,
            })
    return rows, dict(skipped), aa_mismatches


COLUMNS = [
    "gene_id", "gene_name", "transcript_id", "transcript_id_versioned",
    "uniprot_acc", "chrom", "exon_start", "exon_end", "strand",
    "exon_number", "exon_id",
    "cds_start", "cds_end", "aa_start", "aa_end",
    "anchor_source", "aa_range_consistent", "uniprot_seq",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (.gtf or .gtf.gz)")
    ap.add_argument("--mane", required=True, help="MANE summary (.txt or .txt.gz)")
    ap.add_argument("--fasta", required=True, help="UniProt canonical FASTA (.fasta or .fasta.gz)")
    ap.add_argument("--xref", required=True, help="Filtered UniProt Ensembl_TRS xref (3-col TSV, gz ok)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    log = lambda *a: print(*a, file=sys.stderr, flush=True)

    log("Loading UniProt FASTA…")
    fasta = load_fasta(args.fasta)
    log(f"  {len(fasta)} canonical sequences")

    log("Loading MANE summary…")
    mane = load_mane(args.mane)
    log(f"  {len(mane)} MANE Select transcripts")

    log("Loading Ensembl_TRS xref…")
    xref = load_uniprot_xref(args.xref, set(fasta.keys()))
    log(f"  {len(xref)} transcripts with canonical UniProt accession")

    log(f"Parsing GTF: {args.gtf}")
    transcripts = load_gtf(args.gtf)
    log(f"  {len(transcripts)} transcripts")

    log("Picking anchor transcripts…")
    anchors = pick_anchor_transcripts(transcripts, mane)
    src_counts = defaultdict(int)
    for _, s in anchors.values():
        src_counts[s] += 1
    log(f"  anchor sources: {dict(src_counts)}")

    log("Building rows…")
    rows, skipped, aa_mismatches = build_rows(transcripts, anchors, xref, fasta)
    log(f"  rows: {len(rows)}  skipped: {skipped}")
    if aa_mismatches:
        log(f"  AA-length mismatches (CDS_nt/3 != protein_len): {len(aa_mismatches)} transcripts")
        for tid, acc, obs, exp in aa_mismatches[:5]:
            log(f"    {tid} {acc}: CDS implies {obs} aa, UniProt has {exp} aa")
        if len(aa_mismatches) > 5:
            log(f"    … +{len(aa_mismatches) - 5} more")

    log(f"Writing {args.out}")
    with open(args.out, "w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")
    log("Done.")


if __name__ == "__main__":
    main()
