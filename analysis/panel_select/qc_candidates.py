#!/usr/bin/env python3
"""Panel QC: one row per candidate genome with its lineage, its QC flags and
whether it enters the selection pool. Analysis only.

Candidates are the 596 genomes with direct-alignment states
(analysis/panel_checks/genomes.tsv: CX333 333, Marin 90, GenBank 173).

Exclusion flags (any one excludes; bad assemblies are dropped, not repaired):
  homopolymer_run   a single-base run of >= --max-run bp (scanned here)
  indel_excess      1 bp indels vs H37Rv, sublineage robust z >= 5
                    (external_assemblies/qc/outliers.indel.tsv)
  private_hp_excess homopolymer indels in no other genome, same test
                    (outliers.hpindel.tsv)
  sr_disagree_*     hybrid vs its own short-read assembly (qc_summary.tsv)
  cx333_qc          CX333-script verdict not "pass" (qc_summary.tsv)
  genbank:<flags>   GenBank scan flags (Ns, BCG, indel/hp excess), and
                    PRJNA994284 (as in check 3)
  same_as_cx333     a Marin assembly of an isolate already in CX333
Notes only (not excluding): SNP excess/deficit within sublineage (may be
biology), is6110_truncated and ppe38_atypical (biology).

Pool: lineage 1-4 or 7 (barcode lineage) and no exclusion flag. Forced in:
H37Rv and the decided in-graph outgroups GCF_022870265 (L5.1.1) and
GCF_000253355 (L6.1.1), if they pass QC.

    qc_candidates.py --out candidates.tsv
"""
import argparse, csv, gzip, re

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
EA = f"{REPO}/analysis/external_assemblies"
POOL_LIN = ("lineage1", "lineage2", "lineage3", "lineage4", "lineage7")
FORCED = {"GCF_000195955": "reference", "GCF_022870265": "outgroup", "GCF_000253355": "outgroup"}


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def longest_run(path):
    op = gzip.open if path.endswith(".gz") else open
    best, last, n, where = 0, "", 0, ""
    with op(path, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                last, n = "", 0; continue
            for m in re.finditer(r"(.)\1*", l.strip().upper()):
                b, k = m.group(1), len(m.group(0))
                if b == last and m.start() == 0:
                    n += k
                else:
                    last, n = b, k
                if n > best:
                    best, where = n, b
    return best, where


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-run", type=int, default=100)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genomes = [l.rstrip("\n").split("\t") for l in open(f"{REPO}/analysis/panel_checks/genomes.tsv")]
    lin, flags, notes = {}, {}, {}
    for g, _, _ in genomes:
        flags[g], notes[g] = [], []
    for r in tsv(f"{EA}/qc/variants/variant_counts.tsv"):
        lin[r["id"]] = r["barcode_lineage"]
    for r in tsv(f"{EA}/panel_candidates.tsv"):
        lin[r["id"]] = r["lineage"]
    for r in tsv(f"{REPO}/analysis/genbank_scan/genbank_typed.tsv"):
        g = r["accession"]
        lin[g] = r["lineage"]
        if g in flags:
            if r["flags"]:
                flags[g].append("genbank:" + r["flags"])
            if r["bioproject"] == "PRJNA994284":
                flags[g].append("genbank:PRJNA994284")
    lin["GCF_000195955"] = "lineage4.9"
    for kind, name in (("indel", "indel_excess"), ("hpindel", "private_hp_excess"), ("snp", None)):
        for r in tsv(f"{EA}/qc/outliers.{kind}.tsv"):
            f = r["flag"].strip()
            if f and r["accession"] in flags:
                (flags if name else notes)[r["accession"]].append(name or f"snp_{f.lower()}")
    for r in tsv(f"{EA}/qc_summary.tsv"):
        g = r["id"]
        if g not in flags:
            continue
        for f in (r["lr_flags"] or "").split(";"):
            if f.startswith("sr_disagree"):
                flags[g].append(f)
            elif f.startswith(("is6110_truncated", "ppe38_atypical")):
                notes[g].append(f)
        if r["cx333_verdict"] and r["cx333_verdict"] != "pass":
            flags[g].append(f"cx333_qc:{r['cx333_verdict']}({r['cx333_reasons']})")
        if r["same_isolate_in_cx333"]:
            flags[g].append(f"same_as_cx333:{r['same_isolate_in_cx333']}")

    with open(a.out, "w") as fo:
        fo.write("id\tsource\tlineage\tlongest_run\tflags\tnotes\tin_pool\n")
        for g, path, src in genomes:
            n, b = longest_run(path)
            if n >= a.max_run:
                flags[g].append(f"homopolymer_run({n} {b})")
            # dedupe flags, keep order
            fl = list(dict.fromkeys(flags[g]))
            L = lin.get(g, "") or ""
            if fl:
                pool = "no"
            elif g in FORCED:
                pool = "forced_" + FORCED[g]
            elif L.startswith(POOL_LIN):
                pool = "yes"
            else:
                pool = "no"
            fo.write(f"{g}\t{src}\t{L}\t{n}{b}\t{';'.join(fl)}\t{';'.join(dict.fromkeys(notes[g]))}\t{pool}\n")


if __name__ == "__main__":
    main()
