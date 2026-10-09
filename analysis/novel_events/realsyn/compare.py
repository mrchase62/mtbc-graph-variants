#!/usr/bin/env python3
"""Real vs synthetic reads for the 25 in-panel controls (run from runroot).

For each control and caller, typed calls in callable sequence (the
uncallable mask of the shared reference) from the real reads (marinB63 P2;
our caller + assembly from phaseB/out/v3b) and from synthetic reads simulated
from the control's complete assembly (marinSyn25 P2; our caller + assembly
from realsyn/out). Two calls are the same event when their types are
compatible (score.py's table) and they are near each other by score.py's rule
(start within 50 bp, or for deletions end within 50 bp or 50% reciprocal
overlap). real_only = calls the real reads make that the ideal reads of the
same genome do not: caused by the real reads (library artifacts such as
fold-back chimeras) or by differences between the reads' culture and the
assembly. syn_only = calls only the ideal reads make.

Small variants: HaplotypeCaller calls with a non-reference genotype, real vs
synthetic, by (position, ref, alt).

The fold-back rate is from ../controls/summary.tsv. Writes out/compare.tsv.
"""
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from score import COMPAT, load_masks, masked, near, read_proto, read_vcf  # noqa: E402

N = "../analysis/novel_events"
FAM = {"DEL": "DEL", "REPL": "DEL", "INS": "INS", "INV": "INV", "REARR": "REARR"}


def same(a, b):
    return (a["type"] in COMPAT.get(b["type"], ()) or b["type"] in COMPAT.get(a["type"], ())) \
        and (near(a, b) or near(b, a))


def dysgu_syn(y):
    p = f"refbias/marinSyn25/run/p2/{y}.dysgu.vcf"
    return p if os.path.exists(p) and os.path.getsize(p) else f"{N}/realsyn/out/dysgu/{y}.dysgu.vcf"


def small(path):
    out = set()
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 9 and c[9].split(":")[0] not in ("0", ".", "0/0", "./."):
                out |= {(c[1], c[3], a) for a in c[4].split(",") if a != "*"}
    return out


def main():
    ctl = [x.strip() for x in open(f"{N}/controls/controls.txt") if x.strip()]
    masks = load_masks("refbias/marinB63/run/p1/refmap.tsv", f"{N}/out/masks")
    fb = {r["sample"]: r for r in csv.DictReader(open(f"{N}/controls/summary.tsv"), delimiter="\t")}
    truth = {}
    for r in csv.DictReader(open(f"{N}/phaseB/out/truth.tsv"), delimiter="\t"):
        truth.setdefault(r["sample"], []).append(dict(type=r["type"], start=int(r["r_start"]), end=int(r["r_end"])))
    rows = []
    for s in ctl:
        y = "syn_" + s
        src = {
            "ours_asm": (read_proto(f"{N}/phaseB/out/v3b/asm/{s}.events.tsv"),
                         read_proto(f"{N}/realsyn/out/asm/{y}.events.tsv")),
            "dysgu_pass": (read_vcf(f"refbias/marinB63/run/p2/{s}.dysgu.vcf", True),
                           read_vcf(dysgu_syn(y), True)),
            "delly_pass": (read_vcf(f"refbias/marinB63/run/p2/{s}.delly.vcf", True),
                           read_vcf(f"refbias/marinSyn25/run/p2/{y}.delly.vcf", True)),
        }
        row = dict(sample=s, foldback_per_1k=fb[s]["inv_chimera_per_1k"], true_events=len(truth.get(s, [])))
        for k, (real, syn) in src.items():
            keep = lambda cs: [c for c in cs if c["type"] != "BND" and not masked(masks[s], c["start"], c["end"])]
            real, syn = keep(real), keep(syn)
            ro = [c for c in real if not any(same(c, d) for d in syn)]
            so = [c for c in syn if not any(same(c, d) for d in real)]
            row.update({f"{k}_real": len(real), f"{k}_syn": len(syn),
                        f"{k}_real_only": len(ro), f"{k}_syn_only": len(so)})
            if k == "ours_asm":
                row["ours_real_only_types"] = ",".join(
                    f"{t}:{sum(1 for c in ro if c['type'] == t)}" for t in sorted({c["type"] for c in ro}))
        r_sv, s_sv = small(f"refbias/marinB63/run/p2/{s}.vcf.gz"), small(f"refbias/marinSyn25/run/p2/{y}.vcf.gz")
        snp = lambda v: len(v[1]) == len(v[2]) == 1
        row.update(snv_real_only=sum(1 for v in r_sv - s_sv if snp(v)),
                   indel_real_only=sum(1 for v in r_sv - s_sv if not snp(v)),
                   snv_syn_only=sum(1 for v in s_sv - r_sv if snp(v)),
                   indel_syn_only=sum(1 for v in s_sv - r_sv if not snp(v)),
                   small_shared=len(r_sv & s_sv))
        rows.append(row)
    rows.sort(key=lambda r: -float(r["foldback_per_1k"]))
    os.makedirs(f"{N}/realsyn/out", exist_ok=True)
    with open(f"{N}/realsyn/out/compare.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
