#!/usr/bin/env python3
"""Typed false calls per in-panel control isolate, in callable sequence.

An in-panel control's own genome is in the panel, so P1 should pick it and
there should be almost no events. A typed call (BND excluded) is false when
it is not at a true event (score.py's at_event) and does not touch
the uncallable mask. One row per isolate: counts per caller, and the call
types and sizes for our caller + assembly.
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from score import at_event, load_masks, masked, read_proto, read_vcf  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="phaseB/out/samples.tsv")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--asm-dir", required=True)
    ap.add_argument("--p2-dir", required=True)
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--mask-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ctl = [r for r in csv.DictReader(open(a.samples), delimiter="\t") if r["group"] == "in_panel"]
    truth = collections.defaultdict(list)
    for r in csv.DictReader(open(a.truth), delimiter="\t"):
        truth[r["sample"]].append(dict(type=r["type"], start=int(r["r_start"]), end=int(r["r_end"]),
                                       size=max(int(r["ref_len"]), int(r["alt_len"]))))
    masks = load_masks(a.refmap, a.mask_dir)
    refmap = {r["sample"]: r["reference"] for r in csv.DictReader(open(a.refmap), delimiter="\t")}
    rows = []
    for r in ctl:
        s = r["sample"]
        calls = dict(
            ours_asm=read_proto(os.path.join(a.asm_dir, f"{s}.events.tsv")),
            dysgu_pass=read_vcf(os.path.join(a.p2_dir, f"{s}.dysgu.vcf"), True),
            delly_pass=read_vcf(os.path.join(a.p2_dir, f"{s}.delly.vcf"), True))
        row = dict(sample=s, panel_genome=r["id"] and r["reference"], picked=refmap.get(s, ""),
                   picked_self=r["picked_self"], snp_distance=r["snp_distance"],
                   true_events=len(truth[s]))
        for k, cs in calls.items():
            false = [c for c in cs if c["type"] != "BND" and not masked(masks.get(s), c["start"], c["end"])
                     and not any(at_event(c, t) for t in truth[s])]
            row[k] = len(false)
            if k == "ours_asm":
                row["ours_types"] = ",".join(f"{t}:{n}" for t, n in
                                             collections.Counter(c["type"] for c in false).most_common())
                sz = sorted(c["end"] - c["start"] + 1 for c in false)
                row["ours_median_span"] = sz[len(sz) // 2] if sz else ""
        rows.append(row)
    rows.sort(key=lambda x: -x["dysgu_pass"])
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
