#!/usr/bin/env python3
"""Per-contig support for accessory panel sequence, from unplaced reads.

Breadth matters more than read count. A handful of reads piling on one repetitive
stretch of a contig is not evidence the sample carries that sequence, whereas
reads spread across most of its length is. T10 scored both and the disagreement
between them is what separated real carriers from mismapping.
"""
import argparse, collections, csv, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", required=True, help="samtools depth -aa output")
    ap.add_argument("--idxstats", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # contig lengths come from idxstats, not from the depth file: a contig with
    # no coverage at all is the case we most need to report, and -aa is what
    # keeps it present. -a would drop it and make an unobserved contig look
    # absent from the panel instead of absent from the sample.
    length, reads = {}, {}
    for line in open(a.idxstats):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3 or f[0] == "*":
            continue
        length[f[0]] = int(f[1]); reads[f[0]] = int(f[2])

    covered = collections.Counter()
    total_depth = collections.Counter()
    for line in open(a.depth):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        d = int(f[2])
        if d > 0:
            covered[f[0]] += 1
        total_depth[f[0]] += d

    if not length:
        print("no contigs in idxstats -- refusing to emit an empty table",
              file=sys.stderr)
        return 1
    rows = []
    for c, L in sorted(length.items()):
        rows.append(dict(contig=c, len=L, reads=reads.get(c, 0),
                         covered_bp=covered[c],
                         breadth=round(covered[c] / L, 4) if L else 0.0,
                         mean_depth=round(total_depth[c] / L, 4) if L else 0.0))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    strong = sum(1 for r in rows if r["breadth"] >= 0.8)
    any_ = sum(1 for r in rows if r["reads"] > 0)
    print(f"  {len(rows)} contigs: {any_} with any read, {strong} with breadth >= 0.8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
