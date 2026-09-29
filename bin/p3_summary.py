#!/usr/bin/env python3
"""Collate P3: accessory presence/absence and copy-number calls across isolates.

No truth for real isolates, so this reports structure rather than accuracy. Two
things worth watching:

  - how far the unplaced-read fraction falls against a matched reference compared
    with H37Rv. T10 measured 39.5% of H37Rv-unmapped reads finding a home on the
    panel; against a matched reference the unplaced set should be far smaller,
    and a large one would mean the reference is not being used as intended.
  - whether copy-number calls concentrate in loci or in samples. A locus called
    gained in every isolate is more likely a reference artefact or a mapping
    hotspot than 23 independent gains.
"""
import argparse, collections, csv, os, sys


def rd(p):
    # every TSV this project writes is CRLF; csv handles it in text mode but be
    # explicit so a future shell-written input does not surprise us
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--dir", default="refbias/p3")
    ap.add_argument("--loci", default=None)
    ap.add_argument("--breadth", type=float, default=0.8)
    ap.add_argument("--out", default="refbias/p3/p3_summary.tsv")
    a = ap.parse_args()

    samples = [r["sample"] for r in rd(a.refmap)]
    rows, missing = [], []
    contig_hits = collections.Counter()
    locus_calls = collections.defaultdict(collections.Counter)
    for s in samples:
        acc = os.path.join(a.dir, f"{s}.accessory.tsv")
        cn = os.path.join(a.dir, f"{s}.copynumber.tsv")
        if not (os.path.exists(acc) and os.path.exists(cn)):
            missing.append(s); continue
        A, C = rd(acc), rd(cn)
        strong = [r for r in A if float(r["breadth"] or 0) >= a.breadth]
        anyread = [r for r in A if int(r["reads"] or 0) > 0]
        for r in strong:
            contig_hits[r["contig"]] += 1
        calls = collections.Counter(r["call"] for r in C)
        for r in C:
            locus_calls[r["locus"]][r["call"]] += 1
        rows.append(dict(
            sample=s, contigs_any=len(anyread), contigs_strong=len(strong),
            loci_total=len(C), projected=len(C) - calls["no_projection"],
            no_projection=calls["no_projection"], absent=calls["absent"],
            single=calls["single_or_reference"], gain_2x=calls["gain_2x"],
            gain_multi=calls["gain_multi"]))

    if not rows:
        print("no P3 output found", file=sys.stderr); return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"  {len(rows)} isolates"
          f"{'; MISSING: ' + ', '.join(missing) if missing else ''}\n")
    print(f"  {'isolate':<17s}{'contigs':>8s}{'strong':>7s}{'proj':>6s}"
          f"{'noproj':>7s}{'absent':>7s}{'single':>7s}{'gain2x':>7s}{'gainN':>6s}")
    for r in rows:
        print(f"  {r['sample']:<17s}{r['contigs_any']:>8d}{r['contigs_strong']:>7d}"
              f"{r['projected']:>6d}{r['no_projection']:>7d}{r['absent']:>7d}"
              f"{r['single']:>7d}{r['gain_2x']:>7d}{r['gain_multi']:>6d}")

    n = len(rows)
    print(f"\n  accessory contigs at breadth >= {a.breadth}: "
          f"{len(contig_hits)} distinct across {n} isolates")
    for c, k in contig_hits.most_common(8):
        print(f"    {c:<16s} in {k}/{n} isolates")
    gains = {L: c["gain_2x"] + c["gain_multi"] for L, c in locus_calls.items()}
    ubiq = sorted(((k, L) for L, k in gains.items() if k >= n - 1), reverse=True)
    print(f"\n  copy-number loci called gained in >= {n-1}/{n} isolates: {len(ubiq)}")
    for k, L in ubiq[:8]:
        print(f"    {L:<14s} gained in {k}/{n}  <- suspect reference artefact")
    som = sum(1 for L, g in gains.items() if 0 < g < n - 1)
    print(f"  loci gained in some but not all isolates: {som}"
          f"   (these are the informative ones)")
    npj = collections.Counter()
    for r in rows:
        npj[r["no_projection"]] += 1
    print(f"\n  no_projection per isolate: min {min(npj)}, max {max(npj)}"
          f"  -- loci with no equivalent position in that isolate's reference")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
