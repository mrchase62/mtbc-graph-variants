#!/usr/bin/env python3
"""Draw the new 1K cohort (user, 2026-10-13: "keep gwas1000 isolates, equal
250 per lineage"; L1-4 only, L7 dropped).

1. Pool: passes select_isolates.py (pass_standard.tsv), lineage 1-4, and the
   fold-back count (split_reads.awk, `inv_chimera_per_1k`, the "old" scale)
   <= 2.5.
2. Strata: lineage x RRDR (rpoB 761082-761162 in the isolate's TB-Profiler
   dr_variants at freq >= --min-freq), 125 per cell, as gwas1000 was
   stratified (96 per cell).
3. gwas1000 isolates in the pool are kept first, with their gwas1000 RRDR
   label.
4. Each cell is filled from the rest of the pool: candidates are shuffled
   (fixed seed) and taken round-robin over sublineages (`deepest`), so no
   sublineage dominates. RRDR is read from TB-Profiler JSON lazily, only for
   candidates looked at. The collection table has no study column, so there
   is no study cap.

    draw_1k.py --pass pass_standard.tsv --rates rates/*.tsv cohorts_in_use_rates.tsv \
        --gwas1000 runroot/refbias/gwas1000.phenotype.tsv --out cohort_1k.tsv
"""
import argparse, collections, csv, json, random, sys

TBP = ("/n/netscratch/sfortune_lab/Lab/pculviner/notebooks/260728_sv_exploration/"
       "pipeline_outputs/{0}/results/{0}.tbprofiler.json")
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]


def rrdr(sample, min_freq):
    try:
        j = json.load(open(TBP.format(sample)))
    except (OSError, ValueError):
        return None
    return int(any((v.get("gene_name") or v.get("gene")) == "rpoB"
                   and 761082 <= int(v["pos"]) <= 761162
                   and float(v.get("freq", 1)) >= min_freq
                   for v in j.get("dr_variants", [])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="pass_", required=True)
    ap.add_argument("--rates", nargs="+", required=True)
    ap.add_argument("--gwas1000", required=True)
    ap.add_argument("--per-cell", type=int, default=125)
    ap.add_argument("--max-foldback", type=float, default=2.5)
    ap.add_argument("--min-freq", type=float, default=0.9)
    ap.add_argument("--seed", type=int, default=20261013)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fb = {}
    for f in a.rates:
        for r in csv.DictReader(open(f), delimiter="\t"):
            fb[r["sample"]] = float(r["inv_chimera_per_1k"])
    pool, why = {}, collections.Counter()
    for r in csv.DictReader(open(a.pass_), delimiter="\t"):
        if r["lineage"] not in LINEAGES:
            continue
        if r["sample"] not in fb:
            why["no fold-back rate"] += 1
            continue
        if fb[r["sample"]] > a.max_foldback:
            why["fold-back > %g" % a.max_foldback] += 1
            continue
        pool[r["sample"]] = r
    print(f"pool: {len(pool)} lineage 1-4 isolates; excluded {dict(why)}", file=sys.stderr)

    g1k = {r["sample"]: r for r in csv.DictReader(open(a.gwas1000), delimiter="\t")
           if r["lineage"] in LINEAGES}
    cells = {(l, x): [] for l in LINEAGES for x in (0, 1)}
    for s, r in g1k.items():
        if s in pool:
            cells[r["lineage"], int(r["rrdr"])].append((s, "gwas1000"))
    kept = sum(map(len, cells.values()))
    print(f"gwas1000: {kept} of {len(g1k)} lineage 1-4 isolates kept", file=sys.stderr)

    rnd = random.Random(a.seed)
    for l in LINEAGES:
        by_sub = collections.defaultdict(list)
        for s, r in pool.items():
            if r["lineage"] == l and s not in g1k:
                by_sub[r["deepest"]].append(s)
        subs = sorted(by_sub)
        for k in subs:
            rnd.shuffle(by_sub[k])
        rnd.shuffle(subs)
        # round-robin over sublineages until both cells are full
        i = 0
        while any(len(cells[l, x]) < a.per_cell for x in (0, 1)) and any(by_sub.values()):
            k = subs[i % len(subs)]; i += 1
            if not by_sub[k]:
                continue
            s = by_sub[k].pop()
            x = rrdr(s, a.min_freq)
            if x is not None and len(cells[l, x]) < a.per_cell:
                cells[l, x].append((s, "new"))
        for x in (0, 1):
            if len(cells[l, x]) > a.per_cell:  # gwas1000 overfill: keep all
                pass
            print(f"  {l} rrdr={x}: {len(cells[l, x])} "
                  f"({sum(1 for _, o in cells[l, x] if o == 'gwas1000')} gwas1000)",
                  file=sys.stderr)

    with open(a.out, "w") as fo:
        fo.write("sample\tlineage\tdeepest\trrdr\tfoldback_per_1k\torigin\n")
        for (l, x), v in cells.items():
            for s, o in sorted(v):
                fo.write(f"{s}\t{l}\t{pool[s]['deepest']}\t{x}\t{fb[s]}\t{o}\n")
    print(f"written {a.out}: {sum(map(len, cells.values()))} isolates", file=sys.stderr)


if __name__ == "__main__":
    main()
