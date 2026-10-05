#!/usr/bin/env python3
"""Collapse a genome set into non-redundant representatives by pairwise SNP distance.

Panels assembled from public complete genomes are heavily clonally redundant --
outbreak sets and serial submissions of the same isolate are common -- and any
analysis that treats genomes as independent observations is inflated by them. On
this project lineage4.9 looked highly conserved for IS6110 until 34 of its 38
members turned out to sit in a single 50-SNP cluster.

Distances are pairwise-complete: a site missing in one genome is dropped for that
comparison only, rather than for every comparison. Clusters are single-linkage,
which is the conservative choice here -- it merges whenever any pair is within the
threshold, so it never overstates independence.

The representative of each cluster is the member with the most insertions >= 50 bp
relative to the reference, i.e. the structurally best-resolved assembly, since
that is what a pangenome graph is built from. Ties break on accession for
reproducibility.
"""
import argparse, collections, csv, subprocess, sys

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--bcftools", required=True)
    ap.add_argument("--thresholds", type=int, nargs="+", default=[5, 20, 50, 100, 200])
    ap.add_argument("--pick-threshold", type=int, default=50)
    ap.add_argument("--prefer", default="",
                    help="comma-separated accessions that MUST be the "
                         "representative of whichever cluster they fall in, "
                         "overriding --rank-by. For strains of standing interest "
                         "in their own right -- W-148, say -- a cluster-mate is a "
                         "fine phylogenetic stand-in but not a substitute for the "
                         "strain itself.")
    ap.add_argument("--rank-by", help="TSV: accession<TAB>score; higher is a better "
                                      "representative (e.g. insertions >=50bp)")
    ap.add_argument("--out")
    a = ap.parse_args()

    want = [l.strip() for l in open(a.accessions) if l.strip()]
    sm = subprocess.run([a.bcftools, "query", "-l", a.vcf],
                        capture_output=True, text=True).stdout.split()
    keep = [s for s in want if s in sm]
    missing = [s for s in want if s not in sm]
    if missing:
        print(f"  WARNING: {len(missing)} requested genomes are not in the VCF and "
              f"cannot be assessed: {', '.join(missing[:6])}"
              + (" ..." if len(missing) > 6 else ""))
    idx = [sm.index(s) for s in keep]
    pr = subprocess.run([a.bcftools, "query", "-f", "[%GT\t]\n", a.vcf],
                        capture_output=True, text=True)
    rows = []
    for line in pr.stdout.splitlines():
        f = line.rstrip("\t").split("\t")
        if len(f) != len(sm):
            continue
        rows.append([-1 if f[i] in (".", "./.", ".|.") else
                     (0 if f[i] == "0" else 1) for i in idx])
    A = np.array(rows, dtype=np.int8).T
    print(f"  {A.shape[0]} genomes x {A.shape[1]:,} SNP sites")

    ok = A >= 0
    n = A.shape[0]
    D = np.zeros((n, n), dtype=np.int32)
    for i in range(n):
        both = ok[i] & ok
        D[i] = ((A[i] != A) & both).sum(1)
    iu = np.triu_indices(n, 1)
    v = np.sort(D[iu])
    print(f"  pairwise distance: median {int(np.median(v))}, "
          f"min {v[0]}, max {v[-1]}")

    def clusters(t):
        par = list(range(n))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for i in range(n):
            for j in range(i + 1, n):
                if D[i, j] <= t:
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        par[ri] = rj
        g = collections.defaultdict(list)
        for i in range(n):
            g[find(i)].append(i)
        return list(g.values())

    print(f"\n  {'threshold':>10}{'clusters':>10}{'redundancy':>12}{'largest':>9}")
    res = {}
    for t in a.thresholds:
        cl = clusters(t)
        res[t] = cl
        print(f"  {t:>10}{len(cl):>10}{n/len(cl):>11.1f}x{max(len(c) for c in cl):>9}")

    score = {}
    if a.rank_by:
        for line in open(a.rank_by):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                try: score[f[0]] = float(f[1])
                except ValueError: pass

    prefer = {x.strip() for x in a.prefer.split(",") if x.strip()}
    if prefer:
        absent = prefer - set(keep)
        if absent:
            print(f"  NOTE: --prefer names {len(absent)} accession(s) not in this "
                  f"set: {', '.join(sorted(absent))}")

    def pick(members):
        forced = [m for m in members if m in prefer]
        if len(forced) > 1:
            raise SystemExit(f"--prefer names {len(forced)} genomes in one cluster "
                             f"({', '.join(forced)}); pick one")
        if forced:
            return forced[0]
        return max(members, key=lambda m: (score.get(m, 0), m))

    cl = res[a.pick_threshold]
    reps, out_rows = [], []
    for c in sorted(cl, key=len, reverse=True):
        members = sorted(keep[i] for i in c)
        rep = pick(members)
        reps.append(rep)
        for m in members:
            out_rows.append(dict(accession=m, cluster_size=len(c),
                                 representative=rep, is_representative=int(m == rep),
                                 rank_score=score.get(m, "")))
    reps.sort()
    print(f"\n  at <= {a.pick_threshold} SNPs: {len(reps)} non-redundant representatives "
          f"from {n} genomes")
    print(f"\n  clusters with more than one member:")
    for c in sorted(cl, key=len, reverse=True):
        if len(c) < 2:
            continue
        members = sorted(keep[i] for i in c)
        rep = pick(members)
        sub = D[np.ix_(c, c)][np.triu_indices(len(c), 1)]
        print(f"    {len(members):>3} genomes, max internal distance {sub.max():>4} SNPs"
              f"   rep {rep}")
        print(f"        {', '.join(members)}")
    single = sum(1 for c in cl if len(c) == 1)
    print(f"\n    plus {single} singleton genomes")
    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(out_rows)
        print(f"\n  wrote {a.out}")
        rp = a.out.rsplit(".", 1)[0] + ".representatives.txt"
        open(rp, "w").write("\n".join(reps) + "\n")
        print(f"  wrote {rp} ({len(reps)} accessions)")


if __name__ == "__main__":
    main()
