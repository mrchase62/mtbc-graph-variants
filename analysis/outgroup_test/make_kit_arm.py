#!/usr/bin/env python3
"""Arm K of the outgroup test: the canettii rows from direct alignment.

Analysis only. Takes arm G's alignment (canettii read from the graph) and
writes the same alignment as it would be if the canettii genomes were NOT
in the graph:

  * each canettii row is rebuilt from that genome's direct alignment to
    H37Rv (analysis/panel_checks/states/<acc>.npz, the "outgroup kit"):
    base == REF -> REF, == ALT -> ALT, anything else or 0 -> N; node-frame
    columns are N (no H37Rv coordinate);
  * columns where only canettii varies are dropped: without canettii in the
    graph they would not exist. A column is kept if the non-canettii taxa
    show both alleles.

    $MTB_PY make_kit_arm.py --aln G.fasta --sites G.sites.tsv \\
        --states ../panel_checks/states --out K.fasta --sites-out K.sites.tsv
"""
import argparse, collections, csv
import numpy as np

CANETTII = ("GCF_035581225", "GCF_000253375")
B = "NACGT"            # genome_states.py: 1-4 = A,C,G,T; 0 = N


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        out[name] = "".join(buf)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--states", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sites-out", required=True)
    a = ap.parse_args()

    seqs = read_fasta(a.aln)
    sites = list(csv.DictReader(open(a.sites), delimiter="\t"))
    n = len(sites)
    assert all(len(s) == n for s in seqs.values()), "alignment and sites differ"
    names = list(seqs)
    can = [t for t in CANETTII if t in seqs]

    rows = {t: list(seqs[t]) for t in names}
    stats = collections.Counter()
    for t in can:
        st = np.load(f"{a.states}/{t}.npz")["state"]
        for i, s in enumerate(sites):
            old = rows[t][i]
            if not s["chrom"].endswith("NC_000962.3"):
                new = "N"
            else:
                b = B[int(st[int(s["pos"]) - 1])]
                new = s["ref"] if b == s["ref"] else s["alt"] if b == s["alt"] else "N"
            rows[t][i] = new
            stats[(t, "same" if old == new else f"{old}->{new}")] += 1

    other = [t for t in names if t not in can]
    keep = []
    for i, s in enumerate(sites):
        seen = {rows[t][i] for t in other} - {"N"}
        if s["ref"] in seen and s["alt"] in seen:
            keep.append(i)
    with open(a.out, "w") as fh:
        for t in names:
            fh.write(f">{t}\n")
            s = "".join(rows[t][i] for i in keep)
            for j in range(0, len(s), 60):
                fh.write(s[j:j + 60] + "\n")
    with open(a.sites_out, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["column"] + list(sites[0])[1:] + ["g_column"])
        for k, i in enumerate(keep):
            w.writerow([k] + [sites[i][c] for c in list(sites[0])[1:]] + [i])
    print(f"  {len(names)} taxa; columns {n:,} -> {len(keep):,} "
          f"({n - len(keep):,} varied only through canettii)")
    for t in can:
        sub = {k[1]: v for k, v in stats.items() if k[0] == t}
        tot = sum(sub.values())
        print(f"  {t}: " + ", ".join(f"{k} {v:,}" for k, v in
                                     sorted(sub.items(), key=lambda kv: -kv[1]))
              + f"  (of {tot:,})")


if __name__ == "__main__":
    main()
