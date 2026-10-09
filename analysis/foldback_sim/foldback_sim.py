#!/usr/bin/env python3
"""Replace a fraction of simulated read pairs with pairs from fold-back
chimeric fragments, to measure what these library artifacts do downstream.

A fold-back fragment is part A, read forward from the genome, joined to part
B, the reverse complement of nearby sequence, so the molecule doubles back on
the opposite strand (Zhang et al. 2024, BMC Genomics 25:227; Haile et al.
2019, NAR 47:e12). Two junction models:

  ir      the published mechanism (pairing of a partial single strand with an
          inverted repeat on the same molecule). The junction is at a k-mer
          (k=10) whose reverse complement occurs 20-300 bp upstream. A ends
          with that k-mer; B continues from the upstream copy on the other
          strand, so the two parts share the k-mer (microhomology). Sites are
          fixed by the genome, so chimeras recur at the same places.
  random  junction anywhere; B starts 0-200 bp upstream of the junction on the
          other strand, with no homology. Chimeras scatter.

Fragment lengths follow the base reads (normal, mean 450, sd 50). The
junction falls 30 bp or more from either end. Half the fragments are
reverse-complemented, so both strands are represented. Read 1 is the
fragment's first 150 bp, read 2 the reverse complement of its last 150 bp;
substitution errors at 0.2% per base, as in the base reads. The replaced
pairs are the first n in the file (wgsim writes pairs in random order), and
their names and quality strings are kept.

Writes the two FASTQs and a junction table (pair index, A end, B start,
strand, whether read 1 / read 2 crosses the junction).
"""
import argparse
import bisect
import random

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def ir_sites(g, k=10, lo=20, hi=300):
    """(i, j): k-mer at i whose reverse complement starts at j, lo <= i-j <= hi."""
    pos = {}
    for j in range(len(g) - k + 1):
        pos.setdefault(g[j:j + k], []).append(j)
    out = []
    for i in range(hi, len(g) - k + 1):
        km = g[i:i + k]
        if "N" in km:
            continue
        p = pos.get(rc(km))
        if not p:
            continue
        a = bisect.bisect_left(p, i - hi)
        b = bisect.bisect_right(p, i - lo)
        if a < b:
            out.append((i, p[b - 1]))  # nearest upstream copy
    return out


def mutate(s, rate, rng):
    s = list(s)
    for x in range(len(s)):
        if rng.random() < rate:
            s[x] = rng.choice([b for b in "ACGT" if b != s[x]])
    return "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--frac", type=float, required=True, help="fraction of pairs to replace")
    ap.add_argument("--model", choices=("ir", "random"), required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--readlen", type=int, default=150)
    ap.add_argument("--err", type=float, default=0.002)
    ap.add_argument("--out1", required=True)
    ap.add_argument("--out2", required=True)
    ap.add_argument("--junctions", required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    g = read_fasta(a.genome)
    L = a.readlen
    n_pairs = sum(1 for _ in open(a.r1)) // 4
    n = round(a.frac * n_pairs)
    sites = ir_sites(g) if a.model == "ir" and n else []
    if a.model == "ir" and n:
        print(f"inverted-repeat sites: {len(sites):,}")
    f1, f2 = open(a.r1), open(a.r2)
    o1, o2 = open(a.out1, "w"), open(a.out2, "w")
    jt = open(a.junctions, "w")
    jt.write("pair\tA_end\tB_start\tstrand\tr1_crosses\tr2_crosses\n")
    for p in range(n_pairs):
        h1, s1, _, q1 = (f1.readline().rstrip("\n") for _ in range(4))
        h2, s2, _, q2 = (f2.readline().rstrip("\n") for _ in range(4))
        if p < n:
            while True:
                F = max(2 * L, int(rng.gauss(450, 50)))
                la = rng.randint(30, F - 30)  # length of part A
                if a.model == "ir":
                    i, j = sites[rng.randrange(len(sites))]
                    a_end, b_start = i + 10, j  # A = g[a_end-la:a_end]; B = rc(g[b_start-lb:b_start])
                else:
                    a_end = rng.randint(1000, len(g) - 1000)
                    b_start = a_end - rng.randint(0, 200)
                lb = F - la
                if a_end - la < 0 or b_start - lb < 0:
                    continue
                frag = g[a_end - la:a_end] + rc(g[b_start - lb:b_start])
                if "N" in frag:
                    continue
                break
            strand = "+"
            if rng.random() < 0.5:
                frag, strand = rc(frag), "-"
                la = F - la
            s1 = mutate(frag[:L], a.err, rng)
            s2 = mutate(rc(frag)[:L], a.err, rng)
            jt.write(f"{p}\t{a_end}\t{b_start}\t{strand}\t{int(la < L)}\t{int(F - la < L)}\n")
        o1.write(f"{h1}\n{s1}\n+\n{q1[:len(s1)]}\n")
        o2.write(f"{h2}\n{s2}\n+\n{q2[:len(s2)]}\n")
    print(f"pairs {n_pairs:,}; replaced {n:,} ({a.frac:.3%}), model {a.model}")


if __name__ == "__main__":
    main()
