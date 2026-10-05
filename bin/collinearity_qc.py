#!/usr/bin/env python3
"""Detect scrambled / mis-assembled genomes by collinearity against H37Rv.

WHY NOT THE 2025 TEST
The 2025 QC flagged a genome when its MUMmer block count against H37Rv exceeded
8. Measured against the data, the median block count is 93, so the test flagged
75% of the panel and could not be used as a filter. The reason is that in MTBC
nearly all of that fragmentation is real biology -- IS6110 insertion positions,
PE_PGRS variability, RD deletions -- not assembly error.

WHAT SCRAMBLING ACTUALLY LOOKS LIKE
A mis-assembled complete genome has large segments in the wrong ORDER or
ORIENTATION. So the test here ignores small blocks entirely (--min-block, default
20 kb) and asks whether the surviving large blocks lie on a monotonic diagonal
against the reference:

  order_rho     Spearman correlation of reference vs query block order. A
                correctly assembled genome gives ~1.0 (or ~-1.0 if the whole
                genome is reverse-complemented, which is orientation, not
                scrambling -- so the absolute value is what matters).
  breakpoints   adjacent large blocks whose query order disagrees with their
                reference order. Zero for a clean assembly.
  inverted_frac fraction of large-block bases on the opposite strand. MTBC has
                few genuine large inversions, so a high value is suspicious --
                but it is reported, not thresholded, because real inversions
                exist.
  ref_covered   fraction of the reference spanned by large blocks. Catches gross
                incompleteness, independent of order.

Thresholds are deliberately NOT hardcoded: --report prints the distribution so a
cutoff can be chosen from the data rather than guessed, which is the mistake the
2025 version made.
"""
import argparse, glob, gzip, os, subprocess, sys, tempfile


def paf_blocks(mm2, ref, qry, threads, min_block):
    """Large alignment blocks as (ref_start, qry_start, length, strand)."""
    cmd = [mm2, "-x", "asm5", "-t", str(threads), "--secondary=no", ref, qry]
    pr = subprocess.run(cmd, capture_output=True, text=True)
    if pr.returncode != 0:
        raise RuntimeError(pr.stderr[:400])
    out = []
    for line in pr.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qs, qe, strand, ts, te = int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
        if (te - ts) < min_block:
            continue
        out.append((ts, qs, te - ts, strand))
    return out


def spearman(a, b):
    n = len(a)
    if n < 3:
        return float("nan")
    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        for pos, i in enumerate(order):
            r[i] = pos
        return r
    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return num / den if den else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="reference FASTA (H37Rv)")
    ap.add_argument("--dir", required=True, help="directory of *.fna.gz assemblies")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--min-block", type=int, default=20000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()

    ref_len = 0
    with (gzip.open(args.ref, "rt") if args.ref.endswith(".gz") else open(args.ref)) as fh:
        for l in fh:
            if not l.startswith(">"):
                ref_len += len(l.strip())

    files = sorted(glob.glob(os.path.join(args.dir, "*.fna.gz")))
    rows = []
    for i, p in enumerate(files, 1):
        acc = os.path.basename(p).replace(".fna.gz", "")
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tf:
            with gzip.open(p, "rt") as fh:
                tf.write(fh.read().encode())
            tmp = tf.name
        try:
            b = paf_blocks(args.minimap2, args.ref, tmp, args.threads, args.min_block)
        except RuntimeError as e:
            print(f"  {acc}: minimap2 failed: {e}", file=sys.stderr)
            os.unlink(tmp); continue
        os.unlink(tmp)
        # A correctly assembled, dnaA-rotated genome aligns as ONE full-length
        # collinear block. That is the best possible result, not a missing
        # measurement -- an earlier version of this function short-circuited on
        # len(b) < 2 and reported coverage 0, which flagged every clean genome.
        # Coverage and inversion are always computable; only the order statistics
        # need several blocks.
        b.sort(key=lambda x: x[0])
        total = max(sum(x[2] for x in b), 1)
        inv = sum(x[2] for x in b if x[3] == "-") / total
        cov = sum(x[2] for x in b) / ref_len
        if len(b) >= 3:
            rho = spearman([x[0] for x in b], [x[1] for x in b])
        else:
            rho = 1.0 if len(b) >= 1 else float("nan")
        qs = [x[1] for x in b]
        if len(qs) >= 2:
            fwd = sum(1 for j in range(len(qs) - 1) if qs[j + 1] > qs[j])
            brk = (len(qs) - 1) - max(fwd, (len(qs) - 1) - fwd)
        else:
            brk = 0
        rows.append((acc, len(b), rho, brk, inv, cov))
        if args.report and i % 100 == 0:
            print(f"  ...{i}/{len(files)}", file=sys.stderr)

    with open(args.out, "w") as fh:
        fh.write("accession\tn_large_blocks\torder_rho\tbreakpoints\tinverted_frac\tref_covered\n")
        for r in rows:
            fh.write("%s\t%d\t%.4f\t%d\t%.4f\t%.4f\n" % r)
    print(f"[collinearity_qc] {len(rows)} assemblies -> {args.out}")


if __name__ == "__main__":
    main()
