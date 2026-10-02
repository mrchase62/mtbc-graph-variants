#!/usr/bin/env python3
"""Rebuild the H37Rv repeat mask by measurement instead of by gene name.

WHY. The mask this replaces is a list of 181 intervals selected by name: 164
genes whose names match PE* or PPE*, plus 16 element copies. That is a proxy for
"short reads cannot be placed here", and it fails in both directions. It missed
Rv1759c/wag22 entirely -- a PE_PGRS-family gene named `wag22`, matching neither
pattern -- and four intervals inside its repeat tract came through the novel-
deletion screen as the four most widely shared candidates in the cohort, deleted
in 92 to 115 of 200 isolates. Measured directly, 26% of wag22's positions start a
9-mer that recurs three or more times within the gene, against 13% for PE_PGRS28
and 0% for core rpoB. A name-based list can only ever contain the cases someone
thought to name.

WHAT IS MEASURED. Two independent things go wrong for short reads, and they need
separate tests because neither implies the other.

  paralogy     the sequence occurs somewhere else in the genome, so reads from
               one copy pile onto another. Tested by k-mer uniqueness: a position
               is flagged when the k-mer starting there occurs more than once in
               the genome counting both strands. k defaults to 50, which is the
               scale of the anchor a 100-150 bp read needs to be placed
               confidently -- shorter k over-flags, longer k under-flags.
  tandem       the sequence repeats against ITSELF locally, so read placement is
               fine but the alignment slides and indel length is ambiguous. A
               PGRS tract is the type case. Tested by local 9-mer recurrence in a
               sliding window, which is the same statistic that exposed wag22.

A position failing either test is masked. The two are reported separately so the
reason for each interval is on the record, and because they warrant different
treatment downstream: paralogy corrupts depth, tandem repetition corrupts indel
coordinates.

WHAT IS NOT MEASURED, and is therefore carried over. The class labels. p4's
region() splits the mask into `pe_ppe` and `masked` and routes calls
differently, and that distinction is biological, not measurable from sequence
composition. So each measured interval is labelled by what it overlaps in the
gene annotation, and the pe_ppe label is preserved wherever a measured interval
lands on a PE/PPE-family gene. Nothing downstream loses a category.
"""
import argparse, bisect, collections, os, sys


COMP = str.maketrans("ACGTN", "TGCAN")


def read_one_fasta(path, want):
    seq, on = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                on = (name == want) or (want in name)
                continue
            if on:
                seq.append(line.strip())
    return "".join(seq).upper()


def nonunique_positions(seq, k):
    """Positions whose k-mer occurs more than once in the genome, both strands.

    Counted with two forward passes -- one over the sequence, one over its
    reverse complement -- rather than canonicalising each k-mer, which would
    mean 4.4 million reverse-complement operations in Python.
    """
    rc = seq.translate(COMP)[::-1]
    cf = collections.Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    cr = collections.Counter(rc[i:i + k] for i in range(len(rc) - k + 1))
    starts = bytearray(len(seq))
    n = 0
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if "N" in km:
            continue
        if cf[km] + cr[km] > 1:
            n += 1
            starts[i] = 1
    # A non-unique k-mer makes all k of its bases unplaceable, not just the
    # first. Dilating by marking k positions per hit is O(N*k); a sliding
    # window over the start flags is O(N) and gives the identical result.
    flag = bytearray(len(seq))
    run = 0
    for j in range(len(seq)):
        run += starts[j]
        if j >= k:
            run -= starts[j - k]
        if run:
            flag[j] = 1
    return flag, n


def rep_frac(s, k=9, minocc=3):
    if len(s) <= k:
        return 0.0
    c = collections.Counter(s[i:i + k] for i in range(len(s) - k + 1))
    return sum(v for v in c.values() if v >= minocc) / len(s)


def tandem_positions(seq, win, step, thresh):
    """Positions inside a window whose local 9-mer recurrence exceeds thresh."""
    flag = bytearray(len(seq))
    n_win = 0
    for t in range(0, len(seq), step):
        lo = max(0, t - (win - step) // 2)
        hi = min(len(seq), lo + win)
        if rep_frac(seq[lo:hi]) >= thresh:
            n_win += 1
            for j in range(t, min(len(seq), t + step)):
                flag[j] = 1
    return flag, n_win


def to_intervals(flag, min_len, close_gap):
    iv, s = [], None
    for i, v in enumerate(flag):
        if v and s is None:
            s = i
        elif not v and s is not None:
            iv.append([s, i])
            s = None
    if s is not None:
        iv.append([s, len(flag)])
    if close_gap:
        out = []
        for a, b in iv:
            if out and a - out[-1][1] <= close_gap:
                out[-1][1] = b
            else:
                out.append([a, b])
        iv = out
    return [x for x in iv if x[1] - x[0] >= min_len]


def load_genes(path):
    g = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 7 and f[4] == "Gene" and f[1].isdigit():
                g.append((int(f[1]), int(f[2]), f[6] or f[7]))
    return sorted(g)


def load_bed(path):
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1].isdigit():
                out.append((int(f[1]), int(f[2]),
                            f[3].split("|")[0] if len(f) > 3 else ""))
    return sorted(out)


def overlaps(items):
    st = [i[0] for i in items]

    def f(s, e):
        i = bisect.bisect_right(st, e)
        hit = []
        for j in range(max(0, i - 8), min(len(items), i + 1)):
            if min(items[j][1], e) > max(items[j][0], s):
                hit.append(items[j])
        return hit
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default=os.environ.get("MTB_H37RV", ""))
    ap.add_argument("--contig", default="NC_000962.3")
    ap.add_argument("--path-name", default="GCF_000195955#1#NC_000962.3",
                    help="the PanSN name written into the BED, matching the "
                         "existing mask so every consumer keeps working")
    ap.add_argument("--k", type=int, default=50)
    ap.add_argument("--tandem-window", type=int, default=300)
    ap.add_argument("--tandem-step", type=int, default=50)
    ap.add_argument("--tandem-thresh", type=float, default=0.10)
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument("--close-gap", type=int, default=50)
    ap.add_argument("--genes", default="data/annotation/H37Rv_snpeff_dump.txt")
    ap.add_argument("--elements", default="data/annotation/H37Rv_IS6110.pansn.gff")
    ap.add_argument("--old-mask",
                    default="data/annotation/H37Rv_repeat_mask.named.bed",
                    help="the ORIGINAL name-based list, kept as its own file so "
                         "that rebuilding the canonical mask is idempotent. "
                         "Pointing this at the output would union the mask with "
                         "itself on every rebuild.")
    ap.add_argument("--rds", default="data/annotation/known_RDs.bed")
    ap.add_argument("--gene-extend", type=float, default=0.25,
                    help="a gene with at least this measured-masked fraction is "
                         "masked in FULL. A tandem tract does not stop being "
                         "ambiguous at the tile where its 9-mer recurrence dips "
                         "below threshold: the whole repeat domain slides. "
                         "Without this, three of the four wag22 intervals were "
                         "caught and the fourth was not. 0 disables it.")
    ap.add_argument("--union-old", action="store_true",
                    help="keep every interval of --old-mask as well. Measurement "
                         "and naming do not subsume each other -- see the note "
                         "printed by the validation block -- so the shipped mask "
                         "is the union and each interval records its provenance.")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if not a.fasta:
        sys.exit("FATAL: pass --fasta or set MTB_H37RV")
    seq = read_one_fasta(a.fasta, a.contig)
    if not seq:
        sys.exit(f"FATAL: no sequence named {a.contig} in {a.fasta}")
    N = len(seq)
    print(f"  {a.contig}: {N:,} bp")

    print(f"  paralogy: counting {a.k}-mers on both strands ...")
    par, n_par = nonunique_positions(seq, a.k)
    print(f"    {sum(par):,} bp ({sum(par)/N:.2%}) lie in a non-unique "
          f"{a.k}-mer  [{n_par:,} such k-mers]")

    print(f"  tandem: 9-mer recurrence in a {a.tandem_window} bp window, "
          f"step {a.tandem_step}, threshold {a.tandem_thresh:.0%} ...")
    tan, n_win = tandem_positions(seq, a.tandem_window, a.tandem_step,
                                  a.tandem_thresh)
    print(f"    {sum(tan):,} bp ({sum(tan)/N:.2%}) in {n_win:,} flagged windows")

    both = bytearray(x | y for x, y in zip(par, tan))
    n_both = sum(1 for x, y in zip(par, tan) if x and y)
    print(f"  union {sum(both):,} bp ({sum(both)/N:.2%}); both tests agree on "
          f"{n_both:,} bp")

    ivs = to_intervals(both, a.min_len, a.close_gap)
    print(f"  -> {len(ivs):,} intervals at >= {a.min_len} bp "
          f"(gaps <= {a.close_gap} bp closed), "
          f"{sum(e-s for s,e in ivs):,} bp ({sum(e-s for s,e in ivs)/N:.2%})")

    # ---- classify, so p4's pe_ppe / masked split survives -----------------
    old_iv = load_bed(a.old_mask)
    gene_ov = overlaps(load_genes(a.genes))
    el = []
    if os.path.exists(a.elements):
        with open(a.elements) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) > 4 and f[3].isdigit():
                    el.append((int(f[3]), int(f[4]), "element"))
    el_ov = overlaps(sorted(el))

    # ---- gene extension --------------------------------------------------
    # Measurement locates the unreliable sequence; the gene gives its extent.
    if a.gene_extend > 0:
        allg = load_genes(a.genes)
        n_ext = 0
        add = []
        for gs_, ge_, gn_ in allg:
            span = ge_ - gs_
            if span < 200:
                continue
            m = sum(both[gs_:ge_])
            if m and m >= a.gene_extend * span:
                add.append((gs_, ge_, gn_))
                n_ext += 1
        for gs_, ge_, _ in add:
            for j in range(gs_, min(N, ge_)):
                both[j] = 1
        print(f"  gene extension: {n_ext} genes at >= {a.gene_extend:.0%} "
              f"measured-masked are masked in full "
              f"({', '.join(g[2] for g in sorted(add, key=lambda t:-(t[1]-t[0]))[:6])}"
              f"{' ...' if n_ext > 6 else ''})")
        ivs = to_intervals(both, a.min_len, a.close_gap)
        print(f"  -> {len(ivs):,} intervals, {sum(e-s for s,e in ivs):,} bp "
              f"({sum(e-s for s,e in ivs)/N:.2%})")

    rows = []
    cls_n = collections.Counter()
    for s, e in ivs:
        gs = [g[2] for g in gene_ov(s, e)]
        pe = [g for g in gs if g.startswith(("PE", "PPE"))]
        if el_ov(s, e):
            cls = "element"
        elif pe:
            cls = "pe_ppe"
        elif any(g.startswith(("rrs", "rrl", "rrf")) for g in gs):
            cls = "rrna"
        else:
            # which test fired decides the label for everything else
            npar = sum(par[s:e])
            ntan = sum(tan[s:e])
            cls = "paralog" if npar >= ntan else "tandem"
        why = []
        if sum(par[s:e]):
            why.append(f"paralog{100*sum(par[s:e])//(e-s)}")
        if sum(tan[s:e]):
            why.append(f"tandem{100*sum(tan[s:e])//(e-s)}")
        label = (pe[0] if pe else (gs[0] if gs else "intergenic"))
        cls_n[cls] += 1
        rows.append((s, e, f"{cls}|{label}|{'+'.join(why)}"))

    # ---- union with the named mask ---------------------------------------
    # A straight replacement would UN-mask 265,741 bp on H37Rv, because the big
    # PPE genes have unique 50-mers and are not locally repetitive: PPE56 spans
    # 11.5 kb and no measurement here flags it. The named list is not a failed
    # attempt at measuring placeability; it encodes a family judgement that
    # calls in these genes are unreliable for reasons sequence composition does
    # not show. So the two are combined and each interval says where it came
    # from.
    if a.union_old and old_iv:
        for s, e, nm in old_iv:
            base = nm.split("|")[0]
            cls = ("element" if base.startswith("IS")
                   else ("pe_ppe" if base.startswith(("PE", "PPE")) else "named"))
            rows.append((s, e, f"{cls}|{base}|named"))
            cls_n[cls + "(named)"] += 1
        rows.sort(key=lambda t: (t[0], t[1]))

    with open(a.out, "w") as fh:
        for s, e, nm in rows:
            fh.write(f"{a.path_name}\t{s}\t{e}\t{nm}\n")
    print(f"  classes: " + "  ".join(f"{k} {v}" for k, v in cls_n.most_common()))

    # ---- validation -------------------------------------------------------
    # EXACT coverage for the validation numbers. A windowed bisect over 214
    # intervals printed 31 masked RDs where the true count is 29; a bitmap
    # cannot miss an overlap, and the mask is small enough that it is free.
    covbm = bytearray(N)
    for s, e, _ in rows:
        for j in range(s, min(N, e)):
            covbm[j] = 1

    def cov_frac(s, e):
        return sum(covbm[s:e]) / max(1, e - s)
    old = old_iv
    print(f"\n  VALIDATION")
    if old:
        covered = sum(1 for s, e, nm in old if cov_frac(s, e) >= 0.5)
        print(f"    of the old mask's {len(old)} named intervals, {covered} "
              f"({covered/len(old):.0%}) are at least half covered")
        missing = [(s, e, nm) for s, e, nm in old if cov_frac(s, e) < 0.5]
        for s, e, nm in sorted(missing, key=lambda t: -(t[1]-t[0]))[:8]:
            print(f"      not recovered: {nm:<16}{s:>9}-{e:<9}{e-s:>6} bp")
    # the case that motivated this
    print(f"    Rv1759c/wag22 (1,989,832-1,992,577): "
          f"{sum(covbm[1989832:1992577]):,} of 2,745 bp masked"
          f"  <- the old mask covered 0")
    rds = load_bed(a.rds)
    if rds:
        oldbm = bytearray(N)
        for s, e, _ in old:
            for j in range(s, min(N, e)):
                oldbm[j] = 1
        hit = sum(1 for s, e, n in rds if cov_frac(s, e) >= 0.5)
        was = sum(1 for s, e, n in rds
                  if sum(oldbm[s:e]) >= 0.5 * (e - s))
        newly = [n for s, e, n in rds
                 if sum(oldbm[s:e]) < 0.5 * (e - s) <= sum(covbm[s:e])]
        print(f"    known RDs at least half masked: {hit} of {len(rds)}, "
              f"against {was} under the old mask")
        print(f"    RDs NEWLY buried by this mask: {len(newly)}"
              + (" -- " + ", ".join(newly[:6]) if newly else
                 "  <- nothing the old mask was not already hiding"))
    # ---- does p4's pe_ppe / masked routing survive the new name format? ---
    # p4_place.py splits the mask on `"PE" in name or "PPE" in name`, so the
    # class has to be readable by that test or every PE/PPE region silently
    # reroutes to `masked` and changes how those calls are placed.
    pe = sum(1 for _, _, nm in rows if "PE" in nm or "PPE" in nm)
    print(f"    p4 routing: {pe} intervals read as pe_ppe, "
          f"{len(rows)-pe} as masked")
    print(f"\n  -> {a.out}")


if __name__ == "__main__":
    main()
