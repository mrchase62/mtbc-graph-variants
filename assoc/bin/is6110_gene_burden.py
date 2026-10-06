#!/usr/bin/env python3
"""Association on insertion sites, collapsed to genes.

At site level the test is not merely underpowered, it is unidentifiable: a
single-origin insertion is collinear with the clade it defines, so a
convergence test cannot run and a mixed model correctly hands the signal to
population structure. On scale200 only 7 of 828 sites reach the two
independent origins the statistic needs.

Collapsing every site inside a gene into one burden term fixes that without
more isolates, because a gene accumulates origins even where each of its sites
arose once: 37 genes reach two origins at the same cohort size, and 21 reach
three. The unit of the test becomes "this gene was hit by the element
independently, in these isolates".

Both nulls are carried over from assoc_scan.py and mean the same thing. The
branch null draws k branches in proportion to edge length. The region null
draws them from the pooled event branches of the OTHER genes, leave-one-out,
so whatever is specific about where element insertions land -- and they do not
land uniformly -- is in the null rather than being discovered by the test.

Sites with no H37Rv coordinate, or falling in no annotated gene, cannot enter a
gene burden and are counted out loud rather than dropped quietly.
"""
import argparse, bisect, collections, csv, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sv_scatter import tree_frame
from assoc_scan import (descendant_matrix, bh, event_key, leave_one_out,
                        load_lineages, load_carriers)


def load_genes(path):
    """(start, end, strand, name), 1-based inclusive. The snpEff dump's start
    is 0-based (rpoB is `759806 763325`, true span 759807-763325), so reading
    it as 1-based gave every gene one extra base on its left."""
    g = []
    for r in csv.reader(open(path), delimiter="\t"):
        if len(r) > 7 and r[4] == "Gene":
            try:
                g.append((int(r[1]) + 1, int(r[2]), int(r[3]), r[6] or r[7]))
            except ValueError:
                pass
    g.sort()
    return g, [x[0] for x in g]


def deleted_span(r):
    """(first, last) deleted base, 1-based, for a deletion record; else None.

    VCF POS is the anchor base BEFORE the deletion, so it is never deleted.
    The span comes from the catalogue ID `svi:DEL:<first>:<len>`, else from a
    symbolic `<DEL:-len>` ALT."""
    f = (r.get("id") or "").split("#")[0].split(":")
    if len(f) == 4 and f[0] == "svi" and f[1] == "DEL" \
            and f[2].isdigit() and f[3].isdigit():
        return int(f[2]), int(f[2]) + int(f[3]) - 1
    alt = r.get("alt") or ""
    if alt.startswith("<DEL:-") and alt.endswith(">") \
            and alt[6:-1].isdigit() and r["pos"].isdigit():
        return int(r["pos"]) + 1, int(r["pos"]) + int(alt[6:-1])
    return None


def small_deleted_span(r):
    """(first, last) deleted base, 1-based, for a left-anchored small
    deletion (REF longer than ALT and ALT is REF's first base); else None.

    The same rule as deleted_span: the anchor base at POS is not deleted, so
    the record removes POS+1 .. POS+len(REF)-1. The `*` allele is not the
    record's ALT. Insertions have no deleted base and complex replacements no
    anchor, so both stay on the unit at POS."""
    ref = (r.get("ref") or "").upper()
    alts = [x for x in (r.get("alt") or "").upper().split(",") if x != "*"]
    if len(alts) != 1 or not r["pos"].isdigit():
        return None
    alt = alts[0]
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt and \
            not set(ref) - set("ACGTN"):
        return int(r["pos"]) + 1, int(r["pos"]) + len(ref) - 1
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="assoc/scale200/events")
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--genes", default="data/annotation/H37Rv_snpeff_dump.txt")
    ap.add_argument("--cls", default="is6110")
    ap.add_argument("--upstream", type=int, default=200,
                    help="an intergenic site within this many bp of a gene's "
                         "5' end becomes an `up:<gene>` promoter unit. Without "
                         "it 42.7%% of the insertion records with an origin are "
                         "dropped for falling in no gene -- and an element in "
                         "a promoter is one of the better documented ways it "
                         "changes a phenotype. Strand-aware: upstream of a "
                         "minus-strand gene is at HIGHER coordinates.")
    ap.add_argument("--min-origins", type=int, default=2)
    ap.add_argument("--permutations", type=int, default=20000)
    ap.add_argument("--lineages", default="",
                    help="TSV with `sample` and `lineage` columns. Adds a third "
                         "null that permutes the PHENOTYPE within lineage, "
                         "holding each lineage's carrier count at its observed "
                         "value, so between-lineage prevalence cannot generate "
                         "a hit. The region null cannot see that confounder: it "
                         "controls where events fall, not who carries the "
                         "phenotype.")
    ap.add_argument("--min-pool", type=int, default=200)
    ap.add_argument("--min-determinacy", type=float, default=0.80,
                    help="the scan's callability floor, applied per record: a "
                         "record whose reconstruction leaves more than this "
                         "share of branches undetermined contributes no "
                         "origins. Without it 6%% (small), 16-22%% (sv) and "
                         "7-9%% (is6110) of burden origins came from records "
                         "the scan treats as untestable.")
    ap.add_argument("--sv-min-overlap", type=int, default=1,
                    help="--cls sv, and small deletions under --cls small: a "
                         "deletion is credited to EVERY gene it "
                         "removes at least this many bp of. 1 credits any "
                         "partial overlap. A deletion that touches no gene "
                         "falls back to the promoter or intergenic unit at "
                         "its first deleted base.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260925)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    T = tree_frame(os.path.join(a.events, "labelled.nwk"))
    li = {l: j for j, l in enumerate(T["leaves"])}
    D = descendant_matrix(T)
    nd = D.sum(axis=1)
    carriers = load_carriers(a.phenotype, li)
    dep = np.zeros(len(T["leaves"]), dtype=bool)
    for s in carriers:
        dep[li[s]] = True
    ov = (D & dep).sum(axis=1) / np.maximum(nd, 1)
    OVP = None
    if not a.lineages:
        print("  WARNING: no --lineages, so no lineage null and nothing can "
              "pass all three")
    if a.lineages:
        lin = load_lineages(a.lineages, T["leaves"])
        groups = collections.defaultdict(list)
        for i in range(len(T["leaves"])):
            groups[lin.get(i, "unassigned")].append(i)
        DEP = np.zeros((len(T["leaves"]), a.permutations), dtype=np.float32)
        for g, idx in groups.items():
            idx = np.asarray(idx)
            k = int(dep[idx].sum())
            if k == 0:
                continue
            for pp in range(a.permutations):
                DEP[rng.choice(idx, size=k, replace=False), pp] = 1.0
        OVP = (D.astype(np.float32) @ DEP) / np.maximum(nd, 1)[:, None]
        print(f"  lineage null on: {len(groups)} strata, carrier counts held")
    el = np.asarray([(n.edge.length or 0.0) if n.edge is not None else 0.0
                     for n in T["nodes"]])
    w = el / el.sum() if el.sum() > 0 else np.full(T["n"], 1 / T["n"])
    print(f"  {int(dep.sum())} carriers of the phenotype over "
          f"{len(T['leaves'])} tips")

    genes, gstart = load_genes(a.genes)

    def gene_at(p):
        # ONE GENE PER POSITION: where genes overlap (17,874 bp of H37Rv) a
        # point record is credited to the earlier-starting gene only. Deletion
        # spans are credited to every gene they touch (genes_over below).
        i = bisect.bisect_right(gstart, p)
        for s, e, st, n in genes[max(0, i - 3):i]:
            if s <= p <= e:
                return n
        return None

    gmaxlen = max((e - s for s, e, st, n in genes), default=0)

    def genes_over(first, last):
        """every gene sharing at least --sv-min-overlap bp with [first, last]"""
        lo = bisect.bisect_left(gstart, first - gmaxlen)
        hi = bisect.bisect_right(gstart, last)
        return [n for s, e, st, n in genes[lo:hi]
                if min(e, last) - max(s, first) + 1 >= a.sv_min_overlap]

    def unit_at(p):
        """gene body, else promoter, else the intergenic gap itself."""
        g = gene_at(p)
        if g:
            return g, "gene"
        best = None
        for s, e, st, n in genes:
            d = (s - p) if st > 0 else (p - e)
            if 0 < d <= a.upstream and (best is None or d < best[0]):
                best = (d, n)
        if best:
            return f"up:{best[1]}", "promoter"
        # the gap between the nearest gene ending before p and the next
        # starting after it, so an intergenic site is still a named unit
        i = bisect.bisect_right(gstart, p)
        left = genes[i - 1][3] if i > 0 else "start"
        right = genes[i][3] if i < len(genes) else "end"
        return f"ig:{left}-{right}", "intergenic"

    from mtbvartools import CallBytestream
    ev = CallBytestream(os.path.join(a.events, "event"))
    cols = np.asarray(ev.calls.col)
    lab_i = {c: i for i, c in enumerate(cols)}
    by_gene = collections.defaultdict(set)
    reg_of = collections.defaultdict(collections.Counter)
    kind_of = {}
    n_kind = collections.Counter()
    n_site = n_nogene = n_noframe = n_uncall = n_multi = 0
    missing = []
    nbranch = T["n"] - 1
    for r in csv.DictReader(open(os.path.join(a.events, "variants.tsv")),
                            delimiter="\t"):
        if r["class"] != a.cls or int(r["n_gain"]) < 1:
            continue
        n_site += 1
        # SV records carry no FRAME tag -- merge_cohort_vcf.py writes CLASS,
        # SVTYPE and SVLEN but not FRAME -- while still sitting on H37Rv
        # coordinates. Requiring frame == h37rv therefore dropped all 4,145 of
        # them and reported zero testable units, which looks like a finding
        # and is a parsing bug. A blank frame with a numeric position on the
        # H37Rv contig is the H37Rv frame.
        if r["frame"] not in ("h37rv", "") or not r["pos"].isdigit():
            n_noframe += 1
            continue
        # THE SCAN'S CALLABILITY FLOOR. A record mostly undetermined on the
        # tree has not been reconstructed, and its few resolved gains are not
        # origins the burden can count. Burden records are H37Rv-frame and
        # unconditional, so the floor is over the whole tree, as in the scan.
        if 1.0 - int(r["n_undet"]) / max(1, nbranch) < a.min_determinacy:
            n_uncall += 1
            continue
        # A DELETION REMOVES A SPAN, NOT ITS ANCHOR BASE. Crediting it to the
        # unit at POS put a deletion of several genes on one of them, and could
        # put it on the gene ENDING at the anchor base, which it does not touch
        # at all (svi:DEL:3348474:59 and Rv2991). 21% of gwas1000's deletions
        # with an origin removed a gene other than the one credited.
        # The same holds for a SMALL deletion: crediting its anchor base put
        # 99 of scale200's 3,862 small deletions with an origin (gwas1000: 175
        # of 9,133) on a unit other than the genes whose bases they remove.
        span = (deleted_span(r) if a.cls == "sv" else
                small_deleted_span(r) if a.cls == "small" else None)
        if span:
            units = [(g, "gene") for g in genes_over(*span)]
            if not units:
                units = [unit_at(span[0])]
            n_multi += len(units) > 1
        else:
            units = [unit_at(int(r["pos"]))]
        for kind in {kd for _, kd in units}:
            n_kind[kind] += 1
        try:
            v = ev.calls.loc[event_key(r["row_key"])]
        except Exception:
            missing.append(r["row_key"])
            continue
        b = {lab_i[cols[j]] for j in np.where(v == 1)[0]}
        for g, kind in units:
            by_gene[g] |= b
            kind_of[g] = kind
            reg_of[g][r.get("region") or "other"] += 1
    ev.close()
    if missing:
        sys.exit(f"FATAL: {len(missing):,} variants.tsv rows have no event-"
                 f"matrix row, e.g. {missing[:3]}")
    print(f"  {n_site:,} {a.cls} records with at least one origin; "
          f"{n_noframe:,} have no H37Rv coordinate; {n_uncall:,} are below "
          f"the {a.min_determinacy:.0%} callability floor. The rest enter the "
          f"burden over {len(by_gene):,} units"
          + (f" ({n_multi:,} deletions span two or more genes)"
             if a.cls in ("sv", "small") else "") + ":")
    for k in ("gene", "promoter", "intergenic"):
        if n_kind[k]:
            nu = sum(1 for u, kk in kind_of.items() if kk == k)
            print(f"    {k:<12}{n_kind[k]:>6,} records over {nu:,} units")

    tested = {g: sorted(b) for g, b in by_gene.items()
              if len(b) >= a.min_origins}
    print(f"  {len(tested):,} genes reach >= {a.min_origins} independent "
          f"origins and are testable; {len(by_gene) - len(tested):,} do not")
    # STRATIFY. A gene unit's null must come from units like it: PE/PPE
    # carries a 13.0% homoplasy rate against core's 2.0%, and promoter and
    # intergenic units differ again from gene bodies, so one pooled null and
    # one BH family let the easy stratum set the threshold for the hard one.
    # This is the same correction the variant-level scan already applies.
    strat = {g: f"{reg_of[g].most_common(1)[0][0]}:{kind_of[g]}"
             for g in by_gene}
    pool = collections.defaultdict(list)
    for g, b in by_gene.items():
        pool[strat[g]].extend(b)
    print(f"  strata and their pooled event branches:")
    for k in sorted(pool):
        nt = sum(1 for g in tested if strat[g] == k)
        print(f"    {k:<22}{len(pool[k]):>8,} branches   {nt:>5,} testable "
              f"units" + ("" if len(pool[k]) >= a.min_pool
                          else f"   (below {a.min_pool}: no null)"))

    P, out = a.permutations, []
    for g, br in sorted(tested.items()):
        br = np.asarray(br)
        obs = float(ov[br].mean())
        k = len(br)
        d = rng.choice(T["n"], size=(P, k), p=w)
        p_branch = max(1, int((ov[d].mean(axis=1) >= obs).sum())) / P
        pl = leave_one_out(pool[strat[g]], br.tolist())
        if len(pl) >= a.min_pool:
            pa = np.asarray(pl)
            d2 = pa[rng.integers(0, len(pa), size=(P, k))]
            p_region = max(1, int((ov[d2].mean(axis=1) >= obs).sum())) / P
        else:
            p_region = float("nan")
        p_lineage = (max(1, int((OVP[br].mean(axis=0) >= obs).sum())) / P
                     if OVP is not None else float("nan"))
        carr = int((D[br].any(axis=0) & dep).sum())
        out.append(dict(gene=g, origins=k, obs=f"{obs:.4f}",
                        kind=kind_of.get(g, ""), stratum=strat[g],
                        carriers_with_phenotype=carr,
                        p_branch=p_branch, p_region=p_region,
                        p_lineage=p_lineage, null_pool=len(pl)))
    # BH within stratum, for the same reason the null is within stratum.
    for f, q in (("p_branch", "q_branch"), ("p_region", "q_region"),
                 ("p_lineage", "q_lineage")):
        for st in {o["stratum"] for o in out}:
            ok = [o for o in out if o["stratum"] == st and o[f] == o[f]]
            for o, v in zip(ok, bh([o[f] for o in ok])):
                o[q] = float(v)
    for o in out:
        o.setdefault("q_branch", float("nan"))
        o.setdefault("q_region", float("nan"))
        o.setdefault("q_lineage", float("nan"))
    out.sort(key=lambda o: (o["p_region"] if o["p_region"] == o["p_region"]
                            else 9, o["p_branch"]))
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "gene", "kind", "stratum", "origins", "obs",
            "carriers_with_phenotype",
            "p_branch", "q_branch", "p_region", "q_region",
            "p_lineage", "q_lineage", "null_pool"])
        wr.writeheader()
        wr.writerows(out)
    nb = sum(1 for o in out if o["q_branch"] < 0.05)
    nr = sum(1 for o in out if o["q_region"] == o["q_region"]
             and o["q_region"] < 0.05)
    havenull = any(o["p_region"] == o["p_region"] for o in out)
    nl = sum(1 for o in out if o["q_lineage"] == o["q_lineage"]
             and o["q_lineage"] < 0.05)
    allthree = [o for o in out if o["q_branch"] < 0.05
                and o["q_region"] == o["q_region"] and o["q_region"] < 0.05
                and o["q_lineage"] == o["q_lineage"] and o["q_lineage"] < 0.05]
    print(f"\n  q < 0.05: {nb} branch, "
          + (f"{nr} region" if havenull else "no region null")
          + f", {nl} lineage, {len(allthree)} under all three")
    print(f"  {'gene':<14}{'kind':<11}{'orig':>5}{'obs':>7}{'p_branch':>10}"
          f"{'p_region':>10}{'p_lineage':>11}")
    for o in (allthree or out[:10]):
        pr = ("%.5f" % o["p_region"]) if o["p_region"] == o["p_region"] else "-"
        plg = ("%.5f" % o["p_lineage"]) if o["p_lineage"] == o["p_lineage"] else "-"
        print(f"  {o['gene']:<14}{o['kind']:<11}{o['origins']:>5}"
              f"{float(o['obs']):>7.2f}{o['p_branch']:>10.5f}{pr:>10}{plg:>11}")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
