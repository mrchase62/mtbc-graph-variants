#!/usr/bin/env python3
"""One accessory catalogue, from the two that each hold half of it.

THE PROBLEM. Two passes built two views of the same objects and neither is
usable alone. `accessory_loci.tsv` has 802 loci with the classification that
matters -- polymorphic against reference_gap, allele counts, carrier fractions,
H37Rv coverage -- and sequence for only 128 of them. `insgt/assets/
insertions.tsv` has 336 sequence clusters with panel carrier sets and the
routing decision, keyed on its own cluster representatives. Joining them on
position showed 296 of the 336 clusters sit within 200 bp of a table locus, so
they describe the same things, and the split cost a real measurement: the
accessory census could only place variant-bearing nodes against loci that had
sequence, so it saw 112 of 336.

THE LOCUS IS THE UNIT, not the cluster. The table's 802 entries are the
authoritative locus definition; a sequence cluster can absorb records from
several positions and then no longer keys on any of them. So this walks the
loci and attaches, for each, what the graph and the other passes know.

SEQUENCE FOR AS MANY AS THE GRAPH SUPPLIES. The collapsed graph VCF is the
source `insertion_contigs.py` already reads. Taken per locus rather than per cluster it
provides the representative allele at that position directly, which is what the
census needs and what the 674 loci without sequence were missing.
"""
import argparse, subprocess, bisect, collections, csv, gzip, os, sys


def main():
    ap = argparse.ArgumentParser()
    # NO BUILD DEFAULTS. --loci and --graph-vcf were build 7713a8d71d8e's
    # and CX333's; the clusters, routing and census were insgt's CX333
    # tables and gwas1000's census, which a build of a new graph must not
    # inherit. P0 step catalogue passes the build's two files and "" for the
    # rest; docs/PANEL_TREE.md gives the CX333 reproduction with all five.
    ap.add_argument("--loci", required=True,
                    help="<build>/assets/accessory_loci.tsv")
    ap.add_argument("--graph-vcf", required=True,
                    help="<build>/assets/graph_collapsed.vcf.gz, the "
                         "COLLAPSED graph VCF (audit PGB-9). The most-"
                         "carried allele is chosen per record, so on the "
                         "decomposed file, where one allele's carriers are "
                         "split over several records, a split allele lost to "
                         "a less common one (26 of 802 CX333 loci)")
    ap.add_argument("--clusters", default="",
                    help="this graph's insgt clusters; empty: none")
    ap.add_argument("--routing", default="",
                    help="this graph's insgt routing; empty: none")
    ap.add_argument("--census", default="",
                    help="a cohort's accessory census; empty: none")
    ap.add_argument("--window", type=int, default=200,
                    help="how close a graph record or cluster must be to count "
                         "as the same locus")
    ap.add_argument("--min-len", type=int, default=50)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--h37rv", default="",
                    help="H37Rv FASTA; with --blastn, writes h37rv_cov95")
    ap.add_argument("--blastn", default="")
    a = ap.parse_args()
    # EMPTY MEANS NONE; A NAMED PATH MUST EXIST (the user's decision D43,
    # 2026-10-07). A mistyped path was read as absent, so the catalogue was
    # written without route or census columns and no error.
    missing = [f"--{k} {v}" for k, v in (("clusters", a.clusters),
               ("routing", a.routing), ("census", a.census))
               if v and not os.path.isfile(v)]
    if missing:
        sys.exit(f"FATAL: {'; '.join(missing)}: no such file. Pass an empty "
                 f"value for none.")

    loci = list(csv.DictReader(open(a.loci, newline=""), delimiter="\t"))
    print(f"  {len(loci):,} accessory loci")
    lpos = sorted((int(r["pos"]), r["locus_id"]) for r in loci)
    keys = [p for p, _ in lpos]

    def locus_at(p):
        i = bisect.bisect_left(keys, p - a.window)
        best, bd = None, None
        while i < len(keys) and keys[i] <= p + a.window:
            d = abs(keys[i] - p)
            if bd is None or d < bd:
                best, bd = lpos[i][1], d
            i += 1
        return best

    # ---- graph insertion records, grouped onto loci -----------------------
    op = gzip.open if a.graph_vcf.endswith(".gz") else open
    panel, byloc = [], collections.defaultdict(list)
    n_rec = 0
    with op(a.graph_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                panel = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if "," in alt or alt.startswith("<"):
                continue
            d = len(alt) - len(ref)
            if d < a.min_len:
                continue
            lid = locus_at(int(f[1]))
            if lid is None:
                continue
            n_rec += 1
            carriers = [panel[i] for i, g in enumerate(f[9:])
                        if g.split(":")[0] == "1"]
            byloc[lid].append((d, alt[len(ref):].upper(), carriers))
    print(f"  {n_rec:,} graph insertion records assigned to "
          f"{len(byloc):,} loci")

    clus = list(csv.DictReader(open(a.clusters, newline=""), delimiter="\t")) \
        if a.clusters else []
    route = {r["contig"]: r for r in csv.DictReader(
        open(a.routing, newline=""), delimiter="\t")} \
        if a.routing else {}
    cl_by_loc = collections.defaultdict(list)
    for r in clus:
        lid = locus_at(int(r["h37rv_pos"]))
        if lid:
            cl_by_loc[lid].append(r)
    cen = {r["locus"]: r for r in csv.DictReader(
        open(a.census, newline=""), delimiter="\t")} \
        if a.census else {}

    rows, n_seq = [], 0
    with open(a.out_fasta, "w") as fa:
        for r in loci:
            lid = r["locus_id"]
            recs = byloc.get(lid, [])
            # THE MOST-CARRIED ALLELE IS THE REPRESENTATIVE, not the longest.
            # Taking the longest was the same over-merge that insertion_contigs.py
            # already had to fix, one level up: at ACC_1761791 the longest allele
            # is 3,511 bp -- TbD1's 2,153 bp plus an adjacent 1,358 bp element
            # insertion -- carried by 2 panel genomes, while the 2,153 bp TbD1
            # form is carried by 93. Genotyping against the 3,511 bp form gave a
            # TbD1 carrier 0.61 coverage, which the 0.80 threshold then called
            # UNCERTAIN, and 2153/3511 = 0.613 is exactly that ratio.
            #
            # Carriers are still unioned across nested records of the SAME
            # length band, because those are one event described at several
            # bubble depths; alleles of substantially different length are
            # different events and the one most of the panel carries is the one
            # worth genotyping.
            seq, glen, gcarr = "", 0, set()
            if recs:
                best = max(recs, key=lambda t: (len(t[2]), t[0]))
                glen, seq = best[0], best[1]
                lo = 0.7 * glen
                for d, _, c in recs:
                    if lo <= d <= glen / 0.7:
                        gcarr |= set(c)
            if seq:
                n_seq += 1
                fa.write(f">{lid} pos={r['pos']} len={len(seq)} "
                         f"graph_carriers={len(gcarr)} class={r['klass']}\n")
                for i in range(0, len(seq), 60):
                    fa.write(seq[i:i + 60] + "\n")
            cls = cl_by_loc.get(lid, [])
            rt = route.get(cls[0]["contig"], {}) if cls else {}
            c = cen.get(lid, {})
            rows.append(dict(
                locus_id=lid, pos=r["pos"], klass=r["klass"],
                rep_len=r["rep_len"], n_alleles=r["n_alleles"],
                carriers_any=r["carriers_any"], carrier_frac=r["carrier_frac"],
                h37rv_cov=r["h37rv_cov"], novelty=r["novelty"],
                graph_records=len(recs), graph_len=glen,
                graph_carriers=len(gcarr),
                # the LIST, not only the count: the level-1 genotyper asks
                # whether a given sample's own reference is among them, which a
                # count cannot answer
                panel_carriers=",".join(sorted(gcarr)),
                has_sequence=int(bool(seq)),
                cluster=(cls[0]["contig"] if cls else ""),
                n_clusters=len(cls),
                route=rt.get("route", ""), h37rv_frac=rt.get("h37rv_frac", ""),
                element_frac=rt.get("element_frac", ""),
                variant_records=c.get("variant_records", ""),
                cells_alt=c.get("cells_alt", ""),
                cells_ref=c.get("cells_ref", "")))
    # H37RV_COV95 (review 2, R2-IS-1): the share of the catalogue sequence --
    # the sequence reads are aligned to -- that H37Rv carries at 95% identity
    # or more. h37rv_cov counts BLAST matches at any identity, so a locus 85%
    # identical to an H37Rv region read as "already in H37Rv" and was made
    # unmeasurable, though its reads are too divergent to place on H37Rv and
    # do reach the read route's pool (ACC_2867346, ACC_0334653: 1.0 by
    # h37rv_cov, 0.0 at 95%). accessory/bin/locus_presence.py judges
    # blindness on this column. BLAST, not minimap2: 30 CX333 catalogue sequences
    # are shorter than 200 bp, and minimap2's assembly presets do not align
    # them. Blank where the locus has no sequence.
    if a.h37rv and a.blastn:
        res = subprocess.run(
            [a.blastn, "-query", a.out_fasta, "-subject", a.h37rv, "-outfmt",
             "6 qseqid qstart qend qlen pident", "-evalue", "1e-10",
             "-dust", "no"], capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit(f"FATAL: blastn failed: {res.stderr.strip()}")
        spans, qlen = collections.defaultdict(list), {}
        for line in res.stdout.splitlines():
            q, qs, qe, ql, pid = line.split("\t")
            qlen[q] = int(ql)
            if float(pid) >= 95.0:
                spans[q].append(tuple(sorted((int(qs), int(qe)))))
        for r in rows:
            q = r["locus_id"]
            if not r["has_sequence"]:
                r["h37rv_cov95"] = ""
                continue
            tot, cs, ce = 0, None, None
            for lo, hi in sorted(spans.get(q, ())):
                if cs is None or lo > ce + 1:
                    if cs is not None:
                        tot += ce - cs + 1
                    cs, ce = lo, hi
                else:
                    ce = max(ce, hi)
            if cs is not None:
                tot += ce - cs + 1
            r["h37rv_cov95"] = round(tot / qlen[q], 4) if q in qlen else 0.0
        n95 = sum(1 for r in rows if r["h37rv_cov95"] != ""
                  and r["h37rv_cov95"] >= 0.9)
        print(f"  h37rv_cov95: {n95:,} loci at 0.9 or more (read route blind)")
    else:
        print("  note: no --h37rv/--blastn, so no h37rv_cov95 column; the "
              "presence step then judges blindness on h37rv_cov (any identity)")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"\n  MERGED CATALOGUE, {len(rows)} loci")
    print(f"    with sequence now             {n_seq:>5}   (was 128)")
    print(f"    with a graph insertion record {sum(1 for r in rows if r['graph_records']):>5}")
    print(f"    matched to a sequence cluster {sum(1 for r in rows if r['cluster']):>5}")
    print(f"    already carrying a census row {sum(1 for r in rows if r['variant_records'] != ''):>5}")
    rt = collections.Counter(r["route"] for r in rows if r["route"])
    print(f"    routed: " + "  ".join(f"{k}={v}" for k, v in rt.most_common()))
    kl = collections.Counter((r["klass"], bool(r["has_sequence"])) for r in rows)
    print(f"\n    {'class':<16}{'with seq':>10}{'without':>9}")
    for k in sorted({r["klass"] for r in rows}):
        print(f"    {k:<16}{kl[(k, True)]:>10}{kl[(k, False)]:>9}")
    print(f"\n  -> {a.out}\n  -> {a.out_fasta}")


if __name__ == "__main__":
    main()
