#!/usr/bin/env python3
"""Place P1i junction sites in the H37Rv frame from their flanks, both required to agree.

WHY NOT THE GRAPH, FOR THESE COORDINATES
`is6110_project_sites.py` placed the same sites with `odgi position` and the
result did not validate: 51% agreement with flanking sequence, 19% of sites more
than 5 kb from where P1h measured them directly, and the method contradicting
itself on 15% of sites out to 554 kb. `P1I_PROJECTION.md` section 4 has the
reason, and it is specific rather than general -- every coordinate in this arm
is an IS6110 junction, so it sits on graph nodes where sixteen copies collapse,
and a position on such a node has sixteen H37Rv images.

The flanks do not have that problem. They are ordinary chromosomal sequence, so
they place by alignment, and the failure mode becomes "this flank is not unique"
which is detectable rather than silent.

BOTH SIDES MUST AGREE, AND THE GAP BETWEEN THEM IS THE ANSWER
Each site gets a window on each side, placed independently. What matters is not
only that both place, but how far apart they land in H37Rv, because that
distance says whether H37Rv carries a copy at this locus:

    gap ~= 2G + 1          H37Rv has no element here -- the two flanks are
                           adjacent, so this is an insertion relative to H37Rv
    gap ~= 2G + 1 + 1355   H37Rv carries a copy here, and the flanks sit on
                           either side of it
    anything else          the two sides disagree. The site is DROPPED rather
                           than written at whichever coordinate one flank
                           happened to give

That third case is the point of the exercise. A site is reported only when two
independent placements are mutually consistent, so the yield is lower than the
graph's 100% and each survivor means something.

NOTE ON WHAT THE RIGHT FLANK IS. If the matched reference carries a copy at this
locus -- which 132 of 179 sites do, since they sit at an excision join -- the
right flank begins after that copy, not after the junction coordinate. The
crossmap records each excised span, so the element's extent in the reference is
read from it rather than assumed to be 1355 bp.

ACCESSORY SITES ARE NOT FORCED. A site in sequence H37Rv does not carry has no
H37Rv coordinate to find, and this script reports it as unplaceable rather than
attaching it to the nearest anchor. Those keep the graph node key that
`is6110_project_sites.py` gave them; see P1I_PROJECTION.md section 5.
"""
import argparse, collections, csv, hashlib, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from is6110_seam import Seams


def read_fasta_one(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                break
            name = line[1:].split()[0]
        else:
            buf.append(line.strip())
    return "".join(buf)


def load_excisions(path):
    """The reference's removed spans, as is6110_seam.Seams.

    The site's reference coordinate is never the span's first base. clean_to_orig
    maps a clean position to its original frame, and a junction sits at the seam,
    so a stack on the left of the seam reports orig_start - 1 and one on the
    right reports orig_end + 1. Keying the lookup on orig_start missed every
    single one, which put the right-hand window inside the element, where it
    matches sixteen copies and fails the mapping-quality floor. That is what
    drove 124 of 179 sites into "one_flank_unique" on the first run.
    """
    return Seams(path)


def cigar_ref_span(cig):
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cig)
    return sum(int(n) for n, o in ops if o in "MDN=X"), ops


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sites", default="is6110/results/p1i_sites_h37rv.tsv")
    ap.add_argument("--refs", default="refbias/build/7713a8d71d8e/refs")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched")
    ap.add_argument("--h37rv", default="refbias/build/7713a8d71d8e/refs/GCF_000195955.fasta")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--window", type=int, default=300)
    ap.add_argument("--gap", type=int, default=5,
                    help="bases of clearance between the window and the junction, "
                         "so a target-site duplication or a base of wobble does "
                         "not put element sequence inside the window")
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--empty-tol", type=int, default=10,
                    help="tolerance on the adjacent-flank gap")
    ap.add_argument("--occupied-tol", type=int, default=60,
                    help="tolerance on the element-sized gap, wider because "
                         "element length varies between copies")
    ap.add_argument("--element-len", type=int, default=1355)
    ap.add_argument("--ismapper-dir", default="refbias/p1f")
    ap.add_argument("--ism-window", type=int, default=50)
    ap.add_argument("--graph-fallback", action="store_true",
                    help="on a one_flank_unique verdict, accept the graph "
                         "projection when the one flank that placed "
                         "corroborates it. OFF BY DEFAULT: over the seven "
                         "truth genomes it takes recall from 0.966 to 1.000 "
                         "and combined precision from 0.848 to 0.784, and the "
                         "three extra placements cannot be separated from the "
                         "one correct one by tier, geometry or read count. "
                         "One of the three lands 1 bp from an annotated H37Rv "
                         "element and is more likely a real shared copy the "
                         "truth set could not place than a false positive, so "
                         "the precision cost is an upper bound.")
    ap.add_argument("--graph-tol", type=int, default=2000,
                    help="how far the graph position may sit from the one "
                         "flank that placed, before the fallback is refused")
    ap.add_argument("--workdir", default="refbias/p1i/flankplace")
    ap.add_argument("--out", default="is6110/results/p1i_sites_flank.tsv")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    sites = list(csv.DictReader(open(a.sites), delimiter="\t"))
    seqs, exc = {}, {}
    fa = os.path.join(a.workdir, "flanks.fa")
    n_written = 0
    # Each window is aligned once, named by a hash of its sequence, never by
    # its row. minimap2 seeds its tie-breaking with the read name, and a window
    # whose best hit is near-tied (a chimeric primary plus supplementary) gets
    # MAPQ 60 under one name and 1 under another. Named by row, adding one site
    # renumbered every later one and flipped the verdicts of unrelated sites in
    # other samples (5 of 11,524 on gwas1000). Named by sequence, a window's
    # verdict depends only on the window, and samples sharing a reference and
    # a site share one alignment.
    win = {}                                   # (row, side) -> window name
    seen = set()
    def put(fh, i, side, seq):
        h = hashlib.sha1(seq.encode()).hexdigest()[:20]
        win[(i, side)] = h
        if h not in seen:
            seen.add(h)
            fh.write(f">{h}\n{seq}\n")
    with open(fa, "w") as fh:
        for i, s in enumerate(sites):
            ref = s["reference"]
            if ref not in seqs:
                seqs[ref] = read_fasta_one(os.path.join(a.refs, f"{ref}.fasta"))
                exc[ref] = load_excisions(
                    os.path.join(a.crossmap_dir, f"{ref}.crossmap.tsv"))
            g = seqs[ref]
            p = int(s["r_pos"])
            # Step the window past the reference's own copy, on whichever side of
            # the seam this site sits. lstart/rend bracket the element; where the
            # reference carries no copy here both collapse to p.
            # The seam is found with the rule the VCF writer uses (is6110_seam),
            # so a site the writer calls ref_shared always has its windows
            # stepped past the reference's copy. Exact matching used to leave a
            # site 1-3 bp off the seam with one window inside the element.
            lstart, rend = p, p
            hit = exc[ref].nearest(p)
            if hit is not None:
                side, s0, e0, _ = hit
                if side == "L":
                    rend = p + (e0 - s0 + 1)  # element to the RIGHT of the site
                else:
                    lstart = p - (e0 - s0 + 1)  # element to the LEFT of the site
            l2, l1 = lstart - a.gap - 1, lstart - a.gap - a.window
            r1, r2 = rend + a.gap + 1, rend + a.gap + a.window
            s["_l1"], s["_r1"] = l1, r1
            if l1 < 1 or r2 > len(g):
                s["_skip"] = "at contig edge"
                continue
            put(fh, i, "L", g[l1-1:l2])
            put(fh, i, "R", g[r1-1:r2])
            n_written += 1
    print(f"  {len(sites)} sites, {n_written} with both windows extractable, "
          f"{len(seen)} distinct windows aligned")

    p = subprocess.run([a.minimap2, "-a", "-x", "sr", "--secondary=no",
                        a.h37rv, fa], capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"minimap2 failed:\n{p.stderr[-2000:]}")

    hit = {}
    for line in p.stdout.split("\n"):
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 0x904:                       # unmapped, secondary, supplementary
            continue
        if int(f[4]) < a.min_mapq:
            continue
        span, _ = cigar_ref_span(f[5])
        start = int(f[3])
        hit[f[0]] = (f[2], start, start + span - 1, bool(flag & 16))
    placed = {k: hit[h] for k, h in win.items() if h in hit}

    ism = collections.defaultdict(list)
    ism_ran = set()
    for s in {x["sample"] for x in sites}:
        tp = os.path.join(a.ismapper_dir, s, s, "IS6110",
                          f"{s}__NC_000962.3_table.txt")
        if os.path.exists(tp):
            ism_ran.add(s)
            for r in csv.DictReader(open(tp), delimiter="\t"):
                try:
                    ism[s].append((min(int(r["x"]), int(r["y"])),
                                   max(int(r["x"]), int(r["y"])), r.get("call", "")))
                except (ValueError, KeyError):
                    pass

    rows, tally = [], collections.Counter()
    expect_empty = 2 * a.gap + 1
    expect_occ = expect_empty + a.element_len
    for i, s in enumerate(sites):
        out = dict(sample=s["sample"], reference=s["reference"], tier=s["tier"],
                   geometry=s["geometry"], r_pos=s["r_pos"],
                   graph_placement=s["placement"], graph_h37rv=s["h37rv_pos"],
                   node=s["node"], node_offset=s["node_offset"])
        L, R = placed.get((i, "L")), placed.get((i, "R"))
        if s.get("_skip"):
            verdict = "window_at_contig_edge"
        elif L is None or R is None:
            verdict = ("neither_flank_unique" if L is None and R is None
                       else "one_flank_unique")
        elif L[3] != R[3]:
            verdict = "flanks_disagree_strand"
        else:
            rev = L[3]
            # in H37Rv, the bases lying BETWEEN the two placed windows
            gap = (L[1] - R[2] - 1) if rev else (R[1] - L[2] - 1)
            out["h37rv_gap"] = gap
            if abs(gap - expect_empty) <= a.empty_tol:
                verdict = "placed_h37rv_empty"
                out["h37rv_pos"] = (L[1] - a.gap - 1) if rev else (L[2] + a.gap + 1)
            elif abs(gap - expect_occ) <= a.occupied_tol:
                verdict = "placed_h37rv_occupied"
                out["h37rv_pos"] = (L[1] - a.gap - 1) if rev else (L[2] + a.gap + 1)
            else:
                verdict = "flanks_disagree_gap"
        # --- fallback: one flank placed, and the graph agrees with it ---------
        # A `one_flank_unique` verdict is not a contradiction, it is missing
        # evidence: one window landed in repetitive context and could not be
        # placed uniquely, so the gap test cannot run. The graph projection is
        # an independent placement that does not depend on either window being
        # unique, and it was measured at 74.5% agreement with ISMapper against
        # 95% for flank placement -- weaker, so it is accepted here only when
        # the one flank that DID place corroborates it, and it is labelled so
        # downstream can weigh it differently.
        #
        # Measured over the seven truth genomes: flank placement discards 27 of
        # 60 reconciled sites, and the graph had already placed 2 of those
        # correctly -- including the only site the whole contig-break arm
        # appeared to contribute, which turned out to be this filter dropping a
        # site the reads had found.
        if (verdict == "one_flank_unique" and a.graph_fallback
                and s["placement"] == "on_path" and s["h37rv_pos"]):
            gh = int(s["h37rv_pos"])
            anchor = L if L is not None else R
            lo, hi = min(anchor[1], anchor[2]), max(anchor[1], anchor[2])
            if lo - a.graph_tol <= gh <= hi + a.graph_tol:
                verdict = "placed_by_graph"
                out["h37rv_pos"] = gh
                out["h37rv_gap"] = ""
            else:
                tally["graph_fallback_rejected"] += 1

        out["verdict"] = verdict
        tally[verdict] += 1
        if verdict.startswith("placed"):
            h = out["h37rv_pos"]
            # AN ISOLATE ISMAPPER NEVER RAN ON IS NOT A DISAGREEMENT. Writing 0
            # for it makes "no comparison available" indistinguishable from
            # "compared and did not match", and over a cohort where ISMapper ran
            # on a minority that drags the headline rate to nonsense: 18% over
            # scale100 against 95% on the isolates that actually have a table.
            # Blank means not comparable, and the rates below skip it.
            out["ismapper"] = (
                int(any((x - a.ism_window) <= h <= (y + a.ism_window)
                        for x, y, _ in ism.get(s["sample"], [])))
                if s["sample"] in ism_ran else "")
        rows.append(out)

    cols = ["sample", "reference", "tier", "geometry", "r_pos", "verdict",
            "h37rv_pos", "h37rv_gap", "ismapper", "graph_placement",
            "graph_h37rv", "node", "node_offset"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    print(f"\n  {'verdict':26s} {'n':>5s}")
    print("  " + "-" * 32)
    for k, v in tally.most_common():
        print(f"  {k:26s} {v:5d}")
    print("  " + "-" * 32)
    ok = [r for r in rows if r["verdict"].startswith("placed")]
    print(f"  {'placed':26s} {len(ok):5d}  of {len(rows)} = {len(ok)/len(rows):.0%}")

    if ok:
        cmp_ = [r for r in ok if r["ismapper"] != ""]
        n_no = len(ok) - len(cmp_)
        if cmp_:
            m = sum(r["ismapper"] for r in cmp_)
            print(f"\n  ISMapper agreement on placed sites: {m}/{len(cmp_)} = "
                  f"{m/len(cmp_):.0%}, over the {len(ism_ran)} isolates ISMapper "
                  f"ran on")
            for t in sorted({r["tier"] for r in cmp_}):
                sub = [r for r in cmp_ if r["tier"] == t]
                y = sum(r["ismapper"] for r in sub)
                print(f"    tier {t}: {y}/{len(sub)} = {y/len(sub):.0%}")
        if n_no:
            print(f"  {n_no} placed sites are NOT comparable: ISMapper did not "
                  f"run on their isolate.\n  They are blank in the table and "
                  f"excluded above, not counted as disagreements.")
        emp = sum(1 for r in ok if r["verdict"] == "placed_h37rv_empty")
        print(f"\n  of the {len(ok)} placed: {emp} with H37Rv EMPTY at the locus "
              f"(an insertion relative to H37Rv),\n  {len(ok)-emp} with H37Rv "
              f"OCCUPIED (a copy H37Rv carries too)")
        g = [r for r in ok if r["graph_h37rv"]
             and abs(int(r["graph_h37rv"]) - r["h37rv_pos"]) <= 50]
        print(f"  graph agreed with the flank placement on {len(g)}/{len(ok)} "
              f"= {len(g)/len(ok):.0%} of these")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
