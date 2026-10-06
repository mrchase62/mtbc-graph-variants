#!/usr/bin/env python3
"""Project P1i junction sites from their matched reference into the H37Rv frame.

WHY THIS IS A PATH LOOKUP AND NOT AN ALIGNMENT
Each of the 18 SNP-matched references is a path in the 333-path pangenome graph,
so a coordinate on one converts to an H37Rv coordinate by walking the graph:
`odgi position -r GCF_000195955#1#NC_000962.3`. bin/p4_place.sh has done this
since the P4 arm, and t16_graph_frame.py records what it replaced -- a
`paftools liftover` route where 38.8% of records landed on an H37Rv base that
was not their REF allele.

Note that T16 ALSO found the graph frame failed, at ~13% wrong, but for a
different operation: composing inherited variant calls. That verdict is not
inherited here and is not contradicted here; this script converts coordinates
and does not compose anything.

SEQUENCE H37Rv DOES NOT CONTAIN IS THE POINT, NOT AN EDGE CASE
P1I_ACCESSORY_SITES.md measured 29 loci carrying element-mated reads that cannot
be placed on H37Rv at all. A projection that silently dropped them would discard
exactly the sites the matched reference was adopted to find. So every site gets
BOTH readings, following p4_place.py's convention exactly:

  odgi position -r <H37Rv path>   H37Rv coordinate, plus dist.to.ref
  odgi position -v                graph node id and offset

and `dist.to.ref` classifies, with p4_place.py's own three-way split:

  dist == 0                  on_path       a real H37Rv coordinate
  0 < dist <= --near-tol     off_path_near the interior of a small insertion in
                             R relative to H37Rv. p4_place.py section 11
                             measured 96.3% of off-path calls this way; they sit
                             in the same H37Rv ORF as their anchor
  dist > --near-tol          off_path_accessory  sequence H37Rv does not carry.
                             The H37Rv column is then an ANCHOR, not a position,
                             and the node key is the identity

KEY FOR FOLDING INTO A VCF. An on-path site keys on its H37Rv coordinate. An
off-path site must key on `node:<id>:<offset>`, because p4_place.py's reason
applies unchanged: two samples matched to different references carrying the same
accessory locus land on the same node exactly, where the H37Rv anchor agrees
only to within tens of bp. Writing an off-path site at its anchor coordinate
would merge distinct loci and place them at a base that is not theirs.

LABELLING NOTE, 2026-09-21.  This script scores on the A/B/C tiers, which have
been replaced by two independent columns -- `evidence` (two_sided/one_sided, a
quality ordering) and `site_class` (ref_shared/ref_lacking, not an ordering).
The rule now lives in is6110_promote_sites.py and is described in
is6110/docs/PROMOTED_TIER.md.  `tierAB` below means `evidence = two_sided`, and
it excludes promoted one-sided calls, so any ratio it forms against a
depth-based estimate reads low.  The numbers this script produced are still the
numbers that were measured; do not read tier A as more confident than tier B.
"""
import argparse, collections, csv, importlib.util, os, subprocess, sys

# THE COORDINATES GO OUT AND COME BACK IN THE PANEL FRAME. The graph was built
# from the dnaA-rotated panel FASTA; the crossmap positions below come from
# alignment to refbias/build/<id>/refs/<R>.fasta. For 9 of the 18 references
# this arm uses those are the same sequence written from a different origin,
# and for one also on the other strand, so a refs coordinate handed to
# `odgi position` unconverted queries a base 336 bp or 952 kb away. That is
# what P1I_PROJECTION.md originally blamed on collapsed repeat nodes; the
# measurement is in graphframe/docs/GRAPH_FRAME_RESOLUTION.md.
_gfs = importlib.util.spec_from_file_location(
    "graph_frame", os.path.join("graphframe", "bin", "graph_frame.py"))
graph_frame = importlib.util.module_from_spec(_gfs)
_gfs.loader.exec_module(graph_frame)
FRAMES = graph_frame.Frames()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "bin"))
import mtb_norm  # noqa: E402


def read_paths(og, odgi):
    p = subprocess.run([odgi, "paths", "-i", og, "-L"], capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"odgi paths failed:\n{p.stderr[-2000:]}")
    out = {}
    for line in p.stdout.split("\n"):
        line = line.strip()
        if line:
            out[line.split("#")[0]] = line
    return out


def run_position(odgi, og, posfile, args, nt):
    cmd = [odgi, "position", "-i", og, "-F", posfile, "-t", str(nt)] + args
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"odgi position failed:\n{p.stderr[-3000:]}")
    return p.stdout


def parse_position(text, want_node):
    """Same field layout p4_place.py reads: source in col 1, target in col 2,
    dist.to.ref in col 3. Keyed on (path, 0-based REFS pos) so several samples
    sharing a reference position collapse to one query and fan back out -- the
    source coordinate odgi echoes back is in the panel frame and is converted
    here, so nothing outside this function sees a panel coordinate."""
    out = {}
    for line in text.split("\n"):
        if not line or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        try:
            sp = f[0].rsplit(",", 2)
            sacc = sp[0].split("#")[0]
            src = (sp[0], FRAMES.to_refs(sacc, int(sp[1])))
            if want_node:
                # the offset runs along the source's walk; main() restates
                # it as the node's forward offset (mtb_norm.forward_offset)
                nid, off, nst = f[1].split(",")
                val = (int(nid), int(off), nst.strip())
            else:
                tp = f[1].rsplit(",", 2)
                tacc = tp[0].split("#")[0]
                # the refs-to-refs relation from the storage flip, and then
                # odgi's own flag (f[3] here: this output has not been
                # through frame_convert.py). Where it is `-` the source walks
                # the node opposite to the target, and odgi reports a target
                # one PANEL base high and complemented: measured at 8,684 of
                # 8,694 such GCF_000193185 positions (bin/p4_place.py's
                # parse_pos_file). The earlier note here, that odgi already
                # accounts for the inverted step, rested on 4 pilot SNPs.
                st = "-" if FRAMES.flipped(sacc) != FRAMES.flipped(tacc) else "+"
                tpos = int(tp[1])
                if len(f) > 3 and f[3].strip() == "-":
                    tpos -= 1
                    st = "-" if st == "+" else "+"
                val = (FRAMES.to_refs(tacc, tpos) + 1,
                       int(f[2]) if len(f) > 2 else 0, st)
        except (ValueError, IndexError):
            continue
        out.setdefault(src, val)
    return out


def node_occurrences(odgi, og, want, nt):
    """{(path name, node id): times that path visits the node} for `want`.

    A node key is an identity only where the node occurs once in the path:
    node 46966 is a 1 bp node every path walks 85 to 453 times, and keying on
    it merged insertions more than 1 Mb apart into one record (audit P3IS-3).
    `odgi paths -H` is the path-by-node visit count matrix; only the columns
    and rows asked for are kept."""
    p = subprocess.Popen([odgi, "paths", "-i", og, "-H", "-t", str(nt)],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    by_path = collections.defaultdict(set)
    for pth, nid in want:
        by_path[pth].add(nid)
    cols, out = None, {}
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if cols is None:
            idx = {h[5:]: i for i, h in enumerate(f) if h.startswith("node.")}
            cols = {nid: idx[str(nid)] for _, nid in want if str(nid) in idx}
            continue
        for nid in by_path.get(f[0], ()):
            if nid in cols:
                out[(f[0], nid)] = int(f[cols[nid]])
    err = p.stderr.read()
    if p.wait() != 0:
        sys.exit(f"odgi paths -H failed:\n{err[-2000:]}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reconcile", default="is6110/results/p1i_reconcile.tsv")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched")
    ap.add_argument("--graph", default=None)
    ap.add_argument("--odgi", default=os.environ.get("MTB_ODGI", "odgi"))
    ap.add_argument("--h37rv-path", default="GCF_000195955#1#NC_000962.3")
    # NO DEFAULT DIRECTORY. The default was the PILOT's refbias/p1f, so a
    # cohort run that did not pass this joined against the pilot isolates'
    # tables and wrote 0 ("no ISMapper region near") for every other isolate,
    # which is "not compared" written as "compared and disagreed". Empty turns
    # the join off; an isolate without a table in the given directory is blank.
    ap.add_argument("--ismapper-dir", default="",
                    help="this cohort's ISMapper output directory "
                         "(<dir>/<sample>/<sample>/IS6110/...); empty, the "
                         "default, turns the ISMapper join off")
    ap.add_argument("--ism-window", type=int, default=50)
    ap.add_argument("--near-tol", type=int, default=50,
                    help="dist.to.ref at or below this is a small-insertion "
                         "interior rather than accessory sequence; p4_place.py's "
                         "default and its reasoning")
    ap.add_argument("--all-stacks", action="store_true",
                    help="project every stack, not only tier A+B")
    ap.add_argument("--workdir", default="refbias/p1i/project")
    ap.add_argument("--out", default="is6110/results/p1i_sites_h37rv.tsv")
    ap.add_argument("--threads", type=int, default=4)
    # The build's node table, for node lengths: an off-path site's key is
    # node:<id>:<forward offset>, and for a carrier that walks the node in
    # reverse the forward offset is L-1-offset. Default from MTB_BUILD_DIR
    # (refbias_run.sh exports it; P0 step nodes writes the table).
    _b = os.environ.get("MTB_BUILD_DIR", "")
    ap.add_argument("--node-lengths",
                    default=os.path.join(_b, "assets", "node_positions.tsv") if _b else "",
                    help="<build>/assets/node_positions.tsv; default from "
                         "MTB_BUILD_DIR")
    a = ap.parse_args()

    # the build's graph, as its stamp records it; the fallback was a glob
    # over graphs/CX333..., which on a new build projected onto the old graph
    og = a.graph
    if og is None and _b and os.path.exists(os.path.join(_b, "build_info.tsv")):
        og = next((l.rstrip("\n").split("\t")[1]
                   for l in open(os.path.join(_b, "build_info.tsv"))
                   if l.startswith("graph\t")), None)
    if og is None:
        sys.exit("no graph: pass --graph, or set MTB_BUILD_DIR")
    os.makedirs(a.workdir, exist_ok=True)

    ref_of = {r["sample"]: r["reference"] for r in
              csv.DictReader(open(a.refmap), delimiter="\t")}
    paths = read_paths(og, a.odgi)

    sites = []
    for r in csv.DictReader(open(a.reconcile), delimiter="\t"):
        tier = ("A" if r["geometry"] == "tsd" else
                "B" if r["geometry"] == "two_sided_wide" and r["chrom_side"] == "agreed"
                else None)
        if tier is None and not a.all_stacks:
            continue
        if not r["geometry"]:
            continue                      # chrom_side_only rows have no stack
        ref = ref_of.get(r["sample"])
        path = paths.get(ref)
        if path is None:
            print(f"  {r['sample']}: {ref} is not a path in the graph, skipped",
                  file=sys.stderr)
            continue
        sites.append(dict(sample=r["sample"], reference=ref, path=path,
                          tier=tier or r["geometry"], geometry=r["geometry"],
                          clean_pos=int(r["clean_pos"]), r_pos=int(r["orig_pos"]),
                          reads=r["reads"], reads_q=r["reads_q"],
                          chrom_side=r["chrom_side"]))
    if not sites:
        sys.exit("no sites selected")

    # one odgi query per distinct (path, position); several samples may share one
    keys = sorted({(s["path"], s["r_pos"] - 1) for s in sites})
    pf = os.path.join(a.workdir, "query.pos")
    with open(pf, "w") as fh:
        for p, q in keys:
            fh.write(f"{p},{FRAMES.to_panel(p.split('#')[0], q)},+\n")
    print(f"  {len(sites)} sites, {len(keys)} distinct reference positions")

    href = parse_position(run_position(a.odgi, og, pf, ["-r", a.h37rv_path],
                                       a.threads), want_node=False)
    node = parse_position(run_position(a.odgi, og, pf, ["-v"], a.threads),
                          want_node=True)
    print(f"  odgi returned {len(href)} H37Rv projections and {len(node)} node "
          f"positions for {len(keys)} queries")
    # how often the carrier's own path visits each node it landed on; the
    # writer keys on a node only where this is 1 (audit P3IS-3)
    occ = node_occurrences(a.odgi, og, {(k[0], v[0]) for k, v in node.items()},
                           a.threads)
    # node lengths for the reverse-walked nodes only; a missing table is
    # refused where one is needed, never read as forward
    want_len = {str(v[0]) for k, v in node.items()
                if v[2] == "-" and href.get(k) and href[k][1] != 0}
    nlen = {}
    if want_len:
        if not a.node_lengths or not os.path.exists(a.node_lengths):
            sys.exit(f"FATAL: {len(want_len)} nodes are walked in reverse and "
                     f"need their length for a forward offset; no node table "
                     f"at '{a.node_lengths}' (pass --node-lengths "
                     f"<build>/assets/node_positions.tsv)")
        nlen = mtb_norm.load_node_lengths(a.node_lengths, want_len)

    ism = collections.defaultdict(list)
    ism_ran = set()
    for s in ({x["sample"] for x in sites} if a.ismapper_dir else ()):
        tp = os.path.join(a.ismapper_dir, s, s, "IS6110",
                          f"{s}__NC_000962.3_table.txt")
        if os.path.exists(tp):
            ism_ran.add(s)
            for r in csv.DictReader(open(tp), delimiter="\t"):
                try:
                    ism[s].append((int(r["x"]), int(r["y"]), r.get("call", "")))
                except (ValueError, KeyError):
                    pass

    rows, tally = [], collections.Counter()
    for s in sites:
        k = (s["path"], s["r_pos"] - 1)
        hp, np_ = href.get(k), node.get(k)
        if hp is None or np_ is None:
            tally["unprojected"] += 1
            s.update(placement="unprojected", h37rv_pos="", dist_to_ref="",
                     node="", node_offset="", node_occ="", frame_strand="",
                     key="", ismapper="", ism_call="")
            rows.append(s); continue
        h, dist, strand = hp
        nid = np_[0]
        # ONE KEY PER BASE: the node's forward offset, not odgi's offset along
        # this carrier's walk, which differs (L-1-off) for a carrier walking
        # the node in reverse (gwas1000: 71 of 19,990 keyed nodes had sites
        # from both directions)
        noff = mtb_norm.forward_offset(np_[1], np_[2], nlen.get(str(nid)))
        if noff is None and dist != 0:
            tally["unprojected"] += 1
            tally["no_node_length"] += 1
            s.update(placement="unprojected", h37rv_pos="", dist_to_ref="",
                     node="", node_offset="", node_occ="", frame_strand="",
                     key="", ismapper="", ism_call="")
            rows.append(s); continue
        if noff is None:
            noff = ""         # on-path: the column is informational; unknown
        if dist == 0:
            place = "on_path"; key = f"h37rv:{h}"
        elif dist <= a.near_tol:
            place = "off_path_near"; key = f"node:{nid}:{noff}"
        else:
            place = "off_path_accessory"; key = f"node:{nid}:{noff}"
        tally[place] += 1
        hit, call = "", ""
        if place == "on_path" and s["sample"] in ism_ran:
            for x, y, c in ism.get(s["sample"], []):
                if abs(h - x) <= a.ism_window or abs(h - y) <= a.ism_window:
                    hit, call = 1, c; break
            if hit == "":
                hit = 0
        s.update(placement=place, h37rv_pos=h, dist_to_ref=dist,
                 node=nid, node_offset=noff,
                 node_occ=occ.get((s["path"], nid), ""),
                 frame_strand=strand, key=key,
                 ismapper=hit, ism_call=call)
        rows.append(s)

    rows.sort(key=lambda r: (r["sample"], r["clean_pos"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"\n  {'placement':22s} {'n':>5s}")
    print("  " + "-" * 28)
    for k in ("on_path", "off_path_near", "off_path_accessory", "unprojected"):
        if tally[k]:
            print(f"  {k:22s} {tally[k]:5d}")
    print("  " + "-" * 28)
    print(f"  {'total':22s} {len(rows):5d}")
    if tally["no_node_length"]:
        print(f"  {tally['no_node_length']} of the unprojected are off-path "
              f"sites on a reverse-walked node missing from {a.node_lengths}")

    onp = [r for r in rows if r["placement"] == "on_path"
           and r["ismapper"] != ""]
    if not a.ismapper_dir:
        print("\n  ISMapper join off (no --ismapper-dir): the ismapper column "
              "is blank")
    elif len(ism_ran) < len({r["sample"] for r in rows}):
        print(f"\n  ISMapper tables for {len(ism_ran)} of "
              f"{len({r['sample'] for r in rows})} isolates in "
              f"{a.ismapper_dir}; the others are blank, not 0")
    if onp:
        m = sum(1 for r in onp if r["ismapper"] == 1)
        print(f"\n  ISMapper join, on-path sites only: {m}/{len(onp)} = "
              f"{m/len(onp):.1%} within +/-{a.ism_window} bp of an ISMapper region")
        by = collections.Counter()
        for r in onp:
            by[(r["tier"], r["ismapper"] == 1)] += 1
        for t in sorted({r["tier"] for r in onp}):
            y, n = by[(t, True)], by[(t, False)]
            print(f"    tier {t}: {y}/{y+n} = {y/(y+n):.0%}")
        print(f"  off-path sites are NOT joined: ISMapper ran against H37Rv and "
              f"cannot propose\n  a site in sequence H37Rv does not carry.")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
