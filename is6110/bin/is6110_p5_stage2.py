#!/usr/bin/env python3
"""Stage 2: earn REF for the samples that did not report an IS6110 site.

Stage 1 (`is6110_p5_merge.py`) leaves every non-carrier NOCALL, which asserts
nothing but leaves 94.7% of the IS6110 block uninformative. This asks each
non-carrier's own element-free alignment whether it looked and found nothing.

    for each key, take a carrier's own (reference, position) as the source
    project it onto the non-carrier's reference path            -> ABSENT if
                                                                   dist != 0
    convert that reference position into the non-carrier's clean coordinate
    read depth there, and ask whether a junction candidate sits within --window

    depth >= --min-dp and no candidate   REF      it looked and found nothing
    a junction candidate is present      reported, NOT silently called ALT
    otherwise                            NOCALL

WHY THE PROJECTION GOES CARRIER -> TARGET AND NOT VIA H37Rv. 68 of the 170 keys
are off the H37Rv path, so a route through H37Rv cannot carry them -- those are
exactly the sites in sequence H37Rv does not have. Projecting straight from a
carrier's path to the target's path keeps them. It is batched by TARGET path,
one `odgi position -r <target>` per distinct reference with every key's source
in a single input file, so this is ~18 graph loads rather than one per pair.

FRAMES. The carrier position comes from a refs-frame alignment and `odgi` reads
panel coordinates; the returned target position is panel and the crossmap is
refs. Both conversions go through graphframe, and for a target the panel stores
reverse complemented the returned coordinate is converted with the strand, which
graph_frame.to_refs does. See graphframe/docs/GRAPH_FRAME_RESOLUTION.md.

DEPTH FLOOR. `--min-dp` defaults to 5, the same floor `p5_states.py` applies to
a SNP REF call, so an IS6110 REF means the same strength of evidence as every
other REF in the matrix.
"""
import argparse, bisect, collections, csv, importlib.util, os, subprocess, sys

def _load(name, path):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

graph_frame = _load("graph_frame", os.path.join("graphframe", "bin", "graph_frame.py"))


def read_crossmap(path):
    """(orig_end, cum_deleted) rows, so an original coordinate converts to the
    element-free one by subtracting the deletions that precede it."""
    out = []
    for r in csv.DictReader(open(path, newline=""), delimiter="\t"):
        out.append((int(r["orig_end"]), int(r["cum_deleted"])))
    out.sort()
    return out


def orig_to_clean(cm, pos):
    i = bisect.bisect_right([c[0] for c in cm], pos)
    return pos - (cm[i - 1][1] if i else 0)


# map the stage-1 key back to the same short form the carrier table uses.
# Stage-1 keys are  h37rv:<pos>:<ref>><alt>  and  node:<id>:<off>:<ref>><alt>
# (mtb_norm), and the carrier table's `node` column is already "<id>:<off>",
# so the node short form MUST keep the offset. It used to drop it
# ("node:<id>:"), which never matched "node:<id>:<off>:", so no node-frame
# key found a carrier and every node-frame NOCALL was written ABSENT --
# 446,005 cells in gwas1000.
def short(k):
    f = k.split(":")
    if f[0] == "h37rv":
        return f"h37rv:{f[1]}:"
    if f[0] == "node" and len(f) >= 3:
        return f"node:{f[1]}:{f[2]}:"
    raise ValueError(f"unrecognised stage-1 key {k!r}")


# ---- stage-1 input, whole-cohort table or one file per sample ----------------
def stage1_rows(a, sample=None):
    """Stage-1 state rows: one sample's from --stage1-dir, or everyone's."""
    if a.stage1_dir:
        names = [sample] if sample else [r["sample"] for r in csv.DictReader(
            open(a.refmap, newline=""), delimiter="\t")]
        out = []
        for s in names:
            out += list(csv.DictReader(open(os.path.join(a.stage1_dir, f"{s}.tsv"),
                                            newline=""), delimiter="\t"))
        return out
    rows = list(csv.DictReader(open(a.stage1_states, newline=""), delimiter="\t"))
    return [r for r in rows if r["sample"] == sample] if sample else rows


def carriers(a):
    """short key -> (carrier reference, carrier position), the FIRST carrier in
    the key table's order, exactly as the single-process version chose."""
    src = {}
    for x in csv.DictReader(open(a.cohort_keys, newline=""), delimiter="\t"):
        if x["frame"] == "h37rv":
            k = f'h37rv:{x["h37rv_pos"]}:'
        else:
            k = f'node:{x["node"]}:'
        src.setdefault((x["sample"], x["reference"], int(x["r_pos"])), []).append(k)
    carrier_of = {}
    for (cs, cref, cpos), ks in src.items():
        for k in ks:
            carrier_of.setdefault(k, (cref, cpos))
    return carrier_of


def proj_path(a, tref):
    return os.path.join(a.workdir, "proj", f"{tref}.tsv")


# ---- mode `project`: one target reference -----------------------------------
def project(a, tref, fr, paths, og, ref_of, carrier_of):
    """Project every short key a sample of `tref` still needs onto tref's path.

    One odgi load per reference, as before; what changed is that each reference
    is its own task. Written atomically to <workdir>/proj/<tref>.tsv as
    short_key, pos (refs frame, 1-based; empty when odgi returned no row) and
    dist. A reference with no path in the graph, or with nothing to project,
    still gets a header-only file, so a sample task can tell "projected,
    nothing here" from "not run".
    """
    need = set()
    for s, R in ref_of.items():
        if R != tref:
            continue
        for r in stage1_rows(a, s):
            if r["state"] == "NOCALL":
                need.add(short(r["key"]))
    os.makedirs(os.path.join(a.workdir, "proj"), exist_ok=True)
    rows = []
    tpath = paths.get(tref)
    if tpath is not None and need:
        qf = os.path.join(a.workdir, "proj", f".{tref}.{os.getpid()}.pos")
        order = []
        with open(qf, "w") as fh:
            for k in sorted(need):
                cref, cpos = carrier_of.get(k, (None, None))
                if cref is None or cref not in paths:
                    continue
                fh.write(f"{paths[cref]},{fr.to_panel(cref, cpos - 1)},+\n")
                order.append((k, cref, cpos))
        if order:
            p = subprocess.run([a.odgi, "position", "-i", og, "-F", qf, "-r",
                                tpath, "-t", str(a.threads)],
                               capture_output=True, text=True)
            # A failed odgi call returns no rows, and every key for this target
            # would then read as unprojected. Stop instead.
            if p.returncode != 0:
                sys.exit(f"FATAL: odgi position failed for {tref} "
                         f"(exit {p.returncode}):\n{p.stderr[-2000:]}")
            got = {}
            for line in p.stdout.split("\n"):
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) < 3:
                    continue
                sp = f[0].rsplit(",", 2)
                tp = f[1].rsplit(",", 2)
                try:
                    got[(sp[0], int(sp[1]))] = (int(tp[1]), int(f[2]))
                except (ValueError, IndexError):
                    continue
            for k, cref, cpos in order:
                v = got.get((paths[cref], fr.to_panel(cref, cpos - 1)))
                if v is None:
                    # odgi answered but returned no row for this query:
                    # nothing was measured, which is not evidence of absence
                    rows.append((k, "", -1))
                else:
                    tp, dist = v
                    rows.append((k, fr.to_refs(tref, tp) + 1, dist))
        os.remove(qf)
    out = proj_path(a, tref)
    with open(out + ".tmp", "w") as fh:
        fh.write("short_key\tpos\tdist\n")
        for k, pos, dist in rows:
            fh.write(f"{k}\t{pos}\t{dist}\n")
    os.replace(out + ".tmp", out)
    print(f"  {tref}: {len(rows)} keys projected onto its path")


def load_proj(a, tref):
    p = proj_path(a, tref)
    if not os.path.exists(p):
        sys.exit(f"FATAL: no projection for {tref} at {p}; run --mode project")
    proj = {}
    for r in csv.DictReader(open(p, newline=""), delimiter="\t"):
        proj[r["short_key"]] = ((int(r["pos"]) if r["pos"] else None), int(r["dist"]))
    return proj


# ---- mode `sample`: one sample ---------------------------------------------
def decide(a, s, R):
    """Stage-2 rows for one sample, from its reference's projection, its own
    element-free alignment and its junction candidates."""
    rows = stage1_rows(a, s)
    proj = load_proj(a, R)

    cmp_ = os.path.join(a.crossmap_dir, f"{R}.crossmap.tsv")
    cm = read_crossmap(cmp_) if os.path.exists(cmp_) else None
    bam = os.path.join(a.isclean_dir, f"{s}.isclean.bam")
    contig = None
    if os.path.exists(bam):
        r = subprocess.run([a.samtools, "idxstats", bam],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"FATAL: samtools idxstats failed on {bam} "
                     f"(exit {r.returncode}):\n{r.stderr[-2000:]}")
        for line in r.stdout.split("\n"):
            f = line.split("\t")
            if len(f) > 1 and f[0].endswith("_isclean"):
                contig = f[0]; break

    # pass one: resolve each row to either an immediate verdict or a position
    pending, verdict, want = {}, {}, set()
    for i, r in enumerate(rows):
        if r["state"] != "NOCALL":
            continue
        pos, dist = proj.get(short(r["key"]), (None, -1))
        # ABSENT only from a projection that was actually made and landed off
        # the target's path. No projection at all -- target or carrier path not
        # in the graph, no carrier for the key, or odgi returned nothing -- is
        # "not measured", which is NOCALL.
        if pos is None:
            verdict[i] = ("NOCALL", "no projection")
        elif dist != 0:
            verdict[i] = ("ABSENT", f"off the path, dist {dist}")
        elif cm is None or contig is None:
            verdict[i] = ("NOCALL", "no crossmap or no element-free alignment")
        else:
            cpos = orig_to_clean(cm, pos)
            pending[i] = cpos
            want.add(cpos)

    # pass two: ONE samtools call for every position this sample needs. It
    # used to be one `samtools depth -r` per position -- ~56,000 calls and 33
    # minutes over scale100, nearly all process start-up. -a reports every
    # position in the regions, zero-depth ones included, so a position with no
    # coverage comes back rather than silently defaulting.
    depth = {}
    if want:
        bed = os.path.join(a.workdir, f"{s}.depth.bed")
        with open(bed, "w") as fh:
            for p in sorted(want):
                fh.write(f"{contig}\t{p - 1}\t{p}\n")
        r = subprocess.run([a.samtools, "depth", "-a", "-b", bed, bam],
                           capture_output=True, text=True)
        # a failed call would otherwise read as depth 0 everywhere -> NOCALL
        # for the whole sample, with nothing said
        if r.returncode != 0:
            sys.exit(f"FATAL: samtools depth failed on {bam} "
                     f"(exit {r.returncode}):\n{r.stderr[-2000:]}")
        for line in r.stdout.split("\n"):
            f = line.split("\t")
            if len(f) == 3:
                try:
                    depth[int(f[1])] = int(f[2])
                except ValueError:
                    pass
        os.remove(bed)

    cand = []
    jp = os.path.join(a.isclean_dir, f"{s}.junctions.tsv")
    if os.path.exists(jp):
        for r in csv.DictReader(open(jp, newline=""), delimiter="\t"):
            try:
                cand.append(int(r["pos"]))
            except (ValueError, KeyError):
                pass
    cand.sort()

    # pass three: decide
    out = []
    for i, r in enumerate(rows):
        if r["state"] != "NOCALL":
            out.append(dict(r, evidence="reported by this sample"))
            continue
        if i in verdict:
            st, ev = verdict[i]
        else:
            cpos = pending[i]
            j = bisect.bisect_left(cand, cpos - a.window)
            near = j < len(cand) and cand[j] <= cpos + a.window
            dp = depth.get(cpos, 0)
            if near:
                st, ev = "NOCALL", f"a junction candidate sits within {a.window} bp"
            elif dp >= a.min_dp:
                st, ev = "REF", f"depth {dp}, no junction candidate"
            else:
                st, ev = "NOCALL", f"depth {dp} below {a.min_dp}"
        out.append(dict(r, state=st,
                        allele=r["allele"] if st == r["state"] else "",
                        evidence=ev))
    return out, (list(rows[0]) + ["evidence"] if rows else None)


def sample_path(a, s):
    return os.path.join(a.workdir, "samples", f"{s}.tsv")


def write_sample(a, s, out, fields):
    os.makedirs(os.path.join(a.workdir, "samples"), exist_ok=True)
    p = sample_path(a, s)
    with open(p + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or ["sample", "key", "state",
                                                     "allele", "evidence"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(out)
    os.replace(p + ".tmp", p)


# ---- mode `merge`: the cohort table ----------------------------------------
def merge(a, samples):
    """Concatenate the per-sample tables in refmap order into --out, the same
    table the single-process version wrote, and report the counts."""
    counts = collections.Counter()
    missed = 0
    hdr = None
    with open(a.out + ".tmp", "w", newline="") as fh:
        for s in samples:
            p = sample_path(a, s)
            if not os.path.exists(p):
                sys.exit(f"FATAL: no stage-2 table for {s} at {p}; did every "
                         f"sample task finish?")
            with open(p, newline="") as sf:
                first = sf.readline()
                if hdr is None:
                    hdr = first
                    fh.write(first)
                elif first != hdr:
                    sys.exit(f"FATAL: {p} has a different header")
                cols = first.rstrip("\n").split("\t")
                si, ei = cols.index("state"), cols.index("evidence")
                for line in sf:
                    fh.write(line)
                    f = line.rstrip("\n").split("\t")
                    counts[f[si]] += 1
                    missed += f[ei].startswith("a junction candidate")
        if hdr is None:
            fh.write("sample\tkey\tstate\tallele\tevidence\n")
    os.replace(a.out + ".tmp", a.out)
    n = sum(counts.values()) or 1
    print(f"\n  {n} cells after stage 2")
    for k in ("ALT", "REF", "ABSENT", "NOCALL"):
        if counts[k]:
            print(f"    {k:8} {counts[k]:6} ({100*counts[k]/n:5.1f}%)")
    print(f"  {missed} cells hold a junction candidate the detector did not "
          f"call a site;\n  they are left NOCALL and listed, not promoted to ALT")
    print(f"  written: {a.out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("all", "project", "sample", "merge"),
                    default="all",
                    help="all: every step in this process (the old behaviour). "
                         "project: one reference (--ref or --index). sample: one "
                         "sample (--sample or --index). merge: the cohort table. "
                         "The three split modes are what lets a 10,000-isolate "
                         "cohort run as arrays instead of one serial job.")
    ap.add_argument("--ref", default="", help="--mode project: the target reference")
    ap.add_argument("--sample", default="", help="--mode sample: the sample")
    ap.add_argument("--index", type=int, default=0,
                    help="1-based: with project, the index-th distinct reference "
                         "of the refmap, in sorted order (a task past the last "
                         "reference exits 0); with sample, the index-th refmap "
                         "row. For SLURM_ARRAY_TASK_ID")
    ap.add_argument("--cohort-keys", default="is6110/results/p1i_cohort_keys.tsv")
    ap.add_argument("--stage1-states", default="is6110/results/p5_is6110_states.tsv")
    ap.add_argument("--stage1-dir", default="",
                    help="per-sample stage-1 tables, <sample>.tsv, from "
                         "is6110_p5_merge.py --out-states-dir; read instead of "
                         "--stage1-states so a task never parses the cohort's")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--paths", default="refbias/build/7713a8d71d8e/assets/paths.txt")
    ap.add_argument("--graph", default=None)
    ap.add_argument("--odgi", default=os.environ.get("MTB_ODGI", "odgi"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--isclean-dir", default="refbias/p1i")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched")
    ap.add_argument("--min-dp", type=int, default=5)
    ap.add_argument("--window", type=int, default=10)
    ap.add_argument("--workdir", default="refbias/p1i/p5stage2")
    ap.add_argument("--out", default="is6110/results/p5_is6110_states.stage2.tsv")
    a = ap.parse_args()

    os.makedirs(a.workdir, exist_ok=True)
    refmap = list(csv.DictReader(open(a.refmap, newline=""), delimiter="\t"))
    ref_of = {r["sample"]: r["reference"] for r in refmap}
    samples = [r["sample"] for r in refmap]
    refs = sorted(set(ref_of.values()))

    tref = a.ref
    if a.mode == "project" and not tref:
        # resolved before anything is loaded: most tasks of the array are past
        # the cohort's last reference and should cost nothing
        if a.index < 1:
            sys.exit("FATAL: --mode project needs --ref or --index")
        if a.index > len(refs):
            print(f"  index {a.index}: only {len(refs)} references; nothing to do")
            return 0
        tref = refs[a.index - 1]

    if a.mode in ("all", "project"):
        og = a.graph
        if og is None:
            import glob
            g = sorted(glob.glob("graphs/CX333.s10k.k23.K15/*.smooth.final.og"))
            og = g[0] if g else sys.exit("no graph; pass --graph")
        fr = graph_frame.Frames()
        paths = {l.split("#")[0]: l.strip() for l in open(a.paths) if l.strip()}
        carrier_of = carriers(a)
        if a.mode == "project":
            project(a, tref, fr, paths, og, ref_of, carrier_of)
            return 0
        for tref in refs:
            project(a, tref, fr, paths, og, ref_of, carrier_of)

    if a.mode == "sample":
        s = a.sample or (samples[a.index - 1] if 1 <= a.index <= len(samples)
                         else sys.exit("FATAL: --mode sample needs --sample or a "
                                       "valid --index"))
        out, fields = decide(a, s, ref_of[s])
        write_sample(a, s, out, fields)
        n = collections.Counter(r["state"] for r in out)
        print(f"  {s}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())))
        return 0

    if a.mode == "all":
        for s in samples:
            out, fields = decide(a, s, ref_of[s])
            write_sample(a, s, out, fields)

    merge(a, samples)
    return 0


if __name__ == "__main__":
    sys.exit(main())
