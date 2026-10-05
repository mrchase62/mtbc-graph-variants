#!/usr/bin/env python3
"""T16: is composition's core excess an artefact of the H37Rv projection?

T14 found composition carries +16,561 core false positives against direct
calling. T15 excluded representation, T5 excluded mismapping and stage 4 excluded
silence. What was never excluded is the PROJECTION itself, because composition as
built goes through two approximate steps:

  the sample's calls        lifted R -> H37Rv by `paftools liftover` over a
                            minimap2 alignment. stage 4 measured 38.8% of records
                            landing where H37Rv's base is not their REF allele.
  R's own differences       recomputed per reference by `paftools call` on that
                            same alignment.

The graph replaces both with exact operations. R is a path in the pangenome graph,
built from R's own sequence, so an R position maps to an H37Rv position by
`odgi position` -- a path lookup, not an alignment. And R's differences from H37Rv
are already in the graph VCF as the records where R's genotype is ALT: no
recomputation at all.

So this scores the same composed callset assembled two ways, against the same
truth, and attributes the difference to the projection. If the core excess largely
disappears, the projection was the cause and the graph should be the frame. If it
survives, the errors are in the read-level calls against R and the frame choice is
neutral on precision.

`dist.to.ref` from odgi is kept: a non-zero value means the R position has no
H37Rv equivalent, which is the exact form of the "unliftable" signal that
`paftools` could only approximate.
"""
import argparse, collections, csv, gzip, importlib.util, os, subprocess, sys

# RE-RUN 2026-09-19. GRAPH_FRAME_RESOLUTION.md section 6 suspected this test of
# the coordinate frame defect that P4 had: feeding `odgi position` refs-frame
# positions when it reads the panel frame the graph was built in. 81 of the 218
# genomes below are matched to a reference the panel stores rotated or reverse
# complemented, which would have been enough to produce a double-digit error
# rate on its own.
#
# IT IS NOT AFFECTED, and the suspicion was wrong. Measured with
# graphframe/bin/frame_detect.py: the stage 2 simulation aligned to the ROTATED
# panel sequence, so `refbias/work/stage2/distance/*.vcf.gz` is already in the
# panel frame and needs no conversion, where P1 and P2 align to
# refbias/build/<id>/refs/ and theirs is in the refs frame. The two name the
# same contig at the same length, so nothing distinguishes them but the
# sequence. The frame is therefore DETECTED per sample below rather than
# assumed in either direction, and a sample whose frame cannot be decided is
# skipped rather than guessed at.
_gfs = importlib.util.spec_from_file_location(
    "graph_frame", os.path.join("graphframe", "bin", "graph_frame.py"))
graph_frame = importlib.util.module_from_spec(_gfs)
_gfs.loader.exec_module(graph_frame)
FRAMES = graph_frame.Frames()
_fds = importlib.util.spec_from_file_location(
    "frame_detect", os.path.join("graphframe", "bin", "frame_detect.py"))
frame_detect = importlib.util.module_from_spec(_fds)
_fds.loader.exec_module(frame_detect)


def open_maybe_gz(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def strip(x):
    return x.replace("-", "").replace("*", "")


def classify(ref, alt):
    r, a = strip(ref), strip(alt)
    return "SNP" if (len(r) == 1 and len(a) == 1) else "INDEL"


def load_truth(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if not f or f[0] != "V":
            continue
        ref, alt = f[6].upper(), f[7].upper()
        if set(ref + alt) - set("ACGTN*-"):
            continue
        out.append((int(f[2]), classify(ref, alt)))
    return out


def load_calls_rcoord(path):
    out = []
    for line in open_maybe_gz(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[6] not in (".", "PASS"):
            continue
        ref = f[3].upper()
        for alt in f[4].upper().split(","):
            if set(ref + alt) - set("ACGTN"):
                continue
            out.append((int(f[1]), classify(ref, alt)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default="refbias/stage2/samples.txt")
    ap.add_argument("--refmap", default="refbias/stage2/refmap.distance.tsv")
    ap.add_argument("--matched-dir", default="refbias/work/stage2/distance")
    ap.add_argument("--h37rv-dir", default="refbias/work/stage2/h37rv")
    ap.add_argument("--og", required=True)
    ap.add_argument("--graph-vcf",
                    default="graphs/CX333.s10k.k23.K15/all_variants.nolab.vcf.gz")
    ap.add_argument("--h37rv-path", default="GCF_000195955#1#NC_000962.3")
    ap.add_argument("--paths",
                    default="refbias/build/7713a8d71d8e/assets/paths.txt")
    ap.add_argument("--odgi", default=os.environ.get("MTB_ODGI", "odgi"))
    ap.add_argument("--pe-ppe", default="data/annotation/H37Rv_repeat_mask.bed")
    ap.add_argument("--refs", default="refbias/build/7713a8d71d8e/refs",
                    help="reference FASTA directory, used only to detect which "
                         "frame each matched VCF is in")
    ap.add_argument("--window", type=int, default=25)
    ap.add_argument("--workdir", default="refbias/work/t16")
    ap.add_argument("--out", default="refbias/T16.graph_frame.tsv")
    ap.add_argument("--per-genome",
                    default="refbias/T16.graph_frame.per_genome.tsv",
                    help="one row per genome and region. The first run wrote "
                         "only the aggregate, so whether its error rate was "
                         "concentrated in the genomes whose reference is stored "
                         "in a different frame could not be answered from its "
                         "output. It can be answered from this file.")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    pe = []
    for line in open(a.pe_ppe):
        f = line.split()
        if len(f) >= 4 and ("PE" in f[3].upper() or "PPE" in f[3].upper()):
            pe.append((int(f[1]), int(f[2])))
    pe.sort()

    def in_pe(p):
        lo, hi = 0, len(pe) - 1
        while lo <= hi:
            m = (lo + hi) // 2
            s, e = pe[m]
            if p < s:
                hi = m - 1
            elif p >= e:
                lo = m + 1
            else:
                return True
        return False

    # path name per reference id
    pathname = {}
    if os.path.exists(a.paths):
        for line in open(a.paths):
            nm = line.strip()
            if nm:
                pathname[nm.split("#")[0]] = nm

    # R's own differences from H37Rv, straight out of the graph VCF -- these are
    # exact and need no recomputation per reference
    # vg deconstruct emits nested bubbles: LV=0 top-level snarls and LV=1
    # children that decompose them. In all_variants.nolab.vcf.gz that is 25,354
    # parents and 65,950 children. Summing both counts the same underlying
    # difference twice -- it put 1,804 records per matched reference into the
    # first run of this test where the true figure is nearer 1,272, and the
    # surplus scores as false positives. So keep the finer LV=1 decomposition,
    # and keep an LV=0 record only where no LV=1 record overlaps it, which is a
    # simple top-level variant with nothing nested inside.
    child = []
    for line in open_maybe_gz(a.graph_vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t", 8)
        info = dict((kv.split("=", 1) + [""])[:2] for kv in f[7].split(";") if kv)
        if info.get("LV") != "0":
            pos = int(f[1])
            child.append((pos, pos + max(len(f[3]), 1)))
    child.sort()
    cstart = [c[0] for c in child]

    def has_child(pos, end):
        import bisect
        i = bisect.bisect_left(cstart, end)
        j = i - 1
        while j >= 0 and child[j][1] > pos:
            if child[j][0] < end:
                return True
            j -= 1
        return False

    hdr = None
    for line in open_maybe_gz(a.graph_vcf):
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t"); break
    gsamples = hdr[9:]
    gidx = {s: i for i, s in enumerate(gsamples)}
    rvh = collections.defaultdict(list)
    rvh_allele = collections.defaultdict(dict)
    n_nested_dropped = 0
    for line in open_maybe_gz(a.graph_vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        pos, ref = int(f[1]), f[3].upper()
        info = dict((kv.split("=", 1) + [""])[:2] for kv in f[7].split(";") if kv)
        if info.get("LV") == "0" and has_child(pos, pos + max(len(ref), 1)):
            n_nested_dropped += 1
            continue
        alts = f[4].upper().split(",")
        gts = f[9:]
        for ai, alt in enumerate(alts, 1):
            if set(ref + alt) - set("ACGTN"):
                continue
            tag = str(ai)
            kind = classify(ref, alt)
            for i, g in enumerate(gts):
                if g == tag:
                    rvh[gsamples[i]].append((pos, kind))
                    rvh_allele[gsamples[i]][pos] = (ref, alt)
    print(f"  graph VCF: dropped {n_nested_dropped} LV=0 parents that have a "
          f"nested LV=1 child, to avoid counting a difference twice")
    print(f"  graph VCF: {len(gsamples)} genomes, "
          f"median {sorted(len(v) for v in rvh.values())[len(rvh)//2]} "
          f"differences from H37Rv per genome")

    refmap = {}
    for line in open(a.refmap):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            refmap[f[0]] = f[1]

    def score(truth, called):
        tb = collections.defaultdict(list)
        for p, k in truth:
            tb[k].append(p)
        cb = collections.defaultdict(list)
        for p, k in called:
            cb[k].append(p)
        agg = collections.defaultdict(collections.Counter)
        for k, tp in tb.items():
            cp = sorted(cb.get(k, []))
            used = set()
            for p in sorted(tp):
                hit = None
                for i, q in enumerate(cp):
                    if i in used or abs(q - p) > a.window:
                        continue
                    hit = i
                    break
                if hit is not None:
                    used.add(hit)
                    agg["pe" if in_pe(p) else "core"]["TP"] += 1
                else:
                    agg["pe" if in_pe(p) else "core"]["FN"] += 1
            for i, q in enumerate(cp):
                if i not in used:
                    agg["pe" if in_pe(q) else "core"]["FP"] += 1
        for k, cp in cb.items():
            if k not in tb:
                for q in cp:
                    agg["pe" if in_pe(q) else "core"]["FP"] += 1
        return agg

    tot = collections.defaultdict(collections.Counter)
    per_genome = []
    frame_of, frame_tally = {}, collections.Counter()
    n = nolift = 0
    for s in [x.strip() for x in open(a.samples) if x.strip()]:
        R = refmap.get(s)
        tp = os.path.join(a.h37rv_dir, f"{s}.truth.var")
        vm = os.path.join(a.matched_dir, f"{s}.vcf.gz")
        if not R or not os.path.exists(tp) or not os.path.exists(vm):
            continue
        rpath = pathname.get(R)
        if not rpath:
            continue
        calls = load_calls_rcoord(vm)
        if not calls:
            continue
        # position -> kind, so the odgi output can be paired back in O(1). A
        # position carrying two alleles keeps the first, which is what the
        # scorer would collapse them to anyway.
        kindof = {}
        for p, k in calls:
            kindof.setdefault(p, k)
        altof = {}
        for line in open_maybe_gz(vm):
            if line.startswith("#"):
                continue
            ff = line.rstrip("\n").split("\t")
            if ff[6] not in (".", "PASS"):
                continue
            aa = ff[4].upper().split(",")[0]
            if not (set(ff[3].upper() + aa) - set("ACGTN")):
                altof.setdefault(int(ff[1]), aa)
        pf = os.path.join(a.workdir, f"{s}.pp.txt")
        racc = rpath.split("#")[0]
        vframe = frame_of.get(s)
        if vframe is None:
            vframe, nd, frr, frp = frame_detect.detect(vm, racc, FRAMES,
                                                       a.refs)
            frame_of[s] = vframe
            frame_tally[vframe] += 1
            if vframe is None:
                print(f"    {s}: cannot decide the frame of {vm} "
                      f"(refs {frr:.2f}, panel {frp:.2f}) -- skipped",
                      file=sys.stderr)
        if vframe is None:
            continue
        # 'either' means the panel stores this accession identically, so the
        # conversion is the identity and it does not matter which branch runs
        to_panel = (lambda q: FRAMES.to_panel(racc, q)) if vframe == "refs" \
            else (lambda q: q)
        to_refs = (lambda q: FRAMES.to_refs(racc, q)) if vframe == "refs" \
            else (lambda q: q)
        with open(pf, "w") as fh:
            for p, k in calls:
                fh.write(f"{rpath},{to_panel(p - 1)},+\n")
        r = subprocess.run([a.odgi, "position", "-i", a.og, "-F", pf,
                            "-r", a.h37rv_path, "-t", "2"],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                           universal_newlines=True)
        os.remove(pf)
        lifted = []
        allele_at = {}
        for line in r.stdout.splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            try:
                src = f[0].rsplit(",", 2)
                hpos = int(f[1].rsplit(",", 2)[-2])
                dist = int(f[2])
            except (ValueError, IndexError):
                continue
            srcpos = to_refs(int(src[-2])) + 1
            kind = kindof.get(srcpos)
            if kind is None:
                continue
            if dist != 0:
                nolift += 1
                continue
            h1 = hpos + 1
            sa = altof.get(srcpos)
            rv = rvh_allele.get(R, {}).get(h1)
            if rv and sa is not None and sa == rv[0]:
                allele_at[h1] = "revert"
            lifted.append((h1, kind))
        # Compose ALLELES, not positions. Where R already differs from H37Rv and
        # the sample differs from R, the composed variant is (H37Rv allele ->
        # sample allele); if the sample has reverted to H37Rv's allele there is
        # no variant at all. Concatenating the two callsets emits one anyway.
        rvh_at = {p: k for p, k in rvh.get(R, [])}
        overridden = {p for p, _ in lifted}
        composed = [x for x in rvh.get(R, []) if x[0] not in overridden]
        for p, k in lifted:
            if p in rvh_at and allele_at.get(p) == "revert":
                continue
            composed.append((p, k))
        truth = load_truth(tp)
        for reg, c in score(truth, composed).items():
            tot[("graph", reg)].update(c)
            tpn = c["TP"] + c["FP"]
            per_genome.append(dict(
                sample=s, reference=R,
                frame=("same" if FRAMES.identity(racc) else
                       ("rc" if FRAMES.flipped(racc) else "rotated")),
                offset=FRAMES.f[racc][1], region=reg,
                TP=c["TP"], FN=c["FN"], FP=c["FP"],
                ppv=round(c["TP"] / tpn, 4) if tpn else ""))
        n += 1
        if n % 40 == 0:
            print(f"    {n} genomes placed", file=sys.stderr)

    print(f"\n  {n} genomes; {nolift} calls had no H37Rv equivalent "
          f"(dist.to.ref != 0)")
    print(f"  matched VCF frame, detected per genome: "
          + ", ".join(f"{k}={v}" for k, v in sorted(frame_tally.items(),
                                                    key=lambda x: str(x[0])))
          + "\n")
    rows = []
    print(f"  {'frame':<8s} {'region':<8s} {'TP':>8s} {'FP':>8s} {'PPV':>8s}")
    for reg in ("pe", "core"):
        c = tot[("graph", reg)]
        tpn = c["TP"] + c["FP"]
        pv = c["TP"] / tpn if tpn else 0
        print(f"  {'graph':<8s} {reg:<8s} {c['TP']:>8d} {c['FP']:>8d} {pv:>8.4f}")
        rows.append(dict(frame="graph", region=reg, genomes=n, TP=c["TP"],
                         FN=c["FN"], FP=c["FP"], ppv=round(pv, 4)))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    if per_genome:
        with open(a.per_genome, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(per_genome[0]), delimiter="\t")
            w.writeheader(); w.writerows(per_genome)
        byframe = collections.defaultdict(collections.Counter)
        seen = collections.defaultdict(set)
        for r in per_genome:
            g = "same frame" if r["frame"] == "same" else "shifted"
            byframe[(g, r["region"])].update(
                {k: r[k] for k in ("TP", "FN", "FP")})
            seen[g].add(r["sample"])
        print(f"\n  split by whether the panel stores the matched reference in "
              f"the same frame as refs/:")
        print(f"  {'frame':<12s} {'region':<8s} {'genomes':>8s} {'TP':>8s} "
              f"{'FP':>8s} {'PPV':>8s}")
        for g in ("same frame", "shifted"):
            for reg in ("pe", "core"):
                c = byframe[(g, reg)]
                d = c["TP"] + c["FP"]
                print(f"  {g:<12s} {reg:<8s} {len(seen[g]):>8d} {c['TP']:>8d} "
                      f"{c['FP']:>8d} {(c['TP']/d if d else 0):>8.4f}")
        print(f"  written: {a.per_genome}")
    print(f"\n  for comparison, T14 on the same genomes, paftools projection:")
    print(f"    composed core FP 31816 (PPV 0.9135), direct core FP 15255 (0.9546)")
    print(f"    composed PE   FP 11562 (PPV 0.8664), direct PE   FP 25398 (0.7217)")
    print(f"  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
