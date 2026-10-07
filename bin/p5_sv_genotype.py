#!/usr/bin/env python3
"""Earn REF for a structural variant the sample did not have called.

WHY THIS IS NOT "NOT CALLED MEANS ABSENT". bin/p5_sv_matrix.py refuses to write
REF for an SV, and its reason is measured, not cautious: stage 6 put delly's
sensitivity at 0.29 for 50-499 bp deletions and 0.01-0.05 for insertions, so a
missing call is mostly a missed call. Reading silence as absence would assert
hundreds of thousands of absences that were never observed.

The consequence, though, is that the SV block of the merged VCF is 1.80% ALT
and 98.20% NOCALL with no REF anywhere. Measured on scale200, that makes every
one of its 6,415 SV records untestable by any method that has to contrast
presence with absence: with no leaf ever in the reference state, a parsimony
reconstruction has nothing to work from and all 6,415 score zero independent
gains.

So absence is MEASURED here instead of inferred from a caller's silence. For a
deletion it is directly observable: if the sample has read depth across the
interval the deletion would remove, the deletion is not there. That is an
observation of the same kind and the same strength as the one behind every
other REF in the matrix -- `p5_states.py` calls a SNP position REF when the
gVCF says it was assessed with depth, and the floor here is the same --min-dp.

WHAT THIS DELIBERATELY DOES NOT DO.

  Insertions stay NOCALL. Depth across an insertion point says nothing: the
  inserted sequence is not in the reference, so a sample carrying it and a
  sample lacking it both show ordinary depth there. Testing insertion absence
  needs spanning reads and their insert sizes from the BAM, which is a
  different instrument and is not this script. `svtype` is recorded on every
  row so the untested classes stay visible rather than blending into NOCALL.

  A called ALT is never overwritten.

ALT IS EARNED THE SAME WAY REF IS, and the reason is the matched reference.
Many deletions are already carried by the panel genome a sample was aligned to
-- absorbed into the reference -- so the aligner sees contiguous sequence and
no caller ever emits a call. Those samples genuinely have the deletion. They
show up here as an interval that does not project into their own reference at
all, which an earlier version of this script recorded as ABSENT, conflating
"the sequence is not there" with "nothing is known". The flank probes separate
them: reads either side of the interval mean the sample really is sitting on a
reference that lacks it, and the cell is ALT with `inherited` as its evidence
-- the same thing the SV matrix's `component` column already calls inherited
for the rows a caller did reach.

The symmetric case is an interval the sample's reference does have and the
sample's reads do not cover, with covered flanks: that is ALT by depth. It is
the same measurement that earns REF when the interval IS covered, so accepting
one and refusing the other would be asymmetric. Every cell records in its
`evidence` column how it was earned, and `--no-promote` turns both off.

PROBES. --probes interior positions evenly spaced across the interval, each
projected into the sample's own reference the same way p5_merge.sh projects key
positions, because the gVCF is in that frame and the SV matrix is in H37Rv's.
A probe whose projection reports dist != 0 has no equivalent in the sample's
reference at all; when most probes are like that the cell is ABSENT, which is a
different statement from either REF or NOCALL and the matrix already carries it.

Two modes:
  --probes-out   write the H37Rv probe positions, for the projection step
  --out          read the projection back and write this sample's states
"""
import argparse, bisect, collections, csv, gzip, os, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p5_states import load_gvcf, load_gvcf_blocks, make_cov, make_in_spans


def gvcf_depth(path, min_dp):
    """covered(t): a gVCF line of DP >= min_dp spans R position t, and t is
    not inside a deletion the same gVCF calls.

    NOT INSIDE A DELETION GATK ITSELF CALLED (audit P4P5-14). GATK starts the
    next reference block after a deletion's POS, so bases the sample's own
    gVCF calls deleted can sit in GT=0 blocks at real DP -- SAMEA1016081's
    84 bp deletion at 843,010 had blocks inside it at DP 12-25 -- and were
    counted as covered. Every other line still counts as depth."""
    cov_any = make_cov(load_gvcf_blocks(path), min_dp)
    _, _, alt_at = load_gvcf(path)
    in_called_del = make_in_spans(sorted(
        (p0 + len(c), p0 + len(rf) - 1) for p0, (rf, c, _) in alt_at.items()
        if c != "<NON_REF>" and len(c) < len(rf)))

    def covered(t):
        return cov_any(t) and not in_called_del(t)
    return covered


def probe_positions(pos, svlen, n):
    """Interior H37Rv positions across a deletion, evenly spaced."""
    length = abs(int(svlen or 0))
    if length < 2:
        return [pos]
    n = max(1, min(n, length - 1))
    return sorted({pos + max(1, round(length * i / (n + 1)))
                   for i in range(1, n + 1)})


def probe_span(start, end, n):
    """Interior positions across an explicit [start, end] interval."""
    length = max(end - start + 1, 1)
    if length < 2:
        return [start]
    n = max(1, min(n, length))
    return sorted({start + max(0, round((length - 1) * i / (n + 1)))
                   for i in range(1, n + 1)})


def flank_span(start, end, flank):
    return [max(1, start - flank), end + flank]


def read_intervals(path):
    rows = list(csv.DictReader(open(path, newline=""), delimiter="\t"))
    for r in rows:
        r["start"], r["end"] = int(r["start"]), int(r["end"])
        r["svlen"] = int(r["svlen"])
        r["carriers"] = set((r.get("ref_carriers") or "").split(",")) - {""}
    return rows


def flank_positions(pos, svlen, flank):
    """Positions just outside the interval, to ask whether the sample has any
    data in the neighbourhood at all.

    Without these, an interval that does not project is indistinguishable
    between the two cases that matter: the sample's reference genuinely lacks
    the sequence -- which IS the deletion -- and the projection simply failing.
    """
    length = abs(int(svlen or 0))
    return [max(1, pos - flank), pos + length + flank]


SV_META = ["key", "svtype", "h37rv_pos", "svlen", "size_band", "n_alt",
           "n_both_callers", "src", "max_qual", "qual_band", "stage11_ppv",
           "qual_from", "component"]


def read_sv_matrix(path, keep_sample=None):
    """Rows of the SV matrix with the per-row columns and, at most, ONE
    sample's column.

    Every per-sample task used to load the whole matrix -- every row with
    every sample's cell, as dicts. At 997 isolates that is 169 MB of text and
    2.2 GB of memory per task; at 10,000 it would be about 9 GB of text, per
    task. A task only ever reads its own sample's cell, so only that column is
    kept, and the probe step (keep_sample=None) keeps none. A sample the matrix
    does not name gets no cell, which callers read as NOCALL, as before.
    """
    with open(path, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr, [])
        samples = [c for c in hdr if c not in SV_META]
        keep = [(i, c) for i, c in enumerate(hdr) if c in SV_META]
        if keep_sample is not None and keep_sample in hdr:
            keep.append((hdr.index(keep_sample), keep_sample))
        rows = []
        for f in rdr:
            rows.append({c: (f[i] if i < len(f) else "") for i, c in keep})
    return rows, samples


def load_projection(path):
    """source H37Rv position -> (target position, distance)."""
    by_src = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            src = int(f[0].rsplit(",", 2)[-2]) + 1
            tgt = int(f[1].rsplit(",", 2)[-2]) + 1
            dist = int(f[2])
        except (ValueError, IndexError):
            continue
        by_src.setdefault(src, (tgt, dist))
    return by_src


def mapq_depth(a, positions):
    """MAPQ-filtered depth at a set of target-frame positions, in one call.

    One `samtools depth` over a regions file rather than a call per position:
    at ~6,600 probes a call apiece would dominate the task.
    """
    if not positions:
        return {}
    if not a.bam or not os.path.exists(a.bam):
        sys.exit(f"FATAL: --mapq-scope {a.mapq_scope} needs --bam, and "
                 f"{a.bam or '(none)'} is not readable")
    hdr = subprocess.run([a.samtools, "view", "-H", a.bam],
                         capture_output=True, text=True, check=True).stdout
    contig = None
    for line in hdr.splitlines():
        if line.startswith("@SQ"):
            for kv in line.split("\t"):
                if kv.startswith("SN:"):
                    contig = kv[3:]
            break
    if contig is None:
        sys.exit(f"FATAL: no @SQ line in {a.bam}")
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as bf:
        for p in sorted(positions):
            bf.write(f"{contig}\t{p - 1}\t{p}\n")
        bed = bf.name
    try:
        r = subprocess.run([a.samtools, "depth", "-a", "-Q", str(a.min_mapq),
                            "-b", bed, a.bam],
                           capture_output=True, text=True, check=True)
    finally:
        os.unlink(bed)
    out = {}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            out[int(f[1])] = int(f[2])
    return out


def genotype_intervals(a, ivs):
    """Every interval, genotyped the same way, for one sample."""
    for need in ("sample", "gvcf", "projected", "out"):
        if not getattr(a, need):
            sys.exit(f"FATAL: --{need} is required when writing states")
    by_src = load_projection(a.projected)
    cov = gvcf_depth(a.gvcf, a.min_dp)

    # Which intervals get the filtered depth, and the target positions they
    # need. Collected first so samtools is called once for the whole sample.
    def wants_mapq(r):
        if a.mapq_scope == "never":
            return False
        if a.mapq_scope == "all":
            return True
        return r.get("is6110_prox") in ("1", 1, True)
    need_pos = set()
    for r in ivs:
        if not wants_mapq(r):
            continue
        for p in (probe_span(r["start"], r["end"], a.probes)
                  + flank_span(r["start"], r["end"], a.flank)):
            t = by_src.get(p)
            if t is not None and t[1] == 0:
                need_pos.add(t[0])
    mq = mapq_depth(a, need_pos) if need_pos else {}
    n_mq = sum(1 for r in ivs if wants_mapq(r))
    if n_mq:
        print(f"  mapping-quality-filtered depth (Q>={a.min_mapq}) on "
              f"{n_mq:,} of {len(ivs):,} intervals, {len(need_pos):,} positions")
    # ---- second-frame evidence, for the inherited class --------------------
    tf = {}
    if a.twoframe and os.path.exists(a.twoframe):
        for q in csv.DictReader(open(a.twoframe, newline=""), delimiter="\t"):
            if q.get("sample") and q["sample"] != a.sample:
                continue
            try:
                tf[q["interval"]] = (q["h_state"], int(q["h_clip_sides"]))
            except (KeyError, ValueError):
                pass
        print(f"  second-frame evidence for {len(tf):,} intervals from "
              f"{a.twoframe}")
    elif not a.trust_inherited:
        print("  NOTE no --twoframe table: inherited candidates cannot be "
              "confirmed and will be NOCALL rather than ALT")

    def inherited_call(kind, flank_ok):
        """Promote an inherited candidate only on positive evidence."""
        if a.trust_inherited:
            if flank_ok and not a.no_promote:
                return "ALT", kind
            return "NOCALL", f"{kind}, flanks not covered"
        if not flank_ok or a.no_promote:
            return "NOCALL", f"{kind}, flanks not covered"
        hit = tf.get(r["interval"])
        if hit is None:
            return "NOCALL", f"{kind}_unconfirmed, no second frame"
        hs, hclips = hit
        # THE H37Rv-FRAME STATE DECIDES; A CLIP CLUSTER ALONE DOES NOT (audit
        # P4P5-3). `clips >= 1` used to promote before the state was looked
        # at, so an interval whose H37Rv-frame depth covered it in full was
        # "confirmed" by one end's clip cluster, which is common in repeats
        # and at IS6110 junctions: 995 ALT cells in scale200, 7,321 in
        # gwas1000. sv_twoframe.py's own ALT already takes clips into account
        # (both ends clipped at ambiguous depth), so its state is the
        # confirmation. Full depth with clips at an end is conflicting
        # evidence and is NOCALL; full depth without clips contradicts.
        if hs == "ALT":
            return "ALT", f"{kind}_confirmed"
        if hs == "REF":
            if hclips >= 1:
                return "NOCALL", f"{kind}_depth_present_but_clipped"
            return "REF", f"{kind}_contradicted"
        return "NOCALL", f"{kind}_unconfirmed"

    counts, evid, out = collections.Counter(), collections.Counter(), []
    for r in ivs:
        filtered = wants_mapq(r)
        # In an element-proximal interval the gVCF's DP is not evidence of
        # anything: it counts mismapped copies of the repeat. Use the filtered
        # depth for BOTH the interval and its flanks, so the confirmation is
        # made of the same material as the measurement.
        covf = ((lambda t: mq.get(t, 0) >= a.min_dp) if filtered else cov)
        probes = probe_span(r["start"], r["end"], a.probes)
        fl = [t for t in (by_src.get(p)
                          for p in flank_span(r["start"], r["end"], a.flank))
              if t is not None and t[1] == 0]
        flank_ok = any(covf(t) for t, _ in fl)
        ok = [t for t in (by_src.get(p) for p in probes)
              if t is not None and t[1] == 0]
        frac = (sum(1 for t, _ in ok if covf(t)) / len(ok)) if ok else 0.0

        if a.reference and a.reference in r["carriers"]:
            # THE REFERENCE CARRIES THE DELETION -- WHICH IS NOT THE SAME AS
            # THE SAMPLE CARRYING IT. This used to promote straight to ALT on
            # the graph's lookup plus a covered flank, reasoning that the
            # lookup is a fact about the panel rather than an inference. The
            # lookup is indeed a fact; the step from it to the SAMPLE is not.
            # A reference is chosen on whole-genome SNP distance, not per
            # locus, so a sample can sit on a carrier reference and not carry
            # the deletion itself, and the flank check cannot tell -- the
            # flanks are covered either way.
            #
            # Measured in sv2frame/TWOFRAME_RESULTS.md: of 13,700 inherited
            # calls over scale200, the H37Rv frame CONTRADICTS 26.5% outright,
            # showing the sequence present at normal depth, against 4.4% for
            # calls the matched frame actually measured. So the class needs
            # positive evidence, and `inherited_confirm` supplies it from the
            # second frame. Without a second frame the honest state is NOCALL,
            # not ALT. `--trust-inherited` restores the old behaviour for
            # reproducing output written before this change.
            state, why = inherited_call("inherited_reference", flank_ok)
        elif len(ok) < a.min_projected * len(probes):
            if not flank_ok or a.no_promote:
                state, why = "ABSENT", (f"{len(ok)}/{len(probes)} probes "
                                        f"project, flanks not covered")
            else:
                state, why = inherited_call("inherited_unrecorded", flank_ok)
        elif frac >= a.present_frac:
            state, why = "REF", "depth_present"
        elif frac <= a.deleted_frac and flank_ok and not a.no_promote:
            state, why = "ALT", "depth_absent"
        elif frac <= a.deleted_frac:
            state, why = "NOCALL", "uncovered, and so are the flanks"
        else:
            state, why = "NOCALL", "ambiguous depth"
        counts[state] += 1
        evid[why] += 1
        out.append(dict(sample=a.sample, interval=r["interval"],
                        svtype=r["svtype"], state=state, evidence=why,
                        covered_frac=f"{frac:.3f}", depth=("mapq" if filtered
                                                           else "gvcf"),
                        is6110_prox=r.get("is6110_prox", "0"),
                        n_projected=len(ok), n_probes=len(probes)))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "sample", "interval", "svtype", "state", "evidence",
            "covered_frac", "depth", "is6110_prox", "n_projected", "n_probes"])
        w.writeheader()
        w.writerows(out)
    n = len(out)
    print(f"  {a.sample} (ref {a.reference or '?'}): {n:,} intervals  "
          + "  ".join(f"{k} {v:,} ({v / n:.1%})" for k, v in counts.most_common()))
    for k, v in evid.most_common(6):
        print(f"    {k[:52]:<54}{v:>8,}")
    print(f"  -> {a.out}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sv-matrix", default="")
    ap.add_argument("--intervals", default="",
                    help="the catalogue from bin/sv_intervals.py. In this "
                         "mode every sample is genotyped at every interval "
                         "and the caller's per-sample breakpoints define "
                         "nothing -- see that script's header for why.")
    ap.add_argument("--reference", default="",
                    help="this sample's matched reference accession, so an "
                         "interval the reference itself carries is a lookup "
                         "rather than an inference")
    ap.add_argument("--sample", default="")
    ap.add_argument("--gvcf", default="")
    ap.add_argument("--projected", default="")
    ap.add_argument("--probes-out", default="")
    ap.add_argument("--h37rv-path", default="",
                    help="the graph's H37Rv path name, for --probes-out")
    ap.add_argument("--out", default="")
    ap.add_argument("--probes", type=int, default=9)
    ap.add_argument("--min-dp", type=int, default=5,
                    help="same floor p5_states.py uses for a SNP REF call")
    ap.add_argument("--present-frac", type=float, default=0.8,
                    help="fraction of interior probes that must be covered "
                         "before the deletion is called absent")
    ap.add_argument("--deleted-frac", type=float, default=0.2,
                    help="at or below this the depth says the interval IS "
                         "deleted; ALT (depth_absent) when a flank is "
                         "covered, unless --no-promote")
    ap.add_argument("--bam", default="",
                    help="the sample's alignment, for mapping-quality-filtered "
                         "depth. Required when --mapq-scope is not `never`.")
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS",
                                                         "samtools"))
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--mapq-scope", choices=("is6110", "all", "never"),
                    default="is6110",
                    help="where to use mapping-quality-filtered depth instead "
                         "of the gVCF's DP. The gVCF counts every read whatever "
                         "its mapping quality, and inside a repeat family with "
                         "up to about 25 copies per genome reads from other "
                         "copies mismap into the interval -- so a genuinely "
                         "deleted interval can read as covered and be called "
                         "REF. The project measured 46.1% of reads at MAPQ 0 in "
                         "element regions, so the default filters exactly the "
                         "intervals bin/sv_intervals.py flags as "
                         "element-proximal.")
    ap.add_argument("--min-projected", type=float, default=0.5,
                    help="below this fraction of interior probes projecting, "
                         "the sample's own reference has no such interval")
    ap.add_argument("--flank", type=int, default=200,
                    help="distance outside each breakpoint for the flank "
                         "probes that establish the sample has data here")
    ap.add_argument("--twoframe", default="",
                    help="a two-frame table from sv2frame/bin/sv_twoframe.py, "
                         "either this sample's or a merged cohort one. It "
                         "supplies the positive evidence an inherited "
                         "candidate now needs before it is promoted to ALT.")
    ap.add_argument("--trust-inherited", action="store_true",
                    help="promote an inherited candidate on the graph lookup "
                         "and a covered flank alone, as this script did before "
                         "the two-frame check measured a 26.5% contradiction "
                         "rate for that class. Kept for reproducing earlier "
                         "output, not for new runs.")
    ap.add_argument("--no-promote", action="store_true",
                    help="never write ALT from evidence, only from a call")
    a = ap.parse_args()

    if not (a.sv_matrix or a.intervals):
        sys.exit("FATAL: one of --sv-matrix or --intervals is required")
    if a.intervals:
        ivs = read_intervals(a.intervals)
        print(f"  {len(ivs):,} intervals")
        if a.probes_out:
            if not a.h37rv_path:
                sys.exit("FATAL: --probes-out needs --h37rv-path")
            seen = set()
            with open(a.probes_out, "w") as fh:
                for r in ivs:
                    for p in (probe_span(r["start"], r["end"], a.probes)
                              + flank_span(r["start"], r["end"], a.flank)):
                        if p >= 1 and p not in seen:
                            seen.add(p)
                            fh.write(f"{a.h37rv_path},{p - 1},+\n")
            print(f"  -> {a.probes_out}  ({len(seen):,} distinct positions)")
            return 0
        return genotype_intervals(a, ivs)

    rows, samples = read_sv_matrix(a.sv_matrix,
                                   keep_sample=None if a.probes_out else a.sample)
    dels = [r for r in rows if r["svtype"] == "DEL"]

    if a.probes_out:
        if not a.h37rv_path:
            sys.exit("FATAL: --probes-out needs --h37rv-path")
        seen = set()
        with open(a.probes_out, "w") as fh:
            for r in dels:
                ps = (probe_positions(int(r["h37rv_pos"]), r["svlen"],
                                      a.probes)
                      + flank_positions(int(r["h37rv_pos"]), r["svlen"],
                                        a.flank))
                for p in ps:
                    if 1 <= p and p not in seen:
                        seen.add(p)
                        fh.write(f"{a.h37rv_path},{p - 1},+\n")
        print(f"  {len(dels):,} DEL rows -> {len(seen):,} distinct probe "
              f"positions\n  -> {a.probes_out}")
        return 0

    for need in ("sample", "gvcf", "projected", "out"):
        if not getattr(a, need.replace("-", "_")):
            sys.exit(f"FATAL: --{need} is required when writing states")

    # Pair by SOURCE position, exactly as p5_states.py does and for the same
    # reason: odgi re-emits its header between records and a silently
    # misaligned coordinate is the worst thing this stage could produce.
    by_src = {}
    for line in open(a.projected):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            src = int(f[0].rsplit(",", 2)[-2]) + 1
            tgt = int(f[1].rsplit(",", 2)[-2]) + 1
            dist = int(f[2])
        except (ValueError, IndexError):
            continue
        by_src.setdefault(src, (tgt, dist))

    cov = gvcf_depth(a.gvcf, a.min_dp)
    counts = {"ALT": 0, "REF": 0, "ABSENT": 0, "NOCALL": 0}
    evid = {}
    out = []
    for r in rows:
        prior = r.get(a.sample, "NOCALL")
        if prior == "ALT":
            state, why = "ALT", "called"
        elif r["svtype"] != "DEL":
            state, why = prior, "not a deletion: depth cannot test absence"
        else:
            pos, svlen = int(r["h37rv_pos"]), r["svlen"]
            probes = probe_positions(pos, svlen, a.probes)
            ok = [t for t in (by_src.get(p) for p in probes)
                  if t is not None and t[1] == 0]
            fl = [t for t in (by_src.get(p)
                              for p in flank_positions(pos, svlen, a.flank))
                  if t is not None and t[1] == 0]
            flank_ok = any(cov(t) for t, _ in fl)
            if len(ok) < a.min_projected * len(probes):
                # THE INTERVAL IS NOT IN THIS SAMPLE'S OWN REFERENCE. For a
                # deletion that is not missing data -- it is the deletion,
                # already carried by the matched reference the sample was
                # aligned to, so no caller ever emits it. The flank probes are
                # what separate the two readings: reads either side mean the
                # sample really is sitting on a reference that lacks the
                # interval, which is the inherited ALT the matrix's `component`
                # column calls `inherited`. Without covered flanks the
                # projection may simply have failed, and that stays ABSENT.
                if flank_ok and not a.no_promote and a.trust_inherited:
                    state, why = "ALT", ("inherited: the interval is absent "
                                         "from this sample's reference and "
                                         "its flanks are covered")
                elif flank_ok and not a.no_promote:
                    # THE SAME RULE AS INTERVAL MODE (review 3.9). Promoting
                    # on the graph lookup plus a covered flank is what the
                    # two-frame check contradicted 26.5% of the time on
                    # scale200; interval mode stopped doing it, and this mode
                    # has no second frame to confirm with, so the inherited
                    # candidate is NOCALL. --trust-inherited restores the old
                    # behaviour for reproducing earlier output.
                    state, why = "NOCALL", ("inherited_unconfirmed: absent from "
                                            "this sample's reference, flanks "
                                            "covered, no second frame")
                else:
                    state, why = "ABSENT", (f"{len(ok)}/{len(probes)} interior "
                                            f"probes project, flanks "
                                            f"{'covered' if flank_ok else 'not covered'}")
            else:
                frac = sum(1 for t, _ in ok if cov(t)) / len(ok)
                if frac >= a.present_frac:
                    state, why = "REF", f"depth over {frac:.0%} of the interval"
                elif frac <= a.deleted_frac and flank_ok and not a.no_promote:
                    # The sample lacks sequence its own reference has, with
                    # reads either side of it. That is the same kind of
                    # measurement that earns REF two lines above, so refusing
                    # it while accepting REF would be asymmetric.
                    state, why = "ALT", (f"depth_absent: only {frac:.0%} of the "
                                         f"interval covered, flanks covered")
                elif frac <= a.deleted_frac:
                    state, why = prior, (f"only {frac:.0%} covered but the "
                                         f"flanks are not, so this may be a "
                                         f"gap in coverage")
                else:
                    state, why = prior, f"ambiguous depth, {frac:.0%} covered"
        counts[state] = counts.get(state, 0) + 1
        evid[why.split(":")[0]] = evid.get(why.split(":")[0], 0) + 1
        out.append(dict(sample=a.sample, key=r["key"], svtype=r["svtype"],
                        state=state, evidence=why))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "key", "svtype", "state",
                                           "evidence"], delimiter="\t")
        w.writeheader()
        w.writerows(out)
    n = len(out)
    print(f"  {a.sample}: {n:,} SV cells  " +
          "  ".join(f"{k} {v:,} ({v / n:.1%})"
                    for k, v in counts.items() if v))
    for k, v in sorted(evid.items(), key=lambda x: -x[1])[:5]:
        print(f"    {k:<46s}{v:>8,}")
    print(f"  -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
