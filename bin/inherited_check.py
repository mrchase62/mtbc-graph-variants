#!/usr/bin/env python3
"""Check what a sample INHERITS from its matched reference against its own reads.

THE PROBLEM (HANDOFF fix list, ABS-1 and ABS-2). Each sample is aligned to the
panel genome R nearest it, and wherever R differs from H37Rv the pipeline used
to report R's difference as the sample's without reading the sample:

  ABSENT   P5 states a key ABSENT when R has a deletion over its H37Rv
           position. Measured on scale200_fix, 10,247 of 85,465 testable
           ABSENT cells sat in sequence the sample covers at normal depth;
           in core sequence 2,279 + 161 were wrong, confirmed by local assembly
           at 48 of 49 sampled (analysis/inherited_absent/README.md). One was
           SAMEA2297133's intact ctpV (analysis/ctpV_check.md).
  inherited SVs
           an SV of R's is reported as the sample's. Deletions got one check,
           in the H37Rv frame only; insertions got none, so every inherited
           insertion was present in every sample matched to its R
           (analysis/inherited_matched_check/RESULTS.md).

This module holds the read-level checks that replace that assumption, one
definition each, used by P5's small-variant states (bin/p5_states.py), its SV
genotyping (bin/p5_inherited_sv.py, bin/p5_sv_genotype.py) and P4b's masked
label (bin/p4b_place_sv.py). The analysis scripts named above are the reference
versions; the rules and thresholds are theirs, unchanged.

THE MASK. A site is MASKED when its span, or either end of R's deletion within
150 bp (a read length), overlaps the build's repeat_mask.bed: 594 kb of PE/PPE,
paralog, tandem and IS-element sequence. In it depth, clips and split reads are
all unreliable -- 73% of the contradicted ABSENT cells were there -- so a masked
small-variant site is NOCALL with the reason "masked" (the user's choice,
2026-10-08: NOCALL with a reason, not a new state).

THE END RULE STAYS. It was dropped on 2026-10-14 and restored the same day
("yes, restore the end rule and merge"). Measured on ctrl25_p200's 25 finished
assemblies, ABS-1 with the span alone was 523 right / 187 wrong, against
472 / 11 with the end rule. The bridging reads at a deletion end inside a
repeat come from other copies: mar_02_R1896's phiRv1 locus (R deletion
1,779,267-1,788,513) had about 100 reads across both ends while the assembly
lacks phiRv1 there, and 31 of its cells were called REF. Inherited DELETIONS
inside the mask are genotyped (p5_sv_genotype.py) and masked inherited
INSERTIONS stay NOCALL; see there.

THE H37Rv-FRAME JUNCTION CHECK (ABS-1; analysis/inherited_absent/junction_check.py).
For an ABSENT site, R's deletion over it is read off the panel table
(bin/panel_differences.py): the largest R deletion overlapping the site. In
the sample's P1 alignment to H37Rv, at MAPQ >= 20:

  bridging   a primary read aligned without a gap from B-10 to B+10, at each
             deletion end B;
  joining    a read joining the base before the deletion to the base after
             it: a CIGAR deletion from one end to the other, or a soft clip
             >= 5 bp at one end whose supplementary alignment lands at the
             other.

  present_here   both ends bridged by >= max(5, 0.25 x the sample's median
                 depth) reads and fewer than 3 joining: the sequence is here,
                 so ABSENT is wrong and the sample's own H37Rv-frame call
                 stands (REF or ALT)
  deleted_here   3 or more joining and too few bridging: the locus is deleted
                 (depth over it, if any, comes from a copy elsewhere); ABSENT
                 stands
  mixed          both: a duplication, a mixture or a repeat. NOCALL
  unresolved     neither. NOCALL
  r_no_genotype  no R deletion over the site, and R has no genotype at an
                 overlapping graph site (its path does not cross it in place;
                 no single junction can test that). NOCALL
  no_r_deletion  nothing in R explains the call. NOCALL

Calibrated on scale200_fix against SPAdes reassembly of 120 regions: every
decided region agreed with the junction verdict, none contradicted it.

THE MATCHED-FRAME CHECK (ABS-2; analysis/inherited_matched_check/check.py).
An inherited SV is invisible in R's frame -- the sample's reads sit on R --
so it is tested THERE, at R's own junction. R's local haplotype (300 bp of
H37Rv either side of R's allele from the graph VCF) is aligned to R's FASTA
with minimap2 (asm10); the event must sit inside the alignment with 100 bp to
spare, which also checks the placement. Then, in the sample's P2 BAM (R's
frame, MAPQ >= 20):

  deletion (R has a junction j)
      carries    >= 5 reads span j-15..j+15 in one aligned block and fewer
                 than 15% are disrupted (a >= 10 bp soft clip ending within
                 10 bp of j, or a >= 10 bp indel within 15 bp of it)
      lacks      >= 3 disrupted and >= 30% of the reads at j
  insertion (R has extra sequence Y)
      carries    median depth across Y >= 0.5 x the flanks' and both of Y's
                 junctions spanned cleanly (each a "carries" by the rule above)
      lacks      Y's depth <= 0.1 x the flanks' at MAPQ >= 20 and <= 0.2 x at
                 any MAPQ, with >= 2 reads clipped at Y's ends
      repeat     MAPQ >= 20 depth <= 0.1 x but reads there at MAPQ 0
      low_depth  flank depth < 10;  unclear  anything else

TANDEM. When at least half the 20-mers of the changed sequence also occur
within 300 bp either side (a VNTR unit or a duplication), short reads aligned
to R look the same whether or not the sample has R's extra unit, so the event
is NOCALL "tandem" until a repeat-length genotyper is chosen (parked).

Measured on truth7 (seven finished assemblies): matched-frame "carries" was
right at 181/191 deletions (0.95, against 0.87 unchecked), 27/27 outside the
mask; for insertions outside the mask and not tandem 24/26 (0.92, against
0.84). The insertion rules were revised after a first look at truth7 (the
clip requirement for "lacks", clean junctions for "carries"), so those figures
are not blind; the deletion rules were not changed.
"""
import bisect, collections, gzip, hashlib, os, statistics, subprocess, sys, tempfile

MIN_MAPQ = 20           # both checks, as measured
TOL = 10                # H37Rv junction: bridging window and clip tolerance
MIN_CLIP = 5            # H37Rv junction: clip length that counts
MIN_JOIN = 3            # joining reads that say "deleted here"
FLANK = 300             # R's local haplotype: H37Rv bases either side
MASK_PAD = 150          # a read length either side of R's deletion ends
H37RV_ACC = "GCF_000195955"

VERDICTS = ("present_here", "deleted_here", "mixed", "unresolved",
            "r_no_genotype", "no_r_deletion", "masked")


# ---- the repeat mask --------------------------------------------------------
def load_mask(path):
    """repeat_mask.bed as merged 1-based closed intervals. FATAL when absent:
    without it every masked site would be called as if it were core."""
    if not path or not os.path.exists(path):
        sys.exit(f"FATAL: no repeat mask at {path!r}; it is the build's "
                 f"assets/repeat_mask.bed")
    iv = []
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3 and f[1].strip().isdigit():
                iv.append((int(f[1]) + 1, int(f[2])))    # BED -> 1-based closed
    if not iv:
        sys.exit(f"FATAL: {path} has no intervals")
    iv.sort()
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


class Mask:
    """overlaps(s, e): does 1-based closed [s, e] touch the mask?"""

    def __init__(self, intervals):
        self.iv = intervals
        self.starts = [s for s, _ in intervals]

    @classmethod
    def from_bed(cls, path):
        return cls(load_mask(path))

    def overlaps(self, s, e):
        if e < s:
            s, e = e, s
        i = bisect.bisect_right(self.starts, e) - 1
        return i >= 0 and self.iv[i][1] >= s       # merged, so one test suffices

    def site(self, s, e, r_del=None, pad=MASK_PAD):
        """The ABS-1 rule: the site's span, or either end of R's deletion
        within `pad` bp, overlaps the mask."""
        if self.overlaps(s, e):
            return True
        if r_del is not None:
            ds, de = r_del[0], r_del[1]
            return self.overlaps(ds - pad, ds + pad) or self.overlaps(de - pad, de + pad)
        return False


# ---- the panel's own differences (P0 step panel_diffs) ----------------------
def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class PanelDiffs:
    """One accession's rows of <build>/assets/panel_differences.tsv.gz.

    dels: sorted (start, end, vcf_pos, allele) of deleted H37Rv bases
    ins:  {(vcf_pos, svlen): allele} of insertions >= 50 bp
    nogt: sorted (start, end) of graph sites where it has no genotype

    `graph_vcf`, when given, must be the VCF the table was read from (its
    sha256 is in the table's header): P4b, P5 and this table read ONE graph
    VCF per run (ABS-2 item 6)."""

    def __init__(self, path, accession, graph_vcf=None):
        if not path or not os.path.exists(path):
            sys.exit(f"FATAL: no panel difference table at {path!r}; run "
                     f"bin/p0_prepare.sh --step panel_diffs")
        self.accession = accession
        self.dels, self.nogt, self.ins = [], [], {}
        hdr, panel = {}, None
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    k, _, v = line[1:].rstrip("\n").partition("\t")
                    hdr[k] = v
                    continue
                if line.startswith("accession\t"):
                    continue
                if not line.startswith(accession + "\t"):
                    continue
                f = line.rstrip("\n").split("\t")
                s, e, ln, vp, al = int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[7])
                if f[1] == "DEL":
                    self.dels.append((s, e, vp, al))
                elif f[1] == "INS":
                    self.ins[(vp, ln)] = al
                elif f[1] == "NOGT":
                    self.nogt.append((s, e))
        panel = set((hdr.get("accessions") or "").split(",")) - {""}
        if not panel:
            sys.exit(f"FATAL: {path} has no #accessions header; regenerate it")
        if accession != H37RV_ACC and accession not in panel:
            sys.exit(f"FATAL: {accession} is not a panel genome of {path}; the "
                     f"table is from another build")
        if graph_vcf is not None:
            want = hdr.get("graph_vcf_sha256", "")
            have = file_sha256(graph_vcf)
            if want != have:
                sys.exit(f"FATAL: {path} was made from a graph VCF with sha256 "
                         f"{want[:12] or '<none>'}, not {graph_vcf} "
                         f"({have[:12]}); rerun bin/p0_prepare.sh --step panel_diffs")
        self.dels.sort()
        self.nogt.sort()
        self._ds = [d[0] for d in self.dels]
        self._dreach = _reach([d[1] for d in self.dels])
        self._ns = [n[0] for n in self.nogt]
        self._nreach = _reach([n[1] for n in self.nogt])

    def deletions_over(self, s, e):
        """R's deletions overlapping [s, e]."""
        return _over(self.dels, self._ds, self._dreach, s, e)

    def deletion_over(self, s, e):
        """The LARGEST of R's deletions overlapping [s, e], or None
        (junction_check.py's rule)."""
        hits = self.deletions_over(s, e)
        return max(hits, key=lambda d: (d[1] - d[0], -d[0])) if hits else None

    def no_genotype_over(self, s, e):
        return bool(_over(self.nogt, self._ns, self._nreach, s, e))


def _reach(ends):
    out, m = [], 0
    for e in ends:
        m = max(m, e)
        out.append(m)
    return out


def _over(spans, starts, reach, s, e):
    j = bisect.bisect_right(starts, e) - 1
    out = []
    while j >= 0 and reach[j] >= s:
        if spans[j][1] >= s:
            out.append(spans[j])
        j -= 1
    return out


def fetch_alleles(graph_vcf, wanted):
    """{(vcf_pos, allele): (REF, ALT)} for the graph records named, read by
    tabix from the same VCF the panel table was made from."""
    import pysam
    if not wanted:
        return {}
    tb = pysam.TabixFile(graph_vcf)
    ctg = tb.contigs[0]
    out = {}
    for vp in sorted({p for p, _ in wanted}):
        for line in tb.fetch(ctg, vp - 1, vp):
            f = line.split("\t", 5)
            if int(f[1]) != vp:
                continue
            alts = f[4].split(",")
            for p, al in wanted:
                if p == vp and 0 < al <= len(alts):
                    out[(p, al)] = (f[3].upper(), alts[al - 1].upper())
    return out


# ---- the H37Rv frame: the sample's P1 alignment (ABS-1) --------------------
def _aligned_across(read, a, b):
    """The read is aligned without a gap over [a, b] (1-based)."""
    for bs, be in read.get_blocks():               # 0-based half-open
        if bs + 1 <= a and be >= b:
            return True
    return False


def junction_reads(bam, ctg, s, e, min_mapq=MIN_MAPQ):
    """(bridging s, bridging e, clipped at s, clipped at e, joining reads)
    for a deletion of H37Rv bases s..e: junction_check.py's
    boundary_evidence(), unchanged."""
    left, right = s - 1, e + 1                     # last kept base, first kept base
    bridge_s = bridge_e = clip_l = clip_r = 0
    joined, seen = set(), set()
    lo, hi = max(0, s - 300), e + 300
    windows = [(lo, s + 300)] if e - s < 600 else [(lo, s + 300), (e - 300, hi)]
    for wlo, whi in windows:
        for r in bam.fetch(ctg, max(0, wlo), whi):
            if (r.is_unmapped or r.is_secondary or r.is_supplementary
                    or r.mapping_quality < min_mapq):
                continue
            key = (r.query_name, r.is_read1)
            if key in seen:
                continue
            seen.add(key)
            if _aligned_across(r, s - TOL, s + TOL):
                bridge_s += 1
            if _aligned_across(r, e - TOL, e + TOL):
                bridge_e += 1
            ct = r.cigartuples
            if not ct:
                continue
            p = r.reference_start
            for op, n in ct:
                if op == 2 and n >= 1 and abs(p - left) <= TOL \
                        and abs(p + n + 1 - right) <= TOL:
                    joined.add(r.query_name)
                if op in (0, 2, 3, 7, 8):
                    p += n
            sa = r.get_tag("SA").rstrip(";").split(";") if r.has_tag("SA") else []
            sa_pos = [int(x.split(",")[1]) for x in sa if x.split(",")[0] == ctg]
            if ct[-1][0] == 4 and ct[-1][1] >= MIN_CLIP \
                    and abs(r.reference_end - left) <= TOL:
                clip_r += 1
                if any(abs(q - right) <= 2 * TOL + ct[-1][1] for q in sa_pos):
                    joined.add(r.query_name)
            if ct[0][0] == 4 and ct[0][1] >= MIN_CLIP \
                    and abs(r.reference_start + 1 - right) <= TOL:
                clip_l += 1
                if any(abs(q - left) <= 2 * TOL + ct[0][1] for q in sa_pos):
                    joined.add(r.query_name)
    return bridge_s, bridge_e, clip_r, clip_l, len(joined)


def junction_verdict(bridge_s, bridge_e, joined, min_bridge):
    bridge = min(bridge_s, bridge_e)
    if joined >= MIN_JOIN and bridge >= min_bridge:
        return "mixed"
    if joined >= MIN_JOIN:
        return "deleted_here"
    if bridge >= min_bridge:
        return "present_here"
    return "unresolved"


class H37rvFrame:
    """The sample's P1 alignment to H37Rv, and its P1 H37Rv-frame calls.

    Depth is samtools depth -a at MAPQ >= 20 over the whole chromosome, read
    once: it gives the sample's median (the bridging floor scales with it) and
    the depth behind an H37Rv-frame REF. `vcf` is P1's HaplotypeCaller VCF,
    variants only, so REF is "no call here and covered"."""

    def __init__(self, bam, samtools, vcf=None, h37=None, min_mapq=MIN_MAPQ):
        import numpy as np
        import pysam
        if not bam or not os.path.exists(bam):
            sys.exit(f"FATAL: no P1 H37Rv alignment at {bam!r}; the ABSENT "
                     f"check reads it (restore it with bin/archive_alignments.sh "
                     f"--restore, or pass the opt-out flag to reproduce older output)")
        self.bam = pysam.AlignmentFile(bam)
        self.ctg = self.bam.references[0]
        self.min_mapq = min_mapq
        L = self.bam.lengths[0]
        p = subprocess.run([samtools, "depth", "-a", "-Q", str(min_mapq),
                            "-r", self.ctg, bam],
                           check=True, capture_output=True, text=True)
        d = np.zeros(L + 2, dtype=np.int32)
        for line in p.stdout.splitlines():
            f = line.split("\t")
            if len(f) >= 3:
                d[int(f[1])] = int(f[2])
        if not d.any():
            sys.exit(f"FATAL: samtools depth returned no depth for {bam}")
        self.depth = d
        self.baseline = float(np.median(d[1:L + 1]))
        self.min_bridge = max(5.0, 0.25 * self.baseline)
        self._junc = {}
        self.calls, self.call_spans = set(), []
        if vcf is not None:
            self._load_calls(vcf, h37)

    def _load_calls(self, vcf, h37):
        if not os.path.exists(vcf):
            sys.exit(f"FATAL: no P1 H37Rv-frame VCF at {vcf}")
        from mtb_norm import normalise, h37rv_key
        op = gzip.open if vcf.endswith(".gz") else open
        spans = []
        with op(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                gt = f[9].split(":", 1)[0].replace("|", "/").split("/")[0]
                if not gt.isdigit() or gt == "0":
                    continue
                alts = f[4].split(",")
                i = int(gt)
                if i > len(alts) or alts[i - 1].startswith("<") or alts[i - 1] == "*":
                    continue
                pos, ref = int(f[1]), f[3].upper()
                np_, nr, na = normalise(h37, pos, ref, alts[i - 1].upper())
                self.calls.add(h37rv_key(np_, nr, na))
                spans.append((pos, pos + len(ref) - 1))
        spans.sort()
        self.call_spans = spans
        self._cs = [s for s, _ in spans]
        self._creach = _reach([e for _, e in spans])

    def junction(self, s, e):
        """(verdict, evidence dict) for R's deletion of H37Rv bases s..e."""
        k = (s, e)
        if k not in self._junc:
            bs, be, cr, cl, jn = junction_reads(self.bam, self.ctg, s, e,
                                                self.min_mapq)
            self._junc[k] = (junction_verdict(bs, be, jn, self.min_bridge),
                             dict(bridge_start=bs, bridge_end=be,
                                  clip_start=cr, clip_end=cl, joined=jn))
        return self._junc[k]

    def covered(self, s, e, min_dp):
        return bool((self.depth[s:e + 1] >= min_dp).all()) if e >= s else False

    def small_call(self, key, s, e, min_dp):
        """The sample's own H37Rv-frame state for an allele key spanning H37Rv
        s..e: ALT when P1 called that allele, NOCALL when it called another
        allele over the span, REF when nothing was called and every base has
        MAPQ >= 20 depth >= min_dp."""
        if key in self.calls:
            return "ALT", "present_here_alt"
        if self.call_spans and _over(self.call_spans, self._cs, self._creach, s, e):
            return "NOCALL", "present_other_allele"
        if self.covered(s, e, min_dp):
            return "REF", "present_here_ref"
        return "NOCALL", "present_low_depth"


# ---- a deletion the sample did not call (2026-10-16) -----------------------
DEL_RATIO = 0.2         # depth inside / flanks below which nothing is inside
DEL_MAX_BRIDGE = 2      # reads aligned across an end that a deletion allows
DEL_REPEAT = 3.0        # MAPQ>=0 depth this many times MAPQ>=20: a repeat
DEL_FLANK_REL = 0.25    # flanks below this x the sample's median: dropout


def check_deletion(h37f, mask, s, e, min_dp):
    """(state, evidence, junction counts) for H37Rv bases s..e, a deletion of
    50 bp or more the sample has no call for and whose anchor its gVCF does
    not cover (p5_states.py, `no_coverage`).

    The rule is the read check that adjudicated the 1K evaluation's
    disagreements (analysis/sv_eval_1k/readcheck.py check_del), in the
    sample's P1 alignment to H37Rv at MAPQ >= 20:

      flanks    mean depth over F = clamp(L, 100, 500) bases either side
      inside    mean depth over s..e, at MAPQ >= 20 and at any MAPQ
      junction  junction_reads(): reads aligned without a gap across each end
                (bridging: the mean of the two ends) and reads joining the
                two ends (joined)

      flanks < min_dp                                NOCALL no_coverage
          the flanks are not covered either, so nothing says the missing
          depth is the deletion
      joined >= 3, joined >= bridging                ALT del_junction_joined
      inside < 0.2 x flanks, reads there only at
      MAPQ 0 (a repeat)                              NOCALL del_repeat_mapq0
      inside < 0.2 x flanks, but the site in the
      repeat mask or the flanks below 0.25 x the
      sample's median depth                          NOCALL del_depth_dropout
      inside < 0.2 x flanks, bridging <= 2           ALT del_no_depth_inside
      both ends bridged by >= the bridging floor,
      joined < 3, inside >= 0.6 x flanks, unmasked   REF del_spanned_ref
      anything else                                  NOCALL del_depth_unclear

    "Not a repeat": the any-MAPQ depth inside is not above 3 x the MAPQ >= 20
    depth (and 0.2 x the flanks), as readcheck.py has it.

    DEPTH ALONE IS NOT A DELETION IN LOW-COVERAGE SEQUENCE. Measured on the
    finished assemblies before this rule shipped: every depth-only ALT the
    readcheck.py rule gave was inside the repeat mask with flanks at about a
    tenth of the sample's median depth -- GC-rich PE_PGRS sequence, where depth
    falls away in every isolate -- and on marin36_p200 all 6 such cells were
    REF in the assembly (ctrl25_p200 had none). The 20 read-checked cells of
    the 1K evaluation were all of that kind too: the 5 it called clear
    deletions had flanks of 7-21x against medians of 87-94x. So depth alone
    decides ALT only outside the mask and with flanks at a normal share of
    the sample's depth; elsewhere reads must join the two ends. REF needs the
    site outside the mask as well, where bridging reads can come from another
    copy (THE END RULE STAYS, above)."""
    import numpy as np
    L = e - s + 1
    F = int(min(500, max(100, L)))
    d = h37f.depth
    lo = d[max(1, s - F):s]
    hi = d[e + 1:min(len(d) - 1, e + F + 1)]
    flank = float(np.concatenate([lo, hi]).mean()) if len(lo) + len(hi) else 0.0
    inside = float(d[s:e + 1].mean()) if e >= s else 0.0
    _v, jev = h37f.junction(s, e)
    jev = dict(jev, flank=round(flank, 1), inside=round(inside, 1))
    if flank < min_dp:
        return "NOCALL", "no_coverage", jev
    joined = jev["joined"]
    bridge = (jev["bridge_start"] + jev["bridge_end"]) / 2   # readcheck's mean
    if joined >= MIN_JOIN and joined >= bridge:
        return "ALT", "del_junction_joined", jev
    if inside < DEL_RATIO * flank:
        inside0 = _med_depth(h37f.bam, h37f.ctg, s - 1, e, 0, mean=True)
        jev["inside_mapq0"] = round(inside0, 1)
        repeat = inside0 > DEL_REPEAT * max(inside, 0.5) and inside0 > DEL_RATIO * flank
        if repeat:
            return "NOCALL", "del_repeat_mapq0", jev
        if mask.overlaps(s, e) or flank < DEL_FLANK_REL * h37f.baseline:
            return "NOCALL", "del_depth_dropout", jev
        if bridge <= DEL_MAX_BRIDGE:
            return "ALT", "del_no_depth_inside", jev
        return "NOCALL", "del_depth_unclear", jev
    if (min(jev["bridge_start"], jev["bridge_end"]) >= h37f.min_bridge
            and joined < MIN_JOIN and inside >= 0.6 * flank
            and not mask.overlaps(s, e)):
        return "REF", "del_spanned_ref", jev
    return "NOCALL", "del_depth_unclear", jev


def check_absent(diffs, mask, h37f, s, e):
    """ABS-1 for a site spanning H37Rv s..e that P5 found absent from R.

    Returns (verdict, R's deletion or None, evidence dict); verdict is one of
    VERDICTS. The caller turns present_here into the sample's own H37Rv-frame
    call, keeps ABSENT on deleted_here and states NOCALL with the verdict as
    its reason for everything else."""
    d = diffs.deletion_over(s, e)
    if mask.site(s, e, d):
        return "masked", d, {}
    if d is None:
        return ("r_no_genotype" if diffs.no_genotype_over(s, e)
                else "no_r_deletion"), None, {}
    v, ev = h37f.junction(d[0], d[1])
    return v, d, ev


# ---- the matched frame: the sample's P2 alignment to R (ABS-2) -------------
def _rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _cigar_ops(cg):
    ops, n = [], ""
    for ch in cg:
        if ch.isdigit():
            n += ch
        else:
            ops.append((int(n), ch))
            n = ""
    return ops


def local_haplotype(h37, vcf_pos, ref, alt):
    """(query, a, b, changed sequence) for R's allele at a graph record:
    FLANK bases of H37Rv, R's allele, FLANK bases of H37Rv. [a, b) is the
    changed part in the query after the shared prefix and suffix are trimmed
    (a == b for a deletion: R's junction)."""
    p = 0
    while p < min(len(ref), len(alt)) and ref[p] == alt[p]:
        p += 1
    s = 0
    while s < min(len(ref), len(alt)) - p and ref[-1 - s] == alt[-1 - s]:
        s += 1
    left = h37[max(0, vcf_pos - 1 - FLANK):vcf_pos - 1]
    right = h37[vcf_pos - 1 + len(ref):vcf_pos - 1 + len(ref) + FLANK]
    q = left + alt + right
    a = len(left) + p
    b = len(left) + len(alt) - s
    chg = q[a:b] if b > a else h37[vcf_pos - 1 + p:vcf_pos - 1 + len(ref) - s]
    return q, a, b, chg


def is_tandem(q, a, b, chg, k=20, frac=0.5):
    """At least `frac` of the changed sequence's k-mers also occur within the
    haplotype's FLANK bases either side (check.py's rule)."""
    fl = q[:a] + "N" + q[b:]
    km = {fl[i:i + k] for i in range(len(fl) - k + 1)}
    ck = [chg[i:i + k] for i in range(max(0, len(chg) - k + 1))]
    return bool(ck) and sum(x in km for x in ck) / len(ck) >= frac


def place_on_r(events, h37, r_fasta, minimap2, threads=1):
    """Place R's alleles on R's own sequence.

    events: {id: (vcf_pos, ref, alt)}. Returns {id: dict(status, contig,
    start, end, strand, mapq, tandem)}; start..end are R coordinates of the
    changed part (start == end at a deletion's junction). check.py's place()."""
    if not os.path.exists(r_fasta):
        sys.exit(f"FATAL: no FASTA for the matched reference at {r_fasta}")
    out, qs = {}, {}
    for k, (vp, ref, alt) in events.items():
        q, a, b, chg = local_haplotype(h37, vp, ref, alt)
        qs[k] = (q, a, b)
        out[k] = dict(status="unplaced", tandem="yes" if is_tandem(q, a, b, chg) else "no")
    if not qs:
        return out
    with tempfile.TemporaryDirectory(prefix="inhchk_") as tmp:
        fa = os.path.join(tmp, "q.fa")
        names = {}
        with open(fa, "w") as fh:
            for i, (k, (q, _a, _b)) in enumerate(qs.items()):
                names[f"q{i}"] = k
                fh.write(f">q{i}\n{q}\n")
        res = subprocess.run([minimap2, "-x", "asm10", "-c", "--secondary=no",
                              "-t", str(threads), r_fasta, fa],
                             capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit(f"FATAL: minimap2 failed on {r_fasta}:\n{res.stderr[-2000:]}")
    best = {}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if f[0] not in best:
            best[f[0]] = f
    for qn, f in best.items():
        k = names[qn]
        q, a, b = qs[k]
        qlen, qs_, qe, strand = int(f[1]), int(f[2]), int(f[3]), f[4]
        ctg, ts, mapq = f[5], int(f[7]), int(f[11])
        if qs_ > a - 100 or qe < b + 100:
            continue                                 # event not inside the alignment
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
        if cg is None:
            continue
        qa, qb = (a, b) if strand == "+" else (qlen - b, qlen - a)
        qpos, tpos, m = (qs_ if strand == "+" else qlen - qe), ts, {}
        for n, op in _cigar_ops(cg):
            if op == "M":
                for x in (qa, qb):
                    if qpos <= x <= qpos + n:
                        m.setdefault(x, tpos + (x - qpos))
                qpos += n
                tpos += n
            elif op == "I":
                for x in (qa, qb):
                    if qpos <= x <= qpos + n:
                        m.setdefault(x, tpos)
                qpos += n
            elif op == "D":
                tpos += n
        if qa not in m or qb not in m:
            continue
        out[k].update(status="placed", contig=ctg, start=min(m[qa], m[qb]),
                      end=max(m[qa], m[qb]), strand=strand, mapq=mapq)
    return out


def _med_depth(bam, ctg, s, e, minq, mean=False):
    """Median (or mean) depth over 0-based [s, e) at MAPQ >= minq."""
    if e <= s:
        return 0.0
    cov = bam.count_coverage(ctg, s, e, quality_threshold=0,
                             read_callback=lambda r: (not r.is_secondary
                                                      and not r.is_duplicate
                                                      and not r.is_supplementary
                                                      and r.mapping_quality >= minq))
    tot = [sum(x) for x in zip(*cov)]
    if not tot:
        return 0.0
    return float(statistics.mean(tot) if mean else statistics.median(tot))


def matched_junction(bam, ctg, j, L):
    """R's deletion junction at R position j: reads spanning it cleanly
    against reads disrupted there (check.py's del_metrics)."""
    span = dis = 0
    for r in bam.fetch(ctg, max(0, j - 15), min(L, j + 15)):
        if (r.mapping_quality < MIN_MAPQ or r.is_secondary or r.is_duplicate
                or r.is_supplementary):
            continue
        ct = r.cigartuples or []
        bad = False
        if ct and ct[0][0] == 4 and ct[0][1] >= 10 and abs(r.reference_start - j) <= 10:
            bad = True
        if ct and ct[-1][0] == 4 and ct[-1][1] >= 10 and abs(r.reference_end - j) <= 10:
            bad = True
        p = r.reference_start
        for op, n in ct:
            if op in (0, 7, 8):
                p += n
            elif op == 2:
                if n >= 10 and p - 15 <= j <= p + n + 15:
                    bad = True
                p += n
            elif op == 1:
                if n >= 10 and abs(p - j) <= 15:
                    bad = True
        if bad:
            dis += 1
        elif any(bs <= j - 15 and be >= j + 15 for bs, be in r.get_blocks()):
            span += 1
    n = span + dis
    if n < 5:
        call = "low_depth"
    elif dis >= 3 and dis / n >= 0.3:
        call = "lacks"
    elif span >= 5 and dis / n < 0.15:
        call = "carries"
    else:
        call = "unclear"
    return dict(span=span, disrupted=dis, call=call)


def matched_insertion(bam, ctg, s, e, L):
    """R's inserted sequence Y at R s..e: depth across Y against its flanks,
    clips at its ends, and both junctions (check.py's ins_metrics, with its
    2026-10-13 revision)."""
    left = _med_depth(bam, ctg, max(0, s - 300), max(0, s - 50), MIN_MAPQ)
    right = _med_depth(bam, ctg, min(L, e + 50), min(L, e + 300), MIN_MAPQ)
    flank = (left + right) / 2
    y20 = _med_depth(bam, ctg, s, e, MIN_MAPQ)
    y0 = _med_depth(bam, ctg, s, e, 0)
    clips = 0
    for x in (s, e):
        for r in bam.fetch(ctg, max(0, x - 20), min(L, x + 20)):
            if r.mapping_quality < MIN_MAPQ or r.is_secondary or r.is_duplicate:
                continue
            ct = r.cigartuples or []
            if ct and ct[0][0] == 4 and ct[0][1] >= 10 and abs(r.reference_start - x) <= 10:
                clips += 1
            if ct and ct[-1][0] == 4 and ct[-1][1] >= 10 and abs(r.reference_end - x) <= 10:
                clips += 1
    js = matched_junction(bam, ctg, s, L)
    je = matched_junction(bam, ctg, e, L)
    ends_clean = all(j["call"] == "carries" for j in (js, je))
    if flank < 10:
        call = "low_depth"
    elif y20 / flank >= 0.5 and ends_clean:
        call = "carries"
    elif y20 / flank <= 0.1 and y0 / flank <= 0.2 and clips >= 2:
        call = "lacks"
    elif y20 / flank <= 0.1 and y0 / flank > 0.2:
        call = "repeat"
    else:
        call = "unclear"
    return dict(flank=flank, y20=y20, y0=y0, clips=clips,
                span=f"{js['span']}/{je['span']}",
                disrupted=f"{js['disrupted']}/{je['disrupted']}", call=call)


MATCHED_STATE = {"carries": "ALT", "lacks": "REF"}


def combine_frames(h_state, m_call):
    """ABS-2's rule for an inherited deletion seen in both frames.

    Both frames agree -> that state; one decides and the other is NOCALL ->
    that state, saying which frame decided; they conflict -> NOCALL
    "frames_conflict". Returns (state, which)."""
    m = MATCHED_STATE.get(m_call, "NOCALL")
    h = h_state if h_state in ("ALT", "REF") else "NOCALL"
    if h == m:
        return h, ("two_frames" if h != "NOCALL" else "neither_frame")
    if h == "NOCALL":
        return m, "matched_only"
    if m == "NOCALL":
        return h, "h37rv_only"
    return "NOCALL", "frames_conflict"
