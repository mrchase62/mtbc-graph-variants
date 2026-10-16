#!/usr/bin/env python3
"""Check an IS6110 stage-2 ABSENT against the isolate's own reads (fix SV-EVAL-1).

THE PROBLEM. Stage 2 (is6110_p5_stage2.py) projects a carrier's position onto
each non-carrier's matched reference R and, until 2026-10-16, wrote ABSENT
whenever odgi's walk did not land on R's path (dist != 0). It read nothing.
Measured on new1k_p200 (analysis/sv_eval_1k/README.md, problem 1):

  - 46 H37Rv-frame sites were ABSENT in >= 90% of isolates (45,024 cells).
    At a `ref_shared` site the carrier's position sits beside the element
    its own reference carries, so on every reference WITHOUT the element the
    walk crosses the element and stops ~1,355 bp later: 4,230,303 had dist
    1,358 in all 997 non-carriers, while SAMD00089643 has 116x MAPQ >= 20
    depth straight across the site. "R lacks the element" is the non-carrier
    state, not a deleted region.
  - a dist of 1-10 (a SNP bubble beside the site) was ABSENT too: 27,147 of
    the 96,916 H37Rv-frame ABSENT cells.
  - 4 OLD ALT cells were NEW ABSENT where the reads show an IS6110 junction
    at the edge of a deletion: the region is gone AND the element is there.

THE CHECK, in the sample's P1 alignment to H37Rv (MAPQ >= 20, no duplicate,
secondary or QC-fail reads), at the key's H37Rv position P (the base after
which the element would sit), within W = 20 bp:

  junction   a read soft- or hard-clipped >= 15 bp with its clip point within
             W of the site (right clip ending at P, left clip starting at
             P + 1) whose clipped side is IS6110: either the 15 clipped bases
             next to the clip point are IS6110 15-mers (either strand; the
             element taken from an H37Rv copy), or an SA segment continues
             into an IS6110 terminus (is6110_junctions.far_sides + attribute,
             within W of a copy's end). Counted once per read name.
  spanning   a primary read aligned without a gap over P-W .. P+W.

  ALT     junction >= 3 and junction / (junction + spanning) >= 0.5
  REF     junction <= 1 and spanning >= max(8, 0.25 x the sample's median
          MAPQ >= 20 depth)
  NOCALL  otherwise: absent_low_depth (fewer reads in all than that floor),
          absent_mixed (>= 3 junction reads but under half the reads),
          absent_unresolved (anything else)
  ABSENT  no reads at the site at all, and ABS-1's test of R's deletion over
          it (reads joining the deletion's ends) says deleted_here: EV
          junction_joined, as for a small variant (use_region_check)

These are the rules of analysis/sv_eval_1k/readcheck.py check_is, which
agreed with the reads at 29/30 NEW IS6110 gains and 8/13 OLD-only calls, with
two changes:
  - the REF floor is ABS-1's depth-scaled one (inherited_check.H37rvFrame.
    min_bridge uses the same 0.25 x median), not a flat 8;
  - the ALT share is 0.5, not 0.2. Measured on ctrl25_p200 against the 25
    finished assemblies: at 888,777, 244 bp from H37Rv's IS6110-1, five
    isolates had 19-50 junction reads against 31-118 spanning (share
    0.18-0.34) and their assemblies have no insertion there. An isolate is
    haploid; a real insertion leaves almost no read across the site (every
    other junction ALT there had 0-3 spanning reads).

THE MASK. A REF verdict inside the build's repeat_mask.bed (the site P..P+1)
is NOCALL `absent_masked`, as ABS-1 does for small variants: spanning reads in
a repeat may come from another copy. A junction ALT is kept in the mask: a
read clipped INTO IS6110 at this position is not produced by a paralog
without the element. The measurement is in the fork report of 2026-10-16.

The sample's median depth is estimated from MAPQ >= 20 depth at 400 evenly
spaced positions, not a whole-chromosome samtools depth: it only scales a
floor, and stage 2 runs once per sample.
"""
import os
import statistics
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import is6110_junctions as ij  # noqa: E402

MIN_MAPQ = 20
WIN = 20                # clip point / spanning window either side of P
MIN_CLIP = 15           # clip length that counts, and the k-mer tested
MIN_JUNC = 3            # junction reads for ALT
MIN_FRAC = 0.5          # ... and their share of junction + spanning
MAX_JUNC_REF = 1        # junction reads tolerated in a REF
MIN_SPAN = 8            # spanning floor, before depth scaling
SPAN_FRAC = 0.25        # x the sample's median depth
N_DEPTH_PROBES = 400

EV_CODES = ("absent_checked_alt", "absent_checked_ref", "absent_masked",
            "absent_low_depth", "absent_mixed", "absent_unresolved",
            "junction_joined", "absent_no_r_deletion", "absent_r_no_genotype")


def _rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta1(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                buf.append(line.strip())
    return name, "".join(buf).upper()


def element_kmers(seq, k=MIN_CLIP):
    return ({seq[i:i + k] for i in range(len(seq) - k + 1)}
            | {_rc(seq)[i:i + k] for i in range(len(seq) - k + 1)})


def verdict(junc, span, floor):
    """(state, ev) from the counts; the mask is applied by the caller."""
    n = junc + span
    if junc >= MIN_JUNC and junc / n >= MIN_FRAC:
        return "ALT", "absent_checked_alt"
    if junc <= MAX_JUNC_REF and span >= floor:
        return "REF", "absent_checked_ref"
    if n < floor:
        return "NOCALL", "absent_low_depth"
    if junc >= MIN_JUNC:
        return "NOCALL", "absent_mixed"
    return "NOCALL", "absent_unresolved"


class H37rvSiteCheck:
    """One sample's P1 H37Rv alignment, read at IS6110 sites."""

    def __init__(self, bam, gff, h37rv_fasta, mask=None):
        import pysam
        if not bam or not os.path.exists(bam):
            sys.exit(f"FATAL: no P1 H37Rv alignment at {bam!r}; the IS6110 "
                     f"ABSENT check reads it (bin/archive_alignments.sh "
                     f"--restore brings back an archived one, or pass "
                     f"--no-absent-check to reproduce older output)")
        for p, what in ((gff, "IS6110 copy GFF"), (h37rv_fasta, "H37Rv FASTA")):
            if not p or not os.path.exists(p):
                sys.exit(f"FATAL: no {what} at {p!r}; the IS6110 ABSENT check needs it")
        self.bam = pysam.AlignmentFile(bam)
        self.ctg = self.bam.references[0]
        self.els = ij.load_elements(gff)                 # 0-based inclusive
        if not self.els:
            sys.exit(f"FATAL: {gff} lists no IS6110 copy")
        _, h37 = read_fasta1(h37rv_fasta)
        s0, e0 = self.els[0]
        self.kmers = element_kmers(h37[s0:e0 + 1])
        self.mask = mask
        L = self.bam.lengths[0]
        d = []
        for i in range(N_DEPTH_PROBES):
            p = (i + 0.5) * L / N_DEPTH_PROBES
            cov = self.bam.count_coverage(self.ctg, int(p), int(p) + 1,
                                          quality_threshold=0,
                                          read_callback=self._ok_primary)
            d.append(sum(x[0] for x in cov))
        self.baseline = float(statistics.median(d))
        self.floor = max(MIN_SPAN, SPAN_FRAC * self.baseline)

    @staticmethod
    def _ok(r):
        return not (r.is_unmapped or r.is_secondary or r.is_duplicate
                    or r.is_qcfail) and r.mapping_quality >= MIN_MAPQ

    @classmethod
    def _ok_primary(cls, r):
        return cls._ok(r) and not r.is_supplementary

    def _clip_is_element(self, r, at_left, n):
        """The clipped side of r (left when at_left, n bases) is IS6110."""
        q = r.query_sequence or ""
        ct = r.cigartuples
        if (ct[0] if at_left else ct[-1])[0] == 4 and len(q) >= n:
            km = q[n - MIN_CLIP:n] if at_left else q[len(q) - n:len(q) - n + MIN_CLIP]
            if km in self.kmers:
                return True
        far = {}
        for f in ij.far_sides(r, at_left):
            far[f] = far.get(f, 0) + 1
        hit, _, _ = ij.attribute(far, self.els, WIN)
        return hit > 0

    def counts(self, pos):
        """(junction reads, spanning reads) at H37Rv position pos (1-based)."""
        junc, span = set(), 0
        for r in self.bam.fetch(self.ctg, max(0, pos - 400), pos + 400):
            if not self._ok(r) or not r.cigartuples:
                continue
            ct = r.cigartuples
            lc = ct[0][1] if ct[0][0] in (4, 5) else 0
            rc = ct[-1][1] if ct[-1][0] in (4, 5) else 0
            first, last = r.reference_start + 1, r.reference_end
            if rc >= MIN_CLIP and abs(last - pos) <= WIN \
                    and self._clip_is_element(r, False, rc):
                junc.add(r.query_name)
            if lc >= MIN_CLIP and abs(first - (pos + 1)) <= WIN \
                    and self._clip_is_element(r, True, lc):
                junc.add(r.query_name)
            if r.is_supplementary:
                continue
            if any(bs + 1 <= pos - WIN and be >= pos + WIN for bs, be in r.get_blocks()):
                span += 1
        return len(junc), span

    def check(self, pos):
        """(state, ev, junction, spanning) for an IS6110 key at pos."""
        j, s = self.counts(pos)
        st, ev = verdict(j, s, self.floor)
        if st == "REF" and self.mask is not None and self.mask.overlaps(pos, pos + 1):
            st, ev = "NOCALL", "absent_masked"
        elif ev == "absent_low_depth" and j <= MAX_JUNC_REF and self.region is not None:
            st, ev = self.region(pos)
        return st, ev, j, s

    def use_region_check(self, panel_diffs, reference, samtools, graph_vcf=None):
        """No reads at the site at all is most often the site's whole region
        missing in this sample (new1k_p200: 0 junction and 0 spanning reads
        at about 1,100 of the 1,185 low-depth cells of 40 isolates, at sites
        like 3,119,970, 3,378,897 and 2,632,359 in deletion regions). That
        is what ABSENT means, and ABS-1 has the test for it: R's deletion
        over the site, from the panel table, and reads JOINING its two ends
        in this same alignment (inherited_check.check_absent). deleted_here
        -> ABSENT junction_joined, as for a small variant; masked -> NOCALL
        absent_masked; no_r_deletion / r_no_genotype / mixed -> NOCALL
        absent_<verdict>; unresolved stays NOCALL absent_low_depth. The
        H37Rv-frame depth that check needs is read only when a site asks."""
        sys.path.insert(0, os.path.join(_HERE, "..", "..", "bin"))
        import inherited_check as ic
        state = {}

        def region(pos):
            if "h" not in state:
                state["d"] = ic.PanelDiffs(panel_diffs, reference, graph_vcf)
                state["h"] = ic.H37rvFrame(self.bam.filename.decode()
                                           if isinstance(self.bam.filename, bytes)
                                           else self.bam.filename, samtools)
            v, _d, _ev = ic.check_absent(state["d"], self.mask or ic.Mask([]),
                                         state["h"], pos, pos + 1)
            if v == "deleted_here":
                return "ABSENT", "junction_joined"
            if v == "masked":
                return "NOCALL", "absent_masked"
            if v in ("no_r_deletion", "r_no_genotype", "mixed"):
                return "NOCALL", f"absent_{v}"
            return "NOCALL", "absent_low_depth"
        self.region = region

    region = None
