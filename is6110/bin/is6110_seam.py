#!/usr/bin/env python3
"""One rule for "does this IS6110 site sit at the reference's own copy".

A matched reference carries some of an isolate's insertions. The build removed
each such copy and recorded it in <ref>.crossmap.tsv, so a site at one of those
seams is the reference's own copy (ref_shared, a REF allele in that frame), and
a site anywhere else is an insertion the reference lacks (ref_lacking, ALT).

THE RULE USED TO DIFFER BETWEEN SCRIPTS. is6110_write_vcf.py and
is6110_place_by_flank.py required the site to be EXACTLY on a seam base, and
is6110_promote_sites.py allowed 25 bp (review finding 4.4). A site 1-3 bp off
the seam was therefore an ALT in the VCF, and flank placement did not step past
the reference's copy, so one window sat inside the element.

THE SLOP IS 3 bp, FROM THE DATA. Distances from every reconciled gwas1000 site
(2026-09-30 run) to its reference's nearest seam:

    0 bp    8,138          4-10 bp   18
    1 bp      271         11-25 bp   32
    2 bp       60         26-100 bp 106
    3 bp        8

The peak at 0-3 bp is the target-site duplication IS6110 makes (3-4 bp): a
stack's peak can sit on either copy of the duplicated bases. From 4 bp out the
count is flat at about 2 per bp, the background of independent insertions near
a reference copy, so 25 bp of slop would call about 50 of those ref_shared.

SEAM COORDINATES. A junction sits AT a seam, so a site on the left of the
removed span reports orig_start - 1 and one on the right reports orig_end + 1.
"""
import bisect, csv, os

SEAM_SLOP = 3


class Seams:
    """The removed spans of one reference, in its original coordinates."""

    def __init__(self, crossmap_path):
        self.spans = []
        if crossmap_path and os.path.exists(crossmap_path):
            for r in csv.DictReader(open(crossmap_path, newline=""), delimiter="\t"):
                self.spans.append((int(r["orig_start"]), int(r["orig_end"])))
        self.spans.sort()
        # (seam position, side, s0, e0). side L: the site is left of the
        # element, which runs to its right; side R: the element is to its left.
        self.points = sorted([(s0 - 1, "L", s0, e0) for s0, e0 in self.spans]
                             + [(e0 + 1, "R", s0, e0) for s0, e0 in self.spans])
        self._pos = [p[0] for p in self.points]

    def known(self):
        return bool(self.spans)

    def nearest(self, pos, slop=SEAM_SLOP):
        """(side, s0, e0, offset) for the closest seam within slop, else None.

        offset = pos - seam. Ties go to the smaller seam position, so the
        answer never depends on input order.
        """
        i = bisect.bisect_left(self._pos, pos - slop)
        best = None
        while i < len(self.points) and self._pos[i] <= pos + slop:
            sp, side, s0, e0 = self.points[i]
            d = pos - sp
            if best is None or abs(d) < abs(best[3]):
                best = (side, s0, e0, d)
            i += 1
        return best

    def shared(self, pos, slop=SEAM_SLOP):
        return self.nearest(pos, slop) is not None


def load(crossmap_dir, ref, cache=None):
    """Seams for one reference, memoised in `cache` if given."""
    if cache is not None and ref in cache:
        return cache[ref]
    s = Seams(os.path.join(crossmap_dir, f"{ref}.crossmap.tsv"))
    if cache is not None:
        cache[ref] = s
    return s
