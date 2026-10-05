#!/usr/bin/env python3
"""P5 step 2: resolve one sample's state at every key in the cohort union.

A matrix cell must distinguish four things, and a plain VCF conflates the last
three into a missing value:

  ALT     the sample carries a non-reference allele here
  REF     the sample was observed here and matches the canonical reference
  ABSENT  this position has no equivalent in the sample's own reference context
  NOCALL  nothing was observed, so nothing can be said

T6 measured 508,039 cells where the REF/ABSENT distinction was being lost, 494 of
them at lineage-tracking positions. Writing ABSENT as REF asserts the sample has
reference sequence it may not have; writing it as NOCALL discards real
information. Both are wrong in ways that propagate into any downstream tree or
association test, so all four states are kept distinct.

How each is decided:

  ALT     the sample has a placed record whose allele differs from the canonical
          REF. Note the composition: a composed-arm record's own REF is its
          MATCHED reference's base, not H37Rv's, so the allele is re-expressed
          against the canonical REF before comparison.
  REF     either the sample has a record that composes to the canonical allele --
          a REVERSION, where the sample differs from its reference precisely by
          matching H37Rv, measured here at 12.5% of composed/called records and
          which would otherwise enter the matrix as a false variant -- or the
          sample has no record and its GVCF shows coverage at the projected
          position.
  ABSENT  the sample's reference R has no base for the key's H37Rv position:
          R's sequence where the position lands is H37Rv's with a gap that
          contains it (deleted_in_ref). dist.to.ref != 0 alone is NOT that.
  NOCALL  nothing observed, or the evidence does not decide.

WHAT dist.to.ref MEANS (audit P4P5-1). `odgi position -r R` walks the graph
from the H37Rv base to the nearest node on R's path; the distance is the bases
walked and the target is where the walk lands in R. Tested on GCF_014900175
(25 R SNPs, 10 R deletions, each with its neighbours): at an R SNP the H37Rv
base sits on a 1 bp node R does not walk, so dist is 1 -- and the target is
R's own base at that site, the homologous position. Inside an R deletion the
walk crosses the deleted node, dist is roughly its length, and every deleted
base lands on the same base beside the gap. So a non-zero distance says the
H37Rv base is on a node R does not walk -- a substitution as often as a
deletion; writing ABSENT for every one stated "deleted" at about 97% of R's
own SNP sites. The FASTAs now decide: matching flanks around the target
(homologous) mean the position exists and is genotyped there; R's sequence at
the target being a clean H37Rv junction around the position means ABSENT;
anything else is NOCALL.

WHAT R ITSELF CARRIES (audit P4P5-2). A gVCF reference block says the sample
matches R, which is REF only where R matches H37Rv at the key. Where R carries
the key's ALT (a core key, which has no inherited record) the same block means
the sample carries the ALT; where R carries a third allele it decides nothing.

GVCF coverage is read from DP, never MIN_DP. Stage 4 found one GVCF block
spanning positions 14 to 15,224 with DP 59 and MIN_DP 3, and using MIN_DP
rejected 464 of 1,131 inherited differences spuriously.
"""
import argparse, bisect, collections, csv, gzip, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtb_norm import normalise, h37rv_key, node_key

import importlib.util as _ilu
_gfs = _ilu.spec_from_file_location(
    "graph_frame", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "graphframe", "bin", "graph_frame.py"))
graph_frame = _ilu.module_from_spec(_gfs); _gfs.loader.exec_module(graph_frame)

COMP = str.maketrans("ACGTN", "TGCAN")


def op(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_gvcf(path):
    """Reference blocks, and the positions where the sample is NOT reference.

    Returns (blocks, alt_spans, alt_at):
      blocks     sorted (start, end, dp) over lines whose GT is 0 -- reference
                 blocks and hom-ref sites. Only these are evidence of REF.
      alt_spans  sorted (start, end) covering every base of the REF allele of a
                 line whose GT is a non-reference allele. A key inside one is
                 not REF, whatever the blocks say: GATK starts the next block
                 after POS, so bases inside a called deletion can sit in a
                 reference block.
      alt_at     pos -> (ref, called allele, dp) for those lines.

    Counting every line as coverage, as this used to, made a sample REF at a
    key for allele G while its own gVCF called allele A at that very base.
    """
    blocks, alt_spans, alt_at = [], [], {}
    for line in op(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        pos = int(f[1])
        end = pos
        for kv in f[7].split(";"):
            if kv.startswith("END="):
                try:
                    end = int(kv[4:])
                except ValueError:
                    pass
        # DP, not MIN_DP: MIN_DP is the block minimum and rejects real coverage
        keys = f[8].split(":")
        vals = f[9].split(":")
        dp = 0
        if "DP" in keys:
            try:
                dp = int(vals[keys.index("DP")])
            except (ValueError, IndexError):
                dp = 0
        gt = vals[0].replace("|", "/").split("/")[0] if vals else "."
        if gt == "0":
            blocks.append((pos, max(end, pos), dp))
        elif gt.isdigit():
            ref = f[3].upper()
            alts = f[4].upper().split(",")
            i = int(gt)
            called = alts[i - 1] if 0 < i <= len(alts) else "<NON_REF>"
            alt_spans.append((pos, pos + max(len(ref), 1) - 1))
            alt_at[pos] = (ref, called, dp)
        # GT "." -- no call -- is neither REF nor ALT evidence
    blocks.sort(); alt_spans.sort()
    return blocks, alt_spans, alt_at


def load_gvcf_blocks(path):
    """Sorted (start, end, dp) over EVERY gVCF line, reference or variant.

    The original reader, kept for p5_sv_genotype.py, which imports it: depth
    across a deletion interval is read coverage whatever the genotype of each
    line, so counting variant lines is right there. It is NOT evidence of REF
    at a single key -- p5_states.py itself uses load_gvcf() for that.
    """
    blocks = []
    for line in op(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        pos = int(f[1])
        end = pos
        for kv in f[7].split(";"):
            if kv.startswith("END="):
                try:
                    end = int(kv[4:])
                except ValueError:
                    pass
        keys, vals = f[8].split(":"), f[9].split(":")
        dp = 0
        if "DP" in keys:
            try:
                dp = int(vals[keys.index("DP")])
            except (ValueError, IndexError):
                dp = 0
        blocks.append((pos, max(end, pos), dp))
    blocks.sort()
    return blocks


def make_cov(blocks, min_dp):
    starts = [b[0] for b in blocks]

    def covered(p):
        i = bisect.bisect_right(starts, p) - 1
        if i < 0:
            return False
        s, e, dp = blocks[i]
        return s <= p <= e and dp >= min_dp
    return covered


def make_in_spans(spans):
    """Membership in a sorted list of possibly overlapping closed spans."""
    starts = [x[0] for x in spans]
    # running max of span ends, so an early long span is not missed
    reach, m = [], 0
    for _, e in spans:
        m = max(m, e); reach.append(m)

    def inside(p):
        # walk back from the last span starting at or before p while some
        # earlier span could still reach p; the running max bounds the walk
        j = bisect.bisect_right(starts, p) - 1
        while j >= 0 and reach[j] >= p:
            if spans[j][1] >= p:
                return True
            j -= 1
        return False
    return inside


def load_fasta1(path):
    """First record of a FASTA, upper case."""
    seq = []
    with op(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def r_read(rseq, t, strand, n, back=0):
    """n bases of R read along H37Rv's strand, starting `back` bases before
    the base homologous to R position t (1-based). "" off either end.

    On strand `-` H37Rv position p+k sits at R position t-k, so the bases come
    from R leftwards and are complemented."""
    if strand == "-":
        lo, hi = t + back - n, t + back
        if lo < 0 or hi > len(rseq):
            return ""
        return rseq[lo:hi].translate(COMP)[::-1]
    lo = t - 1 - back
    if lo < 0 or lo + n > len(rseq):
        return ""
    return rseq[lo:lo + n]


def _mis(x, y):
    return sum(1 for u, v in zip(x, y) if u != v)


def homologous(h37, rseq, p, t, strand, flank=15, max_mis=2):
    """Does R position t hold the base homologous to H37Rv position p?

    Both flanks, `flank` bases each, must agree gaplessly to within `max_mis`
    mismatches; the base itself is not compared, since it is the one that may
    differ. A deletion fails it: one side of the target is the far side of the
    gap. 2 of 15 tolerates a neighbouring SNP or a short MNP."""
    if p - 1 - flank < 0:
        return False
    hw = h37[p - 1 - flank:p + flank]
    rw = r_read(rseq, t, strand, 2 * flank + 1, back=flank)
    if len(hw) != 2 * flank + 1 or len(rw) != len(hw):
        return False
    return (_mis(hw[:flank], rw[:flank]) <= max_mis
            and _mis(hw[flank + 1:], rw[flank + 1:]) <= max_mis)


def _anchor(h37, kmer, p, lo, hi, side):
    """0-based start of `kmer` in H37Rv nearest p and wholly on one side of
    it ("before": ends at or before p-1; "after": starts at or after p+1), or
    -1."""
    if side == "before":
        return h37.rfind(kmer, lo, p - 1)
    return h37.find(kmer, p, hi)


def deleted_in_ref(h37, rseq, p, t, strand, dist=0, k=12, max_del=50000,
                   reach=200, junction=120):
    """Is R's sequence near t H37Rv with a gap that contains position p?

    odgi lands every base of a deletion on the first R node its walk meets,
    beside the gap or some bases past it, and a deletion's junction is often
    untidy -- a base or two of substitution or shift. So:

      1. Place a k-mer of R at or near t in H37Rv, wholly on one side of p.
      2. Walk from it towards p along R and H37Rv together, stepping over a
         lone mismatch when the next 6 bases agree, to where they part. If
         the walk reaches p, p is in R.
      3. Past that point, find the nearest R k-mer (within `junction` bases)
         that sits in H37Rv wholly on the OTHER side of p.
      4. Between the two anchored runs R has a middle Rm and H37Rv a middle
         Hm holding p. It is a deletion of p only if Rm is shorter than Hm and
         is Hm's prefix plus Hm's suffix, to within a few mismatches, with p
         in the part left out.

    The anchor is sought within `reach` of t; the walk may go as far as
    odgi's own distance further, since the target can sit that far from the
    junction (a 3 kb deletion landed 1.1 kb short of it).

    A divergent or rearranged target, or one with no local homolog, fails,
    and the caller states NOCALL rather than ABSENT."""
    W = reach + min(max(dist, 0), 20000)
    ow = r_read(rseq, t, strand, 2 * W + 1, back=W)
    if len(ow) != 2 * W + 1:
        return False
    lo, hi = max(0, p - 1 - max_del), min(len(h37), p + max_del)
    # 1. an anchor near t
    hit = None
    for off in [0] + [s * o for o in range(4, reach - k, 4) for s in (1, -1)]:
        i = W + off
        if i < 0 or i + k > len(ow):
            continue
        km = ow[i:i + k]
        a = _anchor(h37, km, p, lo, hi, "before")
        b = _anchor(h37, km, p, lo, hi, "after")
        if a < 0 and b < 0:
            continue
        if b >= 0 and (a < 0 or b - (p - 1) <= (p - 1) - a):
            hit = (i, b, "after")
        else:
            hit = (i + k - 1, a + k - 1, "before")     # its last base
        break
    if hit is None:
        return False
    i, h, side = hit
    # 2. walk towards p
    step = -1 if side == "after" else 1
    while 0 < i < len(ow) - 1 and 0 < h < len(h37) - 1:
        ni, nh = i + step, h + step
        if ow[ni] != h37[nh]:
            ahead = [(ni + step * j, nh + step * j) for j in range(1, 7)]
            if not all(0 <= x < len(ow) and 0 <= y < len(h37)
                       and ow[x] == h37[y]
                       for x, y in ahead):
                break
        i, h = ni, nh
    if (side == "after" and h <= p - 1) or (side == "before" and h >= p - 1):
        return False                 # the run reaches p: p is in R
    # 3. the other side's anchor, past the parting point
    other = "before" if side == "after" else "after"
    found = None
    for j in range(1, junction):
        if side == "after":
            e = i - j                # R k-mer ending at e
            if e - k + 1 < 0:
                break
            a = _anchor(h37, ow[e - k + 1:e + 1], p, lo, hi, other)
            if a >= 0:
                found = (e, a + k - 1)
                break
        else:
            s = i + j                # R k-mer starting at s
            if s + k > len(ow):
                break
            b = _anchor(h37, ow[s:s + k], p, lo, hi, other)
            if b >= 0:
                found = (s, b)
                break
    if found is None:
        return False
    # 4. the middles, and the deletion model
    if side == "after":
        rm, hm = ow[found[0] + 1:i], h37[found[1] + 1:h]
        poff = (p - 1) - (found[1] + 1)
    else:
        rm, hm = ow[i + 1:found[0]], h37[h + 1:found[1]]
        poff = (p - 1) - (h + 1)
    if not 0 <= poff < len(hm) or len(rm) >= len(hm):
        return False
    n, m = len(rm), len(hm)
    best = None
    for c in range(n + 1):           # Rm = Hm[:c] + Hm[m-(n-c):]
        mis = _mis(rm[:c], hm[:c]) + _mis(rm[c:], hm[m - (n - c):])
        if best is None or mis < best[0]:
            best = (mis, c)
    mis, c = best
    if mis > max(2, n // 5):
        return False
    return c <= poff < m - (n - c)


def r_allele(h37, rseq, p, t, strand, canon, kalt, flank=10, ext=1000):
    """Which allele of an H37Rv-frame key R itself carries at target t:
    "REF", "ALT", "OTHER" (neither, or a shifted anchor) or "UNK".

    A SNP is R's base. Anything longer is compared as haplotypes: the key's
    REF and ALT, each followed by H37Rv's sequence after the REF allele, read
    up to the first base where they differ plus `flank` -- in a tandem repeat
    that first difference is at the end of the array, not at the key -- and R
    read from t must match one exactly to there, within one mismatch after.
    R's left flank must also match H37Rv's, so a target shifted by a repeat
    unit is not read as an indel."""
    if not canon or not kalt:
        return "UNK"
    if len(canon) == 1 and len(kalt) == 1:
        rb = r_read(rseq, t, strand, 1)
        return ("REF" if rb == canon else "ALT" if rb == kalt else
                "OTHER" if rb else "UNK")
    tail = h37[p - 1 + len(canon):p - 1 + len(canon) + ext]
    refh, alth = canon + tail, kalt + tail
    n = min(len(refh), len(alth))
    d = next((i for i in range(n) if refh[i] != alth[i]), None)
    if d is None:
        return "UNK"
    W = min(d + 1 + flank, n)
    rh = r_read(rseq, t, strand, W)
    if len(rh) < W:
        return "UNK"
    lf = h37[max(0, p - 1 - flank):p - 1]
    if lf and _mis(lf, r_read(rseq, t, strand, len(lf), back=len(lf))) > 1:
        return "OTHER"

    def fits(h):
        return rh[:d + 1] == h[:d + 1] and _mis(rh[d + 1:W], h[d + 1:W]) <= 1
    fa, fr = fits(alth), fits(refh)
    if fa and not fr:
        return "ALT"
    if fr and not fa:
        return "REF"
    return "OTHER"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--node-positions",
                    default="accessory/assets/node_positions.tsv",
                    help="node offsets along each panel accession, from the "
                         "same GFA pass as --node-paths. With it a node-frame "
                         "key whose reference carries the sequence is read in "
                         "the sample's own frame and gets REF when covered.")
    ap.add_argument("--node-paths",
                    default="accessory/assets/node_paths.tsv",
                    help="node to panel-path membership from "
                         "accessory/bin/node_path_membership.py. Without it a "
                         "node-frame key with no record is NOCALL; with it, one "
                         "whose sequence is absent from this sample's own "
                         "reference is ABSENT, which is a measurement.")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--keys", required=True)
    ap.add_argument("--placed", required=True)
    ap.add_argument("--gvcf", required=True)
    ap.add_argument("--projected", required=True,
                    help="odgi position output for the keys, H37Rv -> R")
    ap.add_argument("--h37rv", required=True,
                    help="needed to recompute the same normalised keys the key "
                         "set used; pairing on position alone is what broke")
    ap.add_argument("--ref-fasta", default="",
                    help="the sample's reference R in the refs frame, which "
                         "decides whether a projected position is homologous, "
                         "deleted, and which allele R carries. Default: "
                         "<directory of --h37rv>/<--reference>.fasta, which is "
                         "${BUILD}/refs where p5_merge.sh takes --h37rv from")
    ap.add_argument("--large-indel", type=int, default=50,
                    help="a deletion R carries at least this long is NOCALL, "
                         "not ALT, from a reference block alone: the block "
                         "spans R's junction whether or not the sample has "
                         "the sequence, and GATK cannot call it that long")
    ap.add_argument("--min-dp", type=int, default=5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dense", action="store_true",
                    help="write the old one-row-per-key form. The default is "
                         "sparse: only the cells that are not REF, keyed by "
                         "the key's INDEX in --keys rather than by its name. "
                         "At 997 isolates a dense file is 10 MB per sample and "
                         "12 GB for the cohort, of which 82.7% is the word REF "
                         "beside a key string that is already in keys.tsv; the "
                         "sparse form is about 240 KB. Every consumer reads "
                         "either, and the header carries the key count and a "
                         "checksum of the key list so a file written against a "
                         "different key set is refused rather than silently "
                         "misaligned.")
    a = ap.parse_args()

    keys = list(csv.DictReader(open(a.keys, newline=""), delimiter="\t"))
    placed = list(csv.DictReader(open(a.placed, newline=""), delimiter="\t"))
    h37 = []
    with open(a.h37rv) as fh:
        for line in fh:
            if line.startswith(">"):
                if h37:
                    break
                continue
            h37.append(line.strip())
    h37 = "".join(h37)
    if len(h37) < 4_000_000:
        print(f"H37Rv FASTA looks wrong: {len(h37)} bp", file=sys.stderr)
        return 1
    h37 = h37.upper()
    rfa = a.ref_fasta or os.path.join(os.path.dirname(os.path.abspath(a.h37rv)),
                                      f"{a.reference}.fasta")
    if not os.path.exists(rfa):
        print(f"  FATAL: no FASTA for {a.reference} at {rfa}; without it a "
              f"non-zero projection distance cannot be told from a deletion",
              file=sys.stderr)
        return 1
    rseq = load_fasta1(rfa)

    # Pair by SOURCE position, not by file order. odgi re-emits its header line
    # between records, so the raw line count is about twice the query count, and
    # while the data lines happen to come back 1:1 and in order, relying on that
    # is an assumption with no upside: the source position is right there in
    # column 1. A silently misaligned coordinate is the worst failure this stage
    # could produce.
    h_keys = [k for k in keys if k["frame"] == "h37rv"]
    by_src = {}
    n_lines = n_multi = n_inv = 0
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
        # relation of R's refs sequence to H37Rv's here (frame_convert.py);
        # needed to compare a called base in R with an H37Rv base
        strand = f[3].strip() if len(f) > 3 and f[3].strip() in "+-" else "+"
        # WHERE R WALKS THE NODE BACKWARDS (odgi's own flag, column 5, `-`)
        # the target is one base off and R reads reverse to H37Rv there.
        # Measured at every dist-0 key of four references, reading R's refs
        # sequence around the target against H37Rv's around the key:
        #   col4 col5   homolog          GCF_000193185  GCF_001870145
        #   +    +      t, same strand          49,926              -
        #   +    -      t-1, complemented       12,121              -
        #   -    +      t, complemented              -         62,379
        #   -    -      t+1, same strand             -              5
        # Comparing R's base at t on column 4 alone read the neighbouring,
        # uncomplemented base at 13,508 of GCF_000193185's 67,901 keys.
        if len(f) > 4 and f[4].strip() == "-":
            tgt += -1 if strand == "+" else 1
            strand = "-" if strand == "+" else "+"
            n_inv += 1
        n_lines += 1
        if src in by_src:
            n_multi += 1          # the target path visits this position twice
            continue
        by_src[src] = (tgt, dist, strand)
    projmap = {}
    for k in h_keys:
        projmap[k["key"]] = by_src.get(int(k["h37rv_pos"]))
    n_unproj = sum(1 for v in projmap.values() if v is None)
    if n_lines == 0:
        print("  FATAL: no projected positions parsed", file=sys.stderr)
        return 1

    # Re-derive each record's normalised key exactly as p5_keys.py did. The
    # record's own `key` field is position-based and no longer what the key set
    # uses.
    by_key = {}
    for x in placed:
        if x["frame"] == "h37rv":
            pos = int(x["h37rv_pos"])
            canon = h37[pos - 1:pos - 1 + len(x["ref"])].upper() or x["ref"]
            npos, nref, nalt = normalise(h37, pos, canon, x["alt"])
            k = h37rv_key(npos, nref, nalt)
        else:
            npos, nref, nalt = normalise("", int(x["r_pos"] or 0), x["ref"],
                                         x["alt"])
            k = node_key(x["node"], x["node_offset"], nref, nalt)
        by_key.setdefault(k, dict(x, alt=nalt))

    # WHERE THE SAMPLE HAS A PLACED NON-REFERENCE ALLELE, BY H37Rv POSITION.
    # A key with no record of its own at such a position is not REF: the sample
    # carries some other allele there. A SNP/MNP marks its bases and a deletion
    # the bases it removes (not the anchor it keeps). An insertion leaves its
    # anchor base unchanged, so it blocks only other INDEL keys at that anchor,
    # never a SNP key there.
    other_at, ins_anchor = set(), set()
    for k2, x in by_key.items():
        if not k2.startswith("h37rv:"):
            continue
        f2 = k2.split(":")
        npos = int(f2[1])
        nref, _, nalt = f2[2].partition(">")
        if nref == nalt:
            continue                                # a reversion states REF
        if len(nref) > len(nalt):
            other_at.update(range(npos + len(nalt), npos + len(nref)))
        elif len(nalt) > len(nref):
            ins_anchor.add(npos)
        else:
            other_at.update(range(npos, npos + len(nref)))

    blocks, alt_spans, alt_at = load_gvcf(a.gvcf)
    cov = make_cov(blocks, a.min_dp)
    in_alt = make_in_spans(alt_spans)
    n_other = collections.Counter()
    n_proj = collections.Counter()

    def span_covered(t, strand, n):
        """n R bases from t along H37Rv's strand all in reference blocks."""
        ps = range(t, t + n) if strand != "-" else range(t - n + 1, t + 1)
        return all(cov(x) and not in_alt(x) for x in ps)

    def ref_state(p, canon, strand="+", hp=None, kalt=""):
        """REF, ALT or NOCALL at R position p for a key whose H37Rv reference
        allele is `canon`. With `hp`, the key's H37Rv position, the allele R
        itself carries there is checked before a reference block is read as
        REF: the block says the sample matches R, not that it matches
        H37Rv."""
        if in_alt(p):
            ra = alt_at.get(p)
            # A SNP whose called base IS the H37Rv base is a reversion: the
            # sample differs from R and agrees with H37Rv, so REF is right.
            if ra and len(ra[0]) == 1 and len(ra[1]) == 1 and len(canon) == 1:
                b = ra[1].translate(COMP) if strand == "-" else ra[1]
                if b == canon and ra[2] >= a.min_dp:
                    return "REF"
            n_other["gvcf_non_ref"] += 1
            return "NOCALL"
        if not cov(p):
            return "NOCALL"
        if hp is None:
            return "REF"
        ra = r_allele(h37, rseq, hp, p, strand, canon, kalt)
        if ra == "REF":
            return "REF"
        if ra == "ALT":
            # The sample matches R and R carries the key's allele. A long
            # deletion in R is the exception: the block over R's junction is
            # there whether or not the sample has the sequence.
            if len(canon) - len(kalt) >= a.large_indel:
                n_other["r_carries_large_deletion"] += 1
                return "NOCALL"
            # R's whole allele must be covered, and for an indel the base
            # after it too, so the reads span the junction
            if not span_covered(p, strand,
                                len(kalt) + (len(kalt) != len(canon))):
                n_other["r_carries_alt_span_uncovered"] += 1
                return "NOCALL"
            n_other["r_carries_alt"] += 1
            return "ALT"
        n_other["r_carries_other_allele"] += 1
        return "NOCALL"

    frames = graph_frame.Frames()

    rows = []
    # node -> the panel accessions whose path traverses it, so a node-frame
    # key can be told "absent from this sample's reference" from "unknown".
    nodepaths = {}
    if a.node_paths and os.path.exists(a.node_paths):
        for q in csv.DictReader(open(a.node_paths, newline=""), delimiter="\t"):
            nodepaths[q["node"]] = set(q["paths"].split(",")) if q["paths"] else set()
        print(f"  node-path membership for {len(nodepaths):,} nodes from "
              f"{a.node_paths}")
    else:
        print("  NOTE no --node-paths table: every node-frame key without a "
              "record stays NOCALL, which is 97% of them")

    # (node, accession) -> (start, strand, occurrences) along that accession's
    # own sequence, so a node-frame key can be read in the sample's own frame.
    nodepos = {}
    if a.node_positions and os.path.exists(a.node_positions):
        for q in csv.DictReader(open(a.node_positions, newline=""), delimiter="\t"):
            nodepos[(q["node"], q["accession"])] = (
                int(q["start"]), q["strand"], int(q.get("n_occurrences") or 1),
                int(q["length"]) if q.get("length") else None)
        print(f"  node positions for {len(nodepos):,} (node, accession) pairs")

    counts = {"ALT": 0, "REF": 0, "ABSENT": 0, "NOCALL": 0}
    reversions = 0
    for k in keys:
        canon = k["canonical_ref"].upper()
        rec = by_key.get(k["key"])
        if rec is not None:
            alt = rec["alt"].upper()
            if alt == canon:
                state, reversions = "REF", reversions + 1
                allele = canon
            else:
                state, allele = "ALT", alt
        else:
            pr = projmap.get(k["key"])
            if k["frame"] == "node":
                # AN OFF-PATH KEY THIS SAMPLE PRODUCED NO RECORD FOR. The
                # earlier version wrote NOCALL here and said why: "its own
                # reference may or may not traverse that node, and this pass has
                # not asked". accessory/ACCESSORY_CALL_CENSUS.md measures what
                # not asking cost -- 91 REF calls in 6.9 million accessory
                # cells, 40 of 63 loci with zero REF anywhere, and TbD1's 424
                # variant records genotyped in 0.2% of cells despite two test
                # isolates covering the sequence at 0.99 and 1.00.
                #
                # Node-to-path membership is a property of the graph alone, so
                # it can be asked without a projection and without a realignment
                # -- accessory/bin/node_path_membership.py computes it once.
                # Where the sample's own reference does NOT traverse the node,
                # the sequence is not in its reference context and the honest
                # state is ABSENT, which is the same statement an H37Rv-frame
                # key gets when its position does not project. Measured on
                # gwas1000: 90.25% of node-frame cells are that case.
                #
                # Where the reference DOES traverse it, 7.07% of cells, REF
                # needs depth at the corresponding position and this pass still
                # cannot say -- so those stay NOCALL, and closing them is the
                # second half of the fix.
                nd = k.get("node") or ""
                paths = nodepaths.get(nd)
                if paths is None:
                    state, allele = "NOCALL", ""
                elif a.reference and a.reference not in paths:
                    state, allele = "ABSENT", ""
                else:
                    # THE REFERENCE DOES CARRY THIS SEQUENCE, so the sample was
                    # looked at here and the only question is whether the reads
                    # cover it. That is the same question the H37Rv-frame keys
                    # answer, asked in the sample's own frame: the node's offset
                    # along its reference plus the key's offset within the node
                    # gives the position, and the gVCF gives the depth.
                    #
                    # Node positions come from the same GFA pass as membership,
                    # so this still needs no projection and no realignment.
                    hit = nodepos.get((nd, a.reference))
                    off = int(k.get("node_offset") or 0)
                    if hit is None:
                        state, allele = "NOCALL", ""
                    elif hit[2] > 1:
                        # The reference visits this node more than once, so its
                        # copies share reads and depth at one of them settles
                        # nothing. 8.8% of pairs.
                        state, allele = "NOCALL", ""
                    else:
                        # `start` is 1-based along the PANEL path (the GFA
                        # walk). The key's node offset, as P4 writes it, is
                        # already counted along the path's direction, so the
                        # base is start + off on EITHER strand. Measured on the
                        # pilot's references, reading the refs base and
                        # comparing it with the key's reference allele:
                        #
                        #   node walked   old rule           this rule
                        #   +             56.9%  (start+off, no frame change)
                        #                                    100.0%
                        #   -             15.3%  (start-off)  99.8%
                        #
                        # The + gain is the PANEL -> REFS conversion: the gVCF
                        # is refs-frame and 110 of 333 accessions are rotated
                        # or flipped between the two frames. The - gain is the
                        # offset direction; start - off left the node whenever
                        # off > 0.
                        start = hit[0]
                        pp = start + off
                        p = frames.to_refs(a.reference, pp - 1) + 1
                        st = ref_state(p, canon)
                        state, allele = (st, canon if st == "REF" else "")
            elif pr is None:
                # odgi returned nothing for this position: unmeasured, not
                # shown to be deleted
                n_proj["no projection"] += 1
                state, allele = "NOCALL", ""
            elif pr[1] != 0 and not homologous(h37, rseq, int(k["h37rv_pos"]),
                                               pr[0], pr[2]):
                # The H37Rv base is on a node R does not walk and the target
                # is not its homolog. ABSENT only where R's sequence there is
                # H37Rv with a gap containing the position; a rearranged,
                # divergent or repeat-shifted target decides nothing.
                if deleted_in_ref(h37, rseq, int(k["h37rv_pos"]), pr[0], pr[2],
                                  pr[1]):
                    n_proj["off R's path, deleted in R: ABSENT"] += 1
                    state, allele = "ABSENT", ""
                else:
                    n_proj["off R's path, unresolved: NOCALL"] += 1
                    state, allele = "NOCALL", ""
            elif (int(k["h37rv_pos"]) in other_at
                  or (int(k["h37rv_pos"]) in ins_anchor
                      and len(canon) != len(k.get("canonical_alt") or canon))):
                # the sample has a placed record here for a different allele
                n_other["placed_other_allele"] += 1
                state, allele = "NOCALL", ""
            else:
                if pr[1] != 0:
                    # a substitution in R: the position exists, at the target
                    n_proj["off R's path, homologous: genotyped"] += 1
                kalt = (k.get("canonical_alt") or "").upper()
                st = ref_state(pr[0], canon, pr[2], int(k["h37rv_pos"]), kalt)
                state, allele = (st, canon if st == "REF" else
                                 kalt if st == "ALT" else "")
        counts[state] += 1
        rows.append(dict(sample=a.sample, key=k["key"], state=state,
                         allele=allele, frame=k["frame"],
                         region=k["region"], kind=k["kind"]))

    if a.dense:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)
    else:
        # THE KEY SET IS THE INDEX, so it has to be pinned. A sparse file is
        # meaningless against a different key order, and a silently misaligned
        # coordinate is the worst failure this stage can produce -- so the
        # header carries the count and a checksum of the key names, and the
        # reader refuses a mismatch instead of trusting it.
        import hashlib
        h = hashlib.sha1("\n".join(k["key"] for k in keys).encode()).hexdigest()
        with open(a.out, "w", newline="") as fh:
            fh.write(f"#format\tsparse-v1\n#sample\t{a.sample}\n"
                     f"#default\tREF\n#n_keys\t{len(keys)}\n"
                     f"#keys_sha1\t{h}\n")
            fh.write("idx\tstate\tallele\n")
            for i, r in enumerate(rows):
                if r["state"] == "REF":
                    continue
                fh.write(f"{i}\t{r['state']}\t{r['allele']}\n")
    n = len(rows)
    print(f"  {a.sample}: {len(h_keys)} H37Rv keys, {n_lines} projection lines, "
          f"{n_multi} positions the reference visits more than once, "
          f"{n_unproj} with no projection (NOCALL), {n_inv} on steps R "
          f"walks backwards (target moved to the homologous base)")
    print(f"  {a.sample}: {n} keys  " + "  ".join(
        f"{k} {v} ({100*v/n:.1f}%)" for k, v in counts.items()))
    print(f"    reversions resolved to REF rather than ALT: {reversions}")
    print(f"    not REF because the sample carries another allele there: "
          f"{n_other['placed_other_allele']} by placed record, "
          f"{n_other['gvcf_non_ref']} by its gVCF")
    for why, v in sorted(n_proj.items()):
        print(f"    projection: {why}: {v}")
    for why in ("r_carries_alt", "r_carries_alt_span_uncovered",
                "r_carries_large_deletion", "r_carries_other_allele"):
        print(f"    {why.replace('_', ' ')}: {n_other[why]}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
