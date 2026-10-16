#!/usr/bin/env python3
"""LEVEL 1: does this sample carry the accessory locus? One character, every sample stated.

THE TWO LEVELS. Michael's framing, and it settles a question the pipeline had
been getting wrong. An accessory locus carries two different kinds of
information, with different state spaces:

  level 1   presence or absence of the insert itself. One biallelic character
            per locus. Every sample has a state, read from its own reads.
  level 2   variation inside the insert, among the samples that carry it.
            Conditional on level 1, and for a non-carrier a variant inside the
            insert is neither reference nor unknown -- it is INAPPLICABLE. You
            cannot have a SNP in sequence you do not have.

Conflating them is what produced 132 reference calls in 178 million accessory
cells. This script does level 1 only, and does it for every sample, because
that is the character the convergence test can actually use.

TWO READ ROUTES, CHOSEN PER LOCUS (fix SV-EVAL-2, 2026-10-16).

  pool route        the reads H37Rv could not place cleanly -- unmapped, the
                    clipped tails of clipped reads, and reads aligned with
                    NM >= 5 -- realigned to H37Rv plus the catalogue together.
                    Breadth over the locus positions that are UNIQUE AGAINST
                    H37Rv (below) is presence, regardless of which reference
                    the sample was matched to. Validated on TbD1 in
                    insgt/TBD1_TEST.md.
  junction route    for a locus with fewer than MIN_UNIQUE such positions --
                    its sequence is already in H37Rv, 411 of the P200 build's
                    564 loci are IS6110 copies and other repeats -- the pool
                    can say nothing, because the locus's reads place on
                    H37Rv's own copy. What is still read is the anchor: in
                    the sample's alignment to H37Rv, reads clipped within
                    JWIN bp of the catalogue `pos` whose clipped bases are
                    the locus's sequence (an insertion there), against reads
                    running across `pos` unbroken (none).

WHY UNIQUE POSITIONS. The pool is a biased subset: only reads that failed on
H37Rv. Reads from a part of the locus H37Rv also carries align to H37Rv and
never reach it, so a present MOSAIC locus covered "0.2-0.8 of its length" was
UNCERTAIN, written NOCALL with no reason: 42,830 of the 99,000 measured cells
of new1k_p200, while the read check found all 20 such cells it looked at
present (breadth 0.88-1.00 over the H37Rv-unique positions;
analysis/sv_eval_1k/README.md, problem 2). The unique positions are those
no 150 bp tile of the locus covers that aligns to H37Rv as a read would --
clipped by less than --min-clip and with NM <= PLACE_MAX_NM -- since such a
read lands on H37Rv and stays out of the pool. They depend only on the
catalogue and H37Rv, so they are computed once per catalogue and cached
(--index-dir, beside the output, never in the build).

WHY READ_ROUTE_BLIND NO LONGER DECIDES. It was the rule "0.9 of the sequence
in H37Rv at 95% identity: unmeasurable", applied by locus and every cell
written UNMEASURABLE -> the merge skipped the record. 465 of 564 loci, so the
VCF carried 99. The junction route measures those loci from the reads at
their anchor instead. The function stays, for presence tables written before
(no `route` column), which the merge still reads that way.

COVERAGE, NOT NORMALISED DEPTH, is the pool call. A fully spanned locus
still shows a fraction of genome depth in the pool; a read straddling the
boundary contributes its clipped portion, not its length. Measured on TbD1:
coverage 0.99 and 1.00 in carriers at normalised depth 0.30 and 0.48.

UNFILTERED DEPTH, for the same reason the insertion test needed it: 620 of 725
reads on the TbD1 contig carry MAPQ 0, because accessory sequences share
material and the aligner cannot choose between relatives (the catalogue holds
several near-identical alleles at one site). Both thresholds are recorded.

EVERY CELL CARRIES ITS EVIDENCE CODE (column `evidence_code`), which the
merged VCF writes as FORMAT/EV on the NOCALL cells:
  unique_breadth           pool route decided (PRESENT or ABSENT)
  partial_unique_coverage  pool route, breadth between the two thresholds
  junction_reads           junction route, insertion junctions at the anchor
  spanning_reads           junction route, reads cross the anchor unbroken
  junction_unclear         junction route, both or too few of either
  junction_other_sequence  junction route, reads break at the anchor into
                           sequence that is not the locus's
  deletion_junction        junction route, reads join the anchor to H37Rv
                           sequence further along: a deletion edge
  junction_low_depth       junction route, under JMIN_REF reads at the anchor
  no_anchor                junction route, the locus has no H37Rv position
"""
import argparse, collections, csv, hashlib, os, re, shutil, subprocess, sys, tempfile

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")
BLIND_H37RV_COV = 0.9

# The code version, recorded in each table's sidecar by locus_presence_array.sh
# so that tables made by older code are made again rather than kept.
PRESENCE_VERSION = "2"

# unique positions
TILE, TILE_STEP = 150, 10
PLACE_MAX_NM = 4            # a read with NM <= 4 on H37Rv stays out of the pool
MIN_UNIQUE = 100            # fewer unique positions: the junction route
# junction route (analysis/sv_eval_1k/readcheck.py check_ins / check_is)
JWIN = 30                   # clip point within this of the anchor
JCLIP = 15                  # clip length that counts, and k for the k-mer test
JMIN_JUNC = 3               # junction reads for PRESENT ...
JMIN_FRAC = 0.3             # ... at this share of junction + spanning reads
JMAX_JUNC_ABSENT = 1        # at most this many for ABSENT ...
JMIN_REF = 8                # ... with this many spanning reads
JMAPQ = 20


def read_route_blind(locus):
    """True for a catalogue locus whose sequence H37Rv already carries: its
    reads align to H37Rv, so the unmapped and clipped pool cannot measure it.
    No longer decides a state here (the junction route measures such loci);
    bin/merge_cohort_vcf.py applies it to presence tables written before
    PRESENCE_VERSION 2, which have no `route` column.

    "Carries" means at 95% identity or more over 0.9 of the sequence
    (h37rv_cov95, written by merge_catalogues.py; review 2, R2-IS-1). The
    any-identity h37rv_cov and novelty=copy_number made loci 85% identical to
    H37Rv unmeasurable, though their reads are too divergent to place on H37Rv:
    on the CX333 catalogue 70 of the 674 loci so judged are measurable, among
    them ACC_2867346, ACC_2165937 and ACC_0334653, which vary in gwas1000. A
    catalogue without the column (made before) falls back to the old rule,
    which can only call more loci blind."""
    v = locus.get("h37rv_cov95")
    if v not in (None, ""):
        try:
            return float(v) >= BLIND_H37RV_COV
        except ValueError:
            return True
    if (locus.get("novelty") or "") == "copy_number":
        return True
    try:
        return float(locus.get("h37rv_cov") or 0) >= BLIND_H37RV_COV
    except ValueError:
        return False


def read_fasta(path):
    """{name: upper-case sequence}; FATAL on an empty or duplicated file."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                if name in seqs:
                    sys.exit(f"FATAL: {path} holds {name} twice")
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    if not seqs:
        sys.exit(f"FATAL: no sequences in {path}")
    return seqs


def _rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _run(cmd, what, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"FATAL: {what} exited {r.returncode}:\n{r.stderr[-2000:]}")
    return r


# ---- the per-catalogue index: H37Rv + loci, and each locus's unique positions
def catalogue_index(fasta, h37rv, bwa, index_dir, min_clip):
    """(combined FASTA path, {locus: sorted list of 0-based unique positions}).

    Built once per (catalogue, H37Rv, rule) under index_dir/<key>/ and reused
    by every sample; array tasks that race each build in a private directory
    and the first rename wins. Without index_dir it is built in a temporary
    directory the caller removes."""
    h = hashlib.sha256()
    for p in (fasta, h37rv):
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    h.update(f"v{PRESENCE_VERSION}:{TILE}:{TILE_STEP}:{PLACE_MAX_NM}:{min_clip}".encode())
    key = h.hexdigest()[:16]
    final = os.path.join(index_dir, key) if index_dir else None
    if final and os.path.exists(os.path.join(final, "unique.tsv")):
        return os.path.join(final, "comb.fa"), _load_unique(os.path.join(final, "unique.tsv"))
    base = index_dir or tempfile.gettempdir()
    os.makedirs(base, exist_ok=True)
    work = tempfile.mkdtemp(prefix=f".{key}.", dir=base)
    loci = read_fasta(fasta)
    h37 = read_fasta(h37rv)
    if len(h37) != 1:
        sys.exit(f"FATAL: {h37rv} holds {len(h37)} sequences; expected H37Rv alone")
    clash = set(loci) & set(h37)
    if clash:
        sys.exit(f"FATAL: catalogue and H37Rv share sequence names {sorted(clash)}")
    comb = os.path.join(work, "comb.fa")
    with open(comb, "w") as fh:
        for n, s in list(h37.items()) + list(loci.items()):
            fh.write(f">{n}\n{s}\n")
    _run([bwa, "index", comb], "bwa index")
    h37name = next(iter(h37))
    # tiles of every locus, every TILE_STEP bp plus one flush with its end
    tq = os.path.join(work, "tiles.fq")
    with open(tq, "w") as fh:
        for n, s in loci.items():
            L = len(s)
            t = min(TILE, L)
            starts = sorted(set(range(0, L - t + 1, TILE_STEP)) | {L - t})
            for i in starts:
                q = s[i:i + t]
                fh.write(f"@{n}|{i}\n{q}\n+\n{'I' * len(q)}\n")
    # -a: every alignment, so a tile whose primary hit is a sibling locus still
    # shows its H37Rv hit
    r = _run([bwa, "mem", "-a", "-t", "2", comb, tq], "bwa mem (tiles)")
    shared = collections.defaultdict(set)
    for line in r.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11 or f[2] != h37name or int(f[1]) & 4:
            continue
        ops = CIGAR.findall(f[5])
        lc = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
        rcl = int(ops[-1][0]) if ops and ops[-1][1] in "SH" else 0
        nm = next((int(x[5:]) for x in f[11:] if x.startswith("NM:i:")), 99)
        # a read from here would align to H37Rv unclipped (by less than the
        # pool's clip floor) and near-identical: it never reaches the pool
        if lc < min_clip and rcl < min_clip and nm <= PLACE_MAX_NM:
            n, i = f[0].rsplit("|", 1)
            i = int(i)
            shared[n].update(range(i, i + min(TILE, len(loci[n]))))
    with open(os.path.join(work, "unique.tsv"), "w") as fh:
        fh.write("locus\tlength\tn_unique\tunique\n")
        for n, s in loci.items():
            u = [p for p in range(len(s)) if p not in shared[n]]
            fh.write(f"{n}\t{len(s)}\t{len(u)}\t{_runs(u)}\n")
    os.remove(tq)
    if final:
        try:
            os.rename(work, final)
        except OSError:                       # another task got there first
            shutil.rmtree(work, ignore_errors=True)
        return os.path.join(final, "comb.fa"), _load_unique(os.path.join(final, "unique.tsv"))
    return comb, _load_unique(os.path.join(work, "unique.tsv"))


def _runs(ps):
    """0-based positions as 'a-b,c-d' (closed) runs."""
    out, st, pr = [], None, None
    for p in ps:
        if st is None:
            st = pr = p
        elif p == pr + 1:
            pr = p
        else:
            out.append(f"{st}-{pr}"); st = pr = p
    if st is not None:
        out.append(f"{st}-{pr}")
    return ",".join(out)


def _load_unique(path):
    u = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            n, _, _, runs = line.rstrip("\n").split("\t")
            ps = []
            for r in filter(None, runs.split(",")):
                a, b = r.split("-")
                ps.extend(range(int(a), int(b) + 1))
            u[n] = ps
    return u


# ---- the junction route -------------------------------------------------------
def is_tandem(seq, h37, pos):
    """The inserted sequence is >= 90% identical to the H37Rv bases beside
    the anchor (a duplication): its far junction is then at pos +/- L."""
    L = len(seq)
    for s in (h37[pos:pos + L], h37[max(0, pos - L):pos]):
        if len(s) == L and sum(a == b for a, b in zip(seq, s)) >= 0.9 * L:
            return True
    return False


DEL_REACH = 20000           # a supplementary this close implies a deletion


def _depth20(bam, ctg, p):
    """Primary reads at MAPQ >= JMAPQ over 1-based H37Rv base p."""
    n = 0
    for r in bam.fetch(ctg, p - 1, p):
        if not (r.is_secondary or r.is_supplementary or r.is_duplicate
                or r.mapping_quality < JMAPQ):
            n += 1
    return n


def joins_a_deletion(bam, r, ctg, lc, rcl, s1, e1, local, cache):
    """Does the read join the anchor to H37Rv sequence further along across
    bases the sample LACKS -- a deletion edge, not an insertion?

    Its supplementary alignment (MAPQ >= JMAPQ, same contig; the SA strand is
    not tested, bwa reporting it in the supplementary record's own
    orientation) lies beyond the clipped end with reference skipped -- left
    clip: the SA ends before the read starts; right clip: it starts after the
    read ends -- within DEL_REACH, AND the skipped bases have under a quarter
    of the anchor's local MAPQ >= 20 depth (five points across them). The
    skipped bases being covered means they are still in the sample, so the
    clipped part is a second copy of the far sequence: an insertion of a
    duplicated segment, which still counts.

    Why it is needed: at the ends of H37Rv's own IS6110 copies the graph
    writes the copy's loss in some panel genomes as an "insertion" of the
    flank beyond it (ACC_3552712's sequence is H37Rv 3,554,071 on, just past
    the copy 3,552,713-3,554,067), so the reads of a sample without the copy
    match the locus sequence at the anchor. On ctrl25_p200 that made
    ACC_3796412, ACC_3552712, ACC_1541951 and ACC_2365413 PRESENT in 19-21 of
    25 controls whose assemblies have no insertion there."""
    if not r.has_tag("SA"):
        return False
    for e in filter(None, r.get_tag("SA").split(";")):
        f = e.split(",")
        if len(f) < 5 or f[0] != ctg or int(f[4]) < JMAPQ:
            continue
        sp = int(f[1])
        se = sp + sum(int(n) for n, op in CIGAR.findall(f[3]) if op in REF_CONSUMING) - 1
        span = None
        if lc >= JCLIP and s1 - DEL_REACH <= se < s1 - 1:
            span = (se + 1, s1 - 1)
        elif rcl >= JCLIP and e1 + 1 < sp <= e1 + DEL_REACH:
            span = (e1 + 1, sp - 1)
        if span is None:
            continue
        # spans of one junction differ by a few bases between reads
        key = (span[0] // 50, span[1] // 50)
        if key not in cache:
            a_, b_ = span
            pts = sorted({a_ + (b_ - a_) * k // 6 for k in range(1, 6)})
            d = sum(_depth20(bam, ctg, x) for x in pts) / len(pts)
            cache[key] = d < 0.25 * max(local, 1)
        if cache[key]:
            return True
    return False


def junction_call(bam, ctg, pos, seq, h37):
    """(state, evidence code, junction reads, spanning reads) for an insertion
    of `seq` after H37Rv base `pos` (1-based), read in an H37Rv alignment.

    junction  a read (primary or supplementary, MAPQ >= JMAPQ) soft-clipped by
              >= JCLIP bp with the clip point within JWIN bp of a junction,
              whose JCLIP clipped bases next to the clip point occur in the
              locus sequence (either strand); or a CIGAR insertion of
              >= min(L/2, 30) bp within JWIN bp of it. Junctions are the
              anchor, and for a tandem duplication also pos +/- L.
    spanning  a primary read crossing pos-20..pos+20 with no clip of >= JCLIP
              within 2 JWIN of pos and no insertion or deletion >= JCLIP."""
    if not seq:
        return "UNCERTAIN", "no_sequence", 0, 0
    km = set()
    for s in (seq, _rc(seq)):
        km.update(s[i:i + JCLIP] for i in range(len(s) - JCLIP + 1))
    # the anchor's local depth, against which a skipped span is "lacking"
    local = (_depth20(bam, ctg, max(1, pos - 200)) + _depth20(bam, ctg, pos + 200)) / 2
    dcache = {}
    L = len(seq)
    pts = [pos] + ([pos + L, pos - L] if is_tandem(seq, h37, pos) else [])
    lo, hi = max(0, min(pts) - 400), max(pts) + 400
    junc, other, dels, sa_named, ref = set(), set(), set(), set(), 0
    # read into a list first: joins_a_deletion fetches from the same handle,
    # and a nested fetch on one pysam handle ends the outer iteration silently
    for r in list(bam.fetch(ctg, lo, hi)):
        if (r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail
                or r.mapping_quality < JMAPQ or not r.cigartuples):
            continue
        ct = r.cigartuples
        q = r.query_sequence or ""
        lc = ct[0][1] if ct[0][0] == 4 else 0
        rcl = ct[-1][1] if ct[-1][0] == 4 else 0
        s1, e1 = r.reference_start + 1, r.reference_end    # 1-based closed
        hit = False
        for p in pts:
            if rcl >= JCLIP and abs(e1 - p) <= JWIN \
                    and q[len(q) - rcl:len(q) - rcl + JCLIP] in km:
                hit = True
            if lc >= JCLIP and abs(s1 - (p + 1)) <= JWIN and q[lc - JCLIP:lc] in km:
                hit = True
        near_pt = any((rcl >= JCLIP and abs(e1 - p) <= JWIN)
                      or (lc >= JCLIP and abs(s1 - (p + 1)) <= JWIN) for p in pts)
        if near_pt and joins_a_deletion(bam, r, ctg, lc, rcl, s1, e1, local, dcache):
            # the read joins the two ends of a DELETION (joins_a_deletion),
            # whether or not its clip happens to match the locus
            hit = False
            dels.add(r.query_name)
        big = []
        rp = r.reference_start
        for op, n in ct:
            if op in (1, 2) and n >= JCLIP:
                big.append((op, rp, n))
            if op in (0, 2, 3, 7, 8):
                rp += n
        if any(op == 1 and n >= min(L / 2, 30) and any(abs(bp - p) <= JWIN for p in pts)
               for op, bp, n in big):
            hit = True
        if hit:
            junc.add(r.query_name)
            if r.has_tag("SA"):
                sa_named.add(r.query_name)
        clipped_near = ((lc >= JCLIP and abs(s1 - pos) < 2 * JWIN)
                        or (rcl >= JCLIP and abs(e1 - pos) < 2 * JWIN))
        if clipped_near and not hit:
            other.add(r.query_name)          # some other sequence joins here
        if r.is_supplementary or hit:
            continue
        if s1 <= pos - 20 and e1 >= pos + 20 and not clipped_near and not big:
            ref += 1
    # A DELETION EDGE AT THE ANCHOR. When reads placed by their supplementary
    # alignment join the anchor to H37Rv sequence further along, the shorter
    # clips there (too short for bwa to place) are the same junction, not the
    # locus: the region is deleted, and the locus's presence is not read
    junc -= dels
    if len(dels) >= JMIN_JUNC - 1 and len(dels) >= sum(
            1 for n in junc if n in sa_named):
        return "UNCERTAIN", "deletion_junction", len(junc), ref
    nj = len(junc)
    if nj >= JMIN_JUNC and nj / (nj + ref) >= JMIN_FRAC:
        return "PRESENT", "junction_reads", nj, ref
    no = len(other)
    if nj <= JMAX_JUNC_ABSENT and ref >= JMIN_REF and no / (no + ref) < JMIN_FRAC:
        return "ABSENT", "spanning_reads", nj, ref
    if no >= JMIN_JUNC and nj < JMIN_JUNC:
        # reads break at the anchor but not into this locus: another event
        # (or another allele of this one) is there; not absence
        return "UNCERTAIN", "junction_other_sequence", nj, ref
    if nj + ref < JMIN_REF:
        return "UNCERTAIN", "junction_low_depth", nj, ref
    return "UNCERTAIN", "junction_unclear", nj, ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", default="",
                    help="the sample's matched panel accession, recorded "
                         "beside the read call (reference_carries)")
    # no defaults: accessory/assets/ held CX333's hand-placed catalogue;
    # locus_presence_one.sh passes the build's (P0 step catalogue)
    ap.add_argument("--catalogue", required=True,
                    help="<build>/assets/accessory_catalogue.tsv")
    ap.add_argument("--fasta", required=True,
                    help="<build>/assets/accessory_catalogue.fasta")
    ap.add_argument("--cram", required=True,
                    help="the sample's reads aligned to H37Rv: the pool's source")
    ap.add_argument("--h37rv-bam", default="",
                    help="the sample's P1 alignment to H37Rv, read by the "
                         "junction route (preprocessed reads, as every other "
                         "arm reads); without it the junction route reads --cram")
    ap.add_argument("--h37rv", required=True)
    ap.add_argument("--index-dir", default="",
                    help="cache for the per-catalogue index and unique "
                         "positions, shared by a cohort's tasks; empty: "
                         "built per run in a temporary directory")
    ap.add_argument("--contig", default="NC_000962.3")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-clip", type=int, default=25)
    ap.add_argument("--present-cov", type=float, default=0.80,
                    help="share of the locus's H37Rv-unique positions covered, "
                         "unfiltered, to call it present")
    ap.add_argument("--absent-cov", type=float, default=0.20,
                    help="below this it is called absent; between the two is "
                         "uncertain and says so")
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--bwa", default=os.environ.get("MTB_BWA", "bwa"),
                    help="bwa binary; defaults to MTB_BWA from config/project_env.sh")
    ap.add_argument("--version", action="version", version=PRESENCE_VERSION)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    for p, what in ((a.catalogue, "--catalogue"), (a.fasta, "--fasta"),
                    (a.cram, "--cram"), (a.h37rv, "--h37rv")):
        if not os.path.exists(p):
            sys.exit(f"FATAL: {what} {p} does not exist")
    if a.h37rv_bam and not os.path.exists(a.h37rv_bam):
        sys.exit(f"FATAL: --h37rv-bam {a.h37rv_bam} does not exist")
    cat = list(csv.DictReader(open(a.catalogue, newline=""), delimiter="\t"))
    if not cat or "locus_id" not in cat[0]:
        sys.exit(f"FATAL: {a.catalogue} has no locus_id rows")
    seqs = read_fasta(a.fasta)
    miss = [r["locus_id"] for r in cat if r["locus_id"] not in seqs]
    if miss:
        sys.exit(f"FATAL: {len(miss)} catalogue loci have no sequence in "
                 f"{a.fasta} ({', '.join(miss[:3])}); the two are not one catalogue")
    print(f"  {len(cat):,} accessory loci")

    # EVERY SUBPROCESS IS CHECKED. None was: a CRAM that failed to decode, or a
    # bwa or samtools failure, gave an empty query pool and zero coverage, and
    # every one of the 802 loci was then written ABSENT -- a measurement of
    # nothing, stated as absence.
    def check(rc, what, err=""):
        if rc != 0:
            sys.exit(f"FATAL: {a.sample}: {what} exited {rc}"
                     + (f":\n{err[-2000:]}" if err else ""))

    tmp = tempfile.mkdtemp()
    comb, uniq = catalogue_index(a.fasta, a.h37rv, a.bwa,
                                 a.index_dir or os.path.join(tmp, "index"), a.min_clip)
    h37 = next(iter(read_fasta(a.h37rv).values()))

    # ---- the pool ---------------------------------------------------------
    # One pass over primary records: unmapped reads whole; the clipped tails
    # (>= --min-clip) of mapped reads, whatever their MAPQ; and reads aligned
    # with NM > PLACE_MAX_NM whole, since a locus position "unique against
    # H37Rv" is one whose reads align there with more differences than that.
    # The clipped pool used to be read at MAPQ >= 20 only, and NM was not
    # read: reads of a divergent homolog therefore never reached it.
    q = os.path.join(tmp, "q.fa")
    errf = open(os.path.join(tmp, "stderr.txt"), "w+")   # a file, not a pipe:
    n_q = 0                                             # no deadlock on volume

    def err_text():
        errf.flush(); errf.seek(0); t = errf.read(); errf.seek(0); errf.truncate()
        return t
    with open(q, "w") as fh:
        p = subprocess.Popen([a.samtools, "view", "--reference", a.h37rv,
                              "-F", "0x900", a.cram], stdout=subprocess.PIPE,
                             stderr=errf, text=True)
        for i, line in enumerate(p.stdout):
            f = line.split("\t", 11)
            if len(f) < 10:
                continue
            if int(f[1]) & 4:
                fh.write(f">u{i}\n{f[9]}\n"); n_q += 1
                continue
            nm = re.search(r"\tNM:i:(\d+)", line)
            if nm and int(nm.group(1)) > PLACE_MAX_NM:
                fh.write(f">n{i}\n{f[9]}\n"); n_q += 1
                continue
            if "S" not in f[5]:
                continue
            parts = CIGAR.findall(f[5])
            if not parts:
                continue
            if parts[0][1] == "S" and int(parts[0][0]) >= a.min_clip:
                fh.write(f">l{i}\n{f[9][:int(parts[0][0])]}\n"); n_q += 1
            if parts[-1][1] == "S" and int(parts[-1][0]) >= a.min_clip:
                fh.write(f">r{i}\n{f[9][-int(parts[-1][0]):]}\n"); n_q += 1
        check(p.wait(), "samtools view (pool)", err_text())
    print(f"  {n_q:,} query reads from the unmapped, clipped and divergent pools")
    if n_q == 0:
        # Every real isolate has some unmapped or clipped reads; none at all
        # means the reads were not read, and every locus would come out ABSENT.
        sys.exit(f"FATAL: {a.sample}: no query reads from {a.cram}; refusing to "
                 f"call every locus ABSENT from an empty pool")

    # realigned to H37Rv AND the loci, so a pool read that is H37Rv's own
    # sequence (a sequencing-error read, a clip of a local rearrangement) goes
    # back to H37Rv instead of onto a locus that resembles it
    bam = os.path.join(tmp, "a.bam")
    m = subprocess.Popen([a.bwa, "mem", "-t", "2", "-k", "19", comb, q],
                         stdout=subprocess.PIPE, stderr=errf)
    r = subprocess.run([a.samtools, "sort", "-o", bam, "-"], stdin=m.stdout,
                       capture_output=True, text=True)
    m.stdout.close()
    check(m.wait(), "bwa mem", err_text())
    check(r.returncode, "samtools sort", r.stderr)
    r = subprocess.run([a.samtools, "index", bam], capture_output=True, text=True)
    check(r.returncode, "samtools index", r.stderr)

    bed = os.path.join(tmp, "loci.bed")
    with open(bed, "w") as fh:
        for n, s in seqs.items():
            fh.write(f"{n}\t0\t{len(s)}\n")
    depth0, depth20 = {}, {}
    for Q, store in ((0, depth0), (a.min_mapq, depth20)):
        r = subprocess.run([a.samtools, "depth", "-a", "-Q", str(Q), "-b", bed, bam],
                           capture_output=True, text=True)
        check(r.returncode, "samtools depth", r.stderr)
        for line in r.stdout.splitlines():
            f = line.split("\t")
            if len(f) >= 3 and int(f[2]) > 0:
                store.setdefault(f[0], set()).add(int(f[1]) - 1)

    def breadth(store, lid, ps):
        c = store.get(lid, ())
        return (sum(1 for x in ps if x in c) / len(ps)) if ps else 0.0

    # ---- the junction route's alignment ----------------------------------
    jb = None
    jdone = 0
    if any(len(uniq.get(r["locus_id"], ())) < MIN_UNIQUE for r in cat):
        import pysam
        src = a.h37rv_bam or a.cram
        jb = pysam.AlignmentFile(src, reference_filename=a.h37rv)
        if a.contig not in jb.references:
            sys.exit(f"FATAL: {src} has no contig {a.contig}")

    # ---- the reference route, recorded beside the read call ---------------
    ref_carries = {}
    for r in cat:
        cs = r.get("panel_carriers") or ""
        ref_carries[r["locus_id"]] = (a.reference in cs.split(",")) if a.reference else None

    rows = []
    tally = collections.Counter()
    for r in cat:
        lid = r["locus_id"]
        L = len(seqs[lid])
        allpos = range(L)
        c0, c20 = breadth(depth0, lid, allpos), breadth(depth20, lid, allpos)
        up = uniq.get(lid, [])
        u0, u20 = breadth(depth0, lid, up), breadth(depth20, lid, up)
        nj = nr = ""
        if len(up) >= MIN_UNIQUE:
            route = "pool"
            if u0 >= a.present_cov:
                state, code = "PRESENT", "unique_breadth"
            elif u0 <= a.absent_cov:
                state, code = "ABSENT", "unique_breadth"
            else:
                state, code = "UNCERTAIN", "partial_unique_coverage"
            why = f"{len(up)} H37Rv-unique positions, covered {u0:.2f}"
        else:
            route = "junction"
            pos = r.get("pos") or ""
            if not pos.isdigit():
                state, code = "UNCERTAIN", "no_anchor"
                why = "no H37Rv anchor for the junction route"
            else:
                state, code, nj, nr = junction_call(jb, a.contig, int(pos),
                                                    seqs[lid], h37)
                jdone += 1
                why = (f"{len(up)} H37Rv-unique positions: anchor {pos} read, "
                       f"{nj} junction / {nr} spanning reads")
        rc = ref_carries.get(lid)
        # The reference route is recorded next to the read route rather than
        # overriding it: where they disagree, the sample differs from its own
        # reference, which is a real finding and not an error to hide.
        agree = "" if (rc is None or state == "UNCERTAIN") else (
            "agree" if (rc and state == "PRESENT") or (not rc and state == "ABSENT")
            else "differs")
        tally[(route, state)] += 1
        rows.append(dict(sample=a.sample, locus=lid, klass=r.get("klass", ""),
                         length=r.get("graph_len") or r.get("rep_len", ""),
                         state=state, evidence=why, evidence_code=code,
                         route=route,
                         cov_unfiltered=f"{c0:.3f}", cov_q20=f"{c20:.3f}",
                         n_unique=len(up), cov_unique=f"{u0:.3f}",
                         cov_unique_q20=f"{u20:.3f}",
                         junction_reads=nj, spanning_reads=nr,
                         reference_carries=("" if rc is None else int(rc)),
                         route_agreement=agree,
                         panel_carrier_frac=r.get("carrier_frac", ""),
                         presence_version=PRESENCE_VERSION))
    # Written to a temporary name and renamed: the array skips any sample whose
    # table exists, so a task killed mid-write must not leave one behind.
    with open(a.out + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    os.replace(a.out + ".tmp", a.out)
    shutil.rmtree(tmp, ignore_errors=True)
    t = sum(tally.values())
    print(f"  {a.sample}: " + "  ".join(f"{k[0]}:{k[1]} {v:,} ({v/t:.1%})"
                                        for k, v in sorted(tally.items())))
    dis = sum(1 for r in rows if r["route_agreement"] == "differs")
    print(f"  the read call differs from the sample's own reference at {dis:,} loci")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
