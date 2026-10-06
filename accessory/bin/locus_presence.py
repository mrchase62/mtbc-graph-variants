#!/usr/bin/env python3
"""LEVEL 1: does this sample carry the accessory locus? One character, every sample stated.

THE TWO LEVELS. Michael's framing, and it settles a question the pipeline had
been getting wrong. An accessory locus carries two different kinds of
information, with different state spaces:

  level 1   presence or absence of the insert itself. One biallelic character
            per locus. Every sample has a state, which is UNMEASURABLE for a
            locus whose sequence H37Rv already carries (see read_route_blind).
  level 2   variation inside the insert, among the samples that carry it.
            Conditional on level 1, and for a non-carrier a variant inside the
            insert is neither reference nor unknown -- it is INAPPLICABLE. You
            cannot have a SNP in sequence you do not have.

Conflating them is what produced 132 reference calls in 178 million accessory
cells. This script does level 1 only, and does it for every sample, because
that is the character the convergence test can actually use.

TWO INSTRUMENTS, BECAUSE ONE IS BLIND EITHER WAY.

  reference route   the sample's own matched reference either carries the locus
                    or not -- a graph property, from the accessory catalogue's
                    panel carrier set. Where it carries it, depth over the
                    locus in the sample's own frame decides. This route cannot
                    see a sample that carries an insert its reference lacks.
  read route        the reads H37Rv could not place -- unmapped, and the
                    clipped tails of clipped reads -- aligned to the locus
                    sequence. Coverage of the sequence is presence regardless
                    of which reference the sample was aligned to. This route is
                    the only one that can see the case above, and it is the one
                    validated on TbD1 in insgt/TBD1_TEST.md.

COVERAGE, NOT NORMALISED DEPTH, is the call. The query pool is a biased subset
-- only the reads that failed on H37Rv -- so a fully spanned locus still shows a
fraction of genome depth; a read straddling the boundary contributes its clipped
portion, not its length. Measured on TbD1: coverage 0.99 and 1.00 in carriers at
normalised depth 0.30 and 0.48. Gating on depth would have called them absent.

UNFILTERED DEPTH, for the same reason the insertion test needed it: 620 of 725
reads on the TbD1 contig carry MAPQ 0, because accessory sequences share
material and the aligner cannot choose between relatives. Both thresholds are
recorded and the pair is the discriminator.
"""
import argparse, collections, csv, os, re, subprocess, sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")
BLIND_H37RV_COV = 0.9


def read_route_blind(locus):
    """True for a catalogue locus whose sequence H37Rv already carries
    (novelty copy_number, or h37rv_cov at or above 0.9): its reads align to
    H37Rv, so the unmapped and clipped pool cannot measure it. Shared with
    bin/merge_cohort_vcf.py, which applies it to tables written before."""
    if (locus.get("novelty") or "") == "copy_number":
        return True
    try:
        return float(locus.get("h37rv_cov") or 0) >= BLIND_H37RV_COV
    except ValueError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", default="",
                    help="the sample's matched panel accession, for the "
                         "reference route")
    # no defaults: accessory/assets/ held CX333's hand-placed catalogue;
    # locus_presence_one.sh passes the build's (P0 step catalogue)
    ap.add_argument("--catalogue", required=True,
                    help="<build>/assets/accessory_catalogue.tsv")
    ap.add_argument("--fasta", required=True,
                    help="<build>/assets/accessory_catalogue.fasta")
    ap.add_argument("--cram", required=True)
    ap.add_argument("--h37rv", required=True)
    ap.add_argument("--contig", default="NC_000962.3")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-clip", type=int, default=25)
    ap.add_argument("--present-cov", type=float, default=0.80,
                    help="share of the locus covered, unfiltered, to call it "
                         "present")
    ap.add_argument("--absent-cov", type=float, default=0.20,
                    help="below this it is called absent; between the two is "
                         "uncertain and says so")
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--bwa", default=os.environ.get("MTB_BWA", "bwa"),
                    help="bwa binary; defaults to MTB_BWA from config/project_env.sh")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cat = list(csv.DictReader(open(a.catalogue, newline=""), delimiter="\t"))
    print(f"  {len(cat):,} accessory loci")

    # EVERY SUBPROCESS IS CHECKED. None was: a CRAM that failed to decode, or a
    # bwa or samtools failure, gave an empty query pool and zero coverage, and
    # every one of the 802 loci was then written ABSENT -- a measurement of
    # nothing, stated as absence.
    def check(rc, what, err=""):
        if rc != 0:
            sys.exit(f"FATAL: {a.sample}: {what} exited {rc}"
                     + (f":\n{err[-2000:]}" if err else ""))

    if not os.path.exists(a.fasta + ".bwt"):
        r = subprocess.run([a.bwa, "index", a.fasta], capture_output=True, text=True)
        check(r.returncode, "bwa index", r.stderr)

    # ---- the read route ---------------------------------------------------
    import tempfile
    tmp = tempfile.mkdtemp()
    q = os.path.join(tmp, "q.fa")
    errf = open(os.path.join(tmp, "stderr.txt"), "w+")   # a file, not a pipe:
    n_q = 0                                             # no deadlock on volume
    def err_text():
        errf.flush(); errf.seek(0); t = errf.read(); errf.seek(0); errf.truncate()
        return t
    with open(q, "w") as fh:
        p = subprocess.Popen([a.samtools, "view", "--reference", a.h37rv,
                              "-f", "4", a.cram], stdout=subprocess.PIPE,
                             stderr=errf, text=True)
        for i, line in enumerate(p.stdout):
            f = line.split("\t", 11)
            if len(f) > 9:
                fh.write(f">u{i}\n{f[9]}\n"); n_q += 1
        check(p.wait(), "samtools view (unmapped pool)", err_text())
        p = subprocess.Popen([a.samtools, "view", "--reference", a.h37rv,
                              "-q", str(a.min_mapq), a.cram],
                             stdout=subprocess.PIPE, stderr=errf,
                             text=True)
        for i, line in enumerate(p.stdout):
            f = line.split("\t", 11)
            if len(f) < 10 or "S" not in f[5]:
                continue
            parts = CIGAR.findall(f[5])
            if not parts:
                continue
            if parts[0][1] == "S" and int(parts[0][0]) >= a.min_clip:
                fh.write(f">l{i}\n{f[9][:int(parts[0][0])]}\n"); n_q += 1
            if parts[-1][1] == "S" and int(parts[-1][0]) >= a.min_clip:
                fh.write(f">r{i}\n{f[9][-int(parts[-1][0]):]}\n"); n_q += 1
        check(p.wait(), "samtools view (clipped pool)", err_text())
    print(f"  {n_q:,} query reads from the unmapped and clipped pools")
    if n_q == 0:
        # Every real isolate has some unmapped or clipped reads; none at all
        # means the reads were not read, and every locus would come out ABSENT.
        sys.exit(f"FATAL: {a.sample}: no query reads from {a.cram}; refusing to "
                 f"call every locus ABSENT from an empty pool")

    bam = os.path.join(tmp, "a.bam")
    m = subprocess.Popen([a.bwa, "mem", "-t", "2", "-k", "19", a.fasta, q],
                         stdout=subprocess.PIPE, stderr=errf)
    r = subprocess.run([a.samtools, "sort", "-o", bam, "-"], stdin=m.stdout,
                       capture_output=True, text=True)
    m.stdout.close()
    check(m.wait(), "bwa mem", err_text())
    check(r.returncode, "samtools sort", r.stderr)
    r = subprocess.run([a.samtools, "index", bam], capture_output=True, text=True)
    check(r.returncode, "samtools index", r.stderr)

    cov0, cov20 = {}, {}
    for Q, store in ((0, cov0), (a.min_mapq, cov20)):
        r = subprocess.run([a.samtools, "depth", "-a", "-Q", str(Q), bam],
                           capture_output=True, text=True)
        check(r.returncode, "samtools depth", r.stderr)
        n, c = collections.Counter(), collections.Counter()
        for line in r.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 3:
                continue
            n[f[0]] += 1
            if int(f[2]) > 0:
                c[f[0]] += 1
        for k in n:
            store[k] = c[k] / n[k]

    # ---- the reference route ----------------------------------------------
    # The catalogue records which panel genomes carry each locus, so whether
    # the sample's own reference carries it is a lookup, not a measurement.
    ref_carries = {}
    for r in cat:
        cs = r.get("panel_carriers") or ""
        ref_carries[r["locus_id"]] = (a.reference in cs.split(",")) if a.reference else None

    rows = []
    tally = collections.Counter()
    for r in cat:
        lid = r["locus_id"]
        c0 = cov0.get(lid, 0.0)
        c20 = cov20.get(lid, 0.0)
        if read_route_blind(r):
            # THE READ ROUTE CANNOT SEE THIS LOCUS (audit P3IS-2). Its sequence
            # is already in H37Rv -- novelty=copy_number, h37rv_cov >= 0.9;
            # 674 of the 802 loci, 481 of them IS6110 copies -- so its reads
            # place on H37Rv and never reach the unmapped and clipped pool.
            # Zero coverage here is not absence, and was written ABSENT in
            # every sample (0 PRESENT of 96,200 cells on the IS6110 loci in
            # scale200). What coverage the pool does give is a divergent
            # homolog of unverified meaning, not this locus. Unmeasured; the
            # coverage is kept in the row.
            state, why = "UNMEASURABLE", (
                f"sequence in H37Rv: the read route is blind to it "
                f"(pool coverage {c0:.2f})")
        elif c0 >= a.present_cov:
            state, why = "PRESENT", "reads cover the locus"
        elif c0 <= a.absent_cov:
            state, why = "ABSENT", "reads do not cover the locus"
        else:
            state, why = "UNCERTAIN", f"partial coverage {c0:.2f}"
        rc = ref_carries.get(lid)
        # The reference route is recorded next to the read route rather than
        # overriding it: where they disagree, the sample differs from its own
        # reference, which is a real finding and not an error to hide.
        agree = "" if (rc is None or state == "UNMEASURABLE") else (
            "agree" if (rc and state == "PRESENT") or (not rc and state == "ABSENT")
            else "differs")
        tally[state] += 1
        rows.append(dict(sample=a.sample, locus=lid, klass=r.get("klass", ""),
                         length=r.get("graph_len") or r.get("rep_len", ""),
                         state=state, evidence=why,
                         cov_unfiltered=f"{c0:.3f}", cov_q20=f"{c20:.3f}",
                         reference_carries=("" if rc is None else int(rc)),
                         route_agreement=agree,
                         panel_carrier_frac=r.get("carrier_frac", "")))
    # Written to a temporary name and renamed: the array skips any sample whose
    # table exists, so a task killed mid-write must not leave one behind.
    with open(a.out + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    os.replace(a.out + ".tmp", a.out)
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    t = sum(tally.values())
    print(f"  {a.sample}: " + "  ".join(f"{k} {v:,} ({v/t:.1%})"
                                        for k, v in sorted(tally.items())))
    dis = sum(1 for r in rows if r["route_agreement"] == "differs")
    print(f"  the two routes differ at {dis:,} loci -- the sample against its "
          f"own reference")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
