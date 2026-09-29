#!/usr/bin/env python3
"""LEVEL 1: does this sample carry the accessory locus? One character, every sample stated.

THE TWO LEVELS. Michael's framing, and it settles a question the pipeline had
been getting wrong. An accessory locus carries two different kinds of
information, with different state spaces:

  level 1   presence or absence of the insert itself. One biallelic character
            per locus. EVERY sample has a state; nothing is missing.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", default="",
                    help="the sample's matched panel accession, for the "
                         "reference route")
    ap.add_argument("--catalogue", default="accessory/assets/accessory_catalogue.tsv")
    ap.add_argument("--fasta", default="accessory/assets/accessory_catalogue.fasta")
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
    ap.add_argument("--bwa", default=os.environ.get(
        "MTB_BWA",
        "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_pangenome_qc/bin/bwa"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cat = list(csv.DictReader(open(a.catalogue, newline=""), delimiter="\t"))
    print(f"  {len(cat):,} accessory loci")
    if not os.path.exists(a.fasta + ".bwt"):
        subprocess.run([a.bwa, "index", a.fasta], capture_output=True)

    # ---- the read route ---------------------------------------------------
    import tempfile
    tmp = tempfile.mkdtemp()
    q = os.path.join(tmp, "q.fa")
    n_q = 0
    with open(q, "w") as fh:
        p = subprocess.Popen([a.samtools, "view", "--reference", a.h37rv,
                              "-f", "4", a.cram], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True)
        for i, line in enumerate(p.stdout):
            f = line.split("\t", 11)
            if len(f) > 9:
                fh.write(f">u{i}\n{f[9]}\n"); n_q += 1
        p.wait()
        p = subprocess.Popen([a.samtools, "view", "--reference", a.h37rv,
                              "-q", str(a.min_mapq), a.cram],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
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
        p.wait()
    print(f"  {n_q:,} query reads from the unmapped and clipped pools")

    bam = os.path.join(tmp, "a.bam")
    m = subprocess.Popen([a.bwa, "mem", "-t", "2", "-k", "19", a.fasta, q],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run([a.samtools, "sort", "-o", bam, "-"], stdin=m.stdout,
                   capture_output=True)
    m.wait()
    subprocess.run([a.samtools, "index", bam], capture_output=True)

    cov0, cov20 = {}, {}
    for Q, store in ((0, cov0), (a.min_mapq, cov20)):
        r = subprocess.run([a.samtools, "depth", "-a", "-Q", str(Q), bam],
                           capture_output=True, text=True)
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
        if c0 >= a.present_cov:
            state, why = "PRESENT", "reads cover the locus"
        elif c0 <= a.absent_cov:
            state, why = "ABSENT", "reads do not cover the locus"
        else:
            state, why = "UNCERTAIN", f"partial coverage {c0:.2f}"
        rc = ref_carries.get(lid)
        # The reference route is recorded next to the read route rather than
        # overriding it: where they disagree, the sample differs from its own
        # reference, which is a real finding and not an error to hide.
        agree = "" if rc is None else (
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
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
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
