#!/usr/bin/env python3
"""Checks that tie a P0 build's assets to the graph they claim to describe.

WHY THIS EXISTS (audit section B). P0 took several assets from fixed paths --
the panel SNP matrix from graphs/CX333..., the IS6110 crossmaps from
is6110/assets, the ancestral tree from data/trees/cx333.* -- so a build of a
NEW graph would have been given the old graph's assets with no error. Each
check here compares an asset with the build's own record of the graph
(paths.txt, accessions.txt, refs/) and exits non-zero on a mismatch, so P0
refuses to mark the step done.

    p0_check.py vcf-samples   --vcf V --accessions A [--ref-acc H37RV]
    p0_check.py fasta-names   --fasta F --paths P
    p0_check.py taxa          --fasta F --accessions A
    p0_check.py catalogue     --tsv T --accessions A
    p0_check.py is6110-intervals --isclean-dir D --refs R --accessions A
                              --skip H37RV --out O
    p0_check.py verify        --build B

`verify` is the check at USE: every manifest row under assets/ and
annotation/ is re-hashed, and an asset that is a link leaving the build is
refused. A build is immutable only if nothing outside it can change it (the
repeat mask changed 181 -> 395 intervals under a recorded checksum, audit
P0P2-2).
"""
import argparse, csv, gzip, hashlib, os, subprocess, sys


def lines(p):
    return [l.strip() for l in open(p) if l.strip()]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fasta_names(p):
    op = gzip.open if p.endswith(".gz") else open
    fai = p + ".fai"
    if os.path.exists(fai) and not p.endswith(".gz"):
        return [l.split("\t", 1)[0] for l in open(fai) if l.strip()]
    with op(p, "rt") as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def report_sets(what, have, want):
    extra, miss = sorted(have - want), sorted(want - have)
    if extra or miss:
        print(f"FATAL: {what}: {len(extra)} not in the build "
              f"({', '.join(extra[:5])}{' ...' if len(extra) > 5 else ''}); "
              f"{len(miss)} of the build's missing "
              f"({', '.join(miss[:5])}{' ...' if len(miss) > 5 else ''})",
              file=sys.stderr)
        return 1
    print(f"  {what}: {len(have)} names, identical to the build's")
    return 0


def cmd_vcf_samples(a):
    r = subprocess.run([a.bcftools, "query", "-l", a.vcf],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FATAL: bcftools query -l {a.vcf}: {r.stderr.strip()}")
    have = set(r.stdout.split())
    # the deconstruct reference path is not a sample column
    if a.ref_acc:
        have.add(a.ref_acc)
    return report_sets(f"samples of {a.vcf}", have, set(lines(a.accessions)))


def cmd_fasta_names(a):
    return report_sets(f"sequences of {a.fasta}", set(fasta_names(a.fasta)),
                       set(lines(a.paths)))


def cmd_taxa(a):
    have = set(fasta_names(a.fasta))
    want = set(lines(a.accessions))
    if a.ref_acc:
        have.add(a.ref_acc)
    return report_sets(f"taxa of {a.fasta}", have, want)


def cmd_catalogue(a):
    want = set(lines(a.accessions))
    seen = set()
    for r in csv.DictReader(open(a.tsv, newline=""), delimiter="\t"):
        seen.update(x for x in (r.get("panel_carriers") or "").split(",") if x)
    bad = sorted(seen - want)
    if bad:
        print(f"FATAL: {a.tsv} names {len(bad)} carriers that are not panel "
              f"genomes of this build ({', '.join(bad[:5])})", file=sys.stderr)
        return 1
    print(f"  {a.tsv}: {len(seen)} carrier genomes, all in the build")
    return 0


def cmd_is6110_intervals(a):
    """Per-reference interval counts for the P1 tie-break, checked against
    the build's own references: a crossmap made from another sequence is
    refused, so counts measured on an old panel's assemblies cannot be
    carried into a new build."""
    out, missing, wrong = [], [], []
    for acc in lines(a.accessions):
        if acc == a.skip:
            continue
        xmap = os.path.join(a.isclean_dir, f"{acc}.crossmap.tsv")
        cfai = os.path.join(a.isclean_dir, f"{acc}.isclean.fasta.fai")
        rfai = os.path.join(a.refs, f"{acc}.fasta.fai")
        if not (os.path.exists(xmap) and os.path.exists(cfai)):
            missing.append(acc); continue
        rows = list(csv.DictReader(open(xmap, newline=""), delimiter="\t"))
        removed = int(rows[-1]["cum_deleted"]) if rows else 0
        cname, clen = open(cfai).readline().split("\t")[:2]
        rname, rlen = open(rfai).readline().split("\t")[:2]
        # p1i_build_matched.sh names the clean chromosome <contig>_isclean
        if cname.endswith("_isclean"):
            cname = cname[:-len("_isclean")]
        # the clean chromosome plus what was excised must be the build's
        # reference, by name and by length
        if cname != rname or int(clen) + removed != int(rlen):
            wrong.append(f"{acc} ({cname}:{int(clen) + removed} vs "
                         f"{rname}:{rlen})")
            continue
        out.append((acc, len(rows), removed, sha256(xmap)))
    if missing or wrong:
        if missing:
            print(f"FATAL: {len(missing)} accessions have no crossmap in "
                  f"{a.isclean_dir} ({', '.join(missing[:5])}); build them "
                  f"for every panel genome (is6110/bin/p1i_build_matched.sh "
                  f"with REFS=<build>/refs)", file=sys.stderr)
        if wrong:
            print(f"FATAL: {len(wrong)} crossmaps were not made from this "
                  f"build's references: {', '.join(wrong[:5])}",
                  file=sys.stderr)
        return 1
    with open(a.out, "w") as fh:
        fh.write("reference\tintervals\tbp_removed\tcrossmap_sha256\n")
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"  {len(out)} references with interval counts -> {a.out}")
    return 0


def cmd_verify(a):
    b = a.build.rstrip("/")
    man = os.path.join(b, "manifest.tsv")
    if not os.path.exists(man):
        print(f"FATAL: no manifest in {b}", file=sys.stderr)
        return 1
    real_b = os.path.realpath(b)
    bad = []
    # links that leave the build: the build is not immutable while they exist
    for sub in ("assets", "annotation"):
        d = os.path.join(b, sub)
        for root, _, files in os.walk(d):
            for f in files:
                p = os.path.join(root, f)
                if os.path.islink(p) and not os.path.realpath(p).startswith(
                        real_b + os.sep):
                    bad.append(f"{p} links outside the build, to "
                               f"{os.path.realpath(p)}")
    n = 0
    for r in csv.DictReader(open(man, newline=""), delimiter="\t"):
        # the recorded path is relative to where P0 ran; locate the file by
        # its place inside the build, so a moved build still verifies
        parts = r["path"].split(os.sep)
        top = next((x for x in ("assets", "annotation") if x in parts), "")
        if not top:
            continue
        p = r["path"]
        loc = os.path.join(b, *parts[parts.index(top):])
        if not os.path.exists(loc):
            bad.append(f"{r['asset']}: recorded, now missing ({p})")
            continue
        n += 1
        if sha256(loc) != r["sha256"]:
            bad.append(f"{r['asset']}: content differs from the manifest "
                       f"({loc})")
    if bad:
        print(f"FATAL: build {b} does not match its manifest:", file=sys.stderr)
        for x in bad[:20]:
            print(f"  {x}", file=sys.stderr)
        if len(bad) > 20:
            print(f"  ... and {len(bad) - 20} more", file=sys.stderr)
        return 1
    print(f"  build {b}: {n} manifest assets verified")
    return 0


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    bcf = os.environ.get("MTB_BCFTOOLS", "bcftools")
    p = sp.add_parser("vcf-samples")
    p.add_argument("--vcf", required=True)
    p.add_argument("--accessions", required=True)
    p.add_argument("--ref-acc", default="")
    p.add_argument("--bcftools", default=bcf)
    p = sp.add_parser("fasta-names")
    p.add_argument("--fasta", required=True)
    p.add_argument("--paths", required=True)
    p = sp.add_parser("taxa")
    p.add_argument("--fasta", required=True)
    p.add_argument("--accessions", required=True)
    p.add_argument("--ref-acc", default="")
    p = sp.add_parser("catalogue")
    p.add_argument("--tsv", required=True)
    p.add_argument("--accessions", required=True)
    p = sp.add_parser("is6110-intervals")
    p.add_argument("--isclean-dir", required=True)
    p.add_argument("--refs", required=True)
    p.add_argument("--accessions", required=True)
    p.add_argument("--skip", default="")
    p.add_argument("--out", required=True)
    p = sp.add_parser("verify")
    p.add_argument("--build", required=True)
    a = ap.parse_args()
    return {"vcf-samples": cmd_vcf_samples, "fasta-names": cmd_fasta_names,
            "taxa": cmd_taxa, "catalogue": cmd_catalogue,
            "is6110-intervals": cmd_is6110_intervals,
            "verify": cmd_verify}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
