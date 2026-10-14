#!/usr/bin/env python3
"""Step 7 of the fold-back vs mixed-call proposal: what supports the mixed
calls, read by read. A diagnostic only; nothing is removed or recalled.

Isolates: up to --per-lineage high (fold-back > 4.8 per 1,000, the chimeric
controls' range) and as many low (<= 0.31, the clean controls' range) per
lineage 1-4, drawn at random (fixed seed) from isolates.tsv.

For every unmasked mixed call (mixed_calls.tsv.gz) of those isolates, the
CRAM is read in place at that position (MAPQ >= 20, base quality >= 20,
duplicates/secondary/QC-fail skipped, as Mutect2 would see them). Each read
is classed:
  foldback  SA tag to the same contig, opposite strand, within 1 kb
            (the split_reads.awk definition)
  split     any other SA tag
  clipped   soft clip >= 20 bp, no SA
  plain     none of these
and counted separately for reads carrying the ALT and the REF allele.

    readlevel.py --isolates out/isolates.tsv --calls out/mixed_calls.tsv.gz \
        --out out/readlevel.tsv [--per-lineage 25] [--threads 8]
"""
import argparse, collections, csv, gzip, json, os, random
from multiprocessing import Pool
import pysam

LINEAGES = ("lineage1", "lineage2", "lineage3", "lineage4")
ROOT = os.environ["MTB_CRAM_ROOT"]
CREF = os.environ["MTB_CRAM_REF"]


def read_class(r):
    if r.has_tag("SA"):
        for sa in r.get_tag("SA").rstrip(";").split(";"):
            ch, pos, strand = sa.split(",")[:3]
            if ch == r.reference_name and (strand == "-") != r.is_reverse \
                    and abs(int(pos) - 1 - r.reference_start) <= 1000:
                return "foldback"
        return "split"
    ct = r.cigartuples or []
    if any(op == 4 and n >= 20 for op, n in ct):
        return "clipped"
    return "plain"


def one(job):
    s, path, calls = job
    out = []
    try:
        cram = pysam.AlignmentFile(f"{ROOT}/{path}", "rc", reference_filename=CREF)
    except (OSError, ValueError) as e:
        return [(s, "error", str(e))]
    chrom = cram.references[0]
    for c in calls:
        pos, ref, alt, typ = int(c["pos"]), c["ref"], c["alt"], c["type"]
        n = {f"{a}_{k}": 0 for a in ("alt", "ref") for k in ("foldback", "split", "clipped", "plain")}
        for col in cram.pileup(chrom, pos - 1, pos, truncate=True, min_mapping_quality=20,
                               min_base_quality=20, ignore_orphans=False, stepper="all",
                               flag_filter=0x704):
            for pr in col.pileups:
                r = pr.alignment
                if typ == "snv":
                    if pr.is_del or pr.is_refskip:
                        continue
                    b = r.query_sequence[pr.query_position]
                    allele = "alt" if b == alt else ("ref" if b == ref else None)
                else:
                    ind = pr.indel
                    want = len(alt) - len(ref)
                    allele = "alt" if ind == want and want != 0 else ("ref" if ind == 0 and not pr.is_del else None)
                if allele:
                    n[f"{allele}_{read_class(r)}"] += 1
        out.append([s, pos, ref, alt, typ, c["pass"], c["af"], c["barcode"]] + [n[k] for k in sorted(n)])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isolates", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-lineage", type=int, default=25)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20261013)
    a = ap.parse_args()
    iso = list(csv.DictReader(open(a.isolates), delimiter="\t"))
    rnd = random.Random(a.seed)
    pick = {}
    for l in LINEAGES:
        for grp, test in (("high", lambda x: x > 4.8), ("low", lambda x: x <= 0.31)):
            pool = sorted(r["sample"] for r in iso if r["lineage"] == l and test(float(r["inv_chimera_per_1k"])))
            for s in rnd.sample(pool, min(a.per_lineage, len(pool))):
                pick[s] = grp
    meta = json.load(gzip.open(f"{ROOT}/metadata/samples.json.gz"))
    calls = collections.defaultdict(list)
    with gzip.open(a.calls, "rt") as fh:
        for c in csv.DictReader(fh, delimiter="\t"):
            if c["sample"] in pick:
                calls[c["sample"]].append(c)
    jobs = [(s, meta[s]["path"], calls[s]) for s in pick if s in meta]
    print(f"isolates: {collections.Counter(pick.values())}; calls: {sum(len(v) for v in calls.values())}", flush=True)
    keys = sorted(f"{x}_{k}" for x in ("alt", "ref") for k in ("foldback", "split", "clipped", "plain"))
    with Pool(a.threads) as P, open(a.out, "w") as fo:
        fo.write("\t".join(["sample", "group", "pos", "ref", "alt", "type", "pass", "af", "barcode"] + keys) + "\n")
        for rows in P.imap_unordered(one, jobs):
            for r in rows:
                if r[1] == "error":
                    print("ERROR", r, flush=True); continue
                fo.write("\t".join(map(str, [r[0], pick[r[0]]] + r[1:])) + "\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
