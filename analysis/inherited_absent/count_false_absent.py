#!/usr/bin/env python3
"""For one sample, test every ABSENT call in a merged cohort VCF against the
sample's own reads on H37Rv.

ABSENT says the sample lacks the region. P5 mostly gets it from the matched
reference: a site off R's path where R has a deletion. Nothing then checks the
sample's own H37Rv alignment, which is how SAMEA2297133 (intact ctpV at
72-107x) came out ABSENT inside ctpV (analysis/ctpV_check.md). This measures
how often that happens.

Depth is samtools depth -a with MAPQ >= 20, over each record's H37Rv span:
the REF allele for small records, POS+1..END for symbolic deletions. The
sample's baseline is its genome-wide median of that depth. Classes:

  contradicted  mean >= 0.5 x baseline and no zero-depth base: reads say the
                region is present
  partial       anything between
  supported     mean < 0.1 x baseline: reads agree it is gone

Records with no H37Rv span (node-frame, IS6110, accessory presence) are
counted as not testable.

Report only; nothing in the pipeline reads this.
"""
import argparse
import gzip
import os
import subprocess
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--samtools", required=True)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--out", required=True, help="per-record table")
    ap.add_argument("--summary", required=True, help="one-row summary")
    a = ap.parse_args()

    # ---- depth over the whole chromosome, MAPQ-filtered
    p = subprocess.run([a.samtools, "depth", "-a", "-Q", str(a.min_mapq), a.bam],
                       check=True, capture_output=True, text=True)
    rows = [l.split("\t") for l in p.stdout.splitlines()]
    chrom = rows[0][0]
    n = max(int(r[1]) for r in rows if r[0] == chrom)
    depth = np.zeros(n + 2, dtype=np.int32)
    for r in rows:
        if r[0] == chrom:
            depth[int(r[1])] = int(r[2])
    del rows, p
    base = float(np.median(depth[1:n + 1]))

    # ---- the sample's ABSENT records
    out, counts = [], {"contradicted": 0, "partial": 0, "supported": 0,
                       "not_testable": 0}
    by_class = {}
    with gzip.open(a.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                col = 9 + line.rstrip("\n").split("\t")[9:].index(a.sample)
                continue
            c = line.rstrip("\n").split("\t")
            fmt = c[8].split(":")
            val = dict(zip(fmt, c[col].split(":")))
            if val.get("ST") != "ABSENT":
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, True)
                        for kv in c[7].split(";"))
            klass = info.get("CLASS", "?")
            pos = int(c[1])
            if info.get("FRAME") == "node" or klass not in ("small", "sv"):
                counts["not_testable"] += 1
                by_class.setdefault((klass, "not_testable"), 0)
                by_class[(klass, "not_testable")] += 1
                continue
            if "END" in info:
                s, e = pos + 1, int(info["END"])
            else:
                s, e = pos, pos + len(c[3]) - 1
            if s > e:
                s, e = pos, pos
            d = depth[s:e + 1]
            mean = float(d.mean())
            zero = int((d == 0).sum())
            if mean >= 0.5 * base and zero == 0:
                cls = "contradicted"
            elif mean < 0.1 * base:
                cls = "supported"
            else:
                cls = "partial"
            counts[cls] += 1
            by_class.setdefault((klass, cls), 0)
            by_class[(klass, cls)] += 1
            out.append((a.sample, c[2], klass, info.get("REGION", ""),
                        s, e, e - s + 1, f"{mean:.1f}", zero, cls))

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fo:
        fo.write("sample\tid\tclass\tregion\tstart\tend\tlen\tmean_depth"
                 "\tzero_bases\tverdict\n")
        for r in out:
            fo.write("\t".join(map(str, r)) + "\n")
    with open(a.summary, "w") as fo:
        keys = [(k, v) for k in ("small", "sv")
                for v in ("contradicted", "partial", "supported")]
        keys.append(("other", "not_testable"))
        by_class[("other", "not_testable")] = sum(
            v for (k, c), v in by_class.items() if c == "not_testable")
        fo.write("sample\tbaseline_depth\tabsent_total\t"
                 + "\t".join(counts) + "\t"
                 + "\t".join(f"{k}:{v}" for k, v in keys) + "\n")
        fo.write(f"{a.sample}\t{base:.0f}\t{sum(counts.values())}\t"
                 + "\t".join(str(v) for v in counts.values()) + "\t"
                 + "\t".join(str(by_class.get(k, 0)) for k in keys) + "\n")
    print(a.sample, base, counts, file=sys.stderr)


if __name__ == "__main__":
    main()
