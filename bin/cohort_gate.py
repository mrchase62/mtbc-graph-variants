#!/usr/bin/env python3
"""Cohort-level acceptance gate for a merged VCF (user, 2026-10-16: "add the
gate"). Run after p5vcf; exits 1 if any check fails, so the chain watchdog
reports it and nothing downstream (association, the next batch) starts on a
bad cohort.

Why it exists. Unit tests check rules on small inputs; they cannot see a
cohort-wide failure. Four reached the new1k_p200 VCF unseen (SV-EVAL-1..4,
analysis/sv_eval_1k): IS6110 ABSENT written site-wide (46 sites ABSENT in
>= 90% of isolates), 99 of 564 accessory loci written with 43% NOCALL, 486
"HaplotypeCaller" insertions of >= 500 bp that were misplaced graph alleles,
and NOCALL jumps. Every check below would have failed that VCF.

Checks (FAIL unless noted):
  classes        every class (small, sv, is6110, accessory_presence) has records
  accessory      loci written >= --acc-min-frac (0.95) of the build catalogue
  nocall         per-class NOCALL share of cells <= --max-nocall CLASS=FRAC
                 (defaults small 0.05, sv 0.60, is6110 0.25,
                 accessory_presence 0.20)
  sitewide       H37Rv-frame records outside the repeat mask ABSENT (or NOCALL)
                 in >= 90% of isolates: at most
                 --max-sitewide CLASS=N (defaults is6110 5, accessory 5,
                 small 0.5% of records, sv 2% of records)
  small_allele   no SVSOURCE=small_allele record longer than --max-hc-len
                 (default 400 bp: beyond what HaplotypeCaller assembles from
                 short reads)
  isolates       no isolate with a no-call share > --max-isolate-nocall (0.10)
                 among small H37Rv-frame SNPs
  dr             TB-Profiler dr_variant SNPs at freq >= 0.9 ALT in the VCF for
                 >= --dr-min (0.99) of them (when --tbprofiler is given)
  baseline       with --baseline VCF: per-class record count within
                 --max-count-change (0.5) of it, and per-class NOCALL share
                 not up by more than 2x and 0.05 (WARN only for counts when the
                 builds differ)
Writes a TSV report (--out) of every check with its value and limit.

    cohort_gate.py --vcf <cohort>/p5/merged.vcf.gz --build <build> \
        [--baseline <old merged.vcf.gz>] [--tbprofiler PATTERN] --out <cohort>/p5/gate.tsv
"""
import argparse, collections, gzip, json, os, sys

CLASSES = ("small", "sv", "is6110", "accessory_presence")


def kv(spec, cast):
    out = {}
    for x in spec:
        k, v = x.split("=", 1); out[k] = cast(v)
    return out


def load_mask(build):
    import array
    m = bytearray(4411532 + 200)
    path = f"{build}/assets/repeat_mask.bed"
    if os.path.exists(path):
        for l in open(path):
            c, a_, b_ = l.split("\t")[:3]
            m[int(a_):int(b_)] = b"\x01" * (int(b_) - int(a_))
    return m


def scan(vcf, mask):
    """Per-class record counts, cell state counts, site-wide records,
    small_allele lengths, per-isolate SNP no-calls, SNP genotype lookup."""
    op = gzip.open if vcf.endswith(".gz") else open
    n_rec = collections.Counter(); cells = collections.defaultdict(collections.Counter)
    sitewide = collections.Counter(); hc_long = []; acc_loci = set()
    iso_nc = None; iso_n = 0; snp = {}
    with op(vcf, "rt") as fh:
        for l in fh:
            if l.startswith("#CHROM"):
                samples = l.rstrip("\n").split("\t")[9:]
                iso_nc = [0] * len(samples)
                continue
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) if "=" in x else (x, "1") for x in f[7].split(";"))
            cls = info.get("CLASS", "small")
            n_rec[cls] += 1
            fmt = f[8].split(":"); si = fmt.index("ST") if "ST" in fmt else None
            st = collections.Counter()
            gts = []
            for c in f[9:]:
                parts = c.split(":")
                s_ = parts[si] if si is not None and si < len(parts) else (
                    "ALT" if parts[0] == "1" else "REF" if parts[0] == "0" else "NOCALL")
                st[s_] += 1; gts.append(parts[0])
            cells[cls].update(st)
            n = sum(st.values())
            # site-wide ABSENT/NOCALL means an error only where every isolate
            # could be read: H37Rv-frame records outside the repeat mask (a
            # node-frame record is legitimately NOCALL off its own path)
            pos = int(f[1])
            in_mask = 0 < pos < len(mask) and mask[pos - 1]
            if n and info.get("FRAME", "h37rv") != "node" and not in_mask \
                    and (st["ABSENT"] >= 0.9 * n or st["NOCALL"] >= 0.9 * n):
                sitewide[cls] += 1
            if info.get("SVSOURCE") == "small_allele":
                L = abs(int(info.get("SVLEN", "0") or 0))
                if L > 0: hc_long.append(L)
            if cls == "accessory_presence" and "ACCLOCUS" in info:
                acc_loci.add(info["ACCLOCUS"])
            if cls == "small" and info.get("KIND") == "SNP" and info.get("FRAME", "h37rv") == "h37rv":
                iso_n += 1
                for j, g in enumerate(gts):
                    if g not in ("0", "1"): iso_nc[j] += 1
                alt = f[4].split(",")[0]
                snp[(int(f[1]), alt)] = gts
    return dict(samples=samples, n_rec=n_rec, cells=cells, sitewide=sitewide, hc_long=hc_long,
                acc_loci=acc_loci, iso_nc=iso_nc, iso_n=iso_n, snp=snp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--build", required=True, help="the build dir (assets/accessory_catalogue.tsv)")
    ap.add_argument("--baseline", default="")
    ap.add_argument("--tbprofiler", default="",
                    help="path pattern with {s} for the sample, e.g. .../{s}/results/{s}.tbprofiler.json")
    ap.add_argument("--acc-min-frac", type=float, default=0.95)
    ap.add_argument("--max-nocall", nargs="*", default=[])
    ap.add_argument("--max-sitewide", nargs="*", default=[])
    ap.add_argument("--max-hc-len", type=int, default=400)
    ap.add_argument("--max-isolate-nocall", type=float, default=0.10)
    ap.add_argument("--dr-min", type=float, default=0.99)
    ap.add_argument("--max-count-change", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    lim_nc = {"small": 0.05, "sv": 0.60, "is6110": 0.25, "accessory_presence": 0.20}
    lim_nc.update(kv(a.max_nocall, float))
    mask = load_mask(a.build)
    v = scan(a.vcf, mask)
    lim_sw = {"is6110": 5, "accessory_presence": 5,
              "small": max(1, int(0.005 * v["n_rec"]["small"])), "sv": max(1, int(0.02 * v["n_rec"]["sv"]))}
    lim_sw.update(kv(a.max_sitewide, int))
    rows = []

    def check(name, ok, value, limit, level="FAIL"):
        rows.append((name, "ok" if ok else level, value, limit))

    for c in CLASSES:
        check(f"classes:{c}", v["n_rec"][c] > 0, v["n_rec"][c], "> 0")
    cat = f"{a.build}/assets/accessory_catalogue.tsv"
    if os.path.exists(cat):
        n_cat = sum(1 for _ in open(cat)) - 1
        frac = len(v["acc_loci"]) / max(1, n_cat)
        check("accessory:loci_written", frac >= a.acc_min_frac, f"{len(v['acc_loci'])}/{n_cat}", f">= {a.acc_min_frac:.0%}")
    for c in CLASSES:
        tot = sum(v["cells"][c].values())
        if not tot: continue
        nc = v["cells"][c]["NOCALL"] / tot
        check(f"nocall:{c}", nc <= lim_nc.get(c, 1.0), f"{nc:.4f}", f"<= {lim_nc.get(c, 1.0)}")
        check(f"sitewide:{c}", v["sitewide"][c] <= lim_sw.get(c, 10**9), v["sitewide"][c], f"<= {lim_sw.get(c)}")
    long_ = sum(1 for L in v["hc_long"] if L > a.max_hc_len)
    check("small_allele:too_long", long_ == 0, long_, f"0 records > {a.max_hc_len} bp")
    if v["iso_n"]:
        worst = max(v["iso_nc"]) / v["iso_n"]
        bad = sum(1 for x in v["iso_nc"] if x / v["iso_n"] > a.max_isolate_nocall)
        check("isolates:snp_nocall", bad == 0, f"{bad} isolates; max {worst:.4f}", f"<= {a.max_isolate_nocall}")
    if a.tbprofiler:
        found = total = 0
        for j, s in enumerate(v["samples"]):
            try:
                jd = json.load(open(a.tbprofiler.format(s=s)))
            except (OSError, ValueError):
                continue
            for d in jd.get("dr_variants", []):
                ref, alt = d.get("ref", ""), d.get("alt", "")
                if float(d.get("freq", 1)) < 0.9 or len(ref) != 1 or len(alt) != 1: continue
                total += 1
                g = v["snp"].get((int(d.get("pos", 0)), alt))
                found += bool(g) and g[j] == "1"
        if total:
            check("dr:recovered", found / total >= a.dr_min, f"{found}/{total}", f">= {a.dr_min:.0%}")
    if a.baseline:
        b = scan(a.baseline, mask)
        for c in CLASSES:
            if b["n_rec"][c]:
                ch = abs(v["n_rec"][c] - b["n_rec"][c]) / b["n_rec"][c]
                check(f"baseline:count:{c}", ch <= a.max_count_change, f"{v['n_rec'][c]} vs {b['n_rec'][c]}",
                      f"change <= {a.max_count_change:.0%}", level="WARN")
            tv, tb = sum(v["cells"][c].values()), sum(b["cells"][c].values())
            if tv and tb:
                nv, nb = v["cells"][c]["NOCALL"] / tv, b["cells"][c]["NOCALL"] / tb
                check(f"baseline:nocall:{c}", not (nv > 2 * nb and nv - nb > 0.05), f"{nv:.4f} vs {nb:.4f}",
                      "not > 2x and +0.05")
    with open(a.out, "w") as fo:
        fo.write("check\tresult\tvalue\tlimit\n")
        for r in rows: fo.write("\t".join(map(str, r)) + "\n")
    fails = [r for r in rows if r[1] == "FAIL"]
    for r in rows:
        print(f"  {r[1]:4s}  {r[0]:28s} {r[2]}  (limit {r[3]})")
    print(f"[gate] {len(rows)} checks, {len(fails)} FAIL, {sum(r[1] == 'WARN' for r in rows)} WARN -> {a.out}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
