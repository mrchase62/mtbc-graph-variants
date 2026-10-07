#!/usr/bin/env python3
"""Collate P2's per-isolate callsets: what was produced, and the burden per arm.

There is no truth for real isolates, so this reports counts and coverage rather
than precision. The number to watch against T8 is the small-variant burden
against the matched reference: T8 measured 88.1% fewer calls than against H37Rv
on the screened cohort, and a burden far from that would mean the reference is
not being used as intended rather than that the biology changed.
"""
import argparse, csv, gzip, os, subprocess, sys


def count_records(path):
    if not path or not os.path.exists(path):
        return None
    op = gzip.open if path.endswith(".gz") else open
    n = 0
    with op(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                n += 1
    return n


def stamped(path):
    """yes / no / unknown -- never conflate "could not check" with "not stamped".

    The first version shelled out to bcftools and returned False on OSError, so
    when project_env.sh set MTB_BCFTOOLS without exporting it, Python fell back to
    a bare `bcftools` that is not on PATH and every one of 23 correctly stamped
    VCFs was reported as unstamped. A provenance check that fails in the alarming
    direction wastes an investigation; failing in the reassuring direction would
    be worse. Read the header directly and let a third state carry the difference.
    """
    if not os.path.exists(path):
        return "missing"
    try:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if line.startswith("##MTB_graph_build="):
                    return "yes"
        return "no"
    except (OSError, EOFError) as e:
        print(f"  cannot read header of {path}: {e}", file=sys.stderr)
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--dir", default="refbias/p2")
    ap.add_argument("--work", default="refbias/work/p2")
    ap.add_argument("--p1-work", default="refbias/work/p1",
                    help="pass-1 H37Rv VCFs, for the burden comparison")
    ap.add_argument("--out", default="refbias/p2/p2_summary.tsv")
    ap.add_argument("--allow-missing", action="store_true",
                    help="exit 0 even when some samples have no output. Off by "
                         "default: a short table silently drops those samples "
                         "from every later pass, because P5 sizes itself from "
                         "the refmap")
    a = ap.parse_args()

    rows, incomplete = [], []
    for r in csv.DictReader(open(a.refmap), delimiter="\t"):
        s = r["sample"]
        small = os.path.join(a.dir, f"{s}.vcf.gz")
        gvcf = os.path.join(a.dir, f"{s}.g.vcf.gz")
        delly = os.path.join(a.dir, f"{s}.delly.vcf")
        dysgu = os.path.join(a.dir, f"{s}.dysgu.vcf")
        h37 = os.path.join(a.p1_work, f"{s}.h37rv.vcf.gz")
        ns, nd, ny = count_records(small), count_records(delly), count_records(dysgu)
        nh = count_records(h37)
        if ns is None or nd is None:
            incomplete.append(s); continue
        rows.append(dict(
            sample=s, reference=r["reference"], snp_distance=r["snp_distance"],
            h37rv_small=nh if nh is not None else "",
            matched_small=ns,
            reduction=(round(1 - ns / nh, 4) if nh else ""),
            delly=nd, dysgu=ny if ny is not None else "",
            gvcf=("yes" if os.path.exists(gvcf) else "no"),
            bam=("yes" if os.path.exists(os.path.join(a.work, f"{s}.bam")) else "no"),
            stamped=stamped(small)))

    if not rows:
        print("no P2 output found", file=sys.stderr); return 1
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    # Say so when the H37Rv side is missing: a blank h37rv_small disables
    # p5_sanity's burden check, which then reports 0 rather than failing.
    no_h37 = sum(1 for r in rows if r["h37rv_small"] == "")
    if no_h37:
        print(f"  WARNING: {no_h37} of {len(rows)} isolates have no P1 H37Rv VCF "
              f"in {a.p1_work}; h37rv_small and reduction are blank for them",
              file=sys.stderr)

    print(f"  {len(rows)} isolates with P2 output"
          f"{'; INCOMPLETE: ' + ', '.join(incomplete) if incomplete else ''}\n")
    print(f"  {'isolate':<17s}{'reference':<16s}{'d':>6s}{'H37Rv':>8s}"
          f"{'matched':>9s}{'cut':>8s}{'delly':>7s}{'dysgu':>7s}  stamp")
    # snp_distance can legitimately be blank. An arm that PINS every isolate to
    # one reference -- the two-reference comparison does exactly this -- has no
    # selection distance to report, and requiring one killed that arm's summary
    # with a ValueError, which afterok then propagated into nine cancelled
    # jobs. Blank sorts last rather than crashing.
    def _dist(x):
        v = (x.get("snp_distance") or "").strip()
        try:
            return int(v)
        except ValueError:
            return 1 << 30
    for r in sorted(rows, key=_dist):
        red = f"{100*r['reduction']:.1f}%" if r["reduction"] != "" else "-"
        print(f"  {r['sample']:<17s}{r['reference']:<16s}{r['snp_distance']:>6s}"
              f"{str(r['h37rv_small']):>8s}{r['matched_small']:>9d}{red:>8s}"
              f"{r['delly']:>7d}{str(r['dysgu']):>7s}  {r['stamped']}")

    tot_h = sum(int(r["h37rv_small"]) for r in rows if r["h37rv_small"] != "")
    tot_m = sum(r["matched_small"] for r in rows if r["h37rv_small"] != "")
    print(f"\n  small-variant burden, pooled: H37Rv {tot_h}, matched {tot_m}"
          + (f", {100*(1-tot_m/tot_h):.1f}% fewer" if tot_h else ""))
    print(f"    T8 on the screened 23: 88.1% fewer")
    print(f"  structural: delly {sum(r['delly'] for r in rows)} records, "
          f"dysgu {sum(r['dysgu'] for r in rows if r['dysgu'] != '')}")
    nb = sum(1 for r in rows if r["bam"] == "yes")
    ng = sum(1 for r in rows if r["gvcf"] == "yes")
    unst = [r["sample"] for r in rows if r["stamped"] != "yes"]
    print(f"  retained: {nb}/{len(rows)} BAMs, {ng}/{len(rows)} GVCFs")
    if unst:
        print(f"  small-variant VCFs not confirmed stamped: {', '.join(unst)}",
              file=sys.stderr)
    else:
        print(f"  all {len(rows)} small-variant VCFs carry the build stamp")
    print(f"\n  written: {a.out}")
    if incomplete and not a.allow_missing:
        print(f"FATAL: {len(incomplete)} samples have no usable output "
              f"({', '.join(incomplete[:10])}{' ...' if len(incomplete) > 10 else ''}); "
              f"rerun them, or pass --allow-missing to accept a partial table",
              file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
