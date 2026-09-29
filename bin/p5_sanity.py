#!/usr/bin/env python3
"""Sanity checks over the merged matrix, run after a complete pass.

Seven checks, each chosen because it can fail in a way inspection would not
reveal. None of them needs variant-level truth.

  1  ubiquitous sites      A site all 23 samples carry is a reference-artefact
                           candidate, not a shared variant. Decisive test: look
                           up the panel allele frequency in the graph VCF. If
                           H37Rv carries the minor allele at these positions,
                           they are artefacts of the frame, not biology.
  2  singletons by region  Singletons should concentrate in core sequence simply
                           because core is most of the genome. Enrichment in
                           PE/PPE instead would mean private variation is being
                           manufactured where mapping is hardest.
  3  ALT counts vs pass 1  The matrix is H37Rv-FRAMED, so its per-sample ALT
                           count should track pass one's H37Rv burden, not P2's
                           burden against the matched reference. Comparing it to
                           P2 is a category error: P2 counts differences from R,
                           which matching deliberately minimises, so the ratio
                           would read 2696 for M. canettii and mean nothing.
  4  state balance         No sample should have a wildly different state
                           profile. An outlier means one sample's GVCF or
                           projection behaved differently from the rest.
  5  ABSENT vs distance    ABSENT means a position with no equivalent in that
                           sample's reference, so it should rise with reference
                           distance. If it does not, the ABSENT call is not
                           measuring what it claims.
  6  NOCALL vs depth       NOCALL should be commoner in lower-depth samples. No
                           relationship would suggest NOCALL is an artefact of
                           the projection rather than of observation.
  7  kind balance          The SNP:INDEL ratio should be similar across samples
                           and plausible for MTBC; a sample far off it had a
                           different calling outcome, not different biology.
"""
import argparse, collections, csv, gzip, os, statistics, sys


def rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def _snp_dist(row):
    """snp_distance, or a value that sorts last when it is blank.

    An arm that PINS every isolate to one reference -- the two-reference
    comparison does exactly that -- has no selection distance to report. Three
    separate summaries assumed the field was always an integer, and each one
    died with the same ValueError and took the rest of the chain with it
    through afterok. Fixed in p2_summary.py first, then here, which is the
    argument for fixing a class of defect rather than its instances.
    """
    v = (row.get("snp_distance") or "").strip()
    try:
        return int(v)
    except ValueError:
        return 1 << 30


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="refbias/p5/matrix.tsv")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--p2-summary", default="refbias/p2/p2_summary.tsv")
    ap.add_argument("--cohort", default="refbias/cohort.pilot.tsv")
    ap.add_argument("--graph-vcf",
                    default="graphs/CX333.s10k.k23.K15/all_variants.nolab.vcf.gz")
    ap.add_argument("--out", default="refbias/p5/sanity.tsv")
    a = ap.parse_args()

    meta = {r["sample"]: r for r in rd(a.refmap)}
    depth = {r["sample"]: float(r["meandepth"]) for r in rd(a.cohort)}
    # h37rv_small, not matched_small: the matrix expresses divergence from H37Rv
    p2 = {}
    if os.path.exists(a.p2_summary):
        for r in rd(a.p2_summary):
            if (r.get("h37rv_small") or "").strip().isdigit():
                p2[r["sample"]] = int(r["h37rv_small"])

    with open(a.matrix, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        fixed = hdr.index("n_nocall") + 1
        # skip the artefact-flag columns: they sit between the counts and the
        # samples, and treating them as samples adds two pseudo-isolates
        while fixed < len(hdr) and hdr[fixed] in ("panel_af", "h37rv_minor"):
            fixed += 1
        samples = hdr[fixed:]
        ix = {k: hdr.index(k) for k in ("region", "kind", "h37rv_pos", "n_alt",
                                        "frame")}
        sites = []
        for row in rdr:
            sites.append(row)

    n = len(samples)
    findings = []

    def say(num, title):
        print(f"\n  [{num}] {title}")

    # --- 1: ubiquitous sites vs panel allele frequency -----------------------
    say(1, f"sites carried by all {n} samples")
    ubiq = [r for r in sites if int(r[ix["n_alt"]]) == n]
    ubpos = {int(r[ix["h37rv_pos"]]) for r in ubiq
             if r[ix["frame"]] == "h37rv" and r[ix["h37rv_pos"]]}
    print(f"      {len(ubiq)} sites, {len(ubpos)} with an H37Rv position")
    af = {}
    if ubpos and os.path.exists(a.graph_vcf):
        for line in gzip.open(a.graph_vcf, "rt"):
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            p = int(f[1])
            if p not in ubpos:
                continue
            for kv in f[7].split(";"):
                if kv.startswith("AF="):
                    try:
                        af[p] = max(float(x) for x in kv[3:].split(","))
                    except ValueError:
                        pass
    if af:
        hi = sum(1 for v in af.values() if v >= 0.5)
        print(f"      {len(af)} found in the graph VCF; panel allele frequency "
              f"median {statistics.median(af.values()):.3f}")
        print(f"      {hi} of {len(af)} have panel AF >= 0.5, i.e. H37Rv carries "
              f"the MINOR allele")
        print(f"      -> {100*hi/len(af):.0f}% are reference-artefact positions, "
              f"not shared variants")
        findings.append(dict(check="ubiquitous_sites", value=len(ubiq),
                             detail=f"{hi}/{len(af)} with panel AF>=0.5"))
    else:
        print(f"      no panel allele frequencies recovered; cannot classify")
        findings.append(dict(check="ubiquitous_sites", value=len(ubiq),
                             detail="AF unavailable"))

    # --- 2: singletons by region ---------------------------------------------
    say(2, "singleton sites by region, against all sites by region")
    allr = collections.Counter(r[ix["region"]] for r in sites)
    singr = collections.Counter(r[ix["region"]] for r in sites
                                if int(r[ix["n_alt"]]) == 1)
    print(f"      {'region':<22s}{'sites':>8s}{'singletons':>12s}{'rate':>8s}")
    for reg, tot in allr.most_common():
        s = singr[reg]
        print(f"      {reg:<22s}{tot:>8d}{s:>12d}{100*s/tot:>7.1f}%")
    findings.append(dict(check="singleton_rate",
                         value=sum(singr.values()),
                         detail=f"{100*sum(singr.values())/len(sites):.1f}% of sites"))

    # --- 3-7: per-sample --------------------------------------------------------
    st = {s: collections.Counter() for s in samples}
    kind = {s: collections.Counter() for s in samples}
    for r in sites:
        for s, v in zip(samples, r[fixed:]):
            st[s][v] += 1
            if v == "ALT":
                kind[s][r[ix["kind"]]] += 1

    say(3, "matrix ALT count against pass one's H37Rv burden (same frame)")
    print(f"      {'sample':<17s}{'matrix ALT':>11s}{'H37Rv':>10s}{'ratio':>8s}")
    ratios = []
    for s in samples:
        m, p = st[s]["ALT"], p2.get(s, 0)
        r_ = m / p if p else 0
        ratios.append(r_)
        print(f"      {s:<17s}{m:>11d}{p:>10d}{r_:>8.2f}")
    if ratios:
        med = statistics.median(ratios)
        print(f"      ratio median {med:.2f}")
        print(f"      Expected near 1: both count differences from H37Rv. Below 1")
        print(f"      means the routing dropped direct-arm calls in masked regions")
        print(f"      and replaced them with fewer composed ones; far from 1 in")
        print(f"      either direction means the merge lost or duplicated calls.")
        findings.append(dict(check="alt_vs_p2",
                             value=round(statistics.median(ratios), 3),
                             detail="median matrix/P2 ratio"))

    say(4, "state balance across samples (outliers indicate one sample behaving differently)")
    for stt in ("ALT", "REF", "ABSENT", "NOCALL"):
        vals = [st[s][stt] for s in samples]
        med = statistics.median(vals)
        out = [(s, st[s][stt]) for s in samples
               if med and abs(st[s][stt] - med) > 3 * med]
        print(f"      {stt:<8s} median {med:>8.0f}  range {min(vals)}-{max(vals)}"
              + (f"  OUTLIERS: {out}" if out else ""))
        findings.append(dict(check=f"state_{stt}", value=int(med),
                             detail=f"{min(vals)}-{max(vals)}"))

    say(5, "ABSENT against reference distance (should rise with distance)")
    pairs = [(_snp_dist(meta[s]), st[s]["ABSENT"]) for s in samples
             if s in meta]
    pairs.sort()
    for d, ab in pairs:
        print(f"      d={d:<6d} ABSENT {ab}")
    if len(pairs) > 2:
        xs = [p[0] for p in pairs]; ys = [p[1] for p in pairs]
        mx, my = statistics.mean(xs), statistics.mean(ys)
        num = sum((x - mx) * (y - my) for x, y in pairs)
        den = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
        r_ = num / den if den else 0
        print(f"      Pearson r = {r_:.3f}")
        findings.append(dict(check="absent_vs_distance", value=round(r_, 3),
                             detail="Pearson r"))

    say(6, "NOCALL against mean depth (should fall with depth)")
    pairs = [(depth.get(s, 0), st[s]["NOCALL"]) for s in samples if s in depth]
    if len(pairs) > 2:
        xs = [p[0] for p in pairs]; ys = [p[1] for p in pairs]
        mx, my = statistics.mean(xs), statistics.mean(ys)
        num = sum((x - mx) * (y - my) for x, y in pairs)
        den = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
        r_ = num / den if den else 0
        print(f"      depth {min(xs):.0f}-{max(xs):.0f}x, NOCALL "
              f"{min(ys)}-{max(ys)}, Pearson r = {r_:.3f}")
        findings.append(dict(check="nocall_vs_depth", value=round(r_, 3),
                             detail="Pearson r"))

    say(7, "SNP:INDEL ratio among ALT calls")
    rr = []
    for s in samples:
        sn, ind = kind[s]["SNP"], kind[s]["INDEL"]
        r_ = sn / ind if ind else 0
        rr.append(r_)
    print(f"      median {statistics.median(rr):.2f}, range "
          f"{min(rr):.2f}-{max(rr):.2f}")
    findings.append(dict(check="snp_indel_ratio",
                         value=round(statistics.median(rr), 2),
                         detail=f"{min(rr):.2f}-{max(rr):.2f}"))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["check", "value", "detail"],
                           delimiter="\t")
        w.writeheader(); w.writerows(findings)
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
