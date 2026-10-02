#!/usr/bin/env python3
"""Reconcile every record class and key space across the chain, and fail loudly.

WHY THIS EXISTS. Three times in one week a downstream stage silently ignored
something an upstream stage had added, and each cost a whole arm with no error
raised anywhere:

  the `svi:` key space   bin/sv_intervals.py wrote 3,658 catalogued deletion
                         intervals keyed `svi:`, and merge_cohort_vcf.py read
                         `sv:` on caller clusters. Zero keys in common, so every
                         improvement to the interval arm landed in a file the
                         VCF, the event matrix and the scan all ignored.
  the accessory class    the 802 CLASS=accessory_presence records reached the
                         VCF and the event matrix, and assoc_scan.py's
                         --classes default listed small,sv,is6110. All 802 were
                         invisible -- not skipped with a count, absent. 92 had
                         two or more independent origins and 65 passed the
                         callability floor.
  the callability floor  a conditional variant's branches outside its carrier
                         clades are undetermined by construction, so an 80%
                         whole-tree floor dropped every accessory locus carried
                         by less than about two thirds of the cohort. Six
                         testable variants went with it, including the one
                         inside TbD1.

None of these is a crash. Each is a number quietly becoming zero. The shape they
share is that a class or a stratum present upstream has NO representative
downstream, while the totals still look plausible because the big classes
dominate them.

WHAT IT CHECKS. At each boundary it cross-tabulates records by (class, region)
and compares the strata present on each side:

  FAIL   a stratum with records upstream and ZERO downstream. That is the
         signature of all three bugs above.
  WARN   a stratum whose drop rate is more than --warn-ratio times the global
         drop rate at that boundary, with at least --warn-min records. A
         stratum can legitimately fare worse than average -- the callability
         floor is meant to bite hardest where data is thinnest -- so this is a
         prompt to look, not a failure.
  note   a stratum that is entirely absent from BOTH sides, which usually means
         a cohort simply has none of that kind.

It is deliberately not a test of correctness. It cannot tell whether a filter is
right, only whether a filter has quietly consumed an entire stratum. That is the
class of mistake that actually happened.
"""
import argparse, collections, csv, gzip, os, re, sys


def vcf_strata(path):
    """(class, region) -> count, read from the merged VCF's INFO."""
    c = collections.Counter()
    ids = collections.Counter()
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 8)
            info = f[7]
            m = re.search(r"CLASS=([A-Za-z_0-9]+)", info)
            cls = m.group(1) if m else "?"
            m = re.search(r"REGION=([A-Za-z_0-9]+)", info)
            reg = m.group(1) if m else ""
            c[(cls, reg)] += 1
            pre = f[2].split(":", 1)[0] if ":" in f[2] else "<no prefix>"
            ids[pre] += 1
    return c, ids


def tsv_strata(path, cls_col, reg_col, id_col=None):
    c = collections.Counter()
    ids = collections.Counter()
    if not os.path.exists(path):
        return None, None
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c[(r.get(cls_col, "?"), r.get(reg_col, "") or "")] += 1
            if id_col and r.get(id_col):
                v = r[id_col]
                ids[v.split(":", 1)[0] if ":" in v else "<no prefix>"] += 1
    return c, ids


def stratum(cls, reg):
    """The stratum identity, mirroring assoc_scan.py's own regkey().

    Two transformations happen downstream and both must be reproduced here or
    every affected stratum reads as a total loss plus a new arrival. The scan
    SUBSTITUTES a region for records that carry none -- `accessory_presence`
    for that class, `is6110` for element records, `other` for everything else
    -- and it REFINES the region it was given, `core` becoming `core:genic`.
    The first version of this audit reproduced neither and reported six
    failures that were all the same non-event.
    """
    base = coarse(reg)
    # assoc_scan.py keys a catalogued deletion by its EVIDENCE TIER --
    # `svi:E1_coherent` and so on -- where the event matrix has no region for
    # it at all and the substitution would have been `other`. Mapping it back
    # keeps the two sides comparable; without this every evidence tier reads as
    # a new stratum and `sv/other` as a total loss.
    if base == "svi":
        return (cls, "other")
    if not base:
        if cls == "accessory_presence":
            return (cls, "accessory_presence")
        return (cls, "is6110" if cls == "is6110" else "other")
    return (cls, base)


def coarse(reg):
    """The region label without the refinements a later stage adds.

    assoc_scan.py REFINES the region it was given -- `core` becomes
    `core:genic` or `core:intergenic`, and a conditional variant becomes
    `off_path_accessory:cond_ge200`. Comparing the raw strings across that
    boundary reports every refined stratum as both a total loss and a new
    arrival, which is how the first version of this audit produced seven
    failures that were all the same non-event. Identity is the part before the
    colon.
    """
    return (reg or "").split(":")[0]


def fold(c):
    out = collections.Counter()
    for (cls, reg), n in c.items():
        out[stratum(cls, reg)] += n
    return out


def reconcile(name, up_rows, down, floor, min_gains, fails, warns):
    """The events -> scan boundary, by arithmetic rather than by inspection.

    Every record the scan did not test must be accounted for by one of its two
    documented filters: fewer than --min-gains independent origins, or
    determinacy below --min-determinacy. If those two do not add up to the
    observed loss, something else consumed the records and that is the bug this
    audit exists to catch.
    """
    up = collections.Counter()
    nogain = collections.Counter()
    uncall = collections.Counter()
    for r in up_rows:
        k = stratum(r["class"], r.get("region", ""))
        up[k] += 1
        if int(r["n_gain"]) < min_gains:
            nogain[k] += 1
        elif r["_det"] < floor:
            uncall[k] += 1
    lines = [f"\n  {name}",
             f"    every record not tested must be explained by "
             f"< {min_gains} origins or determinacy < {floor:.0%}",
             f"    {'class':<20}{'region':<22}{'up':>8}{'no orig':>9}"
             f"{'uncall':>8}{'expect':>8}{'scanned':>9}"]
    for k in sorted(set(up) | set(down)):
        u = up.get(k, 0)
        d = down.get(k, 0)
        exp = u - nogain.get(k, 0) - uncall.get(k, 0)
        tag = ""
        if exp != d:
            tag = f"  <== FAIL: {exp - d:+,} unexplained"
            fails.append(f"{name}: {k[0]}/{k[1] or '-'} -- {u:,} records, "
                         f"{nogain.get(k,0):,} without independent origins, "
                         f"{uncall.get(k,0):,} uncallable, so {exp:,} should "
                         f"have been scanned and {d:,} were")
        elif u > 0 and d == 0:
            why = ("no record has two independent origins"
                   if nogain.get(k, 0) == u else
                   "every testable record is uncallable")
            tag = f"  (all accounted for: {why})"
        lines.append(f"    {k[0]:<20}{(k[1] or '-'):<22}{u:>8,}"
                     f"{nogain.get(k,0):>9,}{uncall.get(k,0):>8,}"
                     f"{exp:>8,}{d:>9,}{tag}")
    return lines


def compare(name, up, down, warn_ratio, warn_min, fails, warns):
    """One boundary. Returns the lines to print."""
    tot_up, tot_down = sum(up.values()), sum(down.values())
    rate = 1.0 - (tot_down / tot_up) if tot_up else 0.0
    lines = [f"\n  {name}",
             f"    {tot_up:,} upstream -> {tot_down:,} downstream "
             f"({rate:.1%} dropped overall)"]
    keys = sorted(set(up) | set(down))
    lines.append(f"    {'class':<22}{'region':<26}{'up':>8}{'down':>8}"
                 f"{'dropped':>9}")
    for k in keys:
        u, d = up.get(k, 0), down.get(k, 0)
        if u == 0 and d == 0:
            continue
        dr = 1.0 - (d / u) if u else 0.0
        tag = ""
        if u > 0 and d == 0:
            tag = "  <== FAIL: whole stratum absent downstream"
            fails.append(f"{name}: {k[0]}/{k[1] or '-'} has {u:,} records "
                         f"upstream and none downstream")
        elif u >= warn_min and rate > 0 and dr > warn_ratio * rate:
            tag = f"  <-- WARN: {dr:.0%} dropped vs {rate:.0%} overall"
            warns.append(f"{name}: {k[0]}/{k[1] or '-'} drops {dr:.0%} "
                         f"against a {rate:.0%} baseline ({u:,} records)")
        elif d > u:
            tag = "  (grew: records split or renamed here)"
        lines.append(f"    {k[0]:<22}{(k[1] or '-'):<26}{u:>8,}{d:>8,}"
                     f"{dr:>8.0%}{tag}")
    return lines


def compare_ids(name, up, down, fails, expect=None):
    """Key spaces. `expect[k]` is how many of k's records SHOULD survive, when
    the caller can compute it; without it a key space is only checked for total
    disappearance, which the `svi:` bug would have caught but which also fires
    on a key space whose records are all legitimately untestable."""
    lines = [f"\n  {name} -- ID key spaces"]
    for k in sorted(set(up) | set(down)):
        u, d = up.get(k, 0), down.get(k, 0)
        tag = ""
        # `expect` is authoritative for every key space seen upstream. Testing
        # `k in expect` instead read a legitimately-zero expectation as absent
        # and fell through to the cruder check, which then failed on `sv:` --
        # whose 7,313 presence-only insertion records genuinely have nothing
        # testable and callable in them.
        if expect is not None:
            e = expect.get(k, 0)
            if e != d:
                tag = f"  <== FAIL: {e - d:+,} unexplained (expected {e:,})"
                fails.append(f"{name}: key space `{k}:` should have carried "
                             f"{e:,} records downstream and carried {d:,}")
            elif u > 0 and d == 0:
                tag = "  (all accounted for: none testable and callable)"
        elif u > 0 and d == 0:
            tag = "  <== FAIL: key space never read downstream"
            fails.append(f"{name}: key space `{k}:` has {u:,} records "
                         f"upstream and none downstream")
        lines.append(f"    {k + ':':<22}{u:>8,}{d:>8,}{tag}")
    return lines


def load_carriers(d, events, assoc_bin):
    """locus -> number of applicable branches, exactly as the scan computes it.

    A node is applicable when at least one of its descendant leaves carries the
    locus; the count excludes the root. It is read off the event matrix's own
    tree with the scan's own tree reader and descendant matrix, so the audit
    and the scan cannot disagree about it. Returns (counts, branches in the
    tree).

    It used to approximate the count by twice the carrier count minus one, the
    branch count of a binary subtree on that many leaves. That leaves out the
    path from the carriers' subtree to the root, so it undercounts most for
    loci with few carriers: on gwas1000, ACC_2163649 (83 carriers) came out at
    165 branches against the tree's 244, and the audit expected a variant the
    scan had correctly dropped at 0.64 determinacy.
    """
    import glob as _g
    import numpy as np
    sys.path.insert(0, os.path.abspath(assoc_bin))
    from sv_scatter import tree_frame
    from assoc_scan import descendant_matrix
    T = tree_frame(os.path.join(events, "labelled.nwk"))
    li = {l: j for j, l in enumerate(T["leaves"])}
    D = descendant_matrix(T)
    per = collections.defaultdict(set)
    for f in sorted(_g.glob(os.path.join(d, "*.presence.tsv"))):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("state") == "PRESENT":
                    per[r["locus"]].add(r["sample"])
    out = {}
    for lid, ss in per.items():
        A = np.zeros(len(T["leaves"]), dtype=bool)
        for s2 in ss:
            if s2 in li:
                A[li[s2]] = True
        if not A.any():
            continue
        nba = int(((D & A).sum(axis=1) > 0).sum()) - 1
        if nba > 0:
            out[lid] = nba
    return out, T["n"] - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--vcf", default="")
    ap.add_argument("--events", default="")
    ap.add_argument("--scan", default="")
    ap.add_argument("--min-gains", type=int, default=2,
                    help="must match assoc_scan.py")
    ap.add_argument("--min-determinacy", type=float, default=0.8,
                    help="must match assoc_scan.py")
    ap.add_argument("--branches", type=int, default=0,
                    help="branches in the tree; inferred when omitted")
    ap.add_argument("--accessory-presence", default="",
                    help="needed to reproduce the scan's conditional "
                         "determinacy correction; without it a conditional "
                         "stratum will read as unexplained")
    ap.add_argument("--assoc-bin", default="assoc/bin",
                    help="where assoc_scan.py and sv_scatter.py live; the "
                         "applicable-branch count uses their tree code")
    ap.add_argument("--warn-ratio", type=float, default=2.0)
    ap.add_argument("--warn-min", type=int, default=10)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    C = a.cohort
    vcf = a.vcf or f"refbias/{C}/p5/merged.vcf.gz"
    events = a.events or f"assoc/{C}/events"
    scan = a.scan or f"assoc/{C}/scan.tsv"

    out = [f"=== chain audit, cohort {C}"]
    fails, warns = [], []

    if not os.path.exists(vcf):
        sys.exit(f"FATAL: no VCF at {vcf}")
    v_str, v_ids = vcf_strata(vcf)
    out.append(f"  VCF     {vcf}  {sum(v_str.values()):,} records, "
               f"{len(v_str)} (class, region) strata")

    e_str, e_ids = tsv_strata(os.path.join(events, "variants.tsv"),
                              "class", "region", "id")
    if e_str is None:
        out.append(f"  no event matrix at {events}; stopping after the VCF")
    else:
        out += compare("VCF -> event matrix", fold(v_str), fold(e_str),
                       a.warn_ratio, a.warn_min, fails, warns)
        out += compare_ids("VCF -> event matrix", v_ids, e_ids, fails)

        s_str, s_ids = tsv_strata(scan, "cls", "region", "id")
        if s_str is None:
            out.append(f"\n  no scan at {scan}; stopping after the event matrix")
        else:
            # the scan's own two filters, recomputed here so the audit does not
            # take the scan's word for what it dropped
            ev_rows = list(csv.DictReader(
                open(os.path.join(events, "variants.tsv")), delimiter="\t"))
            carriers, tree_branches = load_carriers(
                a.accessory_presence, events, a.assoc_bin) \
                if a.accessory_presence else ({}, 0)
            nbranch = a.branches or tree_branches
            if not nbranch:
                nbranch = max((int(r["n_gain"]) + int(r["n_loss"])
                               + int(r["n_undet"]) for r in ev_rows),
                              default=1)
            for r in ev_rows:
                und = int(r["n_undet"])
                nb = nbranch
                lid = r.get("acc_locus") or ""
                if lid and lid in carriers:
                    # the same correction the scan applies: branches outside
                    # the carrier clades are inapplicable, not undetermined
                    nb = carriers[lid]
                    und = max(0, und - (nbranch - nb))
                r["_det"] = 1.0 - und / max(1, nb)
            out += reconcile("event matrix -> scan", ev_rows, fold(s_str),
                             a.min_determinacy, a.min_gains, fails, warns)
            exp_ids = collections.Counter()
            for r in ev_rows:
                if int(r["n_gain"]) >= a.min_gains and r["_det"] >= a.min_determinacy:
                    v = r.get("id", "")
                    exp_ids[v.split(":", 1)[0] if ":" in v else "<no prefix>"] += 1
            out += compare_ids("event matrix -> scan", e_ids, s_ids, fails,
                               expect=exp_ids)

    out.append("")
    if fails:
        out.append(f"  {len(fails)} FAILURE(S) -- a stratum or key space present "
                   f"upstream has no representative downstream:")
        out += [f"    - {x}" for x in fails]
    if warns:
        out.append(f"  {len(warns)} warning(s) -- worth a look, not a failure:")
        out += [f"    - {x}" for x in warns]
    if not fails and not warns:
        out.append("  every class, region and key space is represented at every "
                   "boundary.")
    elif not fails:
        out.append("  no stratum was lost entirely.")

    text = "\n".join(out)
    print(text)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text + "\n")
        print(f"\n  -> {a.out}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
