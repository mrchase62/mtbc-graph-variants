#!/usr/bin/env python3
"""Old vs new scale200: the run on the 2026-10-01 code (refbias/scale200,
assoc/scale200, build 7713a8d71d8e) against the rerun on all audit fixes
(refbias/scale200_fix, assoc/scale200_fix, build 7713a8d71d8e-fix1). Same
200 isolates, graph and phenotype file. Report only: reads existing outputs.

Run from runroot. Writes a TSV per section into --out and prints a summary.
"""
import argparse
import collections
import gzip
import os

import pandas as pd

OLD, NEW = "scale200", "scale200_fix"


def read_vcf(path):
    recs, samples = {}, None
    for line in gzip.open(path, "rt"):
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            samples = line.rstrip("\n").split("\t")[9:]
            continue
        c = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in c[7].split(";"))
        fmt = c[8].split(":")
        si = fmt.index("ST") if "ST" in fmt else None
        code = {"ALT": "A", "ABSENT": "B", "REF": "R", "NOCALL": "N"}
        st = "".join((code.get(x.split(":")[si], "?") if si is not None and len(x.split(":")) > si else "?")
                     for x in c[9:])
        recs[c[2]] = dict(pos=int(c[1]), cls=info.get("CLASS", "?"),
                          region=info.get("REGION", info.get("SVTIER", "")),
                          aa=info.get("AA", ""), aa_inv="AA_INVERTED" in info,
                          acclocus=info.get("ACCLOCUS", ""), st=st)
    return recs, samples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rep = []

    def say(s=""):
        print(s)
        rep.append(s)

    # 1. reference choice
    ro = pd.read_csv(f"refbias/{OLD}/p1/refmap.tsv", sep="\t")
    rn = pd.read_csv(f"refbias/{NEW}/p1/refmap.tsv", sep="\t")
    m = ro[["sample", "reference", "snp_distance"]].merge(
        rn[["sample", "reference", "snp_distance"]], on="sample", suffixes=("_old", "_new"))
    ch = m[m.reference_old != m.reference_new]
    ch.to_csv(f"{a.out}/refmap_changes.tsv", sep="\t", index=False)
    say(f"## 1. Matched reference (P1)\n")
    say(f"- {len(m)} samples; reference changed for {len(ch)}.")
    if len(ch):
        d = (ch.snp_distance_new - ch.snp_distance_old)
        say(f"- SNP distance to the reference for those: old median {ch.snp_distance_old.median():.0f}, "
            f"new median {ch.snp_distance_new.median():.0f}; new closer in {int((d < 0).sum())}, "
            f"farther in {int((d > 0).sum())}, equal in {int((d == 0).sum())}.")

    # 2. records and states
    vo, so = read_vcf(f"refbias/{OLD}/p5/merged.vcf.gz")
    vn, sn = read_vcf(f"refbias/{NEW}/p5/merged.vcf.gz")
    assert so == sn, "sample order differs"
    keys = sorted(set(vo) | set(vn))
    rows = []
    for k in keys:
        r = vn.get(k) or vo.get(k)
        rows.append(dict(id=k, cls=r["cls"], in_old=k in vo, in_new=k in vn))
    kt = pd.DataFrame(rows)
    say(f"\n## 2. Records in the merged VCF\n")
    say(f"- old {len(vo):,}, new {len(vn):,}; shared IDs {int((kt.in_old & kt.in_new).sum()):,}.\n")
    t = kt.groupby("cls").agg(old=("in_old", "sum"), new=("in_new", "sum"),
                              shared=("in_old", lambda s: int((s & kt.loc[s.index, "in_new"]).sum())))
    t["old_only"] = t.old - t.shared
    t["new_only"] = t.new - t.shared
    t.to_csv(f"{a.out}/records_by_class.tsv", sep="\t")
    say("| class | old | new | shared | old only | new only |\n|---|---:|---:|---:|---:|---:|")
    for c, r in t.iterrows():
        say(f"| {c} | {r.old:,} | {r.new:,} | {r.shared:,} | {r.old_only:,} | {r.new_only:,} |")

    # states per class, summed over samples
    say(f"\n## 3. Calls (cells = records x {len(so)} samples)\n")
    stc = []
    for tag, v in (("old", vo), ("new", vn)):
        cnt = collections.Counter()
        for r in v.values():
            for ch_ in r["st"]:
                cnt[(r["cls"], ch_)] += 1
        for (c, s), n in cnt.items():
            stc.append(dict(run=tag, cls=c, state={"A": "ALT", "B": "ABSENT", "R": "REF", "N": "NOCALL"}.get(s, s), n=n))
    sd = pd.DataFrame(stc)
    pv = sd.pivot_table(index=["cls", "state"], columns="run", values="n", aggfunc="sum").fillna(0)
    say("| class | state | old | new | change |\n|---|---|---:|---:|---:|")
    for (c, s), r in pv.iterrows():
        say(f"| {c} | {s} | {int(r.old):,} | {int(r.new):,} | {int(r.new - r.old):+,} |")
    say("")
    sd.to_csv(f"{a.out}/states_by_class.tsv", sep="\t", index=False)

    # cell agreement on shared records
    agree = collections.Counter()
    by_cls = collections.defaultdict(collections.Counter)
    for k in set(vo) & set(vn):
        x, y = vo[k]["st"], vn[k]["st"]
        for p, q in zip(x, y):
            agree[(p, q)] += 1
            by_cls[vn[k]["cls"]][(p, q)] += 1
    names = {"A": "ALT", "B": "ABSENT", "R": "REF", "N": "NOCALL"}
    say("Cell transitions on shared records:\n")
    tot = sum(agree.values())
    say("| old -> new | cells | % |\n|---|---:|---:|")
    for (p, q), n in agree.most_common(12):
        say(f"| {names.get(p, p)} -> {names.get(q, q)} | {n:,} | {100 * n / tot:.2f} |")
    pd.DataFrame([dict(cls=c, old=p, new=q, n=n) for c, cc in by_cls.items() for (p, q), n in cc.items()]
                 ).to_csv(f"{a.out}/cell_transitions_by_class.tsv", sep="\t", index=False)

    # 4. ancestral alleles
    say(f"\n## 4. Ancestral alleles (AA) on small records\n")
    so_ = [r for r in vo.values() if r["cls"] == "small"]
    sn_ = [r for r in vn.values() if r["cls"] == "small"]
    say(f"- with AA: old {sum(1 for r in so_ if r['aa']):,} of {len(so_):,}; "
        f"new {sum(1 for r in sn_ if r['aa']):,} of {len(sn_):,}.")
    say(f"- AA_INVERTED (reference carries the derived allele): old {sum(r['aa_inv'] for r in so_):,}, "
        f"new {sum(r['aa_inv'] for r in sn_):,}.")
    shared_small = [k for k in set(vo) & set(vn) if vo[k]["cls"] == "small"]
    diff = [k for k in shared_small if vo[k]["aa"] != vn[k]["aa"]]
    gained = sum(1 for k in diff if not vo[k]["aa"] and vn[k]["aa"])
    lost = sum(1 for k in diff if vo[k]["aa"] and not vn[k]["aa"])
    say(f"- on {len(shared_small):,} shared small records: AA differs on {len(diff):,} "
        f"(gained {gained:,}, lost {lost:,}, changed base {len(diff) - gained - lost:,}).")

    # 5. association
    say(f"\n## 5. Association (same phenotype file)\n")
    qo = pd.read_csv(f"assoc/{OLD}/scan.tsv", sep="\t")
    qn = pd.read_csv(f"assoc/{NEW}/scan.tsv", sep="\t")
    for tag, q in (("old", qo), ("new", qn)):
        surv = (q.q_branch < 0.05) & (q.q_region < 0.05) & (q.q_lineage < 0.05)
        say(f"- {tag}: tested {len(q):,}, q_branch < 0.05 {int((q.q_branch < 0.05).sum())}, "
            f"survivors {int(surv.sum())}; min p_branch {q.p_branch.min():.2g}.")
    j = qo[["id", "gains", "carriers", "p_branch", "q_branch", "region"]].merge(
        qn[["id", "gains", "carriers", "p_branch", "q_branch"]], on="id", how="outer",
        suffixes=("_old", "_new"), indicator=True)
    # IDs that change between runs for reasons other than a change in the
    # event: SV interval IDs (svi:TYPE:start:len) follow the rebuilt interval
    # catalogue, and indel IDs move with left-alignment (D22). An old-only
    # row is matched to a new-only row of the same kind within 20 bp whose
    # length agrees within 10% (indels: inserted/deleted length).
    def parse(i):
        p = i.split(":")
        if p[0] in ("svi", "sv"):
            return (p[0] + p[1], int(p[2]), abs(int(p[3])))
        if p[0] == "h37rv" and ">" in p[-1]:
            r, al = p[-1].split(">")
            return ("indel" if len(r) != len(al) else "snp", int(p[1]), len(al) - len(r))
        return None
    lo = j[j._merge == "left_only"]
    rn_ = j[j._merge == "right_only"]
    cand = collections.defaultdict(list)
    for _, r in rn_.iterrows():
        k = parse(r.id)
        if k:
            cand[k[0]].append((k[1], k[2], r))
    rematch = []
    for _, r in lo.iterrows():
        k = parse(r.id)
        if not k:
            continue
        best = None
        for pos, ln, rr in cand.get(k[0], []):
            if abs(pos - k[1]) <= 20 and abs(ln - k[2]) <= max(2, 0.1 * abs(k[2])):
                if best is None or abs(pos - k[1]) < abs(best[0] - k[1]):
                    best = (pos, rr)
        if best is not None:
            rr = best[1]
            rematch.append(dict(id_old=r.id, id_new=rr.id, gains_old=r.gains_old, gains_new=rr.gains_new,
                                q_branch_old=r.q_branch_old, q_branch_new=rr.q_branch_new))
    rm = pd.DataFrame(rematch)
    rm.to_csv(f"{a.out}/scan_rematched.tsv", sep="\t", index=False)
    if len(rm):
        say(f"- old-only IDs rematched by coordinates to a new ID: {len(rm)} "
            f"(of {len(lo)} old-only); of the rematched, q_branch < 0.05 old {int((rm.q_branch_old < 0.05).sum())}, "
            f"new {int((rm.q_branch_new < 0.05).sum())}, both {int(((rm.q_branch_old < 0.05) & (rm.q_branch_new < 0.05)).sum())}.")
    j.to_csv(f"{a.out}/scan_joined.tsv", sep="\t", index=False)
    po = j.q_branch_old < 0.05
    pn = j.q_branch_new < 0.05
    say(f"- q_branch passes: both {int((po & pn).sum())}, old only {int((po & ~pn).sum())}, "
        f"new only {int((~po & pn).sum())}.")
    lost = j[po & ~pn]
    why = collections.Counter()
    for _, r in lost.iterrows():
        if r._merge == "left_only":
            why["not tested in new"] += 1
        elif r.gains_new < r.gains_old:
            why["fewer gains in new"] += 1
        elif r.gains_new > r.gains_old:
            why["more gains in new"] += 1
        else:
            why["same gains, p changed"] += 1
    say(f"- the {len(lost)} old-only passes: " + ", ".join(f"{k} {v}" for k, v in why.most_common()))
    lost.to_csv(f"{a.out}/qbranch_lost.tsv", sep="\t", index=False)
    for b in ("small_gene", "sv_gene", "is6110_gene"):
        go = pd.read_csv(f"assoc/{OLD}/{b}.tsv", sep="\t")
        gn = pd.read_csv(f"assoc/{NEW}/{b}.tsv", sep="\t")
        qc = [x for x in gn.columns if x.startswith("q_")]
        so2 = set(go[(go[qc] < 0.05).all(axis=1)].iloc[:, 0])
        sn2 = set(gn[(gn[qc] < 0.05).all(axis=1)].iloc[:, 0])
        say(f"- {b}: units old {len(go):,} new {len(gn):,}; survivors old {sorted(so2)} new {sorted(sn2)}")
    open(f"{a.out}/summary.md", "w").write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
