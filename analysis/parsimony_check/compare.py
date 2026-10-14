#!/usr/bin/env python3
"""Validate bin/panel_ancestral.py on the CX333 build (7713a8d71d8e-fix1)
against the two tables it replaces. Analysis only; writes under
analysis/parsimony_check/out/ and a scratch directory, never into the build.
Run from the repository root with $MTB_PY after sourcing
config/project_env.sh:

    $MTB_PY analysis/parsimony_check/compare.py --work <scratch dir>

Four runs of the new script:
  A  full CX333 tree and VCF, outgroups = the two canettii: the same ancestral
     node production's AA used (the MTBC ancestor). Compared with the build's
     assets/ancestral.tsv (SNP AA) and assets/panel_polarity.tsv (canettii
     ET1291's own genotype, every class).
  C  full CX333, ancestral node = MRCA of the lineage 1-4,7 genomes (every
     other leaf an outgroup): the decided node, with the whole panel.
  B  the P200 design emulated on CX333: tree and VCF pruned to the lineage
     1-4,7 genomes + L5.1.1 GCF_022870265 + L6.1.1 GCF_000253355 (in graph)
     + L8 GCF_012923765 and canettii GCF_035581225 as KIT rows (SNP states
     from analysis/panel_checks/states/*.npz; no indel/SV rows). Compared
     with C (does the reduced panel give the full panel's answer) and with
     what the old code would give on that panel: AA by nearest-first, and a
     single-outgroup read of MTB_OUTGROUP = GCF_022870265.
  B2 as B with the outgroup list reversed: must be byte-identical.
"""
import argparse
import collections
import csv
import filecmp
import os
import subprocess
import sys

sys.path.insert(0, "bin")
import panel_ancestral as P  # noqa: E402

RUN = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/runroot"
B = f"{RUN}/refbias/build/7713a8d71d8e-fix1"
STATES = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/panel_checks/states"
QC = "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/qc/qc_master_table.tsv"
OUT = "analysis/parsimony_check/out"
KITEV = "analysis/outgroup_kit/out"     # route A, bin/outgroup_kit.py
H37 = "GCF_000195955"
L56 = ["GCF_022870265", "GCF_000253355"]
KIT = ["GCF_012923765", "GCF_035581225"]
ELEMENT = 2352065
ELEMENT_ALT = "CCTCGCCTGGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTG"


def prune(children, root, keep):
    out = {}

    def rec(n):
        if n not in children:
            return n if n in keep else None
        kids = [k for k in (rec(c) for c in children[n]) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        out[n] = kids
        return n
    return out, rec(root)


def newick(children, root):
    def rec(n):
        return n if n not in children else "(" + ",".join(rec(c) for c in children[n]) + ")"
    return rec(root) + ";\n"


def rows(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def key(r):
    return (r["chrom"], int(r["pos"]), r["ref"], r["alt"])


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FATAL: {' '.join(cmd[:3])}...\n{r.stdout}\n{r.stderr}")
    return r.stdout


def panel_run(tag, tree, vcf, outgroups, kit, work, events=False):
    d = os.path.join(work, tag)
    os.makedirs(d, exist_ok=True)
    cmd = [sys.executable, "bin/panel_ancestral.py", "--tree", tree, "--panel-vcf", vcf,
           "--h37rv", f"{B}/refs/{H37}.fasta", "--ref-sample", H37,
           "--outgroups", ",".join(outgroups), "--bcftools", os.environ["MTB_BCFTOOLS"],
           "--out-ancestral", f"{d}/ancestral.tsv", "--out-polarity", f"{d}/panel_polarity.tsv"]
    for k in kit:
        cmd += ["--kit", f"{k}={STATES}/{k}.npz"]
        if events:
            cmd += ["--kit-events", f"{k}={KITEV}/{k}"]
    log = run(cmd)
    with open(f"{d}/log.txt", "w") as fh:
        fh.write(log)
    return d, log


def kind(ref, alt):
    return P.kind_of(ref, alt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True, help="scratch directory for the runs")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(a.work, exist_ok=True)
    bcf = os.environ["MTB_BCFTOOLS"]
    rep = []

    def say(s=""):
        print(s)
        rep.append(s)

    children, root, leaves = P.parse_newick(open(f"{B}/assets/ancestral_inputs/tree.nwk").read())
    lin = {r["accession"]: r["lineage"] for r in csv.DictReader(open(QC), delimiter="\t")
           if r["in_panel"] == "yes"}
    target = {g for g in leaves if lin.get(g) in
              ("lineage1", "lineage2", "lineage3", "lineage4", "lineage7")}
    target.add(H37)
    others = sorted(leaves - target)
    say(f"CX333 tree: {len(leaves)} leaves; lineage 1-4,7 (with H37Rv): {len(target)}; "
        f"outside: {len(others)}")

    tree_full = f"{B}/assets/ancestral_inputs/tree.nwk"
    vcf_full = f"{B}/assets/graph_collapsed.vcf.gz"
    dA, logA = panel_run("A", tree_full, vcf_full, ["GCF_035581225", "GCF_000253375"], [], a.work)
    dC, logC = panel_run("C", tree_full, vcf_full, others, [], a.work)

    keep = target | set(L56) | set(KIT)
    ch, rt = prune(children, root, keep)
    tree_b = os.path.join(a.work, "p200_emul.nwk")
    open(tree_b, "w").write(newick(ch, rt))
    samples = [s for s in sorted(keep - set(KIT)) if s != H37]
    vcf_b = os.path.join(a.work, "p200_emul.vcf.gz")
    if not os.path.exists(vcf_b):
        run([bcf, "view", "-s", ",".join(samples), "-Oz", "-o", vcf_b, vcf_full])
        run([bcf, "index", "-t", "-f", vcf_b])
    ogs = L56 + KIT
    dB, logB = panel_run("B", tree_b, vcf_b, ogs, KIT, a.work)
    dB2, _ = panel_run("B2", tree_b, vcf_b, list(reversed(ogs)), KIT, a.work)
    same = all(filecmp.cmp(f"{dB}/{f}", f"{dB2}/{f}", shallow=False)
               for f in ("ancestral.tsv", "panel_polarity.tsv"))
    say(f"B vs B2 (outgroup list reversed): outputs identical: {same}")

    # ---------------- A vs the build's tables
    old_aa = {key(r): r for r in rows(f"{B}/assets/ancestral.tsv")}
    old_pp = {key(r): r for r in rows(f"{B}/assets/panel_polarity.tsv")}
    newA_aa = {key(r): r for r in rows(f"{dA}/ancestral.tsv")}
    newA_pp = {key(r): r for r in rows(f"{dA}/panel_polarity.tsv")}

    say("\n== A (same node as production AA: MTBC ancestor, canettii outgroups) vs build assets")
    c = collections.Counter()
    for k, r in old_aa.items():
        n = newA_aa.get(k)
        if n is None:
            c["old site, no new row"] += 1
        elif r["AA"] == n["AA"]:
            c["same AA" + (" (both unresolved)" if r["AA"] == "." else "")] += 1
        elif r["AA"] == ".":
            c["old unresolved -> new resolved"] += 1
        elif n["AA"] == ".":
            c["old resolved -> new unresolved"] += 1
        else:
            c["AA differs"] += 1
    c["new SNP keys not in old table"] = len(set(newA_aa) - set(old_aa))
    c["  of which resolved"] = sum(1 for k in set(newA_aa) - set(old_aa) if newA_aa[k]["AA"] != ".")
    say(f"SNP AA: old table {len(old_aa):,} sites, new {len(newA_aa):,} SNP keys")
    for k_, v in c.items():
        say(f"  {k_:<40}{v:>9,}")
    diffs = [(k, old_aa[k]["AA"], newA_aa[k]["AA"], newA_aa[k]["flag"]) for k in old_aa
             if k in newA_aa and old_aa[k]["AA"] != newA_aa[k]["AA"]]

    # effective direction the event writer used: AA where resolved, else the
    # polarity table's (canettii's own read)
    def eff_old(k):
        r = old_aa.get(k)
        if r and r["AA"] != ".":
            return "ALT" if r["AA"] == k[3] else "REF"
        p = old_pp.get(k)
        return p["ancestral"] if p else None

    def compare(new_pp, eff, label):
        t = collections.defaultdict(collections.Counter)
        allk = set(new_pp) | {k for k in old_pp} | {k for k in old_aa}
        for k in allk:
            o = eff(k)
            n = new_pp[k]["ancestral"] if k in new_pp else None
            cl = kind(k[2], k[3])
            if o is None and n is None:
                t[cl]["both unpolarised"] += 1
            elif o is None:
                t[cl]["gains a direction"] += 1
            elif n is None:
                t[cl]["loses it (new: tied)"] += 1
            elif o == n:
                t[cl]["same"] += 1
            else:
                t[cl][f"flips {o}->{n}"] += 1
        cols = ["same", "gains a direction", "flips REF->ALT", "flips ALT->REF",
                "loses it (new: tied)"]
        say(f"{label}\n  {'class':<7}" + "".join(f"{x:>22}" for x in cols))
        for cl in ("SNP", "MNP", "indel", "SV"):
            say(f"  {cl:<7}" + "".join(f"{t[cl][x]:>22,}" for x in cols))
        return t

    say("")
    compare(newA_pp, eff_old, "every key, direction the event writer would read "
            "(old: AA else canettii's genotype; new: run A)")

    # what the event writer finally uses as the DERIVED allele: AA where it
    # resolved; else the polarity table, an ALT-ancestral row only at
    # panel_af >= 5%; else ALT (unpolarised default)
    def derived(aa_tab, pp, k):
        r = aa_tab.get(k)
        if r and r["AA"] != ".":
            return "REF" if r["AA"] == k[3] else "ALT", "AA"
        p = pp.get(k)
        if p is None:
            return "ALT", "default"
        if p["ancestral"] == "ALT":
            if float(p["panel_af"]) >= 0.05:
                return "REF", "table"
            return "ALT", "default"
        return "ALT", "table"

    def derived_cmp(aa_new, pp_new, label, aa_old=None, pp_old=None):
        aa_old = old_aa if aa_old is None else aa_old
        pp_old = old_pp if pp_old is None else pp_old
        t = collections.defaultdict(collections.Counter)
        for k in set(pp_new) | set(pp_old) | set(aa_old) | set(aa_new):
            o, os_ = derived(aa_old, pp_old, k)
            n, ns = derived(aa_new, pp_new, k)
            cl = kind(k[2], k[3])
            if o != n:
                t[cl][f"derived {o}->{n}"] += 1
            elif os_ == "default" and ns != "default":
                t[cl]["now measured, same"] += 1
            elif ns == "default" and os_ != "default":
                t[cl]["now default, same"] += 1
            else:
                t[cl]["unchanged"] += 1
        cols = ["unchanged", "now measured, same", "now default, same",
                "derived ALT->REF", "derived REF->ALT"]
        say(f"{label}\n  {'class':<7}" + "".join(f"{x:>22}" for x in cols))
        for cl in ("SNP", "MNP", "indel", "SV"):
            say(f"  {cl:<7}" + "".join(f"{t[cl][x]:>22,}" for x in cols))
    self_derived = derived_cmp
    say("")
    self_derived(newA_aa, newA_pp, "derived allele the event writer ends with "
                 "(AA first, then the table with the 5% guard, else ALT): build vs A")

    # ---------------- C and B
    newC_pp = {key(r): r for r in rows(f"{dC}/panel_polarity.tsv")}
    newB_pp = {key(r): r for r in rows(f"{dB}/panel_polarity.tsv")}
    newB_aa = {key(r): r for r in rows(f"{dB}/ancestral.tsv")}

    def direction_table(pp):
        return lambda k: pp[k]["ancestral"] if k in pp else None

    newC_aa = {key(r): r for r in rows(f"{dC}/ancestral.tsv")}
    say("\n== C (full CX333, node = MRCA of L1-4,7) vs build assets")
    compare(newC_pp, eff_old, "every key: old effective direction vs C")
    self_derived(newC_aa, newC_pp, "derived allele the event writer ends with: build vs C")
    say("\n== B (P200 design emulated: L1-4,7 + L5.1 + L6.1, kit L8 + ET1291 SNP rows)")
    say("\n".join("  " + l for l in logB.splitlines() if l.strip()))
    # B vs C on the keys B has
    t = collections.defaultdict(collections.Counter)
    for k in set(newB_pp) | set(newC_pp):
        b, cc = newB_pp.get(k), newC_pp.get(k)
        cl = kind(k[2], k[3])
        if b and cc:
            t[cl]["agree" if b["ancestral"] == cc["ancestral"] else "differ"] += 1
        elif b:
            t[cl]["B only (C tied)"] += 1
        elif cc:
            t[cl]["C only (B tied or key absent from B)"] += 1
    cols = ["agree", "differ", "B only (C tied)", "C only (B tied or key absent from B)"]
    say("\nB vs C (reduced panel vs full panel, same node)\n  " + f"{'class':<7}"
        + "".join(f"{x:>38}" for x in cols))
    for cl in ("SNP", "MNP", "indel", "SV"):
        say(f"  {cl:<7}" + "".join(f"{t[cl][x]:>38,}" for x in cols))

    # old code on the same reduced panel: single-outgroup read of L5.1
    og = "GCF_022870265"
    smp = run([bcf, "query", "-l", vcf_b]).split()
    oi = smp.index(og)
    seq = "".join(l.strip() for l in open(f"{B}/refs/{H37}.fasta") if not l.startswith(">")).upper()
    _, keys, copies, _ = P.read_panel(vcf_b, bcf, seq, "NC_000962.3")
    l5 = {}
    for k, cp in zip(keys, copies):
        g = {c[oi] for c in cp}
        if "1" in g:
            l5[k] = "ALT"
        elif "0" in g:
            l5[k] = "REF"
    say("")
    compare_tab = collections.defaultdict(collections.Counter)
    for k, cp in zip(keys, copies):
        o, n = l5.get(k), (newB_pp[k]["ancestral"] if k in newB_pp else None)
        cl = kind(k[2], k[3])
        if o is None and n is None:
            compare_tab[cl]["both unpolarised"] += 1
        elif o is None:
            compare_tab[cl]["gains a direction"] += 1
        elif n is None:
            compare_tab[cl]["loses it (new: tied)"] += 1
        elif o == n:
            compare_tab[cl]["same"] += 1
        else:
            compare_tab[cl][f"flips {o}->{n}"] += 1
    cols = ["same", "gains a direction", "flips REF->ALT", "flips ALT->REF",
            "loses it (new: tied)", "both unpolarised"]
    say("B vs the old single-outgroup read of MTB_OUTGROUP=GCF_022870265 on the same panel "
        "(what panel_polarity.py would have written)\n  " + f"{'class':<7}"
        + "".join(f"{x:>22}" for x in cols))
    for cl in ("SNP", "MNP", "indel", "SV"):
        say(f"  {cl:<7}" + "".join(f"{compare_tab[cl][x]:>22,}" for x in cols))
    # the old code's full answer on this panel: AA by ancestral_alleles.py's
    # nearest-first rule is not rerun here (it needs the alignment); the old
    # table is the L5.1 read with the panel AF over this VCF
    old_b_pp = {}
    for k, cp in zip(keys, copies):
        if k in l5:
            col = P.pooled_column(cp)
            called = sum(x not in P.MISSING for x in col)
            af = sum(x == "1" for x in col) / called if called else 0.0
            old_b_pp[k] = {"ancestral": l5[k], "panel_af": f"{af:.4f}"}
    say("")
    self_derived(newB_aa, newB_pp, "derived allele the event writer ends with, P200 "
                 "emulation, non-SNP keys: old = L5.1 read (AA as new, SNPs not compared) "
                 "vs B", aa_old=newB_aa, pp_old=old_b_pp)
    agreeAA = sum(1 for k, r in newB_aa.items() if k in newC_pp and r["AA"] != "."
                  and (r["AA"] == k[3]) == (newC_pp[k]["ancestral"] == "ALT"))
    say(f"\nB SNP AA resolved: {sum(r['AA'] != '.' for r in newB_aa.values()):,} of "
        f"{len(newB_aa):,}; agree with C at {agreeAA:,}")

    # ---------------- option A: kit route A event rows (user, 2026-10-14)
    say("\n== option A: L8 and ET1291 indel/SV states from outgroup kit route A")
    # 1. concordance of the kit's event states with the same genomes' graph
    # genotypes, on the full CX333 VCF (both genomes are in its graph)
    _, fkeys, fcopies, _ = P.read_panel(vcf_full, bcf, seq, "NC_000962.3")
    fsmp = run([bcf, "query", "-l", vcf_full]).split()
    fkinds = [kind(k[2], k[3]) for k in fkeys]
    say("kit state (route A) vs the same genome's graph genotype, CX333, non-SNP keys")
    say(f"  {'genome':<15}{'class':<7}{'keys':>8}{'kit called':>12}{'both called':>13}"
        f"{'agree':>8}{'kit1/graph0':>13}{'kit0/graph1':>13}{'concordance':>13}")
    for acc in KIT:
        ev = P.kit_event_states(fkeys, fkinds, seq, "NC_000962.3",
                                P.read_kit_events(f"{KITEV}/{acc}", seq, "NC_000962.3"))
        gi = fsmp.index(acc)
        t = collections.defaultdict(collections.Counter)
        for j, cp in enumerate(fcopies):
            cl = fkinds[j]
            if cl == "SNP":
                continue
            g = {c[gi] for c in cp}
            gg = 1 if "1" in g else 0 if "0" in g else -1
            k_ = int(ev[j])
            t[cl]["keys"] += 1
            t[cl]["kit"] += k_ >= 0
            if k_ >= 0 and gg >= 0:
                t[cl]["both"] += 1
                t[cl]["agree"] += k_ == gg
                t[cl]["10"] += (k_, gg) == (1, 0)
                t[cl]["01"] += (k_, gg) == (0, 1)
        for cl in ("MNP", "indel", "SV"):
            x = t[cl]
            say(f"  {acc:<15}{cl:<7}{x['keys']:>8,}{x['kit']:>12,}{x['both']:>13,}"
                f"{x['agree']:>8,}{x['10']:>13,}{x['01']:>13,}"
                f"{(x['agree'] / x['both'] if x['both'] else 0):>13.4f}")

    # 2. the P200 emulation with the kit's events
    dE, logE = panel_run("BE", tree_b, vcf_b, ogs, KIT, a.work, events=True)
    dE2, _ = panel_run("BE2", tree_b, vcf_b, list(reversed(ogs)), KIT, a.work, events=True)
    sameE = all(filecmp.cmp(f"{dE}/{f}", f"{dE2}/{f}", shallow=False)
                for f in ("ancestral.tsv", "panel_polarity.tsv"))
    say(f"\nBE (B + kit events) vs BE2 (outgroup list reversed): identical: {sameE}")
    say("\n".join("  " + l for l in logE.splitlines() if "kit" in l or l.strip().startswith(
        ("class", "SNP", "MNP", "indel", "SV"))))
    newE_pp = {key(r): r for r in rows(f"{dE}/panel_polarity.tsv")}
    t = collections.defaultdict(collections.Counter)
    for k in set(newE_pp) | set(newB_pp) | set(newC_pp):
        cl = kind(k[2], k[3])
        e_, b_, c_ = (x[k]["ancestral"] if k in x else None for x in (newE_pp, newB_pp, newC_pp))
        if b_ is None and e_ is not None:
            t[cl]["tied in B, resolved in BE"] += 1
            t[cl]["  of which = C" if e_ == c_ else "  of which != C (or C tied)"] += 1
        if b_ is not None and e_ is None:
            t[cl]["resolved in B, tied in BE"] += 1
        if b_ is not None and e_ is not None and b_ != e_:
            t[cl]["B and BE differ"] += 1
        if e_ is not None and c_ is not None:
            t[cl]["BE agrees with C" if e_ == c_ else "BE differs from C"] += 1
        if e_ is None and c_ is not None:
            t[cl]["C resolved, BE tied"] += 1
    cols = ["tied in B, resolved in BE", "  of which = C", "  of which != C (or C tied)",
            "resolved in B, tied in BE", "B and BE differ", "BE agrees with C",
            "BE differs from C", "C resolved, BE tied"]
    say("\nBE vs B (without events) and vs C (full CX333 panel)")
    for cl in ("SNP", "MNP", "indel", "SV"):
        say(f"  {cl:<6}" + "; ".join(f"{c.strip()} {t[cl][c]:,}" for c in cols))

    # ---------------- the element
    say("\ntatC-helY element (68 bp insertion at 2,352,065)")
    for tag, pp in (("A", newA_pp), ("C", newC_pp), ("B", newB_pp), ("BE", newE_pp)):
        hit = [(k, r["ancestral"], r["method"]) for k, r in pp.items()
               if k[1] == ELEMENT and k[3] == ELEMENT_ALT]
        say(f"  run {tag}: " + (", ".join(f"ancestral={x[1]} ({x[2]})" for x in hit) or "TIED/absent"))
    oldel = [(r["ancestral"]) for k, r in old_pp.items() if k[1] == ELEMENT and k[3] == ELEMENT_ALT]
    say(f"  build panel_polarity.tsv (canettii read): {oldel}")

    with open(f"{OUT}/report.txt", "w") as fh:
        fh.write("\n".join(rep) + "\n")
    with open(f"{OUT}/snp_aa_differences_A.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\told_AA\tnew_AA\tnew_flag\n")
        for k, o, n, f in sorted(diffs, key=lambda x: x[0][1]):
            fh.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{o}\t{n}\t{f}\n")
    print(f"-> {OUT}/report.txt")


if __name__ == "__main__":
    main()
