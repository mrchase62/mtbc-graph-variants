#!/usr/bin/env python3
"""Check a cohort against refbias/io_contract.tsv before anything is submitted.

WHY BEFORE. Every pass already validates its own inputs and names the file it
could not find, but that happens inside the job -- after the array is
submitted and after afterok has committed every downstream pass to waiting on
it. A missing input therefore costs a submission and cancels the rest of the
chain. Three times today: p2_summary dying on a blank field cancelled nine
jobs, and P4 missing the P1 H37Rv-frame calls cancelled four more after all
100 of its own tasks failed. This moves the same knowledge in front of the
submission, where the cost of being wrong is a message.

It also answers the question the per-cohort inventory cannot. That table
records what each pass WRITES; this one records what each pass READS, which is
what decides whether a pass can run at all.

Two modes:
  --check   what is missing that a pass needs. The default, and the one to run
            before submitting.
  --report  the full contract resolved against this cohort, inputs and
            outputs, present or not -- the readable I/O map.
"""
import argparse, csv, glob, os, re, sys

def load_registry(path, name):
    for line in open(path):
        if line.startswith("#") or line.startswith("cohort\t"):
            continue
        f = line.rstrip("\n").split("\t")
        if f and f[0] == name:
            return dict(cohort=f[0], table=f[1], crams=f[2], outroot=f[3],
                        passes=f[4], note=f[5] if len(f) > 5 else "",
                        workpfx=f[6] if len(f) > 6 else "")
    raise SystemExit(f"cohort '{name}' is not in {path}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort")
    ap.add_argument("--contract", default="refbias/io_contract.tsv")
    ap.add_argument("--registry", default="refbias/cohorts.tsv")
    ap.add_argument("--build", default="refbias/build/7713a8d71d8e")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--pass", dest="only", default="", help="one pass only")
    ap.add_argument("--will-run", default="",
                    help="comma-separated passes about to be submitted. An "
                         "input a LATER pass needs is satisfied if an EARLIER "
                         "pass in this list produces it, so a fresh cohort is "
                         "not asked to already contain its own outputs.")
    a = ap.parse_args()

    reg = load_registry(a.registry, a.cohort)
    samples = [l.split("\t")[0] for i, l in enumerate(open(reg["table"])) if i]
    sub = {"OUTROOT": reg["outroot"], "WORKPFX": reg["workpfx"],
           "BUILD": a.build, "COHORT": reg["table"], "CRAMS": reg["crams"],
           "REFMAP": f"{reg['outroot']}/p1/refmap.tsv",
           # NAME is the cohort's own tag, for the IS6110 arm's result files,
           # which are named after the cohort rather than placed under OUTROOT.
           "NAME": reg["cohort"]}

    # The IS6110 result tables are not all "<cohort>_...": the pilot's are
    # unprefixed and scale100's carry "scale_", which is what p1i_vcf.sh,
    # p1i_p5states.sh and merge_cohort_vcf.py all use. Resolving the contract's
    # ${NAME}_ literally sent the check looking for pilot_p1i_cohort_keys.tsv
    # and refused a correct submission.
    is_pfx = {"pilot": "", "pilot_rerun": "", "scale100": "scale_"}.get(
        reg["cohort"], reg["cohort"] + "_")

    def resolve(p):
        p = p.replace("is6110/results/${NAME}_", "is6110/results/" + is_pfx)
        for k, v in sub.items():
            p = p.replace("${" + k + "}", v)
        return p

    # ORDER matters: an input is satisfied either because it exists on disk or
    # because a pass running EARLIER in this chain writes it. Without that, a
    # brand-new cohort fails the check on nine inputs that the very submission
    # being checked would create -- which is what happened on the first run of
    # gwas1000, and would have trained everyone to pass IO_CHECK=0.
    ORDER = ["p1", "p2", "p3", "p4", "p4b", "p5", "p1g", "p1i", "p1iv",
             "p1is", "p5svgt", "p5vcf"]
    will = [p for p in ORDER if p in
            {x.strip() for x in a.will_run.split(",") if x.strip()}]

    rows = [r for r in csv.DictReader(
        (l for l in open(a.contract) if not l.startswith("#")), delimiter="\t")]
    registered = set(reg["passes"].split(","))

    # For each pass in the chain, the resolved paths that a pass running
    # earlier in the same chain writes.
    upstream, seen = {}, set()
    for q in will:
        upstream[q] = frozenset(seen)
        seen |= {resolve(r["path"]) for r in rows
                 if r["pass"] == q and r["when"] == "post"}

    problems = []
    print(f"  cohort {a.cohort}: {len(samples)} isolates, outroot {reg['outroot']}")
    print(f"  registered passes: {reg['passes']}\n")

    for r in rows:
        p = r["pass"]
        if a.only and p != a.only:
            continue
        if p not in registered:
            continue
        path = resolve(r["path"])
        if r["scope"] == "sample":
            have = sum(1 for s in samples
                       if os.path.exists(path.replace("${S}", s)))
            want = len(samples)
            ok = have == want
            state = f"{have}/{want}"
        else:
            if path.endswith("/"):
                ok = os.path.isdir(path) and bool(os.listdir(path))
                state = "dir" if ok else "-"
            elif "${S}" in path:
                have = sum(1 for s in samples
                           if os.path.exists(path.replace("${S}", s)))
                ok = have == len(samples)
                state = f"{have}/{len(samples)}"
            else:
                ok = os.path.exists(path)
                state = "yes" if ok else "-"
        # Compare RESOLVED paths, not templates: p2 asks for ${REFMAP} while
        # p1 declares it as ${OUTROOT}/p1/refmap.tsv, and those are the same
        # file under different spellings.
        produced_earlier = path in upstream.get(p, frozenset())
        bad = ((not ok) and r["required"] == "yes" and r["when"] == "pre"
               and not produced_earlier)
        if bad:
            problems.append((p, path, state))
        if a.report:
            mark = ("ok " if ok else
                    "upstream" if produced_earlier else
                    "MISSING" if r["when"] == "pre" else "not yet")
            print(f"  {p:5s} {r['when']:4s} {r['scope']:6s} {mark:8s} {state:>9s}  {path}")

    if a.report:
        print()
    if problems:
        print(f"  {len(problems)} required INPUT(S) missing -- do not submit:")
        for p, path, state in problems:
            print(f"    {p:5s} {state:>9s}  {path}")
        return 1
    print("  every required input for the registered passes is present")
    return 0

if __name__ == "__main__":
    sys.exit(main())
