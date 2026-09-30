#!/usr/bin/env python3
"""Merge the IS6110 arm into P5's key space and per-sample states.

PLAN AND ITS REVIEW: is6110/docs/P5_MERGE_PLAN.md. Stage 1 here; stage 2 adds
the evidence that lets a non-carrier be called REF instead of NOCALL.

WHY A NON-CARRIER IS NOT REF BY DEFAULT. P5 can call a SNP position REF because
the gVCF says the position was assessed and had depth. An IS6110 detector emits
nothing at all for a site it did not call, and that silence has four causes:
the insertion is genuinely absent, the position has no equivalent in this
sample's reference context, there was not enough depth, or the detector did not
look. Defaulting to REF would assert roughly 3,700 unmeasured "no insertion
here" calls over the pilot. `p5_states.py` already takes the conservative line
in the closest case -- an off-path key a sample produced no record for is
NOCALL, "its own reference may or may not traverse that node, and this pass has
not asked" -- and 68 of the 170 IS6110 keys are node-framed.

KEYS GO THROUGH mtb_norm, NOT THROUGH STRING FORMATTING. The arm's own table
writes `h37rv:<pos>` and `node:<id>:<offset>`, while P5's key space is
`h37rv:<pos>:<REF>><ALT>` and `node:<id>:<offset>:<REF>><ALT>`. Mixing the two
is what silently disabled P6's second annotation tier for the whole project, so
these keys are built with the same functions P5 uses.

THE ALT IS SYMBOLIC. `<INS>` -- the representation P1I_VCF_DESIGN.md settled on
by precedent, since dysgu and delly already write symbolic ALTs into the same
directory. It is deliberately NOT passed through `normalise`, which trims shared
prefixes and suffixes and has no meaning for a symbolic allele.
"""
import argparse, collections, csv, importlib.util, os, sys

_s = importlib.util.spec_from_file_location(
    "mtb_norm", os.path.join("bin", "mtb_norm.py"))
mtb_norm = importlib.util.module_from_spec(_s); _s.loader.exec_module(mtb_norm)

ALT = "<INS>"


def load_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort-keys",
                    default="is6110/results/p1i_cohort_keys.tsv")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--h37rv",
                    default="refbias/build/7713a8d71d8e/refs/GCF_000195955.fasta")
    ap.add_argument("--out-keys", default="is6110/results/p5_is6110_keys.tsv")
    ap.add_argument("--out-states",
                    default="is6110/results/p5_is6110_states.tsv",
                    help="the whole-cohort table; empty to skip it")
    ap.add_argument("--out-states-dir", default="",
                    help="also write one table per sample, <dir>/<sample>.tsv, "
                         "for stage 2 to run per sample without reading the "
                         "cohort's (keys x samples rows: 50 million at 10,000)")
    a = ap.parse_args()

    h37 = load_seq(a.h37rv)
    samples = [r["sample"] for r in
               csv.DictReader(open(a.refmap, newline=""), delimiter="\t")]

    rows = list(csv.DictReader(open(a.cohort_keys, newline=""), delimiter="\t"))
    keys, observed = {}, {}
    skipped = collections.Counter()
    for x in rows:
        if x["frame"] == "h37rv":
            if not x["h37rv_pos"]:
                skipped["h37rv frame with no position"] += 1
                continue
            pos = int(x["h37rv_pos"])
            ref = h37[pos - 1]
            k = mtb_norm.h37rv_key(pos, ref, ALT)
            meta = dict(key=k, frame="h37rv", h37rv_pos=pos, node="",
                        node_offset="", canonical_ref=ref, canonical_alt=ALT,
                        region="is6110", kind="IS6110", acc_locus="")
        else:
            nid = x["node"]
            if ":" not in nid:
                skipped["node key without an offset"] += 1
                continue
            node, off = nid.split(":", 1)
            # no H37Rv base to canonicalise against off the path, which is the
            # same position p4_place.py takes for off-path records
            k = mtb_norm.node_key(node, off, "N", ALT)
            meta = dict(key=k, frame="node", h37rv_pos=x["h37rv_pos"] or "",
                        node=node, node_offset=off, canonical_ref="N",
                        canonical_alt=ALT, region="is6110", kind="IS6110",
                        acc_locus="")
        keys.setdefault(k, meta)
        # a carrier's own state: REF where its matched reference already holds a
        # copy at that locus (the site sits at an excision join), ALT otherwise
        observed[(x["sample"], k)] = x["state"]

    with open(a.out_keys, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(next(iter(keys.values()))),
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for k in sorted(keys):
            w.writerow(keys[k])

    counts = collections.Counter()
    FIELDS = ["sample", "key", "state", "allele", "frame", "region", "kind"]
    if a.out_states_dir:
        os.makedirs(a.out_states_dir, exist_ok=True)
    cohort_fh = open(a.out_states, "w", newline="") if a.out_states else None
    w = None
    if cohort_fh:
        w = csv.DictWriter(cohort_fh, fieldnames=FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
    skeys = sorted(keys)
    for s in samples:
        srows = []
        for k in skeys:
            st = observed.get((s, k))
            if st is None:
                st, allele = "NOCALL", ""
            else:
                allele = ALT if st == "ALT" else keys[k]["canonical_ref"]
            counts[st] += 1
            srows.append(dict(sample=s, key=k, state=st, allele=allele,
                              frame=keys[k]["frame"], region="is6110",
                              kind="IS6110"))
        if w:
            w.writerows(srows)
        if a.out_states_dir:
            p = os.path.join(a.out_states_dir, f"{s}.tsv")
            with open(p + ".tmp", "w", newline="") as sf:
                sw = csv.DictWriter(sf, fieldnames=FIELDS, delimiter="\t",
                                    lineterminator="\n")
                sw.writeheader(); sw.writerows(srows)
            os.replace(p + ".tmp", p)
    if cohort_fh:
        cohort_fh.close()

    n = sum(counts.values())
    print(f"  {len(keys)} IS6110 keys x {len(samples)} samples = {n} cells")
    for k in ("ALT", "REF", "ABSENT", "NOCALL"):
        if counts[k]:
            print(f"    {k:8} {counts[k]:6} ({100*counts[k]/n:5.1f}%)")
    for k, v in skipped.items():
        print(f"    skipped: {k}: {v}")
    print(f"  NOCALL is every sample that did not report the site. Stage 2 of "
          f"is6110/docs/P5_MERGE_PLAN.md\n  is what turns the earned ones into "
          f"REF; nothing here asserts an absence.")
    print(f"  written: {a.out_keys}"
          + (f"\n  written: {a.out_states}" if a.out_states else "")
          + (f"\n  written: {a.out_states_dir}/<sample>.tsv" if a.out_states_dir else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
