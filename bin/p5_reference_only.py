#!/usr/bin/env python3
"""Per sample: which calls in the merged VCF rest only on the matched reference.

WHY (HANDOFF fix list, ABS-1 item 4). Every sample is read through the panel
genome R it was matched to, and some of its calls are R's rather than its own:
a node-frame key stated ABSENT because R's path does not traverse the node, an
inherited SV written without a read check. ABS-1 and ABS-2 replaced most of
those with checks against the sample's reads, and every cell that still rests
on R alone -- or that the checks moved -- now says so in FORMAT/EV. This
counts them per sample, so a sample whose calls lean on its reference more
than the others' stands out, and flags the case measured to make that worse:

  ref_deeper   R's sublineage is DEEPER than the sample's own lineage call
               (R is lineage1.2.1.2.1, the sample only lineage1.2.1). On
               scale200_fix those 7 samples had a median of 72 ABSENT cells
               contradicted by their reads, against 24 for the rest
               (analysis/inherited_absent/README.md); SAMEA2297133, whose
               intact ctpV was stated ABSENT, is one of them.

Columns:
  sample, reference, sample_sublineage, ref_sublineage
  ref_deeper              yes / no / unknown (a sublineage call missing)
  reference_only          cells whose EV says they rest on R alone:
                          node_not_on_r, unchecked, inherited_unchecked
  absent_confirmed        ABSENT kept because reads join R's deletion ends
  absent_to_sample_call   ABSENT replaced by the sample's own H37Rv-frame call
                          (reads run across R's deletion)
  masked                  cells NOCALL because the site is in the repeat mask
  frames_conflict         inherited deletions the two frames disagree on
  evidence                every EV code with its count, `code=n|...`

    p5_reference_only.py --vcf p5/merged.vcf.gz --refmap p1/refmap.tsv \\
        --out p5/reference_only.tsv
"""
import argparse, collections, csv, gzip, os, sys

REFERENCE_ONLY = {"node_not_on_r", "unchecked", "inherited_unchecked"}
CONFIRMED = {"junction_joined", "absent_junction_joined"}
TO_SAMPLE = {"present_here_ref", "present_here_alt", "absent_checked_present_ref",
             "absent_checked_present_alt"}


def deeper(sample_sub, ref_sub):
    if not sample_sub or not ref_sub:
        return "unknown"
    return "yes" if ref_sub.startswith(sample_sub + ".") else "no"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    for p in (a.vcf, a.refmap):
        if not os.path.exists(p):
            sys.exit(f"FATAL: {p} does not exist")
    meta = {r["sample"]: r for r in csv.DictReader(open(a.refmap, newline=""),
                                                   delimiter="\t")}
    for col in ("reference", "sample_sublineage", "ref_sublineage"):
        if meta and col not in next(iter(meta.values())):
            sys.exit(f"FATAL: {a.refmap} has no `{col}` column")

    samples, ev = None, None
    op = gzip.open if a.vcf.endswith(".gz") else open
    with op(a.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                ev = [collections.Counter() for _ in samples]
                continue
            f = line.rstrip("\n").split("\t")
            fmt = f[8].split(":")
            if "EV" not in fmt:
                continue                   # no cell of this record has evidence
            k = fmt.index("EV")
            for j, c in enumerate(f[9:]):
                x = c.split(":")
                if len(x) > k and x[k] and x[k] != ".":
                    ev[j][x[k]] += 1
    if samples is None:
        sys.exit(f"FATAL: {a.vcf} has no #CHROM line")
    missing = [s for s in samples if s not in meta]
    if missing:
        sys.exit(f"FATAL: {len(missing)} VCF samples are not in {a.refmap} "
                 f"(e.g. {missing[0]})")

    cols = ["sample", "reference", "sample_sublineage", "ref_sublineage",
            "ref_deeper", "reference_only", "absent_confirmed",
            "absent_to_sample_call", "masked", "frames_conflict", "evidence"]
    n_deep = 0
    with open(a.out + ".tmp", "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for s, c in zip(samples, ev):
            m = meta[s]
            d = deeper(m.get("sample_sublineage", ""), m.get("ref_sublineage", ""))
            n_deep += d == "yes"
            w.writerow([s, m.get("reference", ""), m.get("sample_sublineage", ""),
                        m.get("ref_sublineage", ""), d,
                        sum(v for k, v in c.items() if k in REFERENCE_ONLY),
                        sum(v for k, v in c.items() if k in CONFIRMED),
                        sum(v for k, v in c.items() if k in TO_SAMPLE),
                        sum(v for k, v in c.items() if k.endswith("masked")),
                        sum(v for k, v in c.items() if k.endswith("frames_conflict")),
                        "|".join(f"{k}={v}" for k, v in sorted(c.items()))])
    os.replace(a.out + ".tmp", a.out)
    ro = sorted(sum(v for k, v in c.items() if k in REFERENCE_ONLY) for c in ev)
    print(f"  {len(samples)} samples; calls resting only on the reference: "
          f"median {ro[len(ro) // 2] if ro else 0}, maximum {ro[-1] if ro else 0}")
    print(f"  {n_deep} samples matched to a reference from a deeper sublineage "
          f"than their own lineage call (ref_deeper=yes)")
    print(f"  -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
