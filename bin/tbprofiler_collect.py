#!/usr/bin/env python3
"""Collect tb-profiler per-assembly CSVs into one table and flag mixed lineages.

A mixed or chimeric assembly shows as two incompatible lineage calls, each at a
high fraction -- e.g. lineage2.2.2 and lineage4.4 both at 100%. That is different
from a nested call (lineage2, lineage2.2, lineage2.2.2 all at 100%), which is just
the barcode hierarchy for one strain and is normal.
"""
import argparse, csv, glob, os, re


def parse(path):
    strain, rows = "", []
    with open(path) as fh:
        section = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("Strain,"):
                strain = line.split(",", 1)[1].strip()
            elif line.startswith("Lineage report"):
                section = "lin"; continue
            elif line.startswith("Lineage,Fraction"):
                continue
            elif section == "lin":
                f = line.split(",")
                if len(f) >= 2 and f[0].startswith(("lineage", "La", "M.")):
                    try: rows.append((f[0], float(f[1])))
                    except ValueError: pass
                elif not line.strip():
                    section = None
    return strain, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-fraction", type=float, default=20.0,
                    help="a second lineage above this fraction makes the call mixed")
    a = ap.parse_args()

    panel = {l.split("\t")[0] for i, l in enumerate(open(a.panel)) if i}
    out = []
    for p in sorted(glob.glob(os.path.join(a.dir, "*.results.csv"))):
        acc = os.path.basename(p).replace(".results.csv", "")
        strain, rows = parse(p)
        # collapse the barcode hierarchy: keep only the deepest call per root
        roots = {}
        for lin, frac in rows:
            m = re.match(r"(lineage\d+|La\d+|M\.\w+)", lin)
            if not m: continue
            r = m.group(1)
            if frac >= a.min_fraction and (r not in roots or len(lin) > len(roots[r][0])):
                roots[r] = (lin, frac)
        mixed = len(roots) > 1
        out.append({"accession": acc, "in_panel": acc in panel, "strain": strain,
                    "n_lineage_roots": len(roots), "mixed": mixed,
                    "roots": ";".join(f"{v[0]}:{v[1]:.0f}" for v in roots.values())})
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(out)

    import collections
    print(f"[collect] {len(out)} assemblies -> {a.out}")
    mx = [r for r in out if r["mixed"]]
    print(f"  MIXED-LINEAGE calls: {len(mx)}")
    for r in mx:
        print(f"    {r['accession']:<16s} {'in panel' if r['in_panel'] else 'excluded':<9s} {r['roots']}")
    inp = [r for r in out if r["in_panel"]]
    c = collections.Counter(re.match(r"(lineage\d+|La\d+|M\.\w+)", r["roots"]).group(1)
                            if re.match(r"(lineage\d+|La\d+|M\.\w+)", r["roots"]) else "?"
                            for r in inp if r["roots"])
    print(f"\n  lineage composition of the {len(inp)} panel genomes:")
    for k, n in c.most_common():
        print(f"    {k:<10s} {n}")
    print(f"    (no call: {sum(1 for r in inp if not r['roots'])})")


if __name__ == "__main__":
    main()
