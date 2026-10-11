#!/usr/bin/env python3
"""Before/after of the harmonization pass (bin/harmonize_records.py) on a
cohort's merged VCF. Report only.

  compare_merge.py <before.vcf.gz> <after.vcf.gz>

Prints:
  - records by CLASS (and IS6110LOSS / LOWCONF) before and after;
  - records removed and added by ID, by class;
  - cells whose GT changed in records kept, by class and change;
  - drug-resistance positive controls: every record in the WHO-catalogue
    gene windows below (H37Rv, gene +-200 bp) must keep its ID and GT
    cells, apart from DUP cells; any difference is listed.
"""
import collections
import gzip
import sys

DR = {"gyrA": (7302, 9818), "rpoB": (759807, 763325), "rpoC": (763370, 767320),
      "rpsL": (781560, 781934), "rrs": (1471846, 1473382), "fabG1-inhA": (1673300, 1675011),
      "katG": (2153889, 2156111), "pncA": (2288681, 2289241), "eis": (2714124, 2715332),
      "embB": (4246514, 4249810), "gid": (4407528, 4408202), "ethA": (4326004, 4327473)}


def load(path):
    out = {}
    for line in gzip.open(path, "rt"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        info = dict(x.partition("=")[::2] for x in f[7].split(";"))
        out[f[2]] = (f[0], int(f[1]), info, [c.split(":")[0] for c in f[9:]],
                     [c.split(":")[1] if ":" in c else "" for c in f[9:]])
    return out


def cls(info):
    c = info.get("CLASS", "")
    if "IS6110LOSS" in info:
        c += "+IS6110LOSS"
    return c


def main():
    b, a = load(sys.argv[1]), load(sys.argv[2])
    print("records by class: before -> after")
    cb = collections.Counter(cls(v[2]) for v in b.values())
    ca = collections.Counter(cls(v[2]) for v in a.values())
    for k in sorted(set(cb) | set(ca)):
        print(f"  {k:24s} {cb[k]:7,} -> {ca[k]:7,}")
    gone = collections.Counter(cls(b[k][2]) for k in b.keys() - a.keys())
    new = collections.Counter(cls(a[k][2]) for k in a.keys() - b.keys())
    print("removed:", dict(gone), " added:", dict(new))
    ch = collections.Counter()
    for k in a.keys() & b.keys():
        for x, y, st in zip(b[k][3], a[k][3], a[k][4]):
            if x != y:
                ch[(cls(a[k][2]), f"{x}->{y}", st)] += 1
    print("GT changes in kept records (class, change, new ST):")
    for k, v in sorted(ch.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:,}")
    print("drug-resistance gene windows (+-200 bp):")
    for g, (s, e) in DR.items():
        ids_b = {k for k, v in b.items() if v[0] == "NC_000962.3" and s - 200 <= v[1] <= e + 200}
        ids_a = {k for k, v in a.items() if v[0] == "NC_000962.3" and s - 200 <= v[1] <= e + 200}
        diff_gt = sum(1 for k in ids_a & ids_b for x, y, st in zip(b[k][3], a[k][3], a[k][4])
                      if x != y and st != "DUP")
        dup = sum(1 for k in ids_a & ids_b for st in a[k][4] if st == "DUP")
        print(f"  {g:11s} records {len(ids_b):4d} -> {len(ids_a):4d}; removed "
              f"{sorted(ids_b - ids_a)[:5]}{' ...' if len(ids_b - ids_a) > 5 else ''}; "
              f"added {len(ids_a - ids_b)}; non-DUP GT changes {diff_gt}; DUP cells {dup}")


if __name__ == "__main__":
    main()
