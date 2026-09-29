#!/usr/bin/env python3
"""P5 step 1: the union of record keys across the cohort, with a canonical REF.

Two samples matched to different references can produce the same variant with
different REF strings, because REF is whatever their own reference carries. That
is section 2.2.2's allele-normalisation problem, and it has to be resolved before
a matrix can have one row per site.

The resolution for H37Rv-framed records: **the canonical REF is H37Rv's own base
at that position, read from the reference FASTA**, and the ALT is the sample's
allele. That is what composition means -- (H37Rv base -> sample allele) -- and it
makes the row identical no matter which reference the sample was called against.
A record whose REF disagrees with H37Rv's actual base is counted and reported
rather than silently coerced, since stage 4 found 38.8% of paftools-lifted records
landing where H37Rv's base was not their REF and that was a real defect, not noise.

For node-framed (off-path) records the canonical REF is left as observed. The key
is already exact -- the same node -- but the local sequence a sample reports
depends on its own reference's traversal, and there is no single "H37Rv base" to
normalise to precisely because the node is off the H37Rv path. Disagreements are
recorded per key so the scale of the problem is visible instead of assumed away.
"""
import argparse, collections, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtb_norm import normalise, h37rv_key, node_key


def rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def load_fasta_one(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line.strip())
    return "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--dir", default="refbias/p4")
    ap.add_argument("--h37rv", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    h37 = load_fasta_one(a.h37rv)
    if len(h37) < 4_000_000:
        print(f"H37Rv FASTA looks wrong: {len(h37)} bp", file=sys.stderr); return 1

    keys = {}
    ref_disagree = 0
    ref_variants = collections.defaultdict(set)
    n_rec = 0
    for r in rd(a.refmap):
        p = os.path.join(a.dir, f"{r['sample']}.placed.tsv")
        if not os.path.exists(p):
            continue
        for x in rd(p):
            n_rec += 1
            if x["frame"] == "h37rv":
                pos = int(x["h37rv_pos"])
                # compose against H37Rv's own base, then normalise; the key is a
                # property of the allele and cannot move when samples are added
                canon = h37[pos - 1:pos - 1 + len(x["ref"])].upper()
                if canon != x["ref"].upper():
                    ref_disagree += 1
                npos, nref, nalt = normalise(h37, pos, canon or x["ref"],
                                             x["alt"])
                k = h37rv_key(npos, nref, nalt)
                ref_variants[k].add(nref)
                if k not in keys:
                    keys[k] = dict(key=k, frame="h37rv", h37rv_pos=npos,
                                   node=x["node"], node_offset=x["node_offset"],
                                   canonical_ref=nref, canonical_alt=nalt,
                                   region=x["region"],
                                   kind=("SNP" if len(nref) == 1 and len(nalt) == 1
                                         else "INDEL"),
                                   acc_locus=x["acc_locus"])
            else:
                # off the H37Rv path there is no H37Rv base to compose against,
                # so the observed alleles are trimmed but not left-shifted
                npos, nref, nalt = normalise("", int(x["r_pos"] or 0),
                                             x["ref"], x["alt"])
                k = node_key(x["node"], x["node_offset"], nref, nalt)
                ref_variants[k].add(nref)
                if k not in keys:
                    keys[k] = dict(key=k, frame="node",
                                   h37rv_pos=x["h37rv_pos"], node=x["node"],
                                   node_offset=x["node_offset"],
                                   canonical_ref=nref, canonical_alt=nalt,
                                   region=x["region"],
                                   kind=("SNP" if len(nref) == 1 and len(nalt) == 1
                                         else "INDEL"),
                                   acc_locus=x["acc_locus"])

    if not keys:
        print("no keys found", file=sys.stderr); return 1
    rows = sorted(keys.values(),
                  key=lambda r: (r["frame"] != "h37rv", int(r["h37rv_pos"] or 0)))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    multi = sum(1 for k, v in ref_variants.items() if len(v) > 1)
    print(f"  keys are allele-specific and normalised, so a key cannot change "
          f"when samples are added")
    byframe = collections.Counter(r["frame"] for r in rows)
    byregion = collections.Counter(r["region"] for r in rows)
    print(f"  {n_rec} records -> {len(rows)} distinct keys")
    for k, n in byframe.most_common():
        print(f"    frame {k:<7s}{n:>7d}")
    for k, n in byregion.most_common():
        print(f"    region {k:<20s}{n:>7d}")
    print(f"\n  keys where samples reported DIFFERENT REF strings: {multi} "
          f"({100*multi/len(rows):.2f}%)")
    print(f"    -- this is the allele-normalisation problem, and its size")
    print(f"  H37Rv-framed records whose REF disagreed with H37Rv's actual base: "
          f"{ref_disagree}")
    print(f"    (stage 4 measured 38.8% for paftools-lifted records)")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
