#!/usr/bin/env python3
"""Allele normalisation and cohort-invariant record keys.

The matrix keyed rows on `h37rv:<pos>`, which collapses distinct alleles at one
position into one row. Such a row needs a single canonical REF chosen from among
them, and `p5_keys.py` chose the first-seen record's -- so the choice depended on
sample iteration order. Adding 77 samples to a 23-sample cohort changed the
canonical REF at 17 positions and flipped 68 existing calls between ALT and REF.
A deterministic rule such as longest-REF does not fix it either, because adding a
sample can introduce a longer REF.

**Cohort-invariance requires the key and the canonical REF to be properties of the
allele alone.** So each record is normalised independently -- right-trimmed,
left-aligned, left-trimmed, the same operations `bcftools norm` performs -- and
the key carries the normalised REF and ALT. Adding samples then only ever adds
keys; it can never alter an existing one.

The cost is that two callers' different representations of the same event stay
separate rows. Left-alignment removes most of that, which is why it is done here
rather than trimming alone: without it, `TA->T` at 100 and `AT->A` at 101 are the
same deletion in two rows.
"""


def normalise(seq, pos, ref, alt):
    """Left-align and trim one allele. `pos` is 1-based, `seq` 0-based H37Rv.

    Returns (pos, ref, alt). Mirrors bcftools norm: trim the shared suffix, shift
    a pure indel as far left as the reference allows, then trim the shared prefix
    down to a single anchor base.
    """
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    if not ref or not alt or ref == alt:
        return pos, ref, alt

    # 1. trim the shared suffix, keeping at least one base on each side
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]

    # 2. left-shift a pure indel while the base before it matches the base that
    #    would be lost from the end. Only applies where one allele is a single
    #    base, which is what "pure indel" means after step 1.
    if seq and (len(ref) == 1) != (len(alt) == 1):
        while pos > 1 and len(ref) > 0 and len(alt) > 0 \
                and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
            prev = seq[pos - 2].upper()
            ref = prev + ref[:-1]
            alt = prev + alt[:-1]
            pos -= 1

    # 3. trim the shared prefix, keeping one anchor base
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def h37rv_key(pos, ref, alt):
    return f"h37rv:{pos}:{ref}>{alt}"


def node_key(node, offset, ref, alt):
    return f"node:{node}:{offset}:{ref}>{alt}"


def _selftest():
    # H37Rv-like context; index 0 is position 1
    seq = "ACGTTTTGCA" * 10
    cases = [
        # a shared-suffix pair collapses to the minimal allele
        ((5, "TTTG", "TTG"), "one base deleted from the T-run"),
        # the same deletion written at a different offset must normalise alike
        ((6, "TTG", "TG"), "same deletion, shifted representation"),
        # a SNP is untouched
        ((3, "G", "A"), "SNP"),
        # a long REF with a shared prefix trims to an anchor
        ((1, "ACGT", "A"), "3 bp deletion"),
    ]
    out = []
    for (p, r, a), label in cases:
        out.append((label, normalise(seq, p, r, a)))
    a = normalise(seq, 5, "TTTG", "TTG")
    b = normalise(seq, 6, "TTG", "TG")
    same = (a == b)
    for label, v in out:
        print(f"    {label:<38s} -> {v}")
    print(f"    two representations of one deletion agree: {same}")
    # keys must be stable and allele-specific
    k1 = h37rv_key(*normalise(seq, 5, "TTTG", "TTG"))
    k2 = h37rv_key(*normalise(seq, 5, "T", "TA"))
    print(f"    distinct alleles at one position get distinct keys: {k1 != k2}")
    print(f"      {k1}\n      {k2}")
    return same and k1 != k2


if __name__ == "__main__":
    import sys
    print("  mtb_norm self-test")
    sys.exit(0 if _selftest() else 1)
