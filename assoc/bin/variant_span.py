#!/usr/bin/env python3
"""Which bases a variant record changes, for crediting it to genes.

Shared by is6110_gene_burden.py (which genes a record is credited to) and
assoc_scan.py (whether a record is genic, for its null stratum), so the two
cannot disagree about a record. The rule, the user's decisions D18 and D39
(2026-10-07): the bases a variant changes decide its genes.
"""


def deleted_span(r):
    """(first, last) deleted base, 1-based, for a deletion record; else None.

    VCF POS is the anchor base BEFORE the deletion, so it is never deleted.
    The span comes from the catalogue ID `svi:DEL:<first>:<len>`, else from a
    symbolic `<DEL:-len>` ALT."""
    f = (r.get("id") or "").split("#")[0].split(":")
    if len(f) == 4 and f[0] == "svi" and f[1] == "DEL" \
            and f[2].isdigit() and f[3].isdigit():
        return int(f[2]), int(f[2]) + int(f[3]) - 1
    alt = r.get("alt") or ""
    if alt.startswith("<DEL:-") and alt.endswith(">") \
            and alt[6:-1].isdigit() and r["pos"].isdigit():
        return int(r["pos"]) + 1, int(r["pos"]) + int(alt[6:-1])
    return None


def small_deleted_span(r):
    """(first, last) deleted base, 1-based, for a left-anchored small
    deletion (REF longer than ALT and ALT is REF's first base); else None.

    The same rule as deleted_span: the anchor base at POS is not deleted, so
    the record removes POS+1 .. POS+len(REF)-1. The `*` allele is not the
    record's ALT. Insertions have no deleted base and complex replacements no
    anchor, so both stay on the unit at POS."""
    ref = (r.get("ref") or "").upper()
    alts = [x for x in (r.get("alt") or "").upper().split(",") if x != "*"]
    if len(alts) != 1 or not r["pos"].isdigit():
        return None
    alt = alts[0]
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt and \
            not set(ref) - set("ACGTN"):
        return int(r["pos"]) + 1, int(r["pos"]) + len(ref) - 1
    return None


def changed_span(r, cls):
    """(first, last, inside) for any record: the bases it changes, 1-based.

    THE BASES A VARIANT CHANGES DECIDE ITS GENES (user's decisions D18, D39,
    2026-10-07). `inside` is False for a span of changed bases, credited to
    every gene it touches, and True for an insertion, which changes no base:
    it sits between `first` and `last` (POS and POS+1, the VCF convention,
    symbolic IS6110 ALTs included) and is inside a gene only where the gene
    holds both. Crediting the anchor base instead put an insertion just after
    a gene's last base on that gene, and an MNP or complex record crossing a
    gene boundary on one gene only. None where POS is not a number."""
    if not r["pos"].isdigit():
        return None
    pos = int(r["pos"])
    span = (deleted_span(r) if cls == "sv" else
            small_deleted_span(r) if cls == "small" else None)
    if span:
        return span[0], span[1], False
    ref = (r.get("ref") or "").upper()
    alts = [x for x in (r.get("alt") or "").upper().split(",") if x != "*"]
    if cls != "small" or not ref or not alts:
        # IS6110 and other symbolic insertions, and SV records without a
        # deleted span: an insertion after POS
        return pos, pos + 1, True
    if len(alts) == 1 and len(ref) == 1 and len(alts[0]) > 1 \
            and alts[0][0] == ref:
        return pos, pos + 1, True              # left-anchored insertion
    # SNP, MNP or complex: REF's bases, less a first base every ALT shares
    k = 1 if len(ref) > 1 and all(len(x) > 1 and x[0] == ref[0]
                                  for x in alts) else 0
    return pos + k, pos + len(ref) - 1, False
