#!/usr/bin/env python3
"""Build an IS6110-clean H37Rv: element copies excised, element kept as its own contig.

WHY
Every IS6110 measurement in this arm fights the same thing: H37Rv carries 16
near-identical copies, so a read carrying element sequence has sixteen places it
can align. The aligner may tuck a junction read wholly inside a reference copy
rather than clipping it at its true site, which loses the junction outright, and
reference copies and novel insertions end up needing two separate and
differently-broken measurements (`is6110/README.md` section 4).

Excise every copy and that disappears. On an IS-free backbone a read carrying
element sequence has nowhere else to go, so it MUST clip at its junction, and
H37Rv's own copies become ordinary detectable insertions like any other. The idea
is taken from ISdetector (`docs/ISDETECTOR_EVAL.md`), which excises rather than
N-masks; masking would keep coordinates stable but leaves the copy occupying its
position, so absence stays undetectable and the thing we want fixed is not fixed.

WHY THE ELEMENT IS KEPT AS A SEPARATE CONTIG
A plain excision would break SA-based attribution, because the clipped part of a
junction read would have nothing to realign to and bwa would emit no SA record at
all. Adding the element back as its own contig restores it and makes it
unambiguous: the supplementary alignment target IS the element, rather than
"within 20 bp of the terminus of one of sixteen copies".

It also makes family copy number readable from the same alignment -- depth on the
element contig against depth on the chromosome -- with no division by sixteen and
no collapsed-copy problem. That is a different estimator from P1d and should be
compared with it, not substituted for it.

COORDINATES
Excision shifts everything downstream, so this writes a crossmap and validates it
by round-tripping. Unlike ISdetector, which discovers copies by BLAST at 95%
identity and 0.90 coverage, the sixteen spans here are annotated and exact, so
the arithmetic is exact and no threshold is involved.

Clean -> original is the only direction that is always well defined. Original ->
clean is undefined for a base inside an excised span, and the caller is told so
rather than given a silently wrong answer.

THE TANDEM PAIR
IS6110-12 (3,551,230-3,552,584) and IS6110-13 (3,552,713-3,554,067) sit 128 bp
apart. Excising both leaves a 128 bp island between two artificial joins, too
short to anchor a read pair. It is reported explicitly at build time so it is a
known property of the reference rather than a strange result discovered later.
"""
import argparse, bisect, csv, sys


def read_fasta(path):
    name, buf, out = None, [], {}
    for line in open(path):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def load_spans(gff):
    """(start, end) 1-based inclusive, sorted, merged if they overlap."""
    spans = []
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 4:
            spans.append((int(f[3]), int(f[4])))
    spans.sort()
    merged = []
    for s, e in spans:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def build(seq, spans):
    """Excise spans. Returns (clean_seq, crossmap rows).

    crossmap row: the junction left by one excision, in both frames.
      orig_start/orig_end  the excised span, 1-based inclusive, original frame
      clean_junction       1-based position in the clean sequence of the LAST base
                           before the excision (0 if the span starts at base 1)
      cum_deleted          bases deleted at and before this span
    """
    out, rows, prev, cum = [], [], 0, 0
    for s, e in spans:
        out.append(seq[prev:s - 1])
        cum += (e - s + 1)
        rows.append(dict(orig_start=s, orig_end=e, deleted_len=e - s + 1,
                         clean_junction=(s - 1) - (cum - (e - s + 1)),
                         cum_deleted=cum))
        prev = e
    out.append(seq[prev:])
    return "".join(out), rows


def clean_to_orig(p, rows):
    """1-based clean position -> 1-based original position."""
    shift = 0
    for r in rows:
        if p > r["clean_junction"]:
            shift = r["cum_deleted"]
        else:
            break
    return p + shift


def orig_to_clean(p, rows):
    """1-based original -> 1-based clean, or None if the base was excised."""
    shift = 0
    for r in rows:
        if r["orig_start"] <= p <= r["orig_end"]:
            return None
        if p > r["orig_end"]:
            shift = r["cum_deleted"]
        else:
            break
    return p - shift


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reference", required=True, help="H37Rv FASTA")
    ap.add_argument("--contig", default="NC_000962.3")
    ap.add_argument("--gff", required=True, help="IS6110 element GFF")
    ap.add_argument("--element", required=True,
                    help="canonical element FASTA, added as its own contig")
    ap.add_argument("--element-name", default="IS6110")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-crossmap", required=True)
    ap.add_argument("--annotation-check", default=None,
                    help="optional BED/GFF whose features are round-tripped")
    a = ap.parse_args()

    ref = read_fasta(a.reference)
    key = a.contig if a.contig in ref else next(
        (k for k in ref if k.endswith(a.contig)), None)
    if key is None:
        sys.exit(f"contig {a.contig} not in {a.reference}")
    seq = ref[key].upper()
    spans = load_spans(a.gff)
    if not spans:
        # A REFERENCE WITH NO COPIES IS LEGITIMATE, NOT AN ERROR. Some MTBC
        # lineages carry a handful of IS6110 insertions and some carry about
        # one, so a genome with none is a real genome, not a failed discovery.
        # There is nothing to excise, so the element-free build is the
        # reference unchanged with the canonical element appended as its own
        # contig, and the crossmap is empty because no coordinate moves. This
        # used to exit, which dropped exactly the references where junction
        # detection against the sample's own reference works best.
        print(f"    {a.gff}: no element intervals; building with nothing "
              f"excised, which is correct for a genome that carries none")
    total = sum(e - s + 1 for s, e in spans)
    clean, rows = build(seq, spans)

    if len(clean) != len(seq) - total:
        sys.exit(f"length check failed: {len(clean)} != {len(seq)} - {total}")

    # every base of the clean sequence must equal its original-frame base
    bad = 0
    step = max(1, len(clean) // 200000)
    for i in range(1, len(clean) + 1, step):
        if clean[i - 1] != seq[clean_to_orig(i, rows) - 1]:
            bad += 1
    if bad:
        sys.exit(f"crossmap round-trip failed at {bad} sampled positions")

    # and the two directions must be inverses away from excised spans
    for i in range(1, len(seq) + 1, step):
        c = orig_to_clean(i, rows)
        if c is not None and clean_to_orig(c, rows) != i:
            sys.exit(f"crossmap is not invertible at original position {i}")

    el = read_fasta(a.element)
    elseq = "".join(el.values()).upper()

    with open(a.out_fasta, "w") as fh:
        fh.write(f">{a.contig}_isclean excised {len(spans)} IS6110 copies, "
                 f"{total} bp removed from {len(seq)}\n")
        for i in range(0, len(clean), 60):
            fh.write(clean[i:i + 60] + "\n")
        fh.write(f">{a.element_name} canonical element, {len(elseq)} bp\n")
        for i in range(0, len(elseq), 60):
            fh.write(elseq[i:i + 60] + "\n")

    with open(a.out_crossmap, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["orig_start", "orig_end", "deleted_len",
                                           "clean_junction", "cum_deleted"],
                           delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"  {len(spans)} element copies excised, {total} bp")
    print(f"  chromosome {len(seq)} -> {len(clean)} bp, plus {a.element_name} "
          f"contig at {len(elseq)} bp")
    print(f"  crossmap round-trip verified at {len(range(1, len(clean)+1, step))} "
          f"positions, both directions")

    # short islands left between adjacent excisions
    print("\n  gaps left between adjacent excisions:")
    for i in range(1, len(spans)):
        gap = spans[i][0] - spans[i - 1][1] - 1
        if gap < 1000:
            cj = rows[i - 1]["clean_junction"]
            print(f"    {gap:6d} bp island between {spans[i-1][0]}-{spans[i-1][1]} "
                  f"and {spans[i][0]}-{spans[i][1]}, at clean {cj}-{cj+gap}"
                  + ("   <-- too short to anchor a pair" if gap < 300 else ""))
    print(f"\n  written: {a.out_fasta}\n           {a.out_crossmap}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
