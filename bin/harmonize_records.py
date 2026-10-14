#!/usr/bin/env python3
"""One event, one record: the cross-arm harmonization pass of
merge_cohort_vcf.py (approved 2026-10-11; analysis/qc_harmonize/PLAN.md,
part B).

The merged VCF is filled by arms that never look at each other: small
variants (HaplotypeCaller, plus P4's composed alleles), structural variants
(the deletion catalogue and the delly/dysgu caller matrix) and IS6110
insertions. On scale200_fix 787 SV records (13%) repeated a small-class
indel, 291 caller insertions (14%) were IS6110 insertions the IS6110 arm
already had, and IS6110 copy losses appeared only as anonymous ~1.36 kb
deletions (analysis/harmonize/overlap.py). Each rule below gives an event to
the arm best placed to call it, on H37Rv-frame records:

  1. Indels under 50 bp belong to HaplotypeCaller. In a sample where a
     caller SV record (COMPONENT called or mixed) is ALT and a small-class
     indel of 10-49 bp within HC_WIN (60 bp) is ALT, the SV cell is set to
     no-call (ST=DUP), when the SV is under 50 bp or within PCTSIZE (70%) of
     the indel's size. bwa clips reads at a 20-49 bp indel and the SV callers
     report it again (Phase B: the rule removed 11 false calls, no true one).
  2. Events of 50 bp or more belong to the SV block. A small-class allele
     with an indel of >= 50 bp that matches an SV record (start within
     HC_WIN, same direction, sizes within PCTSIZE, ALT carriers overlapping
     by Jaccard >= 0.5) is dropped; its key goes to the SV record's
     SMALLALLELE, and the SV record's no-call cells (not those rules 1 and 3
     set to DUP) take the small record's state. One with no SV match becomes CLASS=sv itself (its exact REF/ALT
     kept, SVSOURCE=small_allele), so every >= 50 bp event is CLASS=sv.
  3. IS6110 insertions belong to the IS6110 arm. In a sample where a caller
     SV insertion of any size is ALT and an IS6110 record within IS_WIN (100
     bp) is ALT, the SV cell is set to no-call (ST=DUP): the callers report
     the element's junctions as truncated insertions (median 234 bp on
     scale200_fix), which the merge's exact-length rule (1,340-1,370 bp
     within 20 bp) let through. Small-class alleles relabelled by rule 2
     carry their exact sequence, so they are not truncated: they are handed
     over only when IS6110-sized (IS_LEN, 1,300-1,420 bp). Shorter
     insertions near an IS6110 site are other events and are kept.
  4. IS6110 copy losses get one labelled record. A deletion record (either
     SV source, or a small-class allele relabelled by rule 2) whose
     first and last deleted bases lie within LOSS_TOL (100 bp) of the ends of
     an H37Rv IS6110 copy becomes CLASS=is6110 with IS6110LOSS and LOWCONF:
     the IS6110 arm cannot state a loss (its allele is an insertion), and
     short reads do not separate real losses from false ones (Phase B: 35
     real, 16 false at copies; depth inside a lost copy stays near normal
     from reads of the other copies). LOSS_TOL was 50 bp, widened to 100
     bp (user, 2026-10-11): on scale200_fix it added 2 deletions, both
     removing the whole copy plus 40-80 bp of flank; no deletion within 300
     bp of the 16 H37Rv copies removes only part of one. SVLEN keeps the
     full size.
  5. Every record lists the arms that saw it in SOURCES: small, sv_caller,
     sv_catalogue, is6110, accessory; an arm whose cells or record were
     absorbed is added to the record that kept the event.
A caller SV record, or a relabelled small allele, left with no ALT cell is
dropped: every carrier's call was absorbed by the arm that owns the event,
which states the other samples too.

Within a sharded merge, each shard harmonizes its own records; a pair split
across a shard boundary (within 100 bp of it) is not seen. Counted in the
log as a limit, not handled.
"""
import bisect
import collections

HC_WIN = 60
IS_WIN = 100
IS_LEN = (1300, 1420)
PCTSIZE = 0.7
LOSS_TOL = 100
DUP_CELL = ".:DUP"


def _info(rec):
    d = {}
    for x in rec[5]:
        k, _, v = x.partition("=")
        d[k] = v
    return d


def _set_info(info, key, value=None):
    """info list with key set (replaced if present); value None = flag."""
    out = [x for x in info if x.partition("=")[0] != key]
    out.append(key if value is None else f"{key}={value}")
    return out


def _cells(rec):
    return rec[6].split("\t")


def _is_alt(cell):
    return cell.startswith("1:")


def _source(info):
    cls = info.get("CLASS", "")
    if cls == "sv":
        return "sv_catalogue" if info.get("COMPONENT") == "interval" else "sv_caller"
    if cls == "accessory_presence":
        return "accessory"
    return cls or "unknown"


def _small_indel(rec):
    """Signed indel length of a small record (+ insertion, - deletion), or 0."""
    alt = rec[4][0] if rec[4] else ""
    if not alt or alt.startswith("<") or alt == "*":
        return 0
    return len(alt) - len(rec[3])


def _sv_len(info):
    try:
        return abs(int(info.get("SVLEN", "0")))
    except ValueError:
        return 0


def _close_size(a, b):
    return a > 0 and b > 0 and min(a, b) >= PCTSIZE * max(a, b)


def read_is6110_copies(gff):
    """H37Rv IS6110 copies, [(first base, last base)] 1-based, from the
    build's is6110_elements.gff."""
    out = []
    for line in open(gff):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 4 and "IS6110" in line:
            out.append((int(f[3]), int(f[4])))
    return sorted(out)


def harmonize(recs, h37_contig, copies=()):
    """Apply rules 1-5 to merge_cohort_vcf.py record tuples
    (chrom, pos, ident, ref, alts, info, cells). Returns (records, counts)."""
    n = collections.Counter()
    recs = list(recs)
    infos = [_info(r) for r in recs]
    sources = [{_source(i)} for i in infos]
    drop = set()
    cells = {}                        # index -> list of cells, once changed

    def get_cells(i):
        if i not in cells:
            cells[i] = _cells(recs[i])
        return cells[i]

    h37 = [i for i, r in enumerate(recs) if r[0] == h37_contig]
    small_ind = sorted((recs[i][1], i) for i in h37
                       if infos[i].get("CLASS") == "small" and abs(_small_indel(recs[i])) >= 10)
    sm_pos = [p for p, _ in small_ind]
    is_ins = sorted((recs[i][1], i) for i in h37 if infos[i].get("CLASS") == "is6110")
    is_pos = [p for p, _ in is_ins]
    sv_idx = [i for i in h37 if infos[i].get("CLASS") == "sv"]
    sv_sorted = sorted((recs[i][1], i) for i in sv_idx)
    sv_pos = [p for p, _ in sv_sorted]

    def window(pos_list, items, p, w):
        lo = bisect.bisect_left(pos_list, p - w)
        hi = bisect.bisect_right(pos_list, p + w)
        return [i for _, i in items[lo:hi]]

    # rule 1: caller SV cells that are HaplotypeCaller's 10-49 bp indels
    for i in sv_idx:
        if infos[i].get("COMPONENT") not in ("called", "mixed"):
            continue
        L = _sv_len(infos[i])
        near = [j for j in window(sm_pos, small_ind, recs[i][1], HC_WIN)
                if 10 <= abs(_small_indel(recs[j])) < 50
                and (L < 50 or _close_size(L, abs(_small_indel(recs[j]))))]
        if not near:
            continue
        c = get_cells(i)
        for j in near:
            cj = get_cells(j)
            for s, x in enumerate(c):
                if _is_alt(x) and _is_alt(cj[s]):
                    c[s] = DUP_CELL
                    sources[j].add("sv_caller")
                    n["rule1 SV cells handed to HaplotypeCaller"] += 1

    def hand_to_is6110(cands, src):
        for i in cands:
            near = window(is_pos, is_ins, recs[i][1], IS_WIN)
            if not near:
                continue
            c = get_cells(i)
            for j in near:
                cj = get_cells(j)
                for s, x in enumerate(c):
                    if not _is_alt(x):
                        continue
                    if _is_alt(cj[s]):
                        c[s] = DUP_CELL
                        sources[j].add(src)
                        n["rule3 SV INS cells handed to the IS6110 arm"] += 1
                    # the state is the SECOND field: a cell may carry a
                    # trailing evidence field (GT:ST:EV, ABS-1)
                    elif cj[s].split(":")[1:2] in (["REF"], ["ABSENT"]):
                        n["rule3 SV INS ALT where the IS6110 arm says REF/ABSENT (kept)"] += 1

    # rule 3: caller SV insertion cells that are IS6110 insertions
    hand_to_is6110([i for i in sv_idx if infos[i].get("SVTYPE") == "INS"
                    and infos[i].get("COMPONENT") != "interval"], "sv_caller")

    # rule 2: small-class alleles of >= 50 bp
    for _, j in small_ind:
        d = _small_indel(recs[j])
        if abs(d) < 50:
            continue
        typ = "INS" if d > 0 else "DEL"
        best = None
        for i in window(sv_pos, sv_sorted, recs[j][1], HC_WIN):
            if i in drop or infos[i].get("SVTYPE") != typ or not _close_size(_sv_len(infos[i]), abs(d)):
                continue
            ci, cj = get_cells(i), get_cells(j)
            a = {s for s, x in enumerate(ci) if _is_alt(x)}
            b = {s for s, x in enumerate(cj) if _is_alt(x)}
            jac = len(a & b) / len(a | b) if a | b else 1.0
            if jac >= 0.5 and (best is None or jac > best[1]):
                best = (i, jac)
        if best is None:
            recs[j] = recs[j][:5] + (
                _set_info(_set_info(_set_info(_set_info(recs[j][5], "CLASS", "sv"),
                                              "SVTYPE", typ), "SVLEN", str(d)),
                          "SVSOURCE", "small_allele"),) + recs[j][6:]
            infos[j] = _info(recs[j])
            n["rule2 small alleles >= 50 bp relabelled CLASS=sv (no SV record)"] += 1
            continue
        i = best[0]
        ci, cj = get_cells(i), get_cells(j)
        for s, x in enumerate(ci):
            if x.startswith(".:") and x != DUP_CELL and not cj[s].startswith(".:"):
                ci[s] = cj[s]
                n["rule2 SV no-call cells filled from the small allele"] += 1
            elif _is_alt(x) != _is_alt(cj[s]) and not cj[s].startswith(".:") and not x.startswith(".:"):
                n["rule2 cells where the SV record and the small allele disagree (SV kept)"] += 1
        prev = infos[i].get("SMALLALLELE")
        recs[i] = recs[i][:5] + (_set_info(recs[i][5], "SMALLALLELE",
                                           f"{prev},{recs[j][2]}" if prev else recs[j][2]),) + recs[i][6:]
        infos[i] = _info(recs[i])
        sources[i] |= sources[j]
        drop.add(j)
        n["rule2 small alleles >= 50 bp merged into an SV record"] += 1

    relabelled = [j for j in h37 if j not in drop and infos[j].get("SVSOURCE") == "small_allele"]
    # rule 3 for small-class alleles relabelled by rule 2: exact sequence, not
    # truncated, so only IS6110-sized ones are the IS6110 arm's insertion
    hand_to_is6110([j for j in relabelled if infos[j].get("SVTYPE") == "INS"
                    and IS_LEN[0] <= _sv_len(infos[j]) <= IS_LEN[1]], "small")

    # rule 4: deletions that are an H37Rv IS6110 copy
    for i in sv_idx + relabelled:
        if i in drop or infos[i].get("SVTYPE") != "DEL":
            continue
        first = recs[i][1] + 1     # symbolic DEL: anchor base before; allele: shared first base
        if infos[i].get("SVSOURCE") == "small_allele":
            last = recs[i][1] + len(recs[i][3]) - 1
        else:
            try:
                last = int(infos[i].get("END", ""))
            except ValueError:
                continue
        if any(abs(first - a) <= LOSS_TOL and abs(last - b) <= LOSS_TOL for a, b in copies):
            info = _set_info(_set_info(_set_info(recs[i][5], "CLASS", "is6110"), "IS6110LOSS"), "LOWCONF")
            recs[i] = recs[i][:5] + (info,) + recs[i][6:]
            infos[i] = _info(recs[i])
            n["rule4 deletions relabelled as IS6110 copy losses"] += 1

    out = []
    for i, r in enumerate(recs):
        if i in drop:
            continue
        if i in cells:
            c = cells[i]
            if infos[i].get("CLASS") == "sv" and (infos[i].get("COMPONENT") in ("called", "mixed", "inherited")
                                                  or infos[i].get("SVSOURCE") == "small_allele") \
                    and not any(_is_alt(x) for x in c):
                n["caller SV records dropped: no ALT cell left"] += 1
                continue
            r = r[:6] + ("\t".join(c),)
        r = r[:5] + (_set_info(r[5], "SOURCES", ",".join(sorted(sources[i]))),) + r[6:]
        out.append(r)
    return out, n


HEADER = [
    '##INFO=<ID=SOURCES,Number=1,Type=String,Description="Arms that reported '
    'this event: small, sv_caller, sv_catalogue, is6110, accessory. Calls of '
    'the same event by another arm are absorbed into this record '
    '(bin/harmonize_records.py)">\n',
    '##INFO=<ID=SMALLALLELE,Number=1,Type=String,Description="Key of the '
    'small-class allele (exact sequence) merged into this SV record">\n',
    '##INFO=<ID=IS6110LOSS,Number=0,Type=Flag,Description="A deletion of an '
    'H37Rv IS6110 copy (ends within 100 bp of the copy\'s ends)">\n',
    '##INFO=<ID=LOWCONF,Number=0,Type=Flag,Description="Short reads cannot '
    'reliably call this event class (IS6110 copy losses: depth inside a lost '
    'copy stays near normal from reads of the other copies); keep out of '
    'association and convergence tests">\n',
]
