#!/usr/bin/env python3
"""Stage 1 of P1i: find IS6110 intervals in an assembly that carries no element annotation.

WHY THIS HAS TO EXIST
is6110_build_isclean.py excises element intervals and writes a crossmap that
shifts every coordinate downstream of each cut, so it needs boundaries to the
base. For H37Rv those are given: RefSeq annotates sixteen mobile_genetic_element
features and the builder simply reads them. The eighteen SNP-matched references
the pilot selects have no such feature -- their GFFs under
refbias/build_cold/<build>/annotation/ carry the transposase as CDS and nothing
else -- so the intervals must be discovered.

WHAT WAS TRIED AND REJECTED, so it is not tried again
  terminal inverted repeats  H37Rv's annotation gives both 28 bp inverted
      repeats verbatim, and anchoring on them would make boundaries exact by
      construction. But two of H37Rv's sixteen copies carry termini ~20
      mismatches from the canonical pair, so a budget tight enough to be
      specific misses them and one loose enough to catch them is not specific.
  target-site duplication    Scored 0/16 as a boundary refiner: the element's
      own inverted repeats contain 6-mers that outrank the real 3-4 bp
      duplication anywhere within the search window. Checked at the annotated
      boundaries directly, only 10 of 16 carry a 3-5 bp direct repeat at all,
      most of them 3-mers, which occur by chance every few dozen bases. The
      duplication is strong evidence in READ data, where clipped stacks say
      where to look and the offset between them confirms the geometry. It is
      not evidence in a bare genome sequence.

WHAT IS USED
Whole-element alignment with minimap2, the canonical 1355 bp copy as query, and
an interval extended to full query length where the alignment stops short. On
H37Rv this returns the sixteen annotated spans with fourteen agreeing base for
base; see --grade.

THE TWO IT DOES NOT GET EXACTLY are the two copies LONGER than the query -- the
annotation gives them as 1357 and 1375 bp against the canonical 1355 -- so a
1355 bp query cannot span them and the interval falls short by 1 bp and by 19
bp. That is a property of defining the element as "the region homologous to the
canonical copy", not a tuning failure, and no threshold here changes it.
"""
import argparse, gzip, os, re, subprocess, sys, tempfile


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf, out = None, [], {}
    for line in op(path, "rt"):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def annotated_elements(gff_gz):
    """The exact intervals, where a genome has them. H37Rv only, in practice."""
    out = []
    op = gzip.open if gff_gz.endswith(".gz") else open
    for line in op(gff_gz, "rt"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "mobile_genetic_element" and "IS6110" in f[8]:
            out.append((f[0], int(f[3]), int(f[4])))
    return sorted(out, key=lambda r: (r[0], r[1]))


def transposase_loci(gff_gz):
    """Distinct transposase loci, for gates 2 and 3.

    The element's transposase is annotated as CDS in two rows per copy because of
    its programmed frameshift, so rows are grouped by their Parent gene. Product
    naming varies across assemblies -- 'IS3-like element IS6110/IS987 family
    transposase' and 'IS3-like element IS987 family transposase' are both used,
    IS987 being a synonym -- so the test is on those two tokens. The generic 'IS3
    family transposase' is deliberately NOT matched: the IS3 family contains
    other elements and counting them would inflate the gate's expectation.
    """
    spans = {}
    op = gzip.open if gff_gz.endswith(".gz") else open
    for line in op(gff_gz, "rt"):
        f = line.rstrip("\n").split("\t")
        if len(f) <= 8 or f[2] != "CDS":
            continue
        if not re.search(r"IS6110|IS987", f[8]):
            continue
        m = re.search(r"Parent=([^;]+)", f[8])
        key = (f[0], m.group(1) if m else f"{f[0]}:{f[3]}")
        s, e = int(f[3]), int(f[4])
        if key in spans:
            spans[key] = (spans[key][0], min(spans[key][1], s), max(spans[key][2], e))
        else:
            spans[key] = (f[0], s, e)
    return sorted(spans.values(), key=lambda r: (r[0], r[1]))


def run_minimap2(mm2, assembly, query, extra):
    cmd = [mm2, "-c", "-x", "asm20", "-N", "200", "-p", "0.02",
           "--secondary=yes"] + extra + [assembly, query]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"minimap2 failed:\n{p.stderr[-2000:]}")
    return p.stdout


def parse_paf(paf, min_ident, min_cov):
    hits = []
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        strand, tname = f[4], f[5]
        tlen, ts, te = int(f[6]), int(f[7]), int(f[8])
        nmatch, blocklen = int(f[9]), int(f[10])
        ident = nmatch / blocklen if blocklen else 0.0
        cov = (qe - qs) / qlen
        if ident < min_ident or cov < min_cov:
            continue
        # EXTEND TO FULL QUERY. An alignment that stops short of the query's end
        # has the rest of the element sitting in the target unaligned, and
        # excising only the aligned part would leave element sequence behind.
        # Which side the unaligned head goes to depends on strand.
        head, tail = qs, qlen - qe
        if strand == "-":
            head, tail = tail, head
        s = max(1, ts + 1 - head)
        e = min(tlen, te + tail)
        hits.append(dict(contig=tname, start=s, end=e, strand=strand,
                         ident=ident, cov=cov, qs=qs + 1, qe=qe, qlen=qlen,
                         ext5=ts + 1 - s, ext3=e - te))
    return hits


def merge(hits, min_overlap=0.5):
    """Collapse intervals that are the same copy; keep adjacent copies separate.

    A plain overlap test is wrong here, because extending an alignment to full
    query length can push a truncated copy into its neighbour. GCF_039770655
    carries exactly that: a copy truncated by its first 136 bases immediately
    5' of a full-length one, which a plain merge reported as a single 2574 bp
    interval and gate 2 caught as one interval too few. Two copies abutting each
    other is a real arrangement, not an artefact, so the rule is reciprocal:
    intervals are one copy only when they overlap by more than half of the
    shorter. Otherwise the overlap is trimmed off whichever side was extended
    further, since that extension is what created it.
    """
    out = []
    for h in sorted(hits, key=lambda x: (x["contig"], x["start"])):
        h = dict(h)
        if not (out and h["contig"] == out[-1]["contig"] and h["start"] <= out[-1]["end"]):
            out.append(h)
            continue
        prev = out[-1]
        ov = min(prev["end"], h["end"]) - max(prev["start"], h["start"]) + 1
        shorter = min(prev["end"] - prev["start"] + 1, h["end"] - h["start"] + 1)
        if shorter and ov / shorter > min_overlap:
            prev["end"] = max(prev["end"], h["end"])
            if h["ident"] > prev["ident"]:
                prev["ident"], prev["strand"] = h["ident"], h["strand"]
        elif prev.get("ext3", 0) >= h.get("ext5", 0):
            prev["end"] = h["start"] - 1
            out.append(h)
        else:
            h["start"] = prev["end"] + 1
            out.append(h)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--name", default=None, help="accession, for the report")
    ap.add_argument("--query", default="is6110/assets/IS6110.query.fasta")
    ap.add_argument("--minimap2", default=os.environ.get("MTB_MINIMAP2", "minimap2"))
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-coverage", type=float, default=0.10,
                    help="fraction of the query an alignment must cover. Set "
                         "from a sweep, not by taste: on H37Rv the interval set "
                         "is identical at 0.80 and at 0.10 -- 16 intervals, no "
                         "spurious ones -- because it is the identity floor, "
                         "not coverage, that excludes the other IS3-family "
                         "elements. Lowering it costs nothing where the answer "
                         "is known and recovers a truncated copy in "
                         "GCF_022870365 that gate 3 was reporting as missing. "
                         "A fragment of element sequence left in the chromosome "
                         "recruits element reads exactly as a whole copy does, "
                         "so excising it is the point rather than a side "
                         "effect.")
    ap.add_argument("--annotation", default=None,
                    help="that assembly's GFF, for gates 2 and 3")
    ap.add_argument("--grade", action="store_true",
                    help="gate 1: compare against annotated element intervals, "
                         "which only H37Rv has")
    ap.add_argument("--out-gff", default=None)
    a = ap.parse_args()

    name = a.name or os.path.basename(a.assembly).split(".")[0]
    hits = merge(parse_paf(run_minimap2(a.minimap2, a.assembly, a.query, []),
                           a.min_identity, a.min_coverage))

    if a.out_gff:
        with open(a.out_gff, "w") as fh:
            fh.write("##gff-version 3\n")
            for i, h in enumerate(hits, 1):
                fh.write(f"{h['contig']}\tis6110_discover\tmobile_genetic_element\t"
                         f"{h['start']}\t{h['end']}\t{h['ident']:.4f}\t{h['strand']}\t.\t"
                         f"ID=IS6110_disc_{i};Note=discovered by whole-element "
                         f"alignment,identity {h['ident']:.4f},query coverage "
                         f"{h['cov']:.3f}\n")

    lens = sorted(h["end"] - h["start"] + 1 for h in hits)
    print(f"{name}: {len(hits)} intervals, lengths {lens[0] if lens else 0}"
          f"-{lens[-1] if lens else 0}, identity "
          f"{min((h['ident'] for h in hits), default=0):.3f}-"
          f"{max((h['ident'] for h in hits), default=0):.3f}")

    rc = 0
    if a.grade and a.annotation:
        ann = annotated_elements(a.annotation)
        print(f"\n  GATE 1  boundaries against {len(ann)} annotated intervals")
        used, exact = set(), 0
        for c, s, e in ann:
            best, bi = None, None
            for i, h in enumerate(hits):
                if h["contig"] != c or i in used:
                    continue
                if h["start"] <= e and h["end"] >= s:
                    d = abs(h["start"] - s) + abs(h["end"] - e)
                    if best is None or d < best:
                        best, bi = d, i
            if bi is None:
                print(f"    {s}-{e}   NOT FOUND")
                continue
            used.add(bi)
            h = hits[bi]
            ds, de = h["start"] - s, h["end"] - e
            exact += (ds == 0 and de == 0)
            flag = "exact" if (ds == 0 and de == 0) else f"off {ds:+d} / {de:+d}"
            print(f"    {s}-{e}  ->  {h['start']}-{h['end']}   {flag}")
        extra = len(hits) - len(used)
        print(f"    {exact}/{len(ann)} base-exact, {len(used)}/{len(ann)} found, "
              f"{extra} discovered interval(s) with no annotated counterpart")
        if exact != len(ann) or extra:
            print("    GATE 1: FAIL as specified (base-exact on all annotated spans)")
            rc = 1
        else:
            print("    GATE 1: pass")

    if a.annotation:
        tl = transposase_loci(a.annotation)
        # GATE 3 tests COVERAGE BY THE UNION of intervals, not containment in a
        # single one. The gate asks whether any element sequence would be left
        # behind by the excision, and a transposase spanning two abutting
        # intervals is not left behind. GCF_039770655 carries a pair that
        # overlap by about 45 bp -- one copy inserted into the end of another --
        # whose transposase runs 4 bp past the boundary between them; a
        # containment test called that a miss when nothing was missing.
        by_contig = {}
        for h in hits:
            by_contig.setdefault(h["contig"], []).append((h["start"], h["end"]))
        covered, split, insets, uncovered = 0, 0, [], []
        for c, s_, e_ in tl:
            iv = sorted(by_contig.get(c, []))
            need = [(s_, e_)]
            hit_n = 0
            for a_, b_ in iv:
                nxt = []
                for x, y in need:
                    if b_ < x or a_ > y:
                        nxt.append((x, y)); continue
                    hit_n += 1
                    if a_ > x: nxt.append((x, a_ - 1))
                    if b_ < y: nxt.append((b_ + 1, y))
                need = nxt
            if not need:
                covered += 1
                split += (hit_n > 1)
                for a_, b_ in iv:
                    if a_ <= s_ and b_ >= e_:
                        insets.append((s_ - a_, b_ - e_)); break
            else:
                uncovered.append((c, s_, e_, sum(y - x + 1 for x, y in need)))
        print(f"\n  GATE 2  {len(hits)} intervals against {len(tl)} annotated "
              f"transposase loci  -> "
              f"{'ok' if len(hits) >= len(tl) else 'FEWER INTERVALS THAN TRANSPOSASES'}")
        if tl:
            msg = f"  GATE 3  {covered}/{len(tl)} transposases fully covered"
            if split:
                msg += f" ({split} spanning two abutting intervals)"
            if insets:
                i5 = sorted(x for x, _ in insets); i3 = sorted(y for _, y in insets)
                msg += (f"; inset 5' {i5[0]}-{i5[-1]}, 3' {i3[0]}-{i3[-1]} "
                        f"(H37Rv reference value 42-51)")
            print(msg)
            for c, s_, e_, n in uncovered:
                print(f"            UNCOVERED {c}:{s_}-{e_}, {n} bp")
            if covered != len(tl) or len(hits) < len(tl):
                rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
