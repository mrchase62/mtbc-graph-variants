#!/usr/bin/env python3
"""P5 step 2: resolve one sample's state at every key in the cohort union.

A matrix cell must distinguish four things, and a plain VCF conflates the last
three into a missing value:

  ALT     the sample carries a non-reference allele here
  REF     the sample was observed here and matches the canonical reference
  ABSENT  this position has no equivalent in the sample's own reference context
  NOCALL  nothing was observed, so nothing can be said

T6 measured 508,039 cells where the REF/ABSENT distinction was being lost, 494 of
them at lineage-tracking positions. Writing ABSENT as REF asserts the sample has
reference sequence it may not have; writing it as NOCALL discards real
information. Both are wrong in ways that propagate into any downstream tree or
association test, so all four states are kept distinct.

How each is decided:

  ALT     the sample has a placed record whose allele differs from the canonical
          REF. Note the composition: a composed-arm record's own REF is its
          MATCHED reference's base, not H37Rv's, so the allele is re-expressed
          against the canonical REF before comparison.
  REF     either the sample has a record that composes to the canonical allele --
          a REVERSION, where the sample differs from its reference precisely by
          matching H37Rv, measured here at 12.5% of composed/called records and
          which would otherwise enter the matrix as a false variant -- or the
          sample has no record and its GVCF shows coverage at the projected
          position.
  ABSENT  the key's H37Rv position does not project into the sample's reference
          (odgi reports dist.to.ref != 0), so the sample's reference has no such
          position and the sequence is not there to be called.
  NOCALL  projects fine, but the GVCF shows no observation.

GVCF coverage is read from DP, never MIN_DP. Stage 4 found one GVCF block
spanning positions 14 to 15,224 with DP 59 and MIN_DP 3, and using MIN_DP
rejected 464 of 1,131 inherited differences spuriously.
"""
import argparse, bisect, csv, gzip, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtb_norm import normalise, h37rv_key, node_key


def op(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_gvcf_blocks(path):
    """Sorted (start, end, dp) over the reference blocks and variant sites."""
    blocks = []
    for line in op(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        pos = int(f[1])
        end = pos
        for kv in f[7].split(";"):
            if kv.startswith("END="):
                try:
                    end = int(kv[4:])
                except ValueError:
                    pass
        # DP, not MIN_DP: MIN_DP is the block minimum and rejects real coverage
        keys = f[8].split(":")
        vals = f[9].split(":")
        dp = 0
        if "DP" in keys:
            try:
                dp = int(vals[keys.index("DP")])
            except (ValueError, IndexError):
                dp = 0
        blocks.append((pos, max(end, pos), dp))
    blocks.sort()
    return blocks


def make_cov(blocks, min_dp):
    starts = [b[0] for b in blocks]

    def covered(p):
        i = bisect.bisect_right(starts, p) - 1
        if i < 0:
            return False
        s, e, dp = blocks[i]
        return s <= p <= e and dp >= min_dp
    return covered


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--node-positions",
                    default="accessory/assets/node_positions.tsv",
                    help="node offsets along each panel accession, from the "
                         "same GFA pass as --node-paths. With it a node-frame "
                         "key whose reference carries the sequence is read in "
                         "the sample's own frame and gets REF when covered.")
    ap.add_argument("--node-paths",
                    default="accessory/assets/node_paths.tsv",
                    help="node to panel-path membership from "
                         "accessory/bin/node_path_membership.py. Without it a "
                         "node-frame key with no record is NOCALL; with it, one "
                         "whose sequence is absent from this sample's own "
                         "reference is ABSENT, which is a measurement.")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--keys", required=True)
    ap.add_argument("--placed", required=True)
    ap.add_argument("--gvcf", required=True)
    ap.add_argument("--projected", required=True,
                    help="odgi position output for the keys, H37Rv -> R")
    ap.add_argument("--h37rv", required=True,
                    help="needed to recompute the same normalised keys the key "
                         "set used; pairing on position alone is what broke")
    ap.add_argument("--min-dp", type=int, default=5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dense", action="store_true",
                    help="write the old one-row-per-key form. The default is "
                         "sparse: only the cells that are not REF, keyed by "
                         "the key's INDEX in --keys rather than by its name. "
                         "At 997 isolates a dense file is 10 MB per sample and "
                         "12 GB for the cohort, of which 82.7% is the word REF "
                         "beside a key string that is already in keys.tsv; the "
                         "sparse form is about 240 KB. Every consumer reads "
                         "either, and the header carries the key count and a "
                         "checksum of the key list so a file written against a "
                         "different key set is refused rather than silently "
                         "misaligned.")
    a = ap.parse_args()

    keys = list(csv.DictReader(open(a.keys, newline=""), delimiter="\t"))
    placed = list(csv.DictReader(open(a.placed, newline=""), delimiter="\t"))
    h37 = []
    with open(a.h37rv) as fh:
        for line in fh:
            if line.startswith(">"):
                if h37:
                    break
                continue
            h37.append(line.strip())
    h37 = "".join(h37)
    if len(h37) < 4_000_000:
        print(f"H37Rv FASTA looks wrong: {len(h37)} bp", file=sys.stderr)
        return 1

    # Pair by SOURCE position, not by file order. odgi re-emits its header line
    # between records, so the raw line count is about twice the query count, and
    # while the data lines happen to come back 1:1 and in order, relying on that
    # is an assumption with no upside: the source position is right there in
    # column 1. A silently misaligned coordinate is the worst failure this stage
    # could produce.
    h_keys = [k for k in keys if k["frame"] == "h37rv"]
    by_src = {}
    n_lines = n_multi = 0
    for line in open(a.projected):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            src = int(f[0].rsplit(",", 2)[-2]) + 1
            tgt = int(f[1].rsplit(",", 2)[-2]) + 1
            dist = int(f[2])
        except (ValueError, IndexError):
            continue
        n_lines += 1
        if src in by_src:
            n_multi += 1          # the target path visits this position twice
            continue
        by_src[src] = (tgt, dist)
    projmap = {}
    for k in h_keys:
        projmap[k["key"]] = by_src.get(int(k["h37rv_pos"]))
    n_unproj = sum(1 for v in projmap.values() if v is None)
    if n_lines == 0:
        print("  FATAL: no projected positions parsed", file=sys.stderr)
        return 1

    # Re-derive each record's normalised key exactly as p5_keys.py did. The
    # record's own `key` field is position-based and no longer what the key set
    # uses.
    by_key = {}
    for x in placed:
        if x["frame"] == "h37rv":
            pos = int(x["h37rv_pos"])
            canon = h37[pos - 1:pos - 1 + len(x["ref"])].upper() or x["ref"]
            npos, nref, nalt = normalise(h37, pos, canon, x["alt"])
            k = h37rv_key(npos, nref, nalt)
        else:
            npos, nref, nalt = normalise("", int(x["r_pos"] or 0), x["ref"],
                                         x["alt"])
            k = node_key(x["node"], x["node_offset"], nref, nalt)
        by_key.setdefault(k, dict(x, alt=nalt))

    cov = make_cov(load_gvcf_blocks(a.gvcf), a.min_dp)

    rows = []
    # node -> the panel accessions whose path traverses it, so a node-frame
    # key can be told "absent from this sample's reference" from "unknown".
    nodepaths = {}
    if a.node_paths and os.path.exists(a.node_paths):
        for q in csv.DictReader(open(a.node_paths, newline=""), delimiter="\t"):
            nodepaths[q["node"]] = set(q["paths"].split(",")) if q["paths"] else set()
        print(f"  node-path membership for {len(nodepaths):,} nodes from "
              f"{a.node_paths}")
    else:
        print("  NOTE no --node-paths table: every node-frame key without a "
              "record stays NOCALL, which is 97% of them")

    # (node, accession) -> (start, strand, occurrences) along that accession's
    # own sequence, so a node-frame key can be read in the sample's own frame.
    nodepos = {}
    if a.node_positions and os.path.exists(a.node_positions):
        for q in csv.DictReader(open(a.node_positions, newline=""), delimiter="\t"):
            nodepos[(q["node"], q["accession"])] = (
                int(q["start"]), q["strand"], int(q.get("n_occurrences") or 1))
        print(f"  node positions for {len(nodepos):,} (node, accession) pairs")

    counts = {"ALT": 0, "REF": 0, "ABSENT": 0, "NOCALL": 0}
    reversions = 0
    for k in keys:
        canon = k["canonical_ref"].upper()
        rec = by_key.get(k["key"])
        if rec is not None:
            alt = rec["alt"].upper()
            if alt == canon:
                state, reversions = "REF", reversions + 1
                allele = canon
            else:
                state, allele = "ALT", alt
        else:
            pr = projmap.get(k["key"])
            if k["frame"] == "node":
                # AN OFF-PATH KEY THIS SAMPLE PRODUCED NO RECORD FOR. The
                # earlier version wrote NOCALL here and said why: "its own
                # reference may or may not traverse that node, and this pass has
                # not asked". accessory/ACCESSORY_CALL_CENSUS.md measures what
                # not asking cost -- 91 REF calls in 6.9 million accessory
                # cells, 40 of 63 loci with zero REF anywhere, and TbD1's 424
                # variant records genotyped in 0.2% of cells despite two test
                # isolates covering the sequence at 0.99 and 1.00.
                #
                # Node-to-path membership is a property of the graph alone, so
                # it can be asked without a projection and without a realignment
                # -- accessory/bin/node_path_membership.py computes it once.
                # Where the sample's own reference does NOT traverse the node,
                # the sequence is not in its reference context and the honest
                # state is ABSENT, which is the same statement an H37Rv-frame
                # key gets when its position does not project. Measured on
                # gwas1000: 90.25% of node-frame cells are that case.
                #
                # Where the reference DOES traverse it, 7.07% of cells, REF
                # needs depth at the corresponding position and this pass still
                # cannot say -- so those stay NOCALL, and closing them is the
                # second half of the fix.
                nd = k.get("node") or ""
                paths = nodepaths.get(nd)
                if paths is None:
                    state, allele = "NOCALL", ""
                elif a.reference and a.reference not in paths:
                    state, allele = "ABSENT", ""
                else:
                    # THE REFERENCE DOES CARRY THIS SEQUENCE, so the sample was
                    # looked at here and the only question is whether the reads
                    # cover it. That is the same question the H37Rv-frame keys
                    # answer, asked in the sample's own frame: the node's offset
                    # along its reference plus the key's offset within the node
                    # gives the position, and the gVCF gives the depth.
                    #
                    # Node positions come from the same GFA pass as membership,
                    # so this still needs no projection and no realignment.
                    hit = nodepos.get((nd, a.reference))
                    off = int(k.get("node_offset") or 0)
                    if hit is None:
                        state, allele = "NOCALL", ""
                    elif hit[2] > 1:
                        # The reference visits this node more than once, so its
                        # copies share reads and depth at one of them settles
                        # nothing. 8.8% of pairs.
                        state, allele = "NOCALL", ""
                    else:
                        start, strand = hit[0], hit[1]
                        p = start + off if strand == "+" else start - off
                        state, allele = (("REF", canon) if cov(p)
                                         else ("NOCALL", ""))
            elif pr is None or pr[1] != 0:
                state, allele = "ABSENT", ""
            elif cov(pr[0]):
                state, allele = "REF", canon
            else:
                state, allele = "NOCALL", ""
        counts[state] += 1
        rows.append(dict(sample=a.sample, key=k["key"], state=state,
                         allele=allele, frame=k["frame"],
                         region=k["region"], kind=k["kind"]))

    if a.dense:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)
    else:
        # THE KEY SET IS THE INDEX, so it has to be pinned. A sparse file is
        # meaningless against a different key order, and a silently misaligned
        # coordinate is the worst failure this stage can produce -- so the
        # header carries the count and a checksum of the key names, and the
        # reader refuses a mismatch instead of trusting it.
        import hashlib
        h = hashlib.sha1("\n".join(k["key"] for k in keys).encode()).hexdigest()
        with open(a.out, "w", newline="") as fh:
            fh.write(f"#format\tsparse-v1\n#sample\t{a.sample}\n"
                     f"#default\tREF\n#n_keys\t{len(keys)}\n"
                     f"#keys_sha1\t{h}\n")
            fh.write("idx\tstate\tallele\n")
            for i, r in enumerate(rows):
                if r["state"] == "REF":
                    continue
                fh.write(f"{i}\t{r['state']}\t{r['allele']}\n")
    n = len(rows)
    print(f"  {a.sample}: {len(h_keys)} H37Rv keys, {n_lines} projection lines, "
          f"{n_multi} positions the reference visits more than once, "
          f"{n_unproj} with no projection")
    print(f"  {a.sample}: {n} keys  " + "  ".join(
        f"{k} {v} ({100*v/n:.1f}%)" for k, v in counts.items()))
    print(f"    reversions resolved to REF rather than ALT: {reversions}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
