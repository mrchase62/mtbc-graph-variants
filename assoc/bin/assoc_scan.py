#!/usr/bin/env python3
"""Convergence scan with a region-matched null beside the branch-length one.

WHY A SECOND NULL. phyoverlap2 tests an event set against branches drawn at
random in proportion to their length. That null is conditioned on the NUMBER of
events, which is right, but it assumes events of every kind fall on branches
the same way. Measured on scale200 they do not: PE/PPE is 14.9% of small-variant
sites and carries a 12.4% homoplasy rate against core's 2.3%, and 1,275 of the
2,459 homoplastic small variants -- 52% -- sit in it. Against a background drawn
from the whole tree, a variant from a region whose own background is five times
higher is being compared with the wrong population, and almost anything in
PE/PPE will look significant.

So each variant is tested twice:

  branch   phyoverlap2's null. k branches drawn in proportion to edge length.
  region   an empirical null. k branches drawn from the pooled event branches
           of the OTHER variants in the same region, leave-one-out, so
           whatever is region-specific about where those events land is
           already in the null rather than being discovered by the test.

Both p-values are reported. Where they disagree the region null is the one to
believe, and the disagreement itself is the useful output: it is the list of
hits that were artefacts of the background.

FALSE DISCOVERY IS CONTROLLED WITHIN REGION for the same reason -- a
Benjamini-Hochberg pass over a pooled list lets the easy region set the
threshold for the hard one.

EQUIVALENT, AND MUCH FASTER. `phyoverlap2.meanEventOverlap` recomputes a set
intersection per event per permutation. The statistic is the mean over the
event branches of (descendants carrying the label) / (descendants), and that
per-branch ratio does not depend on the variant, so it is computed once for the
whole tree and every permutation becomes the mean of k sampled values. The
scan below checks itself against phyoverlap2 on a sample of variants and
reports the agreement rather than asserting it.

ONE TEST PER PATTERN, NOT PER RECORD (2026-10-16). Records whose state rows
are identical at every node, with the same applicability, are one test unit:
tested once and counted once in BH, every member still written with the
unit's result (collapse_units; --no-collapse-identical turns it off). And
--exclude-acc-is6110 removes the accessory records that are IS6110 copies,
for the sensitivity scan the tail runs beside the full one (read_acc_is6110,
scan_compare.py). Both are recorded in <out stem>.run.tsv; unit membership
is in <out stem>.units.tsv.
"""
import argparse, collections, csv, glob, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sv_scatter import tree_frame
from write_event_matrix import fitch, UNK, A_BIT, D_BIT


def descendant_matrix(T):
    """(n_nodes x n_leaves) boolean: which leaves descend from each node."""
    n_n, leaves = T["n"], T["leaves"]
    D = np.zeros((n_n, len(leaves)), dtype=bool)
    for i in T["post"]:
        if T["is_leaf"][i]:
            D[i, T["leaf_index"][i]] = True
        else:
            for c in T["children"][i]:
                D[i] |= D[c]
    return D


def event_key(row_key):
    """The bytestream key for a variants.tsv row_key, EVERY field of it.

    write_event_matrix.py --dedupe suffix renames the second record at a
    repeated (pos, ref, alt) to (pos, ref, alt, '#2'). Rebuilding the key from
    the first three fields only looked the `#2` record up as the FIRST record
    and tested it on that record's branches -- svi:DEL:836371:1010 in gwas1000
    was scored on the 80 gains of svi:DEL:836371:135 instead of its own 32.
    """
    k = row_key.split("|")
    return (int(k[0]) if k[0].isdigit() else k[0],) + tuple(k[1:])


def leave_one_out(pool, own):
    """`pool` with ONE occurrence of each of `own` removed, as a multiset.

    The previous one-liner, `not (own[x] and own.__setitem__(...))`, kept every
    element -- `__setitem__` returns None, so the test was always true -- and
    the region null drew from a pool still holding the variant's own branches.
    """
    own = collections.Counter(own)
    cand = []
    for x in pool:
        if own[x]:
            own[x] -= 1
        else:
            cand.append(x)
    return cand


def load_lineages(path, leaves):
    """{leaf index: lineage} from a cohort table, refusing a table that cannot
    stratify this tree.

    A missing path, a table without `sample` and `lineage` columns, or one
    whose samples are not the tree's tips all used to leave every tip in one
    `unassigned` stratum, which quietly turns the lineage null into a plain
    permutation. H37Rv and the outgroup are tips with no cohort row, so a few
    unassigned tips are expected; more than MIN_LINEAGE_COVERAGE is refused.
    """
    if not os.path.exists(path):
        sys.exit(f"FATAL: --lineages {path} does not exist; the lineage null "
                 f"would be skipped and written blank")
    li = {l: j for j, l in enumerate(leaves)}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if not {"sample", "lineage"} <= set(rd.fieldnames or []):
            sys.exit(f"FATAL: --lineages {path} has no `sample` and `lineage` "
                     f"columns (has {rd.fieldnames})")
        lin = {}
        for r in rd:
            if r.get("sample") in li:
                lin[li[r["sample"]]] = r.get("lineage") or "unknown"
    cov = len(lin) / max(1, len(leaves))
    print(f"  lineage table: {len(lin):,} of {len(leaves):,} tips assigned "
          f"({cov:.1%})")
    if cov < MIN_LINEAGE_COVERAGE:
        sys.exit(f"FATAL: --lineages {path} assigns only {len(lin):,} of "
                 f"{len(leaves):,} tips; below {MIN_LINEAGE_COVERAGE:.0%} the "
                 f"lineage null is mostly one `unassigned` stratum")
    return lin


# share of tips a lineage table must assign (the rest are H37Rv, the outgroup)
MIN_LINEAGE_COVERAGE = 0.95


def load_carriers(path, li):
    """The phenotype carriers, refusing names that are not tips of the tree.

    A carrier missing from the tree used to be dropped without a word, which
    reads exactly like a non-carrier."""
    carriers = {l.strip() for l in open(path) if l.strip()}
    miss = sorted(s for s in carriers if s not in li)
    if miss:
        sys.exit(f"FATAL: {len(miss):,} of {len(carriers):,} phenotype "
                 f"carriers in {path} are not tips of the tree, e.g. "
                 f"{miss[:5]}")
    return carriers


# ---- ADAPTIVE PERMUTATIONS (review 2, R2-TREES-3 / ASSOC-11; the user's
# option b, 2026-10-06). A p-value is max(1, c) / P, so with P = 20,000 none can
# fall below 5e-5, and in a BH family as large as the small burden's 3,548
# genes a lone true hit could reach at best q = 0.177: it could never be
# significant. Every row is first tested at P; a null whose count c is at or
# below REFINE_BELOW -- p <= 5e-4, where the floor and Monte Carlo error decide
# significance -- is re-tested with REFINE_P fresh permutations, drawn in chunks
# so memory stays bounded, and that p-value replaces the first. BH then runs on
# the refined values. Nothing above the threshold changes.
REFINE_P = 1_000_000
REFINE_BELOW = 10
REFINE_CHUNK = 20_000


def needs_refine(p, P, below=REFINE_BELOW):
    return p == p and p * P <= below + 1e-9


def refine_branch(rng, ov, w, n, br, obs, P2, chunk=REFINE_CHUNK):
    """The branch null re-drawn P2 times: k branches with edge-length weights."""
    k, c, done = len(br), 0, 0
    while done < P2:
        m = min(chunk, P2 - done)
        d = rng.choice(n, size=(m, k), p=w)
        c += int((ov[d].mean(axis=1) >= obs).sum())
        done += m
    return max(1, c) / P2


def refine_region(rng, ov, pa, k, obs, P2, chunk=REFINE_CHUNK):
    """The region null re-drawn P2 times from the stratum's branch pool."""
    c, done = 0, 0
    while done < P2:
        m = min(chunk, P2 - done)
        d2 = pa[rng.integers(0, len(pa), size=(m, k))]
        c += int((ov[d2].mean(axis=1) >= obs).sum())
        done += m
    return max(1, c) / P2


def strat_perm_chunk(rng, strata, n_leaves, m):
    """(n_leaves, m) float32: m permutations of the phenotype, each holding the
    carrier count k fixed within every stratum (idx array, k)."""
    M = np.zeros((n_leaves, m), dtype=np.float32)
    for idx, k in strata:
        if k <= 0:
            continue
        if k >= len(idx):
            M[idx, :] = 1.0
            continue
        pick = np.argpartition(rng.random((m, len(idx))), k - 1, axis=1)[:, :k]
        M[idx[pick], np.arange(m)[:, None]] = 1.0
    return M


def refine_strat(rng, strata, n_leaves, items, P2, chunk=REFINE_CHUNK):
    """The lineage (or level-2) null re-drawn P2 times for several rows at
    once: each chunk of stratified permutations is built once and applied to
    every row. items: (Dsub float32 (k, n_leaves), denominators (k,), obs).
    Returns one p-value per item."""
    counts = [0] * len(items)
    done = 0
    while done < P2:
        m = min(chunk, P2 - done)
        M = strat_perm_chunk(rng, strata, n_leaves, m)
        for i, (Ds, den, obs) in enumerate(items):
            ovp = (Ds @ M) / np.maximum(den, 1)[:, None]
            counts[i] += int((ovp.mean(axis=0) >= obs).sum())
        done += m
    return [max(1, c) / P2 for c in counts]


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    run = 1.0
    for rank in range(n - 1, -1, -1):
        run = min(run, p[o[rank]] * n / (rank + 1))
        q[o[rank]] = run
    return q


# ---- IDENTICAL-PATTERN COLLAPSE (user, 2026-10-16: "add both options to the
# association step"). The same biological event can reach the scan as two or
# more records: an IS6110 insertion written once as CLASS=is6110 and again as
# the CLASS=accessory_presence record of the catalogue locus that IS the copy
# (INFO/ACCIS6110), a split SV event and a small-allele call of it, a catalogued
# deletion and its caller record. Each was tested and each was counted in the
# Benjamini-Hochberg family, so one event spent two or more tests and, where it
# was real, reported itself twice.
#
# WHAT "IDENTICAL" MEANS HERE, precisely. Two tested records form one test unit
# when
#   1. their rows of the event matrix's STATE stream (<events>/ancestor) are
#      identical at every node -- every leaf AND every internal node, each
#      0 ancestral, 1 derived or 2 unknown. The leaves are the carrier vector
#      after the writer's own missing handling (--absent, the reference and
#      outgroup tips, inapplicable isolates of a conditional record); the
#      internal nodes are its reconstruction. The event branches are a
#      function of the state row, so identical rows are identical tests.
#      Requiring the internal nodes too (not just the leaves) keeps apart two
#      records with the same carriers whose reconstructions differ -- the MTBC
#      pin applies only where AA resolved, so one can be pinned and the other
#      not -- because those are not the same test.
#   2. they have the same APPLICABILITY: both unconditional, or both
#      conditional on the same carrier set (the level-2 presence set of their
#      accessory locus; the locus id when no presence tables are given). An
#      unknown leaf that is unknown because the isolate lacks the insert is a
#      different missingness treatment from a no-call, and the level-2 null
#      reads it differently.
# Nothing else -- not class, not region, not position -- enters the key.
#
# WHAT THE UNIT DOES. It is tested once, on its REPRESENTATIVE's ladder of
# region pools and its level-2 carrier set, and it is ONE entry in the BH
# family of the representative's fine region. Its branches enter each region
# pool once, not once per member. Every member keeps its own row in the
# output, carrying the unit's p- and q-values, so no feature disappears and
# the chain audit still finds every tested record; `unit_rep` marks the one
# row per unit to count. The representative is the member of lowest class
# rank -- small, sv, is6110, accessory_presence -- then the first in
# variants.tsv order, so an IS6110 insertion and its ACCIS6110 catalogue copy
# are tested as the insertion. A unit whose members sit in different regions
# is flagged in its rows (`unit_regions`) because its region null is the
# representative's.
#
# WHAT IT DOES NOT DO. It does not merge features, change any record's
# carriers, or touch the merged VCF; hierarchical cases (an IS6110 insertion
# inside an accessory locus, whose carrier sets differ) are NOT identical and
# stay separate tests. --no-collapse-identical restores one test per record,
# bit for bit the scan before this option.
CLASS_RANK = {"small": 0, "sv": 1, "is6110": 2, "accessory_presence": 3}


def collapse_units(rows, sig_of):
    """Group `rows` into test units by `sig_of(row)` (a hashable pattern key,
    or None to keep the row alone). Returns a list of units, each a list of
    indices into `rows` with the representative first; units are ordered by
    their first member's index, so with no collapse the order is `rows`'s."""
    groups, order = {}, []
    for i, r in enumerate(rows):
        k = sig_of(r)
        if k is None:
            order.append([i])
            continue
        if k not in groups:
            groups[k] = []
            order.append(groups[k])
        groups[k].append(i)
    units = []
    for g in order:
        rep = min(g, key=lambda i: (CLASS_RANK.get(rows[i]["class"], 9), i))
        units.append([rep] + [i for i in g if i != rep])
    return units


# ---- THE IS6110 SENSITIVITY SET (same instruction). The catalogue loci whose
# sequence is an IS6110 copy are written by merge_cohort_vcf.py as
# accessory_presence records flagged ACCIS6110, beside the IS6110 arm's own
# insertion records at the same sites. How a copy should be counted against an
# insertion that sits ON an accessory locus is not decided (HANDOFF SV-EVAL,
# decision a), so the tail runs the scan twice -- once on everything, once
# with every flagged record removed -- and compares the two; an association
# that appears in only one of them is the one the decision would move.
def read_acc_is6110(path):
    """The IDs of the accessory_presence records flagged ACCIS6110.

    `path` is the merged VCF (.vcf or .vcf.gz), whose INFO/ACCIS6110 flag is
    read, or a text file of IDs one per line (`#` comments), for a cohort
    whose VCF predates the flag. A VCF without the flag's ##INFO header cannot
    say which loci are copies and is refused rather than read as `none`."""
    if not os.path.exists(path):
        sys.exit(f"FATAL: --exclude-acc-is6110 {path} does not exist")
    ids = set()
    if path.endswith((".vcf", ".vcf.gz")):
        import gzip
        op = gzip.open if path.endswith(".gz") else open
        header = False
        with op(path, "rt") as fh:
            for line in fh:
                if line.startswith("##"):
                    if line.startswith("##INFO=<ID=ACCIS6110,"):
                        header = True
                    continue
                if line.startswith("#"):
                    if not header:
                        sys.exit(f"FATAL: --exclude-acc-is6110 {path} has no "
                                 f"##INFO ACCIS6110 header; it predates the "
                                 f"flag and cannot say which accessory loci "
                                 f"are IS6110 copies")
                    continue
                f = line.rstrip("\n").split("\t", 8)
                info = f[7].split(";")
                if "ACCIS6110" in info and "CLASS=accessory_presence" in info:
                    ids.add(f[2])
        return ids
    for line in open(path):
        s = line.split("#", 1)[0].strip()
        if s:
            ids.add(s)
    return ids


SCAN_FIELDS = [
    "row_key", "id", "cls", "region", "region_pool", "region_level",
    "pos", "gains", "carriers",
    "obs", "obs_cond", "acc_locus", "cond_carriers",
    "cond_carriers_pheno", "p_cond", "q_cond",
    "p_branch", "q_branch", "p_region", "q_region",
    "p_lineage", "q_lineage", "null_pool", "refined",
    # the identical-pattern collapse: the unit this row was tested in, its
    # size, 1 on the one row per unit to count, the BH family it was tested
    # in, the members' fine regions where they differ, and its members' IDs
    "unit_id", "unit_size", "unit_rep", "unit_region", "unit_regions",
    "unit_classes", "unit_members"]


def survives(o):
    """THE SURVIVOR RULE (see the note at its use in main). `o` holds floats,
    or the strings a scan.tsv row carries ('nan' for a null not run)."""
    o = {k: (float(o[k]) if o.get(k) not in (None, "") else float("nan"))
         for k in ("q_lineage", "q_cond", "q_region", "q_branch")}
    lineage_ok = o["q_lineage"] == o["q_lineage"] and o["q_lineage"] < 0.05
    if o["q_cond"] == o["q_cond"]:
        return (o["q_cond"] < 0.05 and o["q_branch"] < 0.05
                and lineage_ok)
    return (o["q_region"] == o["q_region"] and o["q_region"] < 0.05
            and o["q_branch"] < 0.05 and lineage_ok)


def sha16(path):
    """16 hex of a file's sha256, read in blocks (a cohort VCF can be large)."""
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:16]


def side_path(out, tag):
    """<out stem>.<tag>.tsv beside a scan output."""
    stem = out[:-4] if out.endswith(".tsv") else out
    return f"{stem}.{tag}.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="assoc/scale200/events")
    ap.add_argument("--phenotype", required=True,
                    help="one sample name per line: the carriers")
    ap.add_argument("--min-gains", type=int, default=2)
    ap.add_argument("--permutations", type=int, default=20000)
    ap.add_argument("--refine-permutations", type=int, default=REFINE_P,
                    help="permutations for a null whose count is at or below "
                         "--refine-below after the first pass; 0: no refinement")
    ap.add_argument("--refine-below", type=int, default=REFINE_BELOW)
    ap.add_argument("--sv-intervals", default="",
                    help="the deletion catalogue carrying `evidence_tier` from "
                         "bin/retier_intervals.py. With it, a catalogued "
                         "deletion is nulled against other intervals of the "
                         "SAME evidence tier. Without it all 806 that reach "
                         "gwas1000 share one pool, of which 77% are scattered "
                         "repeat-context intervals, so a lineage-coherent "
                         "deletion is measured against a background of "
                         "unreliable ones. REQUIRED when svi: records are "
                         "scanned, and it must tier every one of them; "
                         "`none` pools them untiered on purpose.")
    ap.add_argument("--accessory-presence", default="",
                    help="directory of <sample>.presence.tsv. Enables the "
                         "LEVEL 2 null: for a variant inside an accessory "
                         "locus, the phenotype is permuted within the isolates "
                         "that CARRY that locus, and the per-branch statistic "
                         "is computed over carriers only. Without it a "
                         "conditional variant has no null of its own, because "
                         "its event branches are confined to the carrier clades "
                         "and no pool of unconditional branches is comparable.")
    ap.add_argument("--min-determinacy", type=float, default=0.80,
                    help="a variant whose reconstruction leaves more than this "
                         "fraction of branches undetermined is not something "
                         "the data can speak to, and is removed from the "
                         "DENOMINATOR rather than counted as a non-hit. A rate "
                         "then means `of the variants this data can speak to`. "
                         "The default is set by the positive control: rpoB "
                         "761155 has 0.863, so a threshold of 0.90 removes the "
                         "one variant known to be true by construction, which "
                         "is the clearest evidence a data-quality cut is too "
                         "aggressive.")
    ap.add_argument("--genes", default="data/annotation/H37Rv_snpeff_dump.txt",
                    help="gene intervals, so the null stratum can be crossed "
                         "with genic/intergenic. Measured on scale200 that "
                         "matters: within core, intergenic homoplasy is 4.2% "
                         "against 2.0% genic, and within the element mask 9.4% "
                         "against 2.7%. Comparing an intergenic variant with a "
                         "background dominated by genic ones is the same bias "
                         "the region null exists to remove, one level down.")
    ap.add_argument("--lineages", default="",
                    help="TSV with `sample` and `lineage` columns -- the "
                         "cohort table serves. Enables the third null, which "
                         "permutes the PHENOTYPE within lineage instead of "
                         "permuting branches. RRDR carriage is 57% in lineage 2 "
                         "and 22% in lineage 4, so a variant whose carriers "
                         "track lineage 2 looks enriched without being "
                         "associated; the region null cannot see that, because "
                         "it controls where events fall and not who carries "
                         "the phenotype. Measured: idsB and Rv1754c pass the "
                         "region null and collapse to the cohort baseline here.")
    ap.add_argument("--min-pool", type=int, default=200,
                    help="a region needs at least this many pooled event "
                         "branches before its own null is usable")
    ap.add_argument("--classes",
                    default="small,sv,is6110,accessory_presence",
                    help="record classes to scan. accessory_presence is the "
                         "LEVEL 1 character -- one biallelic call per accessory "
                         "locus, every isolate stated. It was added to the VCF "
                         "and the event matrix before this default was updated, "
                         "so gwas1000's first scan after the rebuild silently "
                         "dropped all 802 of them: 92 had two or more "
                         "independent gains and 65 passed the callability "
                         "floor, and none was tested. A class this list does "
                         "not name is invisible here, not reported as skipped.")
    ap.add_argument("--check", type=int, default=12,
                    help="variants to re-test with phyoverlap2 itself")
    ap.add_argument("--phyoverlap2", default="")
    ap.add_argument("--no-collapse-identical", action="store_true",
                    help="test every record on its own, as before 2026-10-16. "
                         "By default records whose state rows in "
                         "<events>/ancestor are identical at every node, with "
                         "the same applicability, are ONE test unit: tested "
                         "once and counted once in BH (see collapse_units). "
                         "Members are listed in <out stem>.units.tsv.")
    ap.add_argument("--exclude-acc-is6110", default="",
                    help="the IS6110 sensitivity run: remove every "
                         "accessory_presence record flagged ACCIS6110 before "
                         "anything else -- from the tests, the region pools "
                         "and the BH families. The merged VCF (its INFO flag "
                         "is read) or a file of IDs. The tail runs this as a "
                         "second scan beside the full one.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260925)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    run_rec = [("events", a.events), ("phenotype", a.phenotype),
               ("classes", a.classes), ("permutations", a.permutations),
               ("refine_permutations", a.refine_permutations),
               ("seed", a.seed),
               ("collapse_identical",
                "off (--no-collapse-identical)" if a.no_collapse_identical
                else "on")]

    T = tree_frame(os.path.join(a.events, "labelled.nwk"))
    li = {l: j for j, l in enumerate(T["leaves"])}
    D = descendant_matrix(T)
    nd = D.sum(axis=1)

    carriers = load_carriers(a.phenotype, li)
    dep = np.zeros(len(T["leaves"]), dtype=bool)
    for s in carriers:
        dep[li[s]] = True
    print(f"  tree {T['n']} nodes, {len(T['leaves'])} leaves; "
          f"{int(dep.sum())} carriers of the phenotype")

    # The per-branch statistic, once. ov[i] is the fraction of node i's
    # descendant leaves that carry the label.
    ov = (D & dep).sum(axis=1) / np.maximum(nd, 1)

    # ---- the lineage null: P permuted phenotypes, each keeping every
    # lineage's carrier count exactly as observed, so between-lineage
    # prevalence cannot generate a hit and only within-lineage association can.
    # OVP[:, p] is the per-branch statistic under permutation p, computed as
    # one matrix product rather than P passes over the tree.
    OVP = None
    if not a.lineages:
        print("  WARNING: no --lineages, so no lineage null; under the "
              "survivor rule nothing can survive without it")
    if a.lineages:
        lin = load_lineages(a.lineages, T["leaves"])
        groups = collections.defaultdict(list)
        for i in range(len(T["leaves"])):
            groups[lin.get(i, "unassigned")].append(i)
        DEP = np.zeros((len(T["leaves"]), a.permutations), dtype=np.float32)
        for g, idx in groups.items():
            idx = np.asarray(idx)
            k = int(dep[idx].sum())
            if k == 0:
                continue
            for p in range(a.permutations):
                DEP[rng.choice(idx, size=k, replace=False), p] = 1.0
        OVP = (D.astype(np.float32) @ DEP) / np.maximum(nd, 1)[:, None]
        eprint = print
        print(f"  lineage null: {len(groups)} strata, carrier counts held at "
              + ", ".join(f"{g}:{int(dep[np.asarray(ix)].sum())}/{len(ix)}"
                          for g, ix in sorted(groups.items()))[:160])

    # ---- the LEVEL 2 null --------------------------------------------------
    # A variant inside an accessory locus is applicable only to the isolates
    # carrying that locus. Two consequences, and the existing nulls handle
    # neither.
    #
    # First, the STATISTIC. ov[i] is the share of node i's descendant leaves
    # carrying the phenotype, counted over ALL descendants. For a conditional
    # variant that mixes isolates the variant cannot apply to into the
    # denominator, so a branch subtending both carriers and non-carriers has
    # its score diluted by isolates that were never at risk of the variant.
    # The conditional statistic restricts both numerator and denominator to
    # carriers of the locus.
    #
    # Second, the NULL. The region null draws from other variants' event
    # branches, and for a conditional variant there is no comparable pool: its
    # events can only land inside the carrier clades, so branches drawn from
    # unconditional variants come from a different amount of tree. Pooling the
    # conditional variants with each other reached only 41 branches on
    # gwas1000, below the floor. So instead of a pool this null permutes the
    # PHENOTYPE within the carrier set, holding the number of phenotype
    # carriers among them fixed -- the same construction as the lineage null,
    # with the locus's carrier set as the stratum instead of a lineage.
    #
    # Only the variant's own event branches are needed, not every node, so the
    # matrix product is over a handful of rows rather than all 1,997.
    cond_carriers = {}
    if a.accessory_presence:
        per = collections.defaultdict(set)
        pf = sorted(glob.glob(os.path.join(a.accessory_presence,
                                           "*.presence.tsv")))
        for f in pf:
            with open(f, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("state") == "PRESENT":
                        per[r["locus"]].add(r["sample"])
        for lid, ss in per.items():
            A = np.zeros(len(T["leaves"]), dtype=bool)
            for s2 in ss:
                if s2 in li:
                    A[li[s2]] = True
            if A.any():
                cond_carriers[lid] = A
        print(f"  level 2 null: {len(pf):,} presence tables, "
              f"{len(cond_carriers):,} loci with at least one carrier on the tree")

    # How many branches can EVER be determined for a variant at this locus.
    # A node is resolvable only if at least one of its descendant leaves is
    # applicable; everything else is inapplicable, not missing. This count is
    # what the callability floor has to be measured against.
    _nbr_cache = {}

    def applicable_branches(lid):
        if lid not in _nbr_cache:
            A = cond_carriers.get(lid)
            if A is None:
                _nbr_cache[lid] = None
            else:
                has = (D & A).sum(axis=1) > 0
                # exclude the root, which carries no branch above it
                _nbr_cache[lid] = int(has.sum()) - 1
        return _nbr_cache[lid]

    _dep_cache = {}

    def cond_null(lid, br, P):
        """(obs, p, n_carriers, n_pheno_carriers) for one conditional variant."""
        A = cond_carriers.get(lid)
        if A is None or not A.any():
            return None
        nA = int(A.sum())
        kA = int((dep & A).sum())
        DA = D[br]                              # (k, n_leaves)
        ndA = (DA & A).sum(axis=1)
        if not ndA.any():
            return None
        obsA = float(((DA & A & dep).sum(axis=1) / np.maximum(ndA, 1)).mean())
        if kA == 0 or kA == nA:
            # the phenotype is constant among carriers, so there is nothing to
            # permute and no association to detect either way
            return obsA, float("nan"), nA, kA
        key = (lid, P)
        if key not in _dep_cache:
            idx = np.flatnonzero(A)
            M = np.zeros((len(T["leaves"]), P), dtype=np.float32)
            for pp in range(P):
                M[rng.choice(idx, size=kA, replace=False), pp] = 1.0
            _dep_cache[key] = M
        M = _dep_cache[key]
        num = DA.astype(np.float32) @ M          # (k, P)
        ovp = num / np.maximum(ndA, 1)[:, None]
        pv = max(1, int((ovp.mean(axis=0) >= obsA).sum())) / P
        return obsA, pv, nA, kA

    # branch-length weights, as phyoverlap2 draws them; the root has none
    el = np.zeros(T["n"])
    for i, n in enumerate(T["nodes"]):
        el[i] = (n.edge.length or 0.0) if n.edge is not None else 0.0
    w = el / el.sum() if el.sum() > 0 else np.full(T["n"], 1 / T["n"])

    want = {c for c in a.classes.split(",") if c}
    everyrow = list(csv.DictReader(
        open(os.path.join(a.events, "variants.tsv")), delimiter="\t"))
    allrows = [r for r in everyrow if r["class"] in want]
    if a.exclude_acc_is6110:
        # Every listed ID must be an accessory_presence record of THIS event
        # matrix. An ID it does not hold means the list (or VCF) is from
        # another cohort or build, and an ID of another class means the list
        # is not what the option says -- either way a silent partial
        # exclusion would be read as the sensitivity result.
        xid = read_acc_is6110(a.exclude_acc_is6110)
        cls_of = {}
        for r in everyrow:
            cls_of.setdefault(r["id"], set()).add(r["class"])
        absent = sorted(x for x in xid if x not in cls_of)
        wrong = sorted(x for x in xid if x in cls_of
                       and cls_of[x] != {"accessory_presence"})
        if absent:
            sys.exit(f"FATAL: {len(absent):,} of {len(xid):,} ACCIS6110 IDs "
                     f"from {a.exclude_acc_is6110} are not in "
                     f"{a.events}/variants.tsv, e.g. {absent[:3]}; the list "
                     f"is not from this cohort's VCF")
        if wrong:
            sys.exit(f"FATAL: {len(wrong):,} IDs from {a.exclude_acc_is6110} "
                     f"are not accessory_presence records, e.g. {wrong[:3]}")
        n0 = len(allrows)
        allrows = [r for r in allrows if r["id"] not in xid]
        print(f"  IS6110 sensitivity: {len(xid):,} accessory_presence records "
              f"flagged ACCIS6110 excluded ({n0 - len(allrows):,} of them in "
              f"the scanned classes)")
        run_rec += [("exclude_acc_is6110", a.exclude_acc_is6110),
                    ("exclude_acc_is6110_sha16",
                     sha16(a.exclude_acc_is6110)),
                    ("acc_is6110_ids", len(xid)),
                    ("acc_is6110_rows_excluded", n0 - len(allrows))]
    else:
        run_rec.append(("exclude_acc_is6110", "none"))
    nbranch = T["n"] - 1
    for r in allrows:
        # THE CALLABILITY FLOOR HAS TO BE RELATIVE TO WHAT IS APPLICABLE.
        # For a variant inside an accessory locus, every branch outside the
        # carrier clades is undetermined BY CONSTRUCTION -- a non-carrier leaf
        # is inapplicable, so Fitch leaves it unknown and every branch above it
        # unresolved. Measured against the whole tree, the undetermined
        # fraction therefore just reports how rare the insert is:
        #
        #     applicable isolates   median undetermined branches
        #     < 100                 1,996 of 1,996   100%
        #     100-299               1,341             67%
        #     300-599                 942             47%
        #     600-899                 331             17%
        #     >= 900                  101              5%
        #
        # So an 80% floor on the whole tree admits only loci carried by roughly
        # two thirds of the cohort or more, and every restricted locus is
        # excluded before any null can see it. On gwas1000 that silently
        # discarded 6 variants with two or more independent origins, and it is
        # why the only loci reaching the scan were carried by 87% and 94% of
        # isolates -- which is exactly the regime where conditioning changes
        # nothing. The floor is now measured against the branches that can ever
        # be determined for this locus.
        und = int(r["n_undet"])
        nb_r = nbranch
        lid = r.get("acc_locus") or ""
        if lid and cond_carriers:
            nba = applicable_branches(lid)
            if nba and nba > 0:
                # branches outside the carrier clades are inapplicable; they
                # are neither in the numerator nor the denominator
                inapplicable = nbranch - nba
                und = max(0, und - inapplicable)
                nb_r = nba
        r["_det"] = 1.0 - und / max(1, nb_r)
        r["_det_branches"] = nb_r
    rows = [r for r in allrows if int(r["n_gain"]) >= a.min_gains]
    # THE DENOMINATOR. A variant whose branches are mostly undetermined has not
    # been tested and found negative -- it has not been tested. Dropping those
    # from the denominator is what makes a rate honest; leaving them in inflates
    # the count of things looked at and deflates every rate computed from it.
    drop = [r for r in rows if r["_det"] < a.min_determinacy]
    rows = [r for r in rows if r["_det"] >= a.min_determinacy]
    print(f"  {len(rows) + len(drop):,} variants with >= {a.min_gains} "
          f"independent gains")
    print(f"  callable denominator: {len(rows):,} at >= "
          f"{a.min_determinacy:.0%} of branches determined; "
          f"{len(drop):,} dropped as uncallable")
    if drop:
        byc = collections.Counter(r["class"] for r in drop)
        print(f"    dropped by class: " +
              ", ".join(f"{k} {v:,}" for k, v in sorted(byc.items())))

    # Event branches per variant, recomputed from the states so the scan does
    # not depend on a separate bytestream staying in step with variants.tsv.
    import gzip
    ev_of = {}
    lab = T["labels"]
    # variants.tsv carries the counts; the branches come from the event matrix
    from mtbvartools import CallBytestream
    ev = CallBytestream(os.path.join(a.events, "event"))
    cols = np.asarray(ev.calls.col)
    col_i = {c: i for i, c in enumerate(cols)}
    keep, missing = [], []
    for r in rows:
        try:
            v = ev.calls.loc[event_key(r["row_key"])]
        except Exception:
            missing.append(r["row_key"])
            continue
        br = np.where(v == 1)[0]
        if len(br) < a.min_gains:
            continue
        ev_of[r["row_key"]] = np.asarray([col_i[cols[j]] for j in br])
        keep.append(r)
    ev.close()
    # A ROW THE EVENT MATRIX DOES NOT HOLD is a mismatch between the two files,
    # not an untestable variant; it used to be skipped without a count.
    if missing:
        sys.exit(f"FATAL: {len(missing):,} variants.tsv rows have no event-"
                 f"matrix row, e.g. {missing[:3]}")
    rows = keep
    print(f"  {len(rows):,} of them resolved to event branches")

    # ---- the identical-pattern collapse (collapse_units, above) -----------
    # The state rows come from the writer's STATE stream, in the tree's node
    # order (looked up by label, not assumed to share the event stream's).
    if a.no_collapse_identical:
        units = [[i] for i in range(len(rows))]
        n_leaf_units = n_gain_units = len(rows)
    else:
        anc_dir = os.path.join(a.events, "ancestor")
        if not os.path.isdir(anc_dir):
            sys.exit(f"FATAL: no {anc_dir}; the identical-pattern collapse "
                     f"reads each record's state row from it. Pass "
                     f"--no-collapse-identical to test every record alone")
        an = CallBytestream(anc_dir)
        acol = {c: i for i, c in enumerate(np.asarray(an.calls.col))}
        nolab = [l for l in T["labels"] if l not in acol]
        if nolab:
            sys.exit(f"FATAL: {len(nolab):,} tree nodes have no column in "
                     f"{anc_dir}, e.g. {nolab[:3]}; it is not this tree's "
                     f"state matrix")
        node_cols = np.asarray([acol[l] for l in T["labels"]])
        leaf_mask = np.asarray(T["is_leaf"], dtype=bool)
        _app_sig = {}

        def applicability(r):
            lid = r.get("acc_locus") or ""
            if not lid:
                return b"unconditional"
            if lid not in _app_sig:
                A = cond_carriers.get(lid)
                _app_sig[lid] = (b"carriers:" + np.packbits(A).tobytes()
                                 if A is not None else b"locus:" + lid.encode())
            return _app_sig[lid]

        sig_full, sig_leaf, sig_gain = {}, {}, {}
        for r in rows:
            try:
                st = np.asarray(an.calls.loc[event_key(r["row_key"])],
                                dtype=np.uint8)[node_cols]
            except Exception:
                sys.exit(f"FATAL: {r['row_key']} has an event row but no state "
                         f"row in {anc_dir}; the two streams are not one "
                         f"matrix")
            app = applicability(r)
            sig_full[r["row_key"]] = (st.tobytes(), app)
            # two diagnostics, NOT used to collapse: how many units a looser
            # key would give (leaves only; the gain branches only)
            sig_leaf[r["row_key"]] = (st[leaf_mask].tobytes(), app)
            sig_gain[r["row_key"]] = (np.sort(ev_of[r["row_key"]]).tobytes(),
                                      app)
        an.close()
        units = collapse_units(rows, lambda r: sig_full[r["row_key"]])
        n_leaf_units = len(set(sig_leaf.values()))
        n_gain_units = len(set(sig_gain.values()))
    multi = [u for u in units if len(u) > 1]
    print(f"  identical-pattern collapse: "
          + ("OFF (--no-collapse-identical); " if a.no_collapse_identical
             else "")
          + f"{len(rows):,} records -> {len(units):,} test units "
          f"({len(rows) - len(units):,} tests removed; {len(multi):,} units "
          f"of >= 2 records)")
    if multi:
        byc = collections.Counter(
            "+".join(sorted({rows[i]["class"] for i in u})) for u in multi)
        print("    multi-record units by classes: " + ", ".join(
            f"{k} {v:,}" for k, v in byc.most_common()))
    if not a.no_collapse_identical:
        print(f"    (diagnostic, not applied: keyed on the leaves alone "
              f"{n_leaf_units:,} units; on the gain branches alone "
              f"{n_gain_units:,})")
    run_rec += [("records_tested", len(rows)), ("test_units", len(units)),
                ("tests_removed_by_collapse", len(rows) - len(units)),
                ("multi_record_units", len(multi)),
                ("units_if_keyed_on_leaves_only", n_leaf_units),
                ("units_if_keyed_on_gain_branches_only", n_gain_units)]
    reps = [rows[u[0]] for u in units]

    import bisect as _bi
    _g = []
    if a.genes and not os.path.exists(a.genes):
        sys.exit(f"FATAL: --genes {a.genes} does not exist; without it the "
                 f"genic/intergenic split is silently dropped. Pass --genes '' "
                 f"to run without it on purpose")
    if a.genes:
        for row in csv.reader(open(a.genes), delimiter="\t"):
            if len(row) > 7 and row[4] == "Gene":
                try:
                    # the snpEff dump's start is 0-based, its end 1-based;
                    # positions here are 1-based VCF positions
                    _g.append((int(row[1]) + 1, int(row[2])))
                except ValueError:
                    pass
        _g.sort()
    _gs = [x[0] for x in _g]

    _gmax = max((e - s for s, e in _g), default=0)

    def _genic(r):
        """genic when the bases the record changes touch a gene, or, for an
        insertion, a gene holds both its flanks (variant_span; D18, D39).
        The anchor base alone called an MNP crossing into a gene intergenic,
        and an insertion just after a gene's last base genic."""
        from variant_span import changed_span
        first, last, inside = changed_span(r, r.get("class") or "small")
        lo = _bi.bisect_left(_gs, first - _gmax)
        hi = _bi.bisect_right(_gs, last)
        if inside:
            return any(s <= first and last <= e for s, e in _g[lo:hi])
        return any(s <= last and first <= e for s, e in _g[lo:hi])

    # evidence tier per catalogued interval, for the SV region key
    ev_tier = {}
    # A CATALOGUED DELETION WITHOUT A TIER is nulled against every other one,
    # and the chain ran that way from 2026-10-02 without a word: it stopped
    # passing this option, and the only table on disk was from the previous
    # catalogue and held 292 of scale200's 467 current IDs. So a scan holding
    # svi: records now needs the table, built from the same catalogue, or an
    # explicit `--sv-intervals none` to run untiered on purpose.
    n_svi = sum(1 for r in rows if r["id"].startswith("svi:"))
    if a.sv_intervals != "none":
        if n_svi and not a.sv_intervals:
            sys.exit(f"FATAL: {n_svi:,} svi: records to scan and no "
                     f"--sv-intervals; pass the retiered catalogue (or "
                     f"`--sv-intervals none` to pool them untiered)")
        if a.sv_intervals and not os.path.exists(a.sv_intervals):
            sys.exit(f"FATAL: --sv-intervals {a.sv_intervals} does not exist")
    if a.sv_intervals and a.sv_intervals != "none":
        with open(a.sv_intervals, newline="") as fh:
            for q in csv.DictReader(fh, delimiter="\t"):
                if q.get("evidence_tier"):
                    ev_tier[q["interval"]] = q["evidence_tier"]
        print(f"  sv evidence tiers: {len(ev_tier):,} intervals, "
              + "  ".join(f"{k} {v:,}" for k, v in
                          sorted(collections.Counter(ev_tier.values()).items())))
        untiered = sorted({r["id"].split("#")[0] for r in rows
                           if r["id"].startswith("svi:")
                           and r["id"].split("#")[0] not in ev_tier})
        if untiered:
            sys.exit(f"FATAL: {len(untiered):,} scanned svi: IDs have no "
                     f"evidence_tier in {a.sv_intervals}, e.g. {untiered[:3]}; "
                     f"the table is not from this cohort's catalogue")

    def regkey(r):
        # LEVEL 1 is its own region. An accessory presence character is one
        # biallelic call per locus with every sample stated, and its background
        # rate of repeated gain has nothing to do with the small variants that
        # happen to share the coordinate -- pooling it into "other" would null
        # it against a class it has no relation to.
        if r["class"] == "accessory_presence":
            return "accessory_presence"
        # A CATALOGUED DELETION IS KEYED BY ITS EVIDENCE TIER, not lumped into
        # `other` with everything else. The tier is measured -- depth support
        # and the lineage coherence of the deleted isolate set -- and it
        # separates sharply on the RD validation set: of known regions of
        # difference with a checkable label, the coherent tier is concordant
        # for 14 of 37 against 1 of 27 for the rest. Nulling a coherent
        # interval against scattered repeat-context ones compares it to a
        # background it has nothing in common with.
        if ev_tier and r["id"].startswith("svi:"):
            t = ev_tier.get(r["id"].split("#")[0])
            if t:
                return f"svi:{t}"
        base = r["region"] or ("is6110" if r["class"] == "is6110" else "other")
        # LEVEL 2 is conditional, so its pool must be too. A variant applicable
        # to 40 carriers has far fewer branches available to it than one
        # applicable to 900, and pooling the two would compare an observed gain
        # count against a background drawn from a different amount of tree. The
        # bin is coarse on purpose -- it only has to stop the two extremes
        # sharing a null, and a finer split would starve the --min-pool floor.
        na = r.get("n_applicable") or ""
        if r.get("acc_locus") and na.isdigit():
            n = int(na)
            b = "lt50" if n < 50 else ("lt200" if n < 200 else "ge200")
            return f"{base}:cond_{b}"
        if not _g or r["frame"] != "h37rv" or not r["pos"].isdigit():
            return base
        return f"{base}:{'genic' if _genic(r) else 'intergenic'}"

    # ---- THE FALLBACK LADDER ----------------------------------------------
    # The fine key above is the right pool when there is enough of it, and on
    # the first gwas1000 and scale200 runs there often was not: is6110 split
    # genic/intergenic into 109 and 145 branches, each below the 200 floor,
    # so 108 variants got no region null at all while their union of 254 would
    # have cleared it comfortably. Seven of twelve strata on gwas1000 and nine
    # of twelve on scale200 were in that state.
    #
    # So each variant now gets a LADDER of progressively coarser pools and uses
    # the finest one that clears the floor, with the level it landed on recorded
    # in the output. The rungs are chosen so that every one of them is still a
    # defensible comparison:
    #
    #   fine     the key above: region, split by genic/intergenic or by the
    #            conditional applicability bin.
    #   coarse   the same region with those splits dropped. This pools the
    #            accessory loci with each other across applicability bins, and
    #            genic with intergenic. It is the rung the user asked for.
    #   family   node-frame records together, element records together. A
    #            node-frame variant's events land on the graph's own structure
    #            rather than on H37Rv coordinates, which is the property the
    #            region null is trying to hold constant.
    #
    # What the ladder deliberately does NOT do is pool a CONDITIONAL variant
    # with an unconditional one. A level-2 variant is applicable only to the
    # isolates carrying its insert, so its events can only land on that
    # carrier clade, while a level-1 presence character's can land anywhere on
    # the tree. Drawing a level-2 null from level-1 branches would compare an
    # observed count against a background from a different amount of tree --
    # the same error the applicability bins were introduced to prevent. Where
    # the conditional pool cannot reach the floor on its own the honest answer
    # is still no region null, and the run says so per stratum.
    def ladder(r):
        fine = regkey(r)
        cond = bool(r.get("acc_locus") and (r.get("n_applicable") or "").isdigit())
        base = fine.split(":")[0]
        rungs = [fine]
        if cond:
            # conditional variants pool only with other conditional variants
            rungs.append(f"{base}:cond_any")
            rungs.append("conditional_any")
        elif fine.startswith("svi:"):
            rungs.append("svi:any")
            rungs.append("other")
        else:
            if fine != base:
                rungs.append(base)
            fam = ("nodeframe" if base.startswith("off_path")
                   else ("element" if base == "is6110"
                         else ("accessory_l1" if base == "accessory_presence"
                               else base)))
            if fam != base:
                rungs.append(fam)
        out2 = []
        for x in rungs:
            if x not in out2:
                out2.append(x)
        return out2

    # One entry per TEST UNIT: a unit's branches enter each pool once, through
    # its representative's ladder (with the collapse off, one per record).
    pool = collections.defaultdict(list)
    n_var = collections.Counter()
    for r in reps:
        ev = ev_of[r["row_key"]].tolist()
        for rung in ladder(r):
            pool[rung].extend(ev)
            n_var[rung] += 1
    print(f"    {'stratum':<34}{'branches':>10}{'variants':>10}  floor {a.min_pool}")
    for k in sorted(pool):
        print(f"    {k:<34}{len(pool[k]):>10,}{n_var[k]:>10,}"
              f"  {'ok' if len(pool[k]) >= a.min_pool else 'below'}")

    P = a.permutations
    out, aux = [], []
    for r in reps:
        br = ev_of[r["row_key"]]
        k = len(br)
        obs = float(ov[br].mean())
        d = rng.choice(T["n"], size=(P, k), p=w)
        p_branch = max(1, int((ov[d].mean(axis=1) >= obs).sum())) / P
        # Walk the ladder and use the finest rung that clears the floor AFTER
        # leave-one-out, since a stratum holding one variant is all its own
        # branches and would otherwise pass the floor on paper.
        rungs = ladder(r)
        rk, pl, rlevel = rungs[0], [], ""
        for i, rung in enumerate(rungs):
            cand = leave_one_out(pool[rung], br.tolist())
            rk, pl = rung, cand
            if len(cand) >= a.min_pool:
                rlevel = ("fine", "coarse", "family")[min(i, 2)]
                break
        if len(pl) >= a.min_pool:
            pa = np.asarray(pl)
            d2 = pa[rng.integers(0, len(pa), size=(P, k))]
            p_region = max(1, int((ov[d2].mean(axis=1) >= obs).sum())) / P
            nulln = len(pl)
        else:
            p_region, nulln, rlevel = float("nan"), len(pl), "none"
        r["_region_pool"] = rk
        r["_region_level"] = rlevel
        # the level-2 null, for a variant that sits inside an accessory locus
        obs_cond = p_cond = float("nan")
        n_carr = n_carr_pheno = ""
        if cond_carriers and r.get("acc_locus"):
            got = cond_null(r["acc_locus"], br, P)
            if got is not None:
                obs_cond, p_cond, n_carr, n_carr_pheno = got
        if OVP is not None:
            p_lineage = max(1, int((OVP[br].mean(axis=0) >= obs).sum())) / P
        else:
            p_lineage = float("nan")
        # `region` remains the FINE key: it is the interpretive stratum and the
        # one BH groups on, and changing that would move the survivor set for
        # reasons unrelated to pooling. Which pool the region null was actually
        # drawn from is recorded beside it.
        out.append(dict(row_key=r["row_key"], id=r["id"], cls=r["class"],
                        region=rungs[0], pos=r["pos"], gains=k,
                        carriers=r["n_derived_leaves"], obs=f"{obs:.4f}",
                        p_branch=p_branch, p_region=p_region,
                        p_lineage=p_lineage, null_pool=nulln,
                        region_pool=rk, region_level=rlevel,
                        acc_locus=r.get("acc_locus", ""),
                        obs_cond=(f"{obs_cond:.4f}"
                                  if obs_cond == obs_cond else ""),
                        p_cond=p_cond,
                        cond_carriers=n_carr,
                        cond_carriers_pheno=n_carr_pheno))
        aux.append((br, np.asarray(pl) if len(pl) >= a.min_pool else None, obs))

    # adaptive permutations: re-test, with more draws, the nulls at the floor
    for o in out:
        o["refined"] = ""
    if a.refine_permutations:
        P2 = a.refine_permutations
        lin_items, lin_rows = [], []
        cond_jobs = collections.defaultdict(list)
        for o, (br, pa, obs) in zip(out, aux):
            done = []
            if needs_refine(o["p_branch"], P, a.refine_below):
                o["p_branch"] = refine_branch(rng, ov, w, T["n"], br, obs, P2)
                done.append("branch")
            if pa is not None and needs_refine(o["p_region"], P, a.refine_below):
                o["p_region"] = refine_region(rng, ov, pa, len(br), obs, P2)
                done.append("region")
            if OVP is not None and needs_refine(o["p_lineage"], P, a.refine_below):
                lin_items.append((D[br].astype(np.float32), nd[br], obs))
                lin_rows.append(o)
                done.append("lineage")
            if o["acc_locus"] and needs_refine(o["p_cond"], P, a.refine_below):
                cond_jobs[o["acc_locus"]].append((o, br))
                done.append("cond")
            o["refined"] = ",".join(done)
        if lin_items:
            strata = [(np.asarray(ix), int(dep[np.asarray(ix)].sum()))
                      for ix in groups.values()]
            for o, pv in zip(lin_rows, refine_strat(
                    rng, strata, len(T["leaves"]), lin_items, P2)):
                o["p_lineage"] = pv
        for lid, lst in cond_jobs.items():
            A = cond_carriers[lid]
            strata = [(np.flatnonzero(A), int((dep & A).sum()))]
            items = [(D[br].astype(np.float32), (D[br] & A).sum(axis=1),
                      float(o["obs_cond"])) for o, br in lst]
            for (o, _), pv in zip(lst, refine_strat(
                    rng, strata, len(T["leaves"]), items, P2)):
                o["p_cond"] = pv
        nref = collections.Counter(x for o in out for x in o["refined"].split(",") if x)
        print(f"  refined with {P2:,} permutations (count <= {a.refine_below} at "
              f"{P:,}): " + (", ".join(f"{k} {v:,}" for k, v in sorted(nref.items()))
                             or "none"))

    # BH within region, on each null separately
    for field, q in (("p_branch", "q_branch"), ("p_region", "q_region"),
                     ("p_lineage", "q_lineage"), ("p_cond", "q_cond")):
        for rk in {o["region"] for o in out}:
            grp = [o for o in out if o["region"] == rk
                   and not (isinstance(o[field], float) and
                            np.isnan(o[field]))]
            if not grp:
                continue
            qq = bh([o[field] for o in grp])
            for o, v in zip(grp, qq):
                o[q] = float(v)
    for o in out:
        o.setdefault("q_branch", float("nan"))
        o.setdefault("q_region", float("nan"))
        o.setdefault("q_lineage", float("nan"))
        o.setdefault("q_cond", float("nan"))

    # ---- ONE ROW PER RECORD, THE UNIT'S RESULT ON EVERY MEMBER ROW. `out`
    # holds one entry per test unit, in `units` order, and everything above --
    # the tests, the pools, BH -- saw only those. Each member gets a copy of
    # its unit's result under its own row_key, id, class, fine region,
    # position and locus, so every tested feature is still in the file (and
    # the chain audit's per-stratum reconciliation still holds); `unit_rep`
    # is 1 on the row to count. A member's own fine region can differ from
    # the unit's, whose BH family it was tested in (`unit_region`).
    members_out = []
    for k, (o, u) in enumerate(zip(out, units)):
        o["_unit"] = k
        rep = rows[u[0]]
        uid = f"u:{rep['row_key']}"
        regs = list(dict.fromkeys(regkey(rows[i]) for i in u))
        o.update(unit_id=uid, unit_size=len(u), unit_rep=1,
                 unit_region=o["region"],
                 unit_regions=",".join(regs) if len(regs) > 1 else "",
                 unit_classes=",".join(sorted({rows[i]["class"] for i in u})),
                 unit_members=";".join(rows[i]["id"] for i in u))
        for i in u[1:]:
            r = rows[i]
            m = dict(o)
            m.update(row_key=r["row_key"], id=r["id"], cls=r["class"],
                     region=regkey(r), pos=r["pos"],
                     carriers=r["n_derived_leaves"],
                     acc_locus=r.get("acc_locus", ""), unit_rep=0)
            members_out.append(m)
    allout = out + members_out
    allout.sort(key=lambda o: (o["q_region"] if o["q_region"] == o["q_region"]
                               else 9, o["p_branch"], o["_unit"],
                               -o["unit_rep"]))
    out.sort(key=lambda o: (o["q_region"] if o["q_region"] == o["q_region"]
                            else 9, o["p_branch"]))
    with open(a.out, "w", newline="") as fh:
        w2 = csv.DictWriter(fh, delimiter="\t", extrasaction="ignore",
                            fieldnames=SCAN_FIELDS)
        w2.writeheader()
        w2.writerows(allout)
    # the membership table: one line per record, every unit
    with open(side_path(a.out, "units"), "w") as fh:
        fh.write("unit_id\tunit_size\tis_rep\trow_key\tid\tclass\tregion\t"
                 "acc_locus\tunit_region\n")
        for u in units:
            uid = f"u:{rows[u[0]]['row_key']}"
            ureg = regkey(rows[u[0]])
            for j, i in enumerate(u):
                r = rows[i]
                fh.write(f"{uid}\t{len(u)}\t{int(j == 0)}\t{r['row_key']}\t"
                         f"{r['id']}\t{r['class']}\t{regkey(r)}\t"
                         f"{r.get('acc_locus', '')}\t{ureg}\n")
    n_mixed = sum(1 for o in out if o["unit_regions"])
    run_rec += [("units_spanning_regions", n_mixed),
                ("survivors_units", sum(1 for o in out if survives(o))),
                ("survivors_rows", sum(1 for o in allout if survives(o)))]
    with open(side_path(a.out, "run"), "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in run_rec:
            fh.write(f"{k}\t{v}\n")
    print(f"  -> {side_path(a.out, 'units')} (unit membership), "
          f"{side_path(a.out, 'run')} (options and counts)")
    if n_mixed:
        print(f"  {n_mixed:,} units hold records of more than one region; "
              f"each was tested on its representative's region null")

    print(f"\n  significant at q < 0.05, by region and by which null:")
    print(f"  {'region':<34}{'tested':>8}{'branch':>9}{'region':>9}"
          f"{'lineage':>9}{'all three':>11}  null drawn from")
    for rk in sorted({o["region"] for o in out}):
        g = [o for o in out if o["region"] == rk]
        b = sum(1 for o in g if o["q_branch"] < 0.05)
        # A region whose pool is too small has NO region null, which is a
        # different statement from having no hits. Reporting both as 0 would
        # be the same conflation this whole scan exists to avoid.
        usable = [o for o in g if o["q_region"] == o["q_region"]]
        # Which rung of the ladder these variants actually used, so a coarsened
        # null is never mistaken for a fine one.
        lv = collections.Counter((o["region_level"], o["region_pool"]) for o in g)
        via = "  ".join(f"{pl} ({lev}, n={len(pool.get(pl, [])):,})"
                        for (lev, pl), _ in lv.most_common(2))
        if not usable:
            print(f"  {rk:<34}{len(g):>8,}{b:>9,}{'no null':>9}"
                  f"{'-':>9}{'-':>11}  every rung below {a.min_pool}: "
                  f"{', '.join(f'{x}={len(pool.get(x,[])):,}' for x in dict.fromkeys(o['region_pool'] for o in g))}")
            continue
        rr = sum(1 for o in usable if o["q_region"] < 0.05)
        ln = sum(1 for o in g if o["q_lineage"] == o["q_lineage"]
                 and o["q_lineage"] < 0.05)
        both = sum(1 for o in g if o["q_region"] == o["q_region"]
                   and o["q_region"] < 0.05 and o["q_lineage"] == o["q_lineage"]
                   and o["q_lineage"] < 0.05 and o["q_branch"] < 0.05)
        print(f"  {rk:<34}{len(g):>8,}{b:>9,}{rr:>9,}{ln:>9,}{both:>11,}"
              f"  {via}")
    # THE SURVIVOR RULE. For a conditional variant the level-2 null stands in
    # for the region null it cannot have: q_cond is the test that its event
    # branches beat a phenotype shuffled within its own carrier set. For every
    # other variant the rule is unchanged.
    # The lineage null applies to a conditional variant exactly as it does to
    # any other: it asks whether the association is within-lineage rather than
    # a lineage being both common and phenotype-rich, and nothing about
    # conditioning on carriers makes that question go away. The first version
    # of this rule substituted q_cond for q_region and dropped the lineage
    # requirement with it, which let three scale200 variants through on
    # q_lineage 0.0798 -- and they were the same two-origin haplotype at
    # ACC_3846790 that gwas1000 rejects at q_branch 0.31.
    #
    # SURVIVING MEANS q < 0.05 UNDER ALL THREE NULLS -- branch, region (or
    # level 2 for a conditional variant) and lineage -- as HANDOFF 0b defines
    # it and as the `all three` column above counts it. The rule used to omit
    # q_branch for unconditional variants and to count a lineage null that was
    # never run (q_lineage blank) as a pass, so it matched `all three` only
    # while every row had a lineage null and every region pass was also a
    # branch pass. A null that did not run is now `not tested`, never a pass.
    # The rule is survives(), at module level, so scan_compare.py applies the
    # same one. Counted in TEST UNITS: a unit of several records is one
    # association, however many rows carry it.
    surv = [o for o in out if survives(o)]
    if surv:
        print(f"\n  the {len(surv)} that survive all three nulls:")
        print(f"  {'pos':>9} {'region':<16}{'gains':>6}{'carr':>6}{'obs':>7}"
              f"{'p_branch':>10}{'p_region':>10}{'p_lineage':>11}")
        for o in surv[:20]:
            pl = (f"{o['p_lineage']:.5f}"
                  if o["p_lineage"] == o["p_lineage"] else "-")
            print(f"  {o['pos']:>9} {o['region']:<16}{o['gains']:>6}"
                  f"{o['carriers']:>6}{float(o['obs']):>7.2f}"
                  f"{o['p_branch']:>10.5f}{o['p_region']:>10.5f}{pl:>11}")
    print(f"  -> {a.out}")

    # ---- self-check against phyoverlap2 itself
    if not (a.check and a.phyoverlap2 and os.path.isdir(a.phyoverlap2)):
        print("  self-check against phyoverlap2 NOT run (no --phyoverlap2)")
    if a.check and a.phyoverlap2 and os.path.isdir(a.phyoverlap2):
        sys.path.insert(0, a.phyoverlap2)
        import phyoverlap2 as po
        tree_dict = po.getChildNodeDict(
            tree_path=os.path.join(a.events, "labelled.nwk"))
        dep_set = {T["leaves"][i] for i in np.where(dep)[0]}
        idx = rng.choice(len(rows), size=min(a.check, len(rows)),
                         replace=False)
        worst = 0.0
        for i in idx:
            r = rows[int(i)]
            labs = [lab[j] for j in ev_of[r["row_key"]]]
            theirs = po.meanEventOverlap(labs, dep_set, tree_dict)["obs_overlap"]
            mine = float(ov[ev_of[r["row_key"]]].mean())
            worst = max(worst, abs(theirs - mine))
        print(f"  self-check against phyoverlap2 on {len(idx)} variants: "
              f"largest difference in the observed statistic {worst:.2e}")


if __name__ == "__main__":
    main()
