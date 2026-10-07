# Review 2: trees, ancestral states and association (`trees_ancestral_assoc`)

2026-10-06. Branch `audit-fixes` at ffd7bb9. This review was read-only. Nothing
was written to the working tree or the repository except this file. No cluster
jobs were run.

**Scratch:** `/tmp/claude-12043/review2_trees/`:

- `src/` is a `git archive audit-fixes` export;
- `old/` is `git archive main`, which is byte-identical to the working-tree
  copies that production ran;
- `rr/` and `rr2/` are scratch runroots;
- `build/` is a scratch P0 build;
- `toy/` is the hand-checked case;
- `scripts/` holds the independent re-derivations;
- `runs/` holds the logs.

**Summary.** I found no new bug in the fixed code that changes a result.
Verified by hand:

- the Fitch reconstruction, root placement, outgroup alleles, polarity and the
  gain/loss counts;
- the fixed `ancestral_alleles.py` against check 2's column M;
- the leave-one-out, the lineage null and BH.

The fixed chain runs end to end on scale200 from a scratch runroot, and every
change against production is attributed to a named fix. The DR positive
controls are intact. Three things need attention before the rerun:

- **R2-TREES-1 (HIGH, process):** the chain still builds a cohort-only tree,
  and it refuses the CX333-combined tree the user requires. Audit TP-5 was not
  addressed.
- **R2-TREES-2 (MEDIUM):** the new provenance guards do not key on code, the
  tree or the outgroup alignment.
- **R2-TREES-3 (MEDIUM):** at P = 20,000, the permutation floor cannot reach
  q < 0.05 in the large burden families (ASSOC-11 not addressed).

---

## What was run (all reproducible from `scripts/` and `runs/`)

| run | code | inputs | purpose |
|---|---|---|---|
| `old_on_prod` | main (= production) | production `events/` | **reproduces production `scan.tsv` byte for byte**, and `sv_gene.tsv` and `is6110_gene.tsv`, so production is deterministic and differences are attributable |
| `noise` | main, `--seed 777` | production `events/` | the Monte-Carlo noise floor: 0 rows with \|Δp\| > 4 SE; survivors 3 → 3; q_branch passes 137 → 139 |
| `fix_on_prod` | audit-fixes scan and burdens | production `events/` | isolates the scan and burden fixes |
| `rr` (fixed chain) | audit-fixes `cohort_assoc_tail.sh`, steps 1-6 | production merged VCF; production tree (with a `.prov` written by hand, see R2-TREES-1); scratch P0 assets | the fixed tail end to end |
| `rr2` | same | as `rr`, with the merged VCF's AA re-joined from the fixed `ancestral.tsv` (`merge_cohort_vcf.py`'s join rule) | adds the Fault B fix, which reaches the association only through a re-merge |

**Scratch P0 assets,** made with the branch's own code:

- `bin/vcf_collapse.sh` on CX333's decomposed VCF gives
  `build/7713a8d71d8e/assets/graph_collapsed.vcf.gz`: 93,214 keys, 0
  duplicates, in 5 minutes;
- `bin/panel_polarity.py` gives `panel_polarity.tsv`: 90,364 keys, 0 pooled;
- `bin/vcf_to_alignment.py` on the collapsed VCF gives the panel alignment;
- `bin/ancestral_alleles.py` gives `ancestral.tsv`.

The fixed chain's own audit (step 6) passed. Its one warning is the existing
small/masked warning. It reported 0 scan rows whose `gains` differ from
`n_gain`.

---

## Findings

### R2-TREES-1 (HIGH, process): the chain's tree step is not the CX333-combined tree the user requires, and refuses it; TP-5 not addressed

**Where:** `assoc/bin/cohort_assoc_tail.sh:149-163` (audit-fixes).
`analysis/combined_tree/{build_alignment.py,prune_for_cohort.py}` are
untracked and unchanged since 2026-10-01/02.

**What is wrong.** Step 3 still submits `bin/build_snp_tree.sh` on the cohort's
own `og.fasta`, which gives a cohort-only tree. Production's trees, and the
user's standing instruction, are the cohort plus the CX333 panel's SNPs, pruned
for the chain.

The fix also added `.prov` guards. A tree made by the combined route has no
`.prov`, so step 3 now stops with `FATAL: data/trees/<C>.rooted.nwk exists but
was not made from ...`. The only way to use it is to write the `.prov` by hand,
which I had to do in `rr/` and `rr2/`.

Two parts of audit TP-5 have no fix and no DECISION:

1. moving the combined route into the chain;
2. `build_alignment.py` still has the Fault A, missing-as-REF and
   upstream-record outgroup faults, and they apply to all 332 panel columns.
   Lines 110-131: `pan[(pos,ref,alt)] = calls`, so the last duplicate wins, and
   "no record" or a missing GT reads as REF.

**Effect.** At the rerun, following the chain gives scale200 and gwas1000 a
different kind of tree from today's. Every before/after comparison would then
move for a reason unrelated to the fixes. Following the hand route instead
uses an alignment with the unfixed outgroup faults for every panel genome.

**Suggested fix.**

- Make step 3 the combined route: `build_alignment.py` reads the build's
  collapsed VCF with `add_outgroup.outgroup_allele` applied per panel genome,
  then `build_snp_tree.sh`, then `prune_for_cohort.py`.
- Give it a `.prov` that records the panel VCF and code checksums.
- At minimum, document the hand route with a supported way to stamp the tree.

### R2-TREES-2 (MEDIUM): the chain's provenance guards ignore code, tree, outgroup alignment and node-locus inputs

**Where:** `cohort_assoc_tail.sh:109-128` (`_prov` records `build_id`,
`vcf_sha` and, for events only, `polarity_sha`), together with
`:131,140,157,174,190`.

**What is wrong.** Steps 1-4 skip when the product's `.prov` matches the VCF
and build. The record does not include:

- the code of `vcf_to_alignment.py`, `add_outgroup.py` or
  `write_event_matrix.py`;
- the tree checksum (events);
- `og.fasta` and `sites.tsv` (tree and events);
- `node_locus.tsv`, or the presence tables (events).

So these cases print "already built" and silently reuse the stale product:

- a later fix to any of those scripts;
- replacing `data/trees/<C>.rooted.nwk`, for example with the combined tree;
- a regenerated `og.fasta`.

TP-4 asked for code hashes. P0's steps do record `code_sha256`
(`p0_prepare.sh:589, 640`); the chain does not.

**Effect.** None for a first run into new folders. It is the same trap as
TP-4 for any re-run of the tail after a code fix under an unchanged VCF.

**Suggested fix.** Add the script sha256 to each step's `.prov`. Add the tree
checksum, `og.fasta`, `sites.tsv`, `node_locus.tsv` and a presence-table digest
to the events record.

### R2-TREES-3 (MEDIUM): the permutation floor cannot reach significance in large BH families (ASSOC-11 not addressed)

**Where:**

- `assoc/bin/assoc_scan.py:625, 630, 650` and `cond_null`;
- `assoc/bin/is6110_gene_burden.py:306, 311, 314`;
- all use `p = max(1, c) / P` with P = 20,000.

**What is wrong.** BH within a stratum of n units gives a single unit at the
floor q = n / 20,000. On the fixed scale200 run:

| family | units | best q for a lone hit at the floor |
|---|---:|---:|
| small burden `core:gene` | 3,548 | **0.177** |
| small burden `core:promoter` | 869 | 0.043 |
| scan `pe_ppe:genic` | 843 | 0.042 |
| scan `core:genic` | 486 | 0.024 |

gwas1000's families are larger.

- A gene in `core:gene` passes only when about 4 or more units tie at the
  floor; the DR genes do, so they pass together.
- Near-threshold units sit within Monte-Carlo error. katG has q_lineage 0.053
  in both production and fixed runs, from 3 exceedances of 20,000.
- `(c + 1) / (P + 1)` is still not used.

**Effect.** Genuine single hits in big strata are not reportable at P = 20,000.
It is conservative, but it can hide true hits; katG and up:embA are examples.

**Suggested fix.**

- Use `(c + 1) / (P + 1)`.
- Re-permute adaptively (10^5 to 10^6) any unit whose q could cross 0.05.
- Or report the attainable-q floor beside each family.

### R2-TREES-4 (LOW): `vcf_to_alignment.py --ref-sample` writes REF for the reference at node-frame columns

**Where:** `bin/vcf_to_alignment.py`, `cols[args.ref_sample].append(ref)`
after the site filters. It is unconditional.

**What is wrong.** This is the tree-alignment analogue of TP-3. The writer now
treats H37Rv as unknown at node-frame records, but the alignment still gives
H37Rv the node's REF base.

**Evidence:** the fixed scale200 alignment has 67 node-frame columns (`add_outgroup`
reports 67 `off_panel_contig`). In all 67, at least one isolate also carries
REF, so H37Rv is never the only REF.

**Effect:** small. It matters only for a tree built from this alignment, the
cohort-only tree.

**Fix:** write N for the ref-sample when the record's ID starts with `node:` or
its contig is not the reference.

### R2-TREES-5 (LOW): the sibling-allele rule in `vcf_to_alignment.py` assumes records at one locus are adjacent

**Where:** `with_siblings` uses `itertools.groupby(key=locus)`.

**What is wrong.** Node-frame records on one `node_<id>` contig all have
POS = 1 and are not sorted by offset.

**Evidence:** in scale200's merged VCF, 7 node loci appear in two or more
separated runs, for example `node:24243:34`. All 7 are indel-only groups, and
an indel's anchor never changes the base, so there is no effect today.

**Fix:** group with a dict per contig, not `groupby`.

### R2-TREES-6 (LOW, open question): panel AA and the cohort reconstruction disagree at the MTBC node for 80 variants

**Where:** design; not introduced by the fixes.

**What I measured.** I took the fixed run with the fixed AA (`rr2`). For each
AA-polarised record, I read the reconstructed state at the cohort tree's MTBC
node (`n00003`, 200 leaves; the root's children are ET1291 and
(`n00003`, `canettii`)):

- the node is **derived** at **80** AA-polarised variants, against 34,222 + 1,008
  that are consistent;
- 26 of the 80 have 2 or more losses and fewer than 2 gains. Their
  convergence sits in the loss column and is never tested. Examples:
  `1096470 G>C`, 19 losses; `2439519 G>A`, 14 losses.

With the production (Fault B) AA, the same count was 5,505.

**Question for the user.** Should the MTBC node be pinned to AA where AA
resolves? `--root ancestral` pins the true root (the canettii split), which is
not the MTBC node. Or should these variants be reported?

### R2-TREES-7 (LOW, cosmetic, pre-existing): `n_derived_leaves` and the scan's `carriers` count the outgroup, `canettii` and H37Rv tips

In the toy, record 200 reports 3 derived leaves: A, B and the outgroup. Read
the column as "derived tips", not "carrier isolates".

---

## Highest-priority checks

### (a) Fitch, root, outgroup polarity and gains/losses, re-derived

**Constructed tree** (`toy/`):

```
(GCF_035581225, ((A, B), (C, (D, H37Rv))))
```

Nine records, one per fixed behaviour:

| behaviour | audit ID |
|---|---|
| any duplicate with GT 1 → ALT | Fault A |
| `X,*` site whose outgroup is ALT through an upstream MNP | TP-2, TP-1 |
| outgroup inside a deletion → N | TP-1 |
| exact missing → N | TP-1 |
| indel polarised by the table | TP-8 |
| `AA_INVERTED` | |
| a sibling allele 700 T>C / T>G (third base → N; carrier of the other allele → N) | TP-7 |
| node-frame record with the H37Rv tip unknown | TP-3 |

The chain ran on it: `vcf_to_alignment.py`, `add_outgroup.py`,
`panel_polarity.py`, `write_event_matrix.py`. I worked the up and down passes,
root state, events, `n_gain`, `n_loss` and `n_undet` by hand for all nine.
**All nine match.**

For example, 100 C>T: the table says ALT is ancestral, so derived = REF.

- Up pass: n5 = {A,D}, n4 = {D}, n2 = {A,D}, root = {A}.
- Down pass: n4 = D (a gain); D-leaf = A (a loss).
- So 1 gain, 1 loss, root ancestral, as written.

For the node-frame record, the root is derived with 0 gains. The old code gave
a gain from the H37Rv tip it fabricated.

**Real data, full tree.** `scripts/indep_fitch.py` is an independent set-based
Fitch. It has:

- its own newick parser;
- its own polarity from INFO, the panel table and min-af;
- its own outgroup and reference-tip rules.

I ran it over the fixed scale200 event matrix (`rr/`, 403 nodes). **All
85,575 comparable records match the writer exactly:** every node state, every
branch event, `n_gain`, `n_loss`, `n_undet` and the polarity. Not compared:

- 308 level-2 conditional records;
- 1,741 long insertions that are ladder candidates.

**A real subtree:** 4247730 G>A (embB G406). The 6 carriers are tips. Each
tip's parent (n00091, n00094, n00116, n00136, n00172, n00190) is ancestral, so
there are 6 tip gains, 0 losses and obs = 1.0. That matches the scan row
(gains 6, obs 1.0000).

**Ancestral alleles on the CX333 panel tree:** the fixed
`ancestral_alleles.py` was run on production `data/trees/cx333.{rooted.nwk,
snps.fasta,sites.tsv}`.

- **It matches check 2's column M at 72,986 of 72,986 sites**, including all
  15 ties. Its ingroup is the MRCA of 331 leaves; its outgroups, nearest first,
  are GCF_000253375, then GCF_035581225.
- On the branch's P0 route (an alignment made with fixed `vcf_to_alignment.py`
  from the regenerated collapsed VCF):
  - 72,938 of the 72,940 shared sites equal M;
  - the 2 differences become TIED, both outside lineages 1-4, from
    sibling-allele N (TP-7);
  - 46 sites drop out because sibling N now counts toward `--max-missing`.
    That is D12, and the 37 + 9 split matches its "keep 37 panel sites";
  - 316 sites are new.

### (b) The fixed tail end to end on scale200, against production

**Event matrix.** 87,624 rows on both sides; 21,644 change. Every change
carries at least one fix flag; **0 are unexplained**:

| rows | fix |
|---:|---|
| 11,234 (+59 with the next) | the node-frame reference tip (TP-3) |
| 7,470 (+306 with polarity) | the `X,*` site now in the alignment, so the outgroup is known (TP-2) |
| 1,394 (+1) | an outgroup allele change (Fault A / TP-1 / D35) |
| 1,176 | a polarity-table change (TP-8, D13, D34) |
| 63 | the regenerated `node_locus.tsv` (no hand-made inputs) |

On the outgroup cells shared with production (46,379):

- REF → ALT 525;
- REF → N 724;
- N → ALT 144;
- N → REF 2.

These are in line with TP-1's count (386 Fault A + 157 upstream MNP).

**Scan, production vs `rr`:**

- 2,602 → 2,622 rows;
- 14 rows only in production, 34 only in the fixed run, 22 with changed gains.

Every one of those 70 is explained by its event flags:

| rows | cause |
|---:|---|
| 47 | `X,*` outgroup |
| 41 | outgroup allele |
| 27 | polarity |
| 2 | node-frame tip |

All other differences are explained by the scan fixes. `fix_on_prod` isolates
them:

- 467 svi rows move from `other` to `svi:<tier>` (ASSOC-4);
- 3 rows switch genic and intergenic (ASSOC-7);
- 4 `#2` rows now use their own branches (ASSOC-3);
- `null_pool` = pool − own gains in all 2,465 fine-level rows (ASSOC-2);
- the remaining p-value differences are within Monte-Carlo noise once the RNG
  stream shifts. The `noise` run shows that scale.

**Scan, `rr` vs `rr2` (Fault B through the AA re-join):**

- 5,914 records change AA, and exactly those 5,914 change events;
- the swaps are exact gain ↔ loss relabels. For example, 2266624 G>T goes from
  3 gains / 29 losses to 29 gains / 3 losses;
- 60 rows enter the scan and 8 leave.

**Survivors.** 3 in production, `rr` and `rr2`: rpoB 761155, embB 4247429 A>G
and embB 4247730 G>A.

**Burdens:**

| burden | production → `fix_on_prod` → `rr` | causes | all three nulls |
|---|---|---|---|
| small | 5,189 → 5,146 → 5,172 | D17, D40, the gene-coordinate fix, then event changes | the same 6 genes |
| sv | 179 → 263 → 263 | D16, as DECISIONS predicts | 0 |
| is6110 | 77 → 73 → 73 | D17 | 0 |

The 6 small-burden genes are rpoB, rpoC, embB, gyrA, pncA and ethA.

**DR positive controls** (scan, by gene ±200 bp, production → fixed):

| gene | q_branch / q_region / q_lineage, production → fixed |
|---|---|
| rpoB | 0.0018/0.0048/0.0060 → 0.0022/0.0081/0.0061 (survives) |
| embB | survives (2 rows) |
| rpoC | 0.0022 / 0.0146 / 0.109 (lineage fails, as before) |
| katG S315T | q_branch 0.0127 → 0.0049; q_region 0.54 (fails, as before) |
| gyrA | q_lineage 0.073 |
| fabG1 promoter | q_region 0.18 |
| rpsL and rrs | branch only |
| pncA, ethA, gid, eis | no row reaches q_branch < 0.05 (as before) |

In the small burden:

- rpoB, rpoC, embB, gyrA, pncA and ethA pass all three nulls;
- katG has q_lineage 0.053 and up:embA 0.058, unchanged (see R2-TREES-3).

**No control is lost.**

gwas1000 was not run: UNVERIFIED.

### (c) Statistics

Re-derived independently with `scripts/stat_check.py` on the fixed scan
(`rr/`):

- **obs:** 0 mismatches over 2,622 rows.
- **Leave-one-out:** `null_pool` equals the stratum total minus the variant's
  own gains in all 2,465 fine-level rows. `leave_one_out` is a true multiset
  subtraction, and the floor is now checked after it.
- **p-values:** for 9 rows (rpoB, embB, rpsL, the fabG1 promoter, gyrA A>G and
  A>C, katG, and two masked rows), independent simulations agree within
  Monte-Carlo error. The simulations used P = 100,000 for the branch and
  region nulls and 20,000 for the lineage null, with a different RNG. The
  largest difference is katG p_branch, 3.6e-4 against 1.5e-4, which is 1.6 SE.
  - The branch null weights branches by edge length; the root has weight 0.
  - The region null draws from the leave-one-out pool.
  - The lineage null holds each lineage's carrier count, with H37Rv and the
    outgroup in an `unassigned` stratum with 0 carriers, and `canettii` in
    `unknown`.
- **BH:** recomputed within region for each null; 0 mismatches. The family is
  the fine region key, even when the pool came from a coarser rung (documented
  in the code).
- **Survivor rule D15:** implemented as decided:
  - q_branch is required for every variant;
  - a NaN q_lineage fails;
  - a conditional variant needs q_cond, q_branch and q_lineage.
- **Lineage table guard:** it needs at least 95% of tips assigned; scale200 has
  200 of 202 (99%).

**Still not changed** (audit LOW, no DECISION):

- `max(1, c) / P` (R2-TREES-3);
- the branch and region nulls sample with replacement while the observed
  branches are distinct;
- the lineage null permutes tips within top-level lineage only (ASSOC-12).

---

## Checked and found sound

- **`add_outgroup.py`:**
  - the allele rule in the toy and on scale200;
  - every record is filed under every site its REF span covers;
  - multiallelic and third-base cases give N;
  - node-frame sites give N;
  - `--panel-vcf` and `--outgroup` are required, and the chain passes the
    build's.
- **`panel_polarity.py`:**
  - pools by (chrom, pos, ref, alt), one row per key;
  - carriers over genomes called;
  - the writer refuses a table with repeated keys or without a `chrom`
    column;
  - the PanSN contig is stripped, so a node-frame POS = 1 cannot match.
- **`write_event_matrix.py`:**
  - `primary_alt` for the outgroup and the table;
  - `ref_tip_gt` (is6110, sv and accessory_presence stay REF, as D14 says);
  - `og_gt[keep]` after dedupe (TP-6);
  - a missing table is fatal;
  - `--node-locus` is required with presence tables;
  - its `summary.txt` records the table.
- **`ancestral_alleles.py`:**
  - the ingroup is the MRCA of non-outgroups;
  - it stops if the MRCA contains an outgroup;
  - nearest-first tie-breaking gives the same result as Fitch's down pass
    through (CIPT, ET1291);
  - P0 step `ancestral` records the tree, alignment, sites and code checksums,
    and `set -e` stops a failed run before the marker is written.
- **`build_snp_tree.sh`:** the outgroup comes from `MTB_OUTGROUP`, and the
  chain passes the build's.
- **`cohort_assoc_tail.sh`:**
  - the runroot-only check;
  - the build stamp check;
  - outgroup from `build_info` with an `MTB_OUTGROUP` fallback;
  - step 4a `node_locus`;
  - step 5 retier, which tiers all 467 svi IDs;
  - `LINTAB` is fatal when missing;
  - the audit reads the burdens and `gains` vs `n_gain`.
- **`retier_intervals.py`:**
  - the BED-against-1-based overlap with a running maximum of ends;
  - missing BED files are fatal.
- **`sv_scatter.tree_frame`:** refuses non-unique labels.
- **The scan's and burden's guards:**
  - `event_key` keeps the `#n` suffix;
  - a row missing from the event matrix is fatal;
  - a missing `--genes` file is fatal;
  - genes are read as 0-based start + 1.
- **The burden's crediting:**
  - D16 (any overlap) and D40 (small-deletion spans);
  - the D17 floor (`n_undet` over n−1 branches, as in the scan).
- **Tests:** `tests/test_audit_polarity.py` and `test_audit_assoc.py` pass, 43
  tests, run from the export.
  - They are behavioural: exact AA, root, `n_gain` and outgroup-alignment
    assertions; a two-outgroup Fault B tree; true multiset leave-one-out.
  - The `EventWriter` and `PanelPolarityTable` tests `skipTest` with a message
    when `MTB_PY_VT` or bcftools is missing. They ran here.
  - No test asserts a whole per-branch event matrix (TP-11's suggestion).
    The toy above could be one.

## Audit findings not addressed

- **TP-5** (MEDIUM): the combined-tree route is not in the chain, and
  `build_alignment.py`'s outgroup faults are unfixed. See R2-TREES-1.
- **ASSOC-11** (LOW): the p-value form and resolution. See R2-TREES-3.
- **ASSOC-12** (LOW): the lineage null remains a within-top-lineage tip
  permutation. The help text describes it, but nothing changed.
- **ASSOC-1:** deferred by the user. `canettii` is still a tip and is sister
  to the MTBC under the root.
- **TP-4** is only partly fixed. See R2-TREES-2.

## Open questions

1. R2-TREES-6: pin or report the 80 variants where AA and the cohort MTBC node
   disagree?
2. Rerun prerequisite, confirmed: the Fault B fix reaches association only
   through a re-merge.
   - The merged VCF's AA must come from the regenerated `ancestral.tsv`.
   - On scale200 that changes 5,914 records, with +60 and −8 scan rows.
   - P0 step `ancestral` must be run with explicit `ANC_TREE`, `ANC_ALN` and
     `ANC_SITES`. Per `docs/PANEL_TREE.md`, a new panel tree is needed for a
     new graph.
3. Losses at known resistance sites (an assoc audit open question): katG S315T
   still has 15 gains and 12 losses. That is parsimony with equal costs
   preferring a gain plus a reversal. If the user wants reversals of DR
   mutations discouraged, that needs a weighted or Dollo-type
   reconstruction, which is a method change and not a fix.
4. gwas1000 was not re-run: UNVERIFIED.
