# Audit: association statistics and their chain (`assoc_tests`)

Date: 2026-10-05. Read-only audit. Scratch scripts are in
`/tmp/claude-12043/audit_assoc_tests/`. They were run against the production
outputs in the working tree (`W` = `/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome`).

**Which copy runs.** The 2026-10-02 production chain was the **repo** copy,
`assoc/bin/cohort_assoc_tail.sh`. It was launched by
`analysis/combined_tree/assoc_tail.sbatch` (absolute path) from the working-tree
root, so every `assoc/bin/*.py` it calls is the **working-tree** copy, and the
audit is the repo's `bin/audit_chain.py`. The burdens in production
(`small_gene.tsv`, `sv_gene.tsv`, `is6110_gene.tsv`, 18:03 to 18:33) came from
`burden_only.sbatch` (slurm jobs 50014892 to 50014895). `scan.tsv` came from
job 49992590 (gwas1000) and 50003176 (scale200). The working tree also holds
older copies of `cohort_assoc_tail.sh` and `bin/audit_chain.py`. They are
stale (see ASSOC-6).

Severity counts: **HIGH 1, MEDIUM 5, LOW 8.**

---

## ASSOC-1 (HIGH): a read-based *M. canettii* isolate is a tip of both association trees and dominates the branch-length null

- **Where:** `W/refbias/cohort.gwas1000.tsv:2` and `W/refbias/cohort.scale200.tsv:23` (sample `canettii`, deepest `M.canetti`).
  - It passes through `analysis/combined_tree/prune_for_cohort.py:45` (repo), which keeps every VCF sample.
  - It then reaches the branch null in `W/assoc/bin/assoc_scan.py:288-292,519-520` and `W/assoc/bin/is6110_gene_burden.py:104-106,207-208`.
  - Production runs these copies.
- **What is wrong.** The cohort tables include the pilot's *M. canettii* read sample, so it is a leaf of `assoc/<C>/events/labelled.nwk`.
  - It hangs off `n00002`, the root's non-outgroup child, with the longest branch in the tree:
    - gwas1000: 0.0314, which is **11.4%** of total tree length;
    - scale200: 0.0528, which is **25.8%**.
  - Its phenotype is 0, so its per-branch statistic is 0.
  - The branch null draws branches in proportion to length, so 11% and 26% of all null draws are a guaranteed 0. That lowers the null mean:
    - gwas1000: E[ov] per draw is 0.379 with the tip, 0.427 without it;
    - scale200: 0.164 with it, 0.221 without it.
  - Gains on its branch also count as independent origins. Every allele shared between *canettii* and one MTBC clade becomes "convergent".
  - The outgroup `GCF_035581225` has edge length 0 and does no harm. The second *canettii* is the problem.
- **Evidence** (`canettii.py`, `canettii2.py`):

  | | gwas1000 | scale200 |
  |---|---|---|
  | scanned rows with a gain on the `canettii` branch | 227 of 7,914 | 120 of 2,602 |
  | of those, fewer than 2 gains without it (testable only because of it) | 174 | 83 |
  | q_branch < 0.05, as run (my rerun, P=5000) | 68 (file: 62) | 138 (file: 137) |
  | q_branch < 0.05, `canettii` removed from null weights and event sets | **38** | **31** |

  The rows that lose branch significance include gyrA 7572, embB 4247730 G>A, 4249583 and 4247574, the eis promoter 2715344 and 2715369, and several svi deletions.
- **Effect on current outputs.**
  - Every `p_branch`/`q_branch` in `scan.tsv` and the three burden files is anti-conservative.
  - The headline counts in HANDOFF 0b are inflated: "11,322 / 3,607 variants with 2 or more independent gains" and "callable 7,914 / 2,602".
  - Survivors are decided mainly on q_region and q_lineage (ASSOC-9). The region pool holds few `canettii`-branch entries (about 1.4% of core:genic). The 28 and 3 survivors all sit at the permutation floor, so they very likely stand. This is UNVERIFIED for near-threshold rows.
  - HANDOFF's "katG and the embA promoter ... pass the branch null" (scale200) rests on this inflated null.
- **Suggested fix.** Drop non-MTBC isolates from the association tree. Either pass a drop list to `prune_for_cohort.py`, or prune via `write_event_matrix.py`'s existing prune path (line 427). Alternatively, take `canettii` out of `cohort.gwas1000.tsv` and `cohort.scale200.tsv` for association use. At minimum, give the branch-null weight of any tip outside the phenotype universe a value of zero, and exclude gains on that branch from `n_gain`.

## ASSOC-2 (MEDIUM): the region null's leave-one-out never removes anything

- **Where:** `W/assoc/bin/assoc_scan.py:527-529` and `W/assoc/bin/is6110_gene_burden.py:209-211`. Production runs these copies.
- **What is wrong.**

  ```python
  own = collections.Counter(br.tolist())
  cand = [x for x in pool[rung] if not (own[x] and own.__setitem__(x, own[x] - 1))]
  ```

  - When `own[x] > 0`, `own[x] and own.__setitem__(...)` evaluates to `None`, so `not None` is True and the element is **kept**.
  - When `own[x] == 0`, it evaluates to `0`, so `not 0` is True and the element is kept.
  - Nothing is ever excluded. With `own = Counter([1])`, the expression on `[1, 1, 2]` returns `[1, 1, 2]`; the intended result is `[1, 2]`.
- **Evidence.**
  - `null_pool` equals the full stratum pool in every row. For example, embB 4247431 G>C (11 gains) and G>T (4 gains) both report 10,699, which is the core:genic pool. My independent re-derivation with true leave-one-out gives 10,688 and 10,695.
  - Rows whose own branches are more than 5% of their pool:
    - scan: 6 in gwas1000 (accessory_presence, up to 11%);
    - burdens: 11 small, 10 sv and 6 is6110 in gwas1000; 8 small and 7 sv in scale200.
  - The comment at line 521 says the floor is checked "AFTER leave-one-out". It is not: **15 scale200 sv-burden units** have a region null that true leave-one-out would leave below the 200 floor.
- **Effect.** The bias is conservative: a variant's own high-ov branches sit in its null. No current survivor is lost on it alone. The largest own-share rows are not near threshold, with q_region of 0.43 or more.
- **Suggested fix.** Use explicit multiset subtraction:

  ```python
  own = collections.Counter(br.tolist()); cand = []
  for x in pool[rung]:
      if own[x]: own[x] -= 1
      else: cand.append(x)
  ```

  Apply the same change in the burden.

## ASSOC-3 (MEDIUM): `--dedupe suffix` records are tested on the other record's branches

- **Where:** `W/assoc/bin/assoc_scan.py:364-367` and `W/assoc/bin/is6110_gene_burden.py:163-166`. Production runs these copies.
- **What is wrong.** The bytestream key is rebuilt from the first three `|` fields only:
  `key = (int(k[0]) ..., k[1], k[2])`.
  - `write_event_matrix.py --dedupe suffix` writes repeated keys as `(pos, ref, alt, '#2')`.
  - So `ev.calls.loc[(pos, ref, alt)]` returns the **first** record's event row for the `#2` record.
  - This is pattern 1 again: same (pos, ref, alt), wrong record.
- **Evidence.** Scan rows whose `gains` differ from `variants.tsv` `n_gain`:
  - gwas1000: **7 rows**, all `svi:...#2`. Example: `836370|G|<DEL>,*|#2` (svi:DEL:836371:1010) has n_gain 32, but the scan used 80 gains and obs 0.3610, identical to svi:DEL:836371:135.
  - scale200: **4 rows**.
  - In the sv burden, the `#2` records' own branches are replaced by the first record's: 10 gwas1000 and 6 scale200 records with at least one gain.
- **Effect.** 11 scan rows and about 16 burden records. None is a survivor. The closest is gwas1000 svi:DEL:3942723:231 at p_branch 1.5e-4 on the wrong branch set.
- **Suggested fix.** Use the full key: `key = (int(k[0]) if k[0].isdigit() else k[0],) + tuple(k[1:])`. Better, keep a row-index map from `variants.tsv` order. Also make `audit_chain.py` compare per-row `gains` with `n_gain`.

## ASSOC-4 (MEDIUM): the chain dropped the SV evidence-tier stratification, and its input is stale

- **Where:** repo `assoc/bin/cohort_assoc_tail.sh:103-106` does not pass `--sv-intervals`. The option is handled at `W/assoc/bin/assoc_scan.py:80-88,396-425`. Production runs the repo chain.
- **What is wrong.** The option's own help says it should be passed. Without it, all catalogued deletions share one `other` pool and one BH family, "of which 77% are scattered repeat-context intervals". This is pattern 3 and pattern 4 together.
  - The pre-rerun scan (`assoc/archive/stale_pre_rerun_20261002/gwas1000/assoc/scan.tsv`) used it, with regions `svi:E1_coherent` (119), `svi:E3_scattered_in_repeat` (621) and so on. `GWAS1000_RESULTS.md:415-424` describes that result.
  - The current scans have no `svi:` region: all 1,079 (gwas1000) and 467 (scale200) deletions are in `other`.
  - The table that carries `evidence_tier`, `refbias/assets/sv_intervals.retiered.tsv` (2026-09-28), is from the **pre-rerun** catalogue. Only 398 of 1,079 current gwas1000 svi IDs and 292 of 467 scale200 IDs appear in it. The current `refbias/<C>/p5/sv_intervals.tsv` has `support_tier` but no `evidence_tier`.
- **Effect.**
  - The region null and BH for every svi deletion now run against a mixed background.
  - Branch-null passes in `other`: 10 (gwas1000) and 100 (scale200). Region-null survivors: 0 either way.
  - The documented method is not what ran.
- **Suggested fix.**
  - Run `bin/retier_intervals.py` on each cohort's current catalogue.
  - Pass `--sv-intervals <that file>` in the chain.
  - Make the scan fail, or at least warn, when svi records exist but no tier table was given.

## ASSOC-5 (MEDIUM): the SV burden assigns a deletion to the unit at its anchor base only

- **Where:** `W/assoc/bin/is6110_gene_burden.py:158-161`, which calls `unit_at(int(r["pos"]))` (lines 119-136). Production runs this copy with `--cls sv`.
- **What is wrong.** A deletion is placed by its VCF POS, the anchor base before the deleted span, and nothing else.
  - A deletion that removes several genes is credited to one gene, or to a promoter or intergenic unit.
  - The unit can even be the gene ending at the anchor base. Example: svi:DEL:3348474:59 is credited to Rv2991 but deletes none of it.
- **Evidence** (`svspan.py`, svi:DEL records with at least one gain):

  | | gwas1000 | scale200 |
  |---|---|---|
  | records | 1,784 | 854 |
  | delete a gene other than the unit credited | **380 (21%)** | 166 (19%) |
  | span 2 or more genes | 328 | 134 |
  | credited to promoter/intergenic but delete genic sequence | 92 | 52 |

- **Effect.** The `sv_gene.tsv` unit composition and origins in both cohorts. There are no survivors either way, but the gene list is not "genes deleted".
- **Suggested fix.** For sv records, parse the span (ID `svi:DEL:<start>:<len>`, or the record's END/SVLEN) and add the record's branches to **every** gene the span overlaps. Fall back to `unit_at` only when the span overlaps no gene.

## ASSOC-6 (MEDIUM): stale copies of the chain and the audit in the working tree, and the repo script's usage line runs the stale one

- **Where.**
  - Repo `assoc/bin/cohort_assoc_tail.sh:4` documents `bash assoc/bin/cohort_assoc_tail.sh gwas1000 [phenotype.txt]`.
  - The chain must be run from the working-tree root (HANDOFF 0b), where that relative path resolves to `W/assoc/bin/cohort_assoc_tail.sh` (2026-09-28).
  - That copy also calls `W/bin/audit_chain.py` (2026-09-28).
- **What is wrong.** Following the documented usage runs the stale script. Compared with the repo copy (`diff` above), the stale script:
  - passes no `--accessory-presence` to the scan;
  - passes no `--lineages` to the IS6110 burden;
  - has no small or SV burden;
  - runs the old audit, whose `2c - 1` applicable-branch estimate false-fails (`W/bin/audit_chain.py:240-256`).

  These are exactly the three omissions HANDOFF 0b records as fixed.
- **Effect.** No current output; the 2026-10-02 runs used the absolute repo path. It is a trap for the next rerun.
- **Suggested fix.**
  - Replace the working-tree copies with a stub that execs the repo copy, or delete them.
  - Change the usage line to the absolute repo path, as `assoc_tail.sbatch` does.

## ASSOC-7 (LOW): gene coordinates read 0-based as 1-based; overlapping genes credited to one gene only

- **Where:** `W/assoc/bin/is6110_gene_burden.py:33-42,112-117,126` and `W/assoc/bin/assoc_scan.py:381-393`. Production runs these copies.
- **What is wrong.** `data/annotation/H37Rv_snpeff_dump.txt` gives 0-based starts and 1-based ends. For example:
  - rpoB is `759806 763325`; the true span is 759807-763325;
  - the Chromosome line starts at 0;
  - Cds IDs read `CDS_Chromosome_759807_...`.

  The code tests `s <= p <= e` with 1-based VCF positions, so each gene gains one extra base on its left side. Promoter distances for plus-strand genes (`d = s - p`) are short by one.

  `gene_at` returns only the first matching gene. A variant in an overlap (17,874 bp of H37Rv) counts toward the earlier-starting gene only.

  I checked that the 3-gene look-back window misses no position.
- **Evidence.** 3,087 genome positions are genic only because of the off-by-one. Records with at least one gain at those positions: 134 in gwas1000 (8 with 2 or more gains, in the scan's genic/intergenic split) and 74 in scale200 (3).
- **Effect.** A handful of records move between gene and promoter units, or between core:genic and core:intergenic. No DR control is affected:
  - fabG1 c-15t gets d=14 but is still `up:fabG1`;
  - embA -12/-16 and eis c-14t are still promoters.
- **Suggested fix.**
  - Load genes as `(int(r[1]) + 1, int(r[2]), ...)`.
  - Credit a variant to every overlapping gene, or document the one-gene rule.

## ASSOC-8 (LOW): the burdens apply no callability floor, unlike the scan

- **Where:** `W/assoc/bin/is6110_gene_burden.py:149` admits every record with `n_gain >= 1` regardless of `n_undet`.
- **Evidence.** Share of burden gains that come from records below the scan's 80% floor:

  | class | gwas1000 | scale200 |
  |---|---|---|
  | small | 6% | 6% |
  | sv | 16% | 22% |
  | is6110 | 7% | 9% |

  Conditional (accessory) records do not reach the burdens; they have no H37Rv frame.
- **Effect.** Units, sv units most of all, are built partly from reconstructions the scan treats as untestable. The "testable units" counts in HANDOFF 0b include them.
- **Suggested fix.** Apply the same `--min-determinacy` per record, or report the floor-passing origins beside the total.

## ASSOC-9 (LOW): the survivor rule does not match "q < 0.05 under all three nulls"

- **Where:** `W/assoc/bin/assoc_scan.py:642-649`.
- **What is wrong.**
  - For a non-conditional variant, `survives()` requires q_region and q_lineage only, **not q_branch**.
  - A NaN q_lineage (no lineage table) counts as a pass.
  - The printed "all three" column (line 625) does require q_branch.
  - HANDOFF 0b defines "survive" as all three.
- **Effect.** None today: survivors equal the "all three" counts (28 = 24 + 4; 3). They would diverge once ASSOC-1 is fixed, and whenever `--lineages` is missing.
- **Suggested fix.** Require q_branch for every variant. Treat a missing lineage null as "not tested", never as a pass.

## ASSOC-10 (LOW): remaining silent defaults (pattern 3)

| where | behaviour |
|---|---|
| repo `cohort_assoc_tail.sh:95-96` | `LINTAB` can resolve to an empty or nonexistent path with no error. `burden_only.sbatch` has a FATAL here; the chain does not. |
| `assoc_scan.py:173`, burden `:86` | A nonexistent `--lineages` path makes `p_lineage` NaN, which ASSOC-9 counts as a pass. |
| `assoc_scan.py:175-180`, burden `:88-93` | There is no check that the lineage table's samples or columns match the tree. A mismatch puts every tip in `unassigned`, which silently turns the lineage null into a plain permutation. |
| `assoc_scan.py:155-159`, burden `:79-83` | Phenotype names that are not tips are dropped silently. Checked: 459/459 and 57/57 match today. |
| `assoc_scan.py:365-369`, burden `:165-168` | `except Exception: continue` on the event lookup. The audit catches this for the scan but never reads the burdens. |
| `assoc_scan.py:381` | A missing `--genes` file silently removes the genic/intergenic split. |
| `assoc_scan.py:664` | `--phyoverlap2` is never passed, so the self-check against phyoverlap2 never runs. `phyoverlap2.py` source is not on disk, only `__pycache__`. |
| `bin/audit_chain.py` | It never reads `small_gene.tsv`, `sv_gene.tsv` or `is6110_gene.tsv`. |

**Fix:** fail on each of these, or print the count and fail when it falls below a threshold.

## ASSOC-11 (LOW): permutation resolution and p-value form

- **Where:** `assoc_scan.py:285,520,537,551` and burden `:208,215,218`.
- **What is wrong.**
  - `p = max(1, c)/P` is used in place of `(c + 1)/(P + 1)`, which is slightly anti-conservative when c ≥ 1.
  - With P = 20,000 and BH over 3,562 to 6,457 units per stratum, near-threshold calls rest on 1 to 3 exceedances. Example: scale200 katG has p_lineage 0.00015 (3 of 20,000) and q_lineage 0.053; up:embA has q_lineage 0.057. HANDOFF reports these as "narrowly miss", but that difference is within Monte Carlo error. A Poisson 95% interval on 3 exceedances is about 0.6 to 8.8.
  - All three nulls share one RNG stream. Adding `--lineages` changed `ig:Rv2813-Rv2814c`'s q_branch from 0.064 (HANDOFF 0b) to 0.076 (current `is6110_gene.tsv`).
- **Suggested fix.**
  - Use `(c + 1)/(P + 1)`.
  - Rerun units with q in [0.01, 0.2] with 10^6 permutations.
  - Give each null its own seeded generator.

## ASSOC-12 (LOW): the lineage null is a tip permutation, not a phylogenetic null

- **Where:** `assoc_scan.py:172-189` and burden `:86-102`.
- **What is wrong.** The phenotype is shuffled among tips within the top-level lineage only (`lineage1` ... `lineage9`). This destroys within-lineage clustering of the phenotype. Under the permutation, internal branches score near the lineage prevalence, so events on internal branches inside resistant clades beat it easily.
- **Evidence.** In gwas1000 the lineage null passes rows that the branch null does not:
  - pe_ppe:genic: 13 lineage passes, 0 branch;
  - masked: 12 lineage passes, 0 branch;
  - other: 2 lineage passes, 0 region.

  In gwas1000, carriers are 96/192 in each of lineages 1 to 4 by design, so there the null is close to an unstratified permutation.
- **Effect.** Interpretation only. Survivors also need the region null.
- **Suggested fix.** Document it as a between-lineage-prevalence control only. Alternatively, stratify by `deepest` sublineage, or permute within clades cut at a depth.

## ASSOC-13 (LOW): chain products made outside the chain, and existence-only guards (patterns 4 and 5)

- **Where:** repo `cohort_assoc_tail.sh:44,52,60,74`.
- **Products in HANDOFF 0b that `cohort_assoc_tail.sh` does not produce:**
  - `data/trees/<C>_cx333.*`, from `build_alignment_cohort.sbatch` and `build_snp_tree.sh`;
  - the pruned chain tree `data/trees/<C>.rooted.nwk`, from `prune_for_cohort.py` in `assoc_tail.sbatch`;
  - the phenotype `assoc/<C>/rrdr_carriers.txt`, made by hand on 2026-09-25 and not regenerated;
  - the iTOL files.
- **Effect of the step-3 fallback.** If `<C>.rooted.nwk` is absent, step 3 submits a **cohort-only** tree. That contradicts the CX333 practice.
- **Effect of the guards.** Steps 1 to 4 skip on existence, not on input freshness. A rerun on a new VCF or tree silently reuses old `snps.fasta`, `og.fasta` and `events/` unless they are archived first.
- **Checked:** today's inputs are all newer than the 2026-10-01 VCFs.
- **Checked:** `rrdr_carriers.txt` matches the current calls:
  - scale200: 57/57 against RRDR ALT calls in the current VCF;
  - gwas1000: equals `phenotype.tsv` `rrdr == 1`, which is the collection label; 6 more isolates carry an RRDR ALT call, consistent with `pipeline_checks`.
- **Suggested fix.** Fold the CX333 alignment, tree and prune, plus a phenotype builder, into the chain. Guard each step on its inputs' mtime or checksum (or write a manifest), not on its output's existence.

## ASSOC-14 (LOW): `tree_frame` does not check that node labels are unique

- **Where:** `W/assoc/bin/sv_scatter.py:35`. `order = {n.label: i ...}` collapses duplicate or `None` labels silently, which corrupts `children` and `parent`.
- **Effect.** None in the chain: it reads only `labelled.nwk`, whose uniqueness `write_event_matrix.py:446` asserts. It is a hazard for anyone who passes a raw IQ-TREE tree.
- **Fix:** `assert len(order) == len(nodes)`.
- **Status of the two diagnostic scripts.** `sv_scatter.py` and `sv_scatter_d.py` are not in the chain; only `tree_frame` and `events` are imported. In `sv_scatter.py:142`, `r.get("reference") or list(r.values())[4]` is a positional fallback. Diagnostic only.

---

## Re-derived numbers

I re-derived these from `events/` with independent code (`rederive.py`, `rederive2.py`, `burden_rederive.py`). That code builds the tree from `nodes.tsv`, maps events by **label** rather than by column index, and reads gene spans in true 1-based coordinates.

| row | gains (mine / file) | obs (mine / file) | other |
|---|---|---|---|
| gwas1000 rpoB 761155 C>T | 105 / 105 | 0.9199 / 0.9199 | p_branch 5e-5 / 5e-5; 262 derived leaves |
| gwas1000 katG 2155168 C>G | 97 / 97 | 0.7645 / 0.7645 | p_branch at floor in both |
| gwas1000 embB 4247429 A>G | 47 / 47 | 0.9416 / 0.9416 | p_branch at floor in both |
| gwas1000 fabG1 promoter 1673425 C>T | 59 / 59 | 0.7642 / 0.7642 | p_branch at floor in both |
| gwas1000 embB 4247431 G>C | 11 / 11 | 0.7424 / 0.7424 | p_branch 0.0026 / 0.0029; p_region 0.0326 / 0.033; p_lineage 0.0269 / 0.0277 |
| gwas1000 embB 4247431 G>T | 4 / 4 | 1.0 / 1.0 | p_branch 0.0052 / 0.0054; p_region 0.0243 / 0.027; p_lineage 0.0604 / 0.0612 |
| burden gwas1000 rpoB, katG, embB, pncA, Rv2000 (random) | origins 352, 195, 288, 208, 48 (equal) | equal | carriers_with_phenotype equal |
| burden scale200 rpoB, katG, embB | origins 54, 28, 66 (equal) | equal | carriers_with_phenotype equal |

- Event bytestream column order equals the tree's preorder labels in both cohorts (1,997 and 403 columns), so the scan's index-based mapping is valid.
- BH q-values recomputed within region from the file's p-values match.

## DR positive controls: where each appears

- **gwas1000, `scan.tsv`: 28 survivors.** All have p at the floor in all three nulls, except 1673432, whose p_lineage is 3e-4.
  - rpoB: 761095, 761109, 761110, 761139 (C>A, G, T), 761155 (S450L), 761161;
  - rpoC: 764817 (T>C, T>G) and 764840;
  - katG: 2155168 (S315T);
  - inhA promoter, `up:fabG1`: 1673425 (c-15t) and 1673432 (t-8c);
  - inhA: 1674481;
  - embB: 4247429 (A>G, A>C; M306), 4247431 G>A, 4247730 (G406) and 4248003 (Q497);
  - **embA promoter: 4243217 and 4243221.** These are present but missing from HANDOFF 0b's list;
  - gyrA: 7570 (A90V) and 7582 (D94G);
  - rpsL: 781687 (K43R) and 781822 (K88R);
  - rrs: 1472362 and 1473246 (a1401g).
- **gwas1000, `small_gene.tsv`: 18 units under all three nulls,** as HANDOFF states.
  - Genes: rpoB, rpoC, rpoA, katG, inhA, embB, gyrA, rpsL, rrs, gid, pncA, ethA, thyA.
  - Promoters: up:fabG1, up:ahpC, up:embA, up:pncA, up:eis.
- **scale200, `scan.tsv`: 3 survivors:** 761155, 4247429 A>G and 4247730 G>A.
  - katG S315T fails the region null (q 0.55; 15 gains, **12 losses**, 59 carriers).
  - gyrA D94G and the fabG1 promoter pass the branch null only. That is expected at 57 phenotype carriers.
- **scale200, `small_gene.tsv`:** rpoB, rpoC, embB, gyrA, pncA and ethA. katG and up:embA are at q_lineage 0.053 and 0.057 (see ASSOC-11).
- **IS6110 and SV burdens:** no DR unit, as expected. 0 survivors in both cohorts.

## Checked and found sound

- **Sample joins.**
  - All tree leaves are in the lineage tables except H37Rv and the outgroup, which go to the `unassigned` stratum with 0 carriers.
  - There are no duplicate carriers.
  - Lineage strata in the logs match the tables: gwas1000 has 96/192 in each of lineages 1 to 4.
- **Event lookup.** The event bytestream column order equals the tree preorder, and gains equal `n_gain` for every non-suffixed row.
- **Descendant matrix and per-branch statistic.** `descendant_matrix`, `ov`, and the observed statistic agree with an independent implementation to 4 decimals.
- **Callability floor.**
  - `n_undet` excludes the root and `nbranch = n - 1`.
  - The conditional correction subtracts exactly the inapplicable branches. Fitch's downpass keeps all-unknown subtrees UNK (`write_event_matrix.py:315-316`), so every inapplicable branch is undetermined.
  - `applicable_branches` excludes the root.
- **Level-2 null.** It permutes within carriers, holds kA fixed, and skips the test when kA is 0 or equal to nA. q_cond is computed within region.
- **Lineage null.** The construction holds each stratum's carrier count, and `OVP` indexing is correct.
- **BH.** The step-up implementation `bh()` is correct. It is applied within region per null, and NaNs are excluded.
- **Region ladder.** Conditional rows never pool with unconditional ones. `region_level` is recorded.
- **Gene lookup and promoter strand logic.** The 3-gene look-back misses no genic position (brute force over 4.41 Mb). Minus-strand promoters sit at higher coordinates, which is correct (up:pncA, up:eis).
- **`audit_chain.py` (repo).**
  - Its applicable-branch count now uses the scan's code.
  - The `expect` logic for zero-expectation key spaces is right.
  - The fallback branch count (max of gain + loss + undet) equals the true count whenever any variant is wholly unresolved, which is true in both cohorts.
- **Current audits.** Both pass. One WARN each, on small/masked.
- **Freshness of inputs.** `data/trees/<C>.*` and `events/` are newer than the 2026-10-01 VCFs, and `rrdr_carriers.txt` agrees with current calls.
- **Missing phenotypes.** None among tips. The phenotype is a carrier list; every non-listed tip, including H37Rv and the outgroup, is a non-carrier by design.

## Open questions

1. **Losses of known resistance mutations.** In gwas1000, rpoB S450L has 26 losses and katG S315T has 41; in scale200, katG S315T has 12. Reversions are implausible, so these look like Fitch placing gains high and then losses. Gains on high branches dilute `obs` with non-carrier descendants. This belongs to the event-matrix and polarity area; worth a check there. It also explains `carriers_with_phenotype` equal to 459 for units such as Rv2000.
2. **Sampling with replacement.** The branch and region nulls sample with replacement, while observed event sets are distinct branches. I could not compare this with phyoverlap2, whose `.py` source is absent.
3. **`sv:` records never reach any test.** In gwas1000, all 17,832 `sv:` caller records have 0 gains (9,319 DEL and 5,729 INS unresolved at the root), so no insertion SV is ever tested. The audit reports this as "accounted for". Is that intended, or an upstream genotype-missingness problem?
4. **Does ASSOC-1 change any survivor through the region null?** It needs a rerun with `canettii` pruned. UNVERIFIED.
5. **Is `canettii` in `cohort.gwas1000.tsv` and `cohort.scale200.tsv` deliberate** (a reference-bias control) for the genotyping chain? If so, it should be excluded only from the association tree.
