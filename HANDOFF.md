# Handoff: one rerun to bring production data up to date

Updated 2026-10-02. Repository `mtbc-graph-variants`, branch `main`
(protected: no force-push, no deletion). PR #1 was merged as `b58fc4f` and
the `review-fixes` branch deleted. The last commit is `58a6c73`. Section 2 is
committed. `analysis/` and the paper PDF are untracked.

**Both production cohorts are current**, and so are their association results
(section 0b).

## The rules

- **All production data must be produced by the current code.** No analysis,
  comparison or figure is valid until it is.
- **One rerun, not a series.** Every known fix and every open decision is
  settled before anything is submitted. As of this update, that is true.
- **scale200 first, then gwas1000.** scale200 is the validating cohort.
  gwas1000 is not submitted until scale200 has passed section 5.
- **Cost is approved before submission.** Fairshare is tight.

## Where things are

| what | path |
|---|---|
| Repository (all code changes go here) | `/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants` |
| Working tree (production data, scratch) | `/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome` |
| Durable mirror of the working tree | `/n/boslfs02/LABS/sfortune_lab/Lab/mchase/MtbPangenome` |
| Run root (repo code + working-tree data) | `mtbc-graph-variants/runroot` (gitignored) |
| Read collection (CRAMs) | set in `config/site.local.sh` (gitignored; copied from the working tree on 2026-10-01) |

Jobs are submitted from `runroot`, so they run the repository's code against
the working tree's data. Do not edit code in the working tree.

**The read collection is a colleague's transient copy on netscratch**, dated
2026-07-28, and netscratch is purged without warning. The rerun needs it for
accessory presence and the SV two-frame step, so do not let it wait.

## 0. Progress

**scale200 is current** (2026-10-01, finished 15:48). Every P3 to P5 file
postdates the fixes, two-frame and interval states exist for all 200
samples, and the VCF passed the gate. One P5 states task hit an NFS
stale-handle race (fixed in `0bc2a10`). After the requeue Slurm still
cancelled its dependents, so p5, p5svgt and p5vcf were resubmitted.

| record class | before | after | why |
|---|---|---|---|
| small variants | 78,727 | 79,429 | P4 and P5 review fixes |
| SV | 5,836 | 6,693 | 1,563 recurrent catalogue intervals plus 5,130 caller rows, 2,947 of them UNCATALOGUED deletions that used to be dropped |
| IS6110 | 944 | 778 | 461 rows moved to a shared key (4.5), seam reclassification (4.4), 5 DR rescues |
| accessory presence | 802 | 802 | unchanged |
| with AA | 34,464 | 34,908 | ancestral fix (3.7) |

Measured cost, using the billing weights from sacct: about 128 billing-hours
and 81 CPU-hours.

| step | billing-hours |
|---|---|
| SV projection, for the new catalogue's probes | 59 |
| p3 | 20 |
| p5states, both submissions | 15 |
| p4b | 14 |
| p4 | 6 |
| everything else | under 4 each |

IS6110 projection came from the store: the longest task took 26 s.

**Requests right-sized for gwas1000 (2026-10-01).** Peaks were measured with
seff and `/usr/bin/time`. The per-sample passes peak at about 0.65 GB, which
is odgi loading the graph, against 8 GB requested. P5 states peaks at 0.49 GB
on a gwas1000 sample, against 32 GB requested. Memory costs a quarter of a
CPU per GB on shared.

| step | new request | old request |
|---|---|---|
| p3, p4b, SV projection | 4 cores, 2 GB | 4 cores, 8 GB |
| p4 | 2 cores, 2 GB | 4 cores, 8 GB |
| p5states | 2 cores, 4 GB | 4 cores, 32 GB |
| p5keys, p5pre | 1 to 2 cores, 8 GB | 4 cores, 32 GB |
| SV genotyping per sample | 1 core, 2 GB | 4 cores, 8 GB |

At these requests, scale200 would have cost about 80 billing-hours instead
of 128. Estimate for gwas1000: roughly 250 to 350 billing-hours, about 100 of
them SV projection for its 150 references.

**DR placement: reverted.** A flank-consensus placement for rescued DR
copies was added after the detettore6110 comparison, then removed at the
user's instruction: detettore6110 was for comparison only, and none of its
method goes into the pipeline. The rescue is back to keying rescued copies at
the array's first base. The comparison's finding stands as a known
limitation. In 19 of gwas1000's 24 rescued isolates, the copy is H37Rv's own
DR copy (3,120,523 to 3,121,897), flanked by spacers the matched reference
lacks, so `h37rv:3119185` with H37Rv "empty" is the wrong key and state for
them. The affected passes (p1iv onward for both cohorts) were rerun on the
reverted code.

**I/O contract.** `refbias/io_contract.tsv` lived only in the working tree.
It still required the dense `matrix.tsv`, which P5 no longer writes, so it
blocked every submission that did not include p5. An updated copy is now in
the repository, and `runroot/refbias/io_contract.tsv` points at it.

**detettore6110 comparison.** On 14 gwas1000 isolates, copy number agrees
within 2. Of 79 H37Rv sites detettore6110 placed, 63 match ours within 10 bp,
and 8 more are in our calls under graph-node keys. See
`analysis/detettore_eval/`.

**gwas1000 is current** (2026-10-02): 200,444 merged VCF records, and the
VCF passed the gate.

## 0b. Association chains (2026-10-02)

Both cohorts' association chains were rerun on the current data.

### Trees: cohort plus CX333

Each cohort's SNPs were aligned together with the 333 CX333 panel genomes'
SNPs:

- **alignment:** `analysis/combined_tree/build_alignment.py`;
- **tree:** IQ-TREE GTR+F+ASC+G4, via `bin/build_snp_tree.sh`, rooted on
  ET1291 (`GCF_035581225`);
- **cut for the chain:** `analysis/combined_tree/prune_for_cohort.py` keeps
  the cohort isolates, H37Rv and the outgroup. The event reconstruction needs
  genotypes the panel genomes lack.

| cohort | combined tree | chain tree |
|---|---|---|
| gwas1000 | 1,330 tips | 999 tips |
| scale200 | 533 tips | 202 tips |

Trees are in `data/trees/<cohort>_cx333.*` and `data/trees/<cohort>.rooted.nwk`.
Tree checks: every lineage is monophyletic (`analysis/combined_tree/tree_checks.md`).

### Three chain omissions, found and fixed

Each one silently dropped part of the results. None raised an error.

| commit | omission | effect |
|---|---|---|
| `a76636c` | `cohort_assoc_tail.sh` never passed `--accessory-presence` to `assoc_scan.py` | no carrier-only (level-2) null for any accessory variant, and the callability floor was measured against the whole tree. gwas1000 lost 5 testable variants. |
| `a76636c` | `audit_chain.py` estimated a locus's applicable branches as twice its carriers minus one | after the scan fix, the audit still flagged one variant (`node:68834:0:A>G`, 83 carriers: 165 estimated branches, 244 in the tree) that the scan had correctly dropped. The audit now reads the tree with the scan's own code. |
| `58a6c73` | the gene burden was never given `--lineages`, and the small-variant and SV burdens (`--cls small`, `--cls sv`) were run by hand and never added to the chain | `p_lineage` was blank in every burden. The rerun archived `small_gene.tsv` and `sv_gene.tsv` and produced neither, so the resistance genes looked missing. |

The chain script and the audit are now in the repo:
`assoc/bin/cohort_assoc_tail.sh` and `bin/audit_chain.py`. The repo chain
script still calls the working tree's `assoc/bin/*.py` (relative to the
current directory), so it is run from the working tree. The launchers are in
`analysis/combined_tree/`:

- `assoc_tail.sbatch`: the whole chain;
- `burden_only.sbatch`: one burden, `CLS=is6110|small|sv`;
- `audit_only.sbatch`: the audit alone.

### Results

The phenotype is RRDR carriage (`assoc/<cohort>/rrdr_carriers.txt`): 459
carriers of 999 tips for gwas1000, and 57 of 202 for scale200. "Survive" means
q < 0.05 under all three nulls: branch, region and lineage.

| | gwas1000 | scale200 |
|---|---|---|
| variants with 2 or more independent gains | 11,322 | 3,607 |
| callable (80% or more of applicable branches determined) | 7,914 | 2,602 |
| variant-level survivors (`scan.tsv`) | 28 | 3 |
| small-variant gene burden, testable units / survivors (`small_gene.tsv`) | 6,457 / 18 | 5,189 / 6 |
| SV gene burden, testable units / survivors (`sv_gene.tsv`) | 340 / 0 | 179 / 0 |
| IS6110 gene burden, testable units / survivors (`is6110_gene.tsv`) | 263 / 0 | 77 / 0 |
| chain audit (`chain_audit.txt`) | pass | pass |

**Positive controls: the resistance genes are found.**

- **gwas1000, variant level.** The 28 survivors are:
  - rpoB: S450L, H445 and D435;
  - rpoC, at 764817 and 764840;
  - katG S315T;
  - inhA, and the inhA promoter at c-15t and t-8c;
  - embB: M306V/I, G406 and Q497;
  - gyrA: A90V and D94G;
  - rpsL: K43R and K88R;
  - rrs, including a1401g.
- **gwas1000, small-variant burden.** The 18 survivors are the same units as
  on 2026-09-28. Genes: rpoB, rpoC, rpoA, katG, inhA, embB, gyrA, rpsL, rrs,
  gid, pncA, ethA, thyA. Promoters: fabG1, ahpC, embA, pncA, eis. This test
  finds the loss-of-function genes (pncA, ethA, gid) that the variant-level
  scan cannot.
- **scale200, small-variant burden.** The 6 survivors are rpoB, rpoC, embB,
  gyrA, pncA and ethA. katG and the embA promoter narrowly miss on the
  lineage null (q 0.053 and 0.057, against 0.040 and 0.015 on 2026-09-28)
  but pass the branch null. This is a power limit at 57 carriers.
- **IS6110 burden.** It is not expected to find resistance genes. Only 3 of
  gwas1000's 2,020 placed insertions fall in a resistance gene, each with one
  origin. Its closest unit is ig:Rv2813-Rv2814c, the DR region, at q_branch
  0.064. Its origins are keyed at the DR array start, which is the known
  limitation.

rpoB is partly circular, because RRDR carriage defines the phenotype. The
other genes are co-resistance in MDR isolates.

**One warning, not a failure, in both audits:** small/masked drops 1% of its
records between the VCF and the event matrix (against 0% overall). It was
present before the rerun.

### Other checks on the current data (`analysis/`)

| check | result | where |
|---|---|---|
| RRDR calls against the reads | sensitivity 1.000, specificity 0.994 | `pipeline_checks/` |
| lineage barcode | 97.5% exact, no cross-lineage calls | `pipeline_checks/` |
| matched reference is the closest panel genome by SNPs | gwas1000: 979 of 997 (98.2%), the rest rank 2 to 5, at most 20 SNPs farther. scale200: 192 of 200 (96.0%), 6 rank 2 to 5, 2 rank 6 (exact ties in P1), at most 19 SNPs farther. No isolate matched across lineages. | `reference_match/` |
| accessory-genome variant audit, gwas1000 | per isolate, including IS6110 | `accessory_audit/` |
| IS6110 simulated-read benchmark against detettore6110 | comparison only; nothing adopted | `is6110_simbench/REPORT.md` |

### Cost

| step | billing-hours |
|---|---|
| gwas1000 + CX333 tree | 49 |
| scale200 + CX333 tree | 10 |
| chains, scan rerun, audit, burdens | 8 |

### Superseded outputs

All in the working tree, synced to durable storage:

- `assoc/archive/stale_pre_rerun_20261002/<cohort>/`: before the rerun;
- `assoc/archive/scan_no_accpres_20261002/gwas1000/`: the scan without the
  presence tables;
- `assoc/archive/burden_no_lineage_20261002/<cohort>/`: the IS6110 burden
  without the lineage null.

The current results are the files directly under `assoc/<cohort>/`. The last
sync was 2026-10-02, after the burdens.

## 1. State before the rerun (2026-09-30)

The review fixes (`deefef7`, `0a3f712`, 2026-09-29) changed P3, P4, P4b, P5
states, SV genotyping, IS6110 stage 2 and the merged VCF. On 2026-09-30 only
p1iv, p1is and p5vcf were rerun, on top of P3 to P5 outputs from 09-17 to
09-28. Today's merged VCFs mix code generations. The earlier "regenerated
scale200" was a test copy that has since been deleted.

| pass | status in all five cohorts |
|---|---|
| p1, p2, p1g, p1i | current (section 3) |
| p3, p4, p4b, p5 | **stale**: review fixes |
| SV genotyping, interval mode | **stale**, and was never part of the chain (section 2.3) |
| accessory presence | current for scale200 and gwas1000; **missing** for pilot, l7, scale100 |
| p1iv, p1is | **stale**: section 2.1 |
| merged VCF | **stale** |

Three more cohort roots, `l49`, `rerun` and `scale_h37rv`, are just as stale.
Regenerate or retire them; do not leave them looking current.

`analysis/is6110_paper_check.md` is marked provisional and must be
regenerated after the rerun.

## 2. What changed since the last commit, all tested

33 regression tests pass (`source config/project_env.sh; $MTB_PY tests/run_tests.py`).

### 2.1 IS6110

- **DR-array rescue** (committed in `5e13986`). It now also handles
  references with two DR clusters. Of gwas1000's 997 isolates, 24 gain a DR
  call and none of the existing calls change.
- **One seam rule (review 4.4, your decision: fix).** Promotion, flank
  placement and the VCF writer now share `is6110/bin/is6110_seam.py`, with a
  3 bp slop. Distances to the nearest seam peak at 0 to 3 bp, the
  target-site duplication, and are flat at about 2 per bp beyond that. On
  gwas1000, 260 sites become ref_shared and 187 more are placed in H37Rv.
- **One key per insertion (review 4.5, your decision: fix).** Keys of one
  insertion placed up to 6 bp apart across isolates now share one canonical
  position. 6 bp is where the excess over background ends: 29 against 12 at
  6 bp, 16 against 12 at 7 bp. Two keys carried by one isolate are never
  merged. On gwas1000, distinct H37Rv keys go from 2,024 to 1,501, and the
  singleton fraction from 66.4% to 59.8% (the paper reports 60%). Each row
  keeps its own placement in `h37rv_pos_placed`.
- **Projection reuse.** The p1is projection step, about 80% of an IS6110
  rerun's cost, now keeps results in `refbias/build/<id>/proj_is6110`.
  p1iv files the previous run's projections there before rewriting the key
  table. On a scale200 copy, all 84,232 keys came from the store, and the
  projection files were identical.

### 2.2 SV genotyping

- **Inherited calls need a second frame in matrix mode too (review 3.9,
  your decision: fix).**
- **SVLEN is the merged span (review 5.7, your decision: fix).** This changes
  `svi:` interval IDs, which is free now because everything reruns.

### 2.3 SV interval arm, now part of the chain

The merged VCF's deletion block comes from interval-mode genotyping. Before
today, three things were run by hand, outside the chain:

- the catalogue, built once from scale200's caller matrix into
  `refbias/assets/sv_intervals.tsv`;
- the two-frame tables, from `sv2frame/bin/sv_twoframe.sh` (working tree
  only);
- the interval-mode genotyping itself.

The result was that **gwas1000 was genotyped against scale200's catalogue, and
the merge dropped every gwas1000 caller deletion that catalogue lacked.**

Pass p5svgt now runs, per cohort, into `<root>/p5/`:

1. `--catalogue`: `sv_intervals.tsv`, built from graph deletions plus the
   cohort's caller deletions called in **2 or more isolates** (your decision,
   `SVCAT_MIN_CARRIERS=2`). For gwas1000 that is 1,132 graph and about 4,400
   caller intervals.
2. `--twoframe`: `twoframe/<s>.2frame.tsv`, from the read collection's H37Rv
   CRAM and the p1g alignment. It takes about 20 s a sample, and a missing
   p1g alignment is now an error.
3. probes, projection per reference, then interval genotyping. Each sample
   requires a two-frame table newer than the catalogue.

The merge now drops only the caller deletions a catalogue interval covers.
The others are kept presence-only and flagged `UNCATALOGUED`. Matrix-mode SV
states are used only when newer than `sv_matrix.tsv`, and interval states only
when newer than the catalogue.

Tested on one scale200 sample against the shared catalogue:

- the interval states are byte-identical to production;
- the two-frame H37Rv-frame columns are identical too.

### 2.4 Accessory presence and ancestral alleles

- **Accessory presence is part of pass p3** (`accessory/bin/locus_presence_array.sh`),
  writing `accessory/<cohort>/`. Existing tables are kept, because their
  output is byte-identical between the code generations; it was tested again
  on one sample. Writes are now atomic.
- **Ancestral alleles at the MTBC ancestor (review 3.7, your decision: fix).**
  4,759 of 5,130 tied sites are resolved, no previously resolved site
  changes, and 371 stay tied. The table is now a build asset
  (`bin/p0_prepare.sh --step ancestral`), and `p5_finish.sh` refuses to run
  without it.
- **Reference selection (review 6.11, your decision: leave).**

### 2.5 Smaller fixes

- `bin/archive_alignments.sh` now keeps and archives the p1g alignment, which
  the two-frame step reads. It used to call it unused and could delete it.
- `config/site.local.sh` was missing from the repository config, so the run
  root had no read collection. It is copied in (gitignored).

## 3. What does not need rerunning

P1, P2, P1g and the P1i scan: the review only made their failures fatal and
their writes atomic. None of the 7,180 P1 and P2 logs shows simulated reads.
FASTQ truncation (review 6.5) cannot be checked, because the FASTQs are
deleted; treat it as a residual risk. Every cohort has its P2 BAM and P1 VCF
for every sample.

The rebuilt frame table is identical to the fallback every run has used, and
the cached P5 projections stay valid.

## 4. The rerun

### 4.1 Once, before any cohort (writes into the build; seconds)

From `runroot`:

```
B=refbias/build/7713a8d71d8e
OG=$(awk -F'\t' '$1=="graph"{print $2}' $B/build_info.tsv)
BUILD_ID=7713a8d71d8e OG=$OG bash bin/p0_prepare.sh --step ancestral
BUILD_ID=7713a8d71d8e OG=$OG bash bin/p0_prepare.sh --step frames
```

### 4.2 Per cohort, scale200 first

The per-sample passes skip samples whose outputs already exist, so move the
stale ones aside first. On scratch, `mv` is instant and frees no space.

```
C=scale200; R=refbias/scale200
A=refbias/archive/stale_pre_rerun_$(date +%Y%m%d)/$C
mkdir -p "$A"
for p in p3 p4 p4b p5; do mv "$R/$p" "$A/$p"; done
bash bin/refbias_run.sh "$C" --from p3 --dry-run     # check every job and dependency
bash bin/refbias_run.sh "$C" --from p3
```

Outroots: pilot is `refbias`, scale100 is `refbias/scale`, l7 is `refbias/l7`,
scale200 is `refbias/scale200` and gwas1000 is `refbias/gwas1000`.

`--from p3` submits:

- p3, with accessory presence;
- p4, p4b and p5;
- p1g and p1i, which keep their current outputs;
- p1iv, which files the old projections, runs the DR rescue and recalls with
  the new seam and key rules;
- p1is, with projection reuse;
- p5svgt, the interval arm;
- p5vcf, which ends with `vcf_gate.sh`.

### 4.3 Cost

Not measured for P3 to P5, and Slurm accounting before 2026-09-30 is not
available. These parts are known:

| part | scale200 | gwas1000 |
|---|---|---|
| IS6110 projection | near zero, from the store | near zero, from the store |
| SV two-frame | about 1 CPU-hour | about 5.5 CPU-hours |
| SV catalogue projection | near zero | about 85 CPU-hours, one-time |
| accessory presence | skipped, current | skipped, current |

Use `sacct` on scale200's run to measure P3 to P5, and scale by 997/200 for
gwas1000 before submitting it.

## 5. Checks before calling a cohort done

- `bash bin/refbias_run.sh <cohort> --status` shows every pass complete.
- Every file under `<root>/p3`, `p4`, `p4b` and `p5` is newer than the
  commit that carries section 2.
- The p5vcf log shows the VCF gate passed, a catalogue of the expected size,
  and the UNCATALOGUED count.
- Diff cell counts against the archived generation, and explain the changes
  using `CODE_REVIEW.md` "Behaviour changes to expect on re-run" and
  section 2 here.
- For scale200 only, before gwas1000: confirm the IS6110 projection step sent
  almost nothing to odgi, and that two-frame and interval states exist for
  all 200 samples.

## 6. After the cohorts are current

1. Done 2026-10-02: the paper comparison was regenerated on the current
   gwas1000 keys (`analysis/is6110_paper_check.md`).
2. Done 2026-10-02: synced to durable storage with `bin/sync_back.sh`.
   Durable storage holds 63 GB, and space is limited, so check sizes before
   the next sync.
3. Archive alignments with `bin/archive_alignments.sh --archive` if scratch
   space is needed.
4. Delete the archived stale generations once the new ones are checked:
   `refbias/archive/pre_review_fixes_20260930/` and the
   `stale_pre_rerun_*` directories.
5. Retire `refbias/assets/sv_intervals.tsv` and `sv2frame/<cohort>/`. Nothing
   reads them now.

## 7. Other open items

- **Untracked:** `analysis/` and the paper PDF. `analysis/dr_test/` is
  throwaway.
- **The association chain's Python modules** (`assoc/bin/*.py`) exist only in
  the working tree. The repo's chain script and audit call them there.
- **Outgroup review:** the output of `check_high_mac_no_panel.py` has not
  been reviewed. The old manual patch, which set 6 outgroup columns to N, was
  not re-applied.
- **The 1% small/masked warning** in both chain audits predates the rerun and
  has not been explained.
- **Not started:** the per-reference coordinate map from the 10k scaling
  plan, and a review of about 15 GB of "possible junk".
- **The 50 (sample, key) pairs** that occur twice in gwas1000's key table
  predate today's changes. They are worth a look, but they do not block.

## 8. Working rules for the next session

- Code changes go only in the repository. The working tree is data, written
  only by approved pipeline runs.
- Conda environments go under the home directory. `gh` is in `~/.conda/envs/gh`.
- Any output-affecting code change makes the affected passes stale for every
  cohort. Record it here, and settle it before the next rerun, not after.
