# Handoff: one rerun to bring production data up to date

Updated 2026-10-05. Repository `mtbc-graph-variants`, branch `main`
(protected: no force-push, no deletion). `analysis/` and the paper PDF are
untracked; `analysis/` lives on durable storage with the repository.

## Start here (2026-10-05)

**Audit (section 0m): 97 findings, and the current association results are
provisional.** Steps 1-3 (code into the repository, then the section A and B
fixes, with tests) were approved on 2026-10-05.

**Production state.**

- **scale200 and gwas1000 are current with the code,** and so are their
  association results (section 0b).
- **Two faults found on 2026-10-04 are open and unfixed** (section 0j). Both
  touch event polarity at about 0.4 to 1% of sites:
  - **A:** `assoc/bin/add_outgroup.py` keeps the last duplicate decomposed
    record, so the canettii outgroup column is wrong at about 500 sites per
    cohort;
  - **B:** `bin/ancestral_alleles.py`'s `AA` node includes the second
    canettii, which changes 132 lineage 1-4 variable sites.
- **The user's decision is pending:**
  - (a) fix at the rebuild;
  - (b) measure the effect first;
  - (c) fix now with one planned rerun, scale200 first.

  **Do not fix or rerun anything without that decision.**

**Forward work (discussion only, nothing built):**

- **the lineage 1-4 panel proposal:** `analysis/strategy/L1_4_PANEL_PROPOSAL.md`,
  sections 0g, 0i, 0j and 0k. The outgroup topology was corrected in 0j:
  lineages 5 and 6 are the nearest outgroups. Section 0k recommends lineage
  5/6 inside the graph and lineage 8 and canettii outside it;
- **the supporting analyses:**
  - external assembly QC (0f);
  - the GenBank scan (0h);
  - the panel checks (0j, `analysis/panel_checks/`).

**State (2026-10-06): section 0o.**

- **Every audit fix is on branch `audit-fixes`** (283 tests), not yet in
  main. The current association results remain provisional.
- **Waiting on the user:** the 44 decisions, and D41.
- **Step 4 is done:** section 0q lists 12 items to fix before any rerun.
  **All 12 are fixed on `audit-fixes`, and D41 (b)** (293 tests). Waiting on the user:
  the decisions, the canettii isolate and R2-TREES-6, before the merge into
  main.
- **Then:** the new panel, graph and build, then scale200 first.

**Before anything else after downtime:**

1. **Check that netscratch survived.** The working tree and the colleague's
   CRAM copy (below) live there.
2. **Compare the working tree with its durable mirror** if anything looks
   missing. The last `bin/sync_back.sh` ran on 2026-10-05, before the
   maintenance.

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
| SNPs called against the matched reference (P2) against H37Rv (P1), same reads | gwas1000: median 302 against 1,636, every isolate fewer, median per-isolate reduction 6.4x. scale200: 212 against 1,095, every isolate fewer, 5.8x. Largest for lineage 2 (12.6x, 13.6x), smallest for lineage 4 (3.8x, 4.5x) and gwas1000's lineage 9 (3.3x). Figure `snp_distance.png`. | `reference_match/` |
| accessory-genome variant audit, gwas1000 | per isolate, including IS6110 | `accessory_audit/` |
| IS6110 simulated-read benchmark against detettore6110 | comparison only; nothing adopted | `is6110_simbench/REPORT.md` |

### RD cross-reference against Behruznia et al., eLife 2024 (`analysis/rd_crossref/`)

This reruns the working tree's `accessory/bin/rd_crossref.py`, unchanged, on
each cohort's current catalogue (`refbias/<cohort>/p5/sv_intervals.tsv`) and
interval states. The 2026-09-28 result was on scale200's pre-rerun catalogue.

| | scale200 09-28 (stale) | scale200 now | gwas1000 now |
|---|---|---|---|
| known RDs with a catalogued interval | 129 / 135 | 123 / 135 | 131 / 135 |
| RDs under 1 kb | 40 / 44 | 37 / 44 | 42 / 44 |
| polymorphic insertion loci that define a lineage | 0 of 57 | 0 of 57 | 0 of 110 |
| RD lineage concordance: agree / unresolvable / mismatch | 33 / 3 / 19 | 41 / 3 / 20 | 37 / 4 / 29 |

- **Both of the paper's predictions hold.** Deletions recover the known RDs,
  and insertions do not define lineages.
- **scale200 lost 6 RDs.** All were caller-only singletons that the
  recurrent-only catalogue (`SVCAT_MIN_CARRIERS=2`) now excludes.
- **The concordance mismatches predate the rerun.** Most share their interval
  with another RD (16 of 29 in gwas1000) or sit near IS6110 (18 of 29). The
  mismatches have not been resolved per RD.

### Stably inherited deletions that are not known RDs (`analysis/rd_crossref/novel_stable.py`)

These are single-origin, never-reverted deletion intervals (at least 2
carriers, 80% or more of branches resolved) with no overlap with the 135
known RDs.

| | gwas1000 | scale200 |
|---|---|---|
| single-origin, never reverted | 167 | 101 |
| not a known RD | 89 | 51 |
| of those, carriers form exactly one clade | 73 | 42 |
| in several panel assemblies | 19 | 30 |

Lineage-wide examples in gwas1000, all in several panel assemblies:

- 135 bp in rho: all of lineage 6 and lineage 9, nothing else;
- 119 bp in Rv0725c/sppA: all of lineage 5;
- 63 bp in ctpG: all of lineage 3, plus 1 lineage 4 (probably a miscall);
- 105 bp in Rv0209: 43 lineage 1.2.1.2;
- 4.75 kb over Rv1353c to Rv1356c: 27 lineage 4.1.2.1.

Caveats:

- Most of the 89 are under 1 kb with 2 to 4 carriers.
- Many are in PE/PPE genes.
- 51 are caller-only.
- "Not a known RD" means not in that 135-RD list.

**Literature check (2026-10-03, `analysis/rd_crossref/LITERATURE.md`).**

Sources: a web search, the Bespiatykh et al. 2021 PDF in the working tree, and
an overlap test against RDscan's `RD.bed`. That file is the paper's full
curated set of 79 lineage-specific RDs, including its new RD301 to RD317 and
RD743.

Already described:

- the 4.75 kb deletion (Rv1353c to Rv1356c, lineage 4.1.2.1) is **HSD3, also
  called RD145**, a Haarlem-specific deletion. It is missing from our 135-RD
  list;
- the nrdZ deletion matches **RD121**;
- the 3.5 kb lineage 2.2.1 deletion lies inside **RDoryx_4**, an independent
  deletion of the same region.

No report found for the five lineage-wide hits, all of them in several panel
assemblies:

- rho (Rv1297): 135 bp, in-frame (codons about 115 to 159), lineages 6 and 9.
  rho is essential, and this probably falls in the mycobacteria-specific
  N-terminal insertion domain (inferred from position);
- ctpG: 63 bp, in-frame, lineage 3;
- Rv0209: 105 bp, in-frame, lineage 1.2.1.2;
- the 3' ends of sppA and Rv0725c: 119 bp, all of lineage 5;
- the 3' ends of arsB2 and Rv3579c: 104 bp plus 65 bp, 41 lineage 5
  isolates.

Bespiatykh's method calls only deletions over 200 bp (low-coverage regions
of 100 bp or more). All five are 63 to 135 bp, so their absence from that map
is expected. The graph and interval approach finds lineage markers that
coverage-based RD scans miss.

The ctpG to Rv1996 block is a deletion hotspot: RD743 (lineage 5), RD174
(lineage 4.3.4), and our lineage 3 in-frame ctpG deletion about 570 bp
upstream of RD743.

**Supplementary tables, added to the working tree 2026-10-03.** These are
Bespiatykh's Tables S1 to S3 (`msphere.00535-21-st00*.xlsx`) and Behruznia's
files 9 to 11 (`584580_file*.xlsx`).

- **Bespiatykh Table S3** (the complete RD list, 213 rows, 187 with
  coordinates):
  - None of the five lineage-wide hits overlaps any RD.
  - 13 of gwas1000's 89 overlap a listed RD, mostly as small deletions
    inside larger historical RDs.
  - The full matches are RD145 (HSD3), RD121, RDoryx_4 and DS20 (2.27 kb,
    PPE37 and metH, lineage 4.1.1.3).
- **Behruznia file 11** (the RDs and genes per lineage, from Panaroo and
  Pangraph, matched by locus tag):
  - scale200's 1,449 bp PE_PGRS53 (Rv3507) deletion in 7 lineage 4.2
    isolates is their new **"New-L4.2"**. It has the same gene and lineage,
    an independent confirmation of a deletion our method found.
  - None of the five lineage-wide hits has a gene in their tables.
- **Conclusion.** The five lineage-wide deletions are not described in either
  paper's full tables, and both methods are blind to them by design:
  - Bespiatykh calls only deletions over 200 bp, and these are 63 to 135 bp;
  - Behruznia scores whole-gene presence, and these are in-frame losses of
    21 to 45 codons, or removals of 3' gene ends.

  They are candidate previously unreported small lineage markers. Next steps:
  confirm them by PCR or against long-read assemblies, and search the
  literature beyond these two papers.

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

## 0c. Scaling to 10,000-40,000 samples (assessment, 2026-10-03)

**10,000 is feasible after a known set of engineering changes. 40,000 is not
feasible as the pipeline stands.** Per-sample compute is no longer the limit.
The limits are the cohort-wide steps, the association tree, job submission and
storage.

Figures are measured on gwas1000 and scale200 unless marked as estimates.
Projections come from the working tree's `refbias/SCALING_MEASURED.md` and
`refbias/WHY_SCALE200_DID_NOT_PREDICT_GWAS1000.md`. Some of their open items
were checked against the current code; parts of those documents predate it.

### What scales well

- **Per-sample passes are linear, with flat memory.** The gwas1000 rerun from
  P3 onward cost 348 billing-hours (330 CPU-hours) for 997 isolates, about 0.35
  billing-hours per sample. It covered P3 to P5, IS6110 stage 2, SV
  genotyping and accessory presence. No per-sample task ran longer than
  15 minutes.

  | step | billing-hours |
  |---|---|
  | SV projection | 158 |
  | p3 | 66 |
  | p4b | 49 |
  | p5states | 19 |
  | p5sv2frame | 16 |
  | p4 | 14 |
  | everything else | under 7 each |

- **P5 states is fixed.** It was 74% of the compute in older measurements
  and grew as n^1.39. With the projection cache it took 19 billing-hours on
  gwas1000, with a longest task of 1.3 minutes.
- **The dense `matrix.tsv` is gone.** P5 writes sparse per-sample states. The
  merge and the VCF took minutes at 1,000 samples (p5vcf 3.6 min).
- **SV projection is per reference** (up to 333), not per sample.

### Estimated cost from raw reads

The per-sample passes cost about 1 to 1.5 CPU-hours per sample, including P1
and P2, which the rerun reused. The docs project about 18,000 CPU-hours at
10k.

| cohort | compute (estimate) | scratch disk |
|---|---|---|
| 10k | 15,000 to 20,000 CPU-hours, about 10,000 billing-hours | about 5.5 TB |
| 40k | about 60,000 to 80,000 CPU-hours | about 22 TB |

The disk figures use gwas1000's 0.55 GB per sample: 282 GB of outputs plus
270 GB of P1 and P2 alignments. Archived as CRAM, alignments take about
90 MB per sample.

### What breaks

1. **Job submission.**
   - The cluster limit is 10,100 submitted jobs, counting array tasks, so one
     10k pass already hits it. The runner must submit in chunks.
   - A single failed task cancels every afterok dependent, and at 40k failures
     are certain. Dependencies need automatic retry first.
   - At the default throttle of 50, 10k takes about 28 days. It needs 200
     to 500.
2. **Records grow as about n^0.5.** They went from 87,702 at 200 samples to
   200,444 at 1,000. Estimates: about 650k at 10k, about 1.3M at 40k.
   - **The states array** (samples × keys, uint8) grows from 0.2 GB to about
     6.5 GB at 10k and about 50 GB at 40k. At 40k it needs memory-mapping or
     a sparse format.
   - **`p5_validate`** (called from `p5_finish.sh`) compares every pair of
     samples, so its cost grows with the square of the cohort. It should
     switch to a sample of pairs.
3. **The association chain is the hardest wall.**
   - The event matrix is variants × branches: 200k × 2k at 1,000 samples,
     about 1.3M × 80k at 40k (about 10^11 cells).
   - The scan's descendant matrix and permutation tables each reach several
     GB at 40k.
   - This needs a sparse or chunked redesign, not a larger memory request.
4. **The tree.**
   - IQ-TREE GTR+F+ASC+G4 with 1,000 UFBoot cost 49 billing-hours for 1,330
     tips and 10 for 533, growing faster than the tip count.
   - At 10k to 40k tips this is impractical. The options are a faster method
     (FastTree or VeryFastTree) or placement onto a fixed backbone (UShER or
     MAPLE).
   - It is a scientific decision, because the tree drives the convergence
     tests.
5. **Storage.**
   - About 22 TB of scratch at 40k, on netscratch, which is purged without
     warning.
   - Durable storage is limited (63 GB used now), so alignments must be
     archived as CRAM batch by batch.
   - The read collection is a colleague's transient copy.
6. **The reference panel, which no cohort size fixes.**
   - The QC-pass pool is 47,815 isolates: lineage 4 24,783, lineage 2 14,335,
     lineage 3 4,950, lineage 1 3,406, lineage 6 176, lineage 5 116.
   - The 333-genome panel is thin for lineage 1, where matched references are
     up to about 900 SNPs away, and for lineages 5 and 6.
   - A 40k cohort would mostly add lineages 2 and 4.

### Recommendation

1. Run a 2,000 to 3,000 cohort as a scaling test before 10k, recording Elapsed
   and MaxRSS for every step and fitting the curves. scale200 was a
   correctness test, not a scaling test: gwas1000, 5 times larger, cost 11 to
   55 times more in three steps and broke four.
2. Before 10k: chunked, retry-aware submission and a higher throttle (item 1);
   P5 at scale (item 2); a decision on the tree method (item 4).
3. Before 40k: also the association chain redesign (item 3) and a storage
   plan (item 5).

### Measured costs and cheaper ways to run the current methods (2026-10-04)

Full analysis: `analysis/scaling/COSTS_AND_IMPROVEMENTS.md`. gwas1000 (997
samples) cost about 2,130 billing-hours, about 2.1 per sample. At current
requests 10k would be about 21,000.

| step | billing-hours | requested | used |
|---|---:|---|---|
| P2 (matched-reference alignment and calling) | 932 | 8 cores, 16 GB | 1.6 cores, 1.3 GB (max 6.3) |
| P1 (H37Rv alignment and reference selection) | 597 | 8 cores, 16 GB | 1.9 cores, 1.4 GB (max 2.4) |
| SV projection (per reference; 150 used) | 158 | 4 cores, 2 GB | right-sized |
| P1g, P1i (IS6110 stage 1) | 101 + 100 | 4 cores, 8 GB | 66% CPU, 1.0 GB |
| P3 to P5, tree, association, the rest | about 270 | mostly right-sized | |

P1 and P2 are 72% of the cost and use about a quarter of their cores and a
tenth of their memory.

**Proposed, in order. Nothing is changed yet.**

| | change | saving | method change? |
|---|---|---|---|
| A | P1 and P2 to 4 cores, 4 GB | about -54% on them, about 8,000 billing-hours at 10k | no |
| B | P1g and P1i to 2 GB | about -25% | no |
| C | serial_requeue for short per-sample tasks | half the billing again | no, but needs E first |
| D | a cheaper P1 reference selector (panel-informative sites from a read subsample, or k-mer sketch distance) | most of P1 | **yes**: inventory everything that reads P1's VCF, then validate against gwas1000's refmap |
| E | several samples per array task, retry before releasing dependents, throttle 200 to 500 | needed at 10k regardless (one 10k pass hits the 10,100-job limit) | no |
| F | precompute per-reference products once for all 333 references | new cohorts pay nothing for them | no |
| G | reuse per-sample results across cohorts, keyed by (sample, build version) | | no; valid only while the code is unchanged |
| H | CRAM-archive alignments per batch | 5.5 TB down to 1.5 TB at 10k | no |
| I | the cohort-level steps (states array, validator, event matrix, tree) | | design work |

**F, measured 2026-10-04: per-reference precompute is not worth it.**
Precomputing the full position set costs about 6.9 billing-hours per reference
(section 0f), about 1,200 for the remaining 173. Keep the lazy store.

**J (idea, a method change, not done): projection from pairwise alignment
instead of `odgi position`.** odgi walks the 540 MB graph for every position
(1.7 h per reference). A whole-genome alignment of each reference to H37Rv
(minimap2 asm5) takes about 1.5 s and gives most of the same coordinate map.

- Validate against the 159 filled references in the store, outside complex
  regions.
- If they agree, P5 key and SV-probe projection costs drop by orders of
  magnitude.
- Needs full validation and approval before it goes near the chain.

**Estimated cost per sample:**

| changes | per sample |
|---|---|
| current | 2.1 |
| A + B | about 1.0 |
| plus C | about 0.5 |
| plus D | about 0.35 |

**Order:** A and B, measured on scale200 first; then E; then F, G, H; then D
after validation; then I before any cohort above about 2,000.

### Option: one linear CX333 reference instead of per-isolate matched references (2026-10-03)

Raised for scaling beyond about 1,000 samples. Discussed, not built.

**The idea.** Build one linear reference for bwa from the CX333 graph:

- **Backbone:** an ancestral sequence reconstructed from the graph, rooted on
  ET1291 (M. canettii). It would be a local, QC'd equivalent of MTBC0.
- **Extra contigs:** the non-redundant accessory sequence above a frequency
  floor (say 5% of genomes or more), as decoy or ALT contigs.
- **Coordinates:** a graph-derived liftover to H37Rv.

**What an imputed (ancestral) genome buys:**

- it sits mid-tree, so it is closer to every lineage than H37Rv (MTBC0 gave a
  17.6% smaller truth set in REFEVAL);
- it restores sequence H37Rv's branch lost (TbD1, RvDs);
- polarity is native;
- one coordinate system, with standard tools.

**What it costs, all seen in MTBC0 (section 0d):**

- input assembly errors averaged in invisibly (PPE38 lost);
- collapsed copy-number regions (one truncated IS6110, a CRISPR arrangement
  no real genome has);
- none of the derived accessory sequence. CX333 holds 366 kb MTBC0 lacks.

A CX333-built backbone avoids the input-quality problem: 333 QC'd assemblies,
IS6110 at real copy number, provenance for every position.

**Expected performance, from the existing measurements:**

- **Core SNPs: no gain.** Every linear reference tested is within about 1.5
  points (REFEVAL, corrected; sensitivity and PPV):

| reference | SNP sensitivity | SNP PPV |
|---|---|---|
| H37Rv | 0.958 | 0.896 |
| MTBC0 | 0.942 | 0.883 |
| H37Rv + accessory contigs | 0.954 | 0.907 |

  MTBC divergence is about 0.05%, so bwa places reads against any of them.
  Sensitivity is limited by PE/PPE repeats and by GATK being blind at 500 bp
  and above (T1). A backbone fixes neither.
- **Accessory: a gain.** Ancestral regions become mappable, and accessory
  contigs give T4's 0.885 / 0.912 in novel sequence.
- **False-positive suppression: smaller than now.** Per-isolate matched
  references cut apparent variant burden by 88% on real reads (T8). A
  mid-tree reference achieves only part of that (MTBC0: 17.6%).

**Why it matters for scaling: corrected 2026-10-04, it mostly does not.** A
single reference would remove every per-reference cost:

- **SV projection:** run per reference, it was 158 of gwas1000's 348
  billing-hours (section 0c);
- **P1 reference selection;**
- **the per-reference P3/P4 coordinate projection** and the per-reference
  IS6110 projection store;
- **one bwa index** instead of one per panel genome.

But those costs are bounded by the panel (at most 333 references), not by the
number of samples, and can be precomputed once and reused across cohorts. They
are not a scaling bottleneck. The 10k to 40k blockers in this section are
independent of the reference choice. Design, assumptions and recommendation:
`analysis/refeval_cx333/DESIGN_NOTES.md`. **Recommendation: keep matched
references.** The CX333 linear reference is optional research only. Its
backbone is built (`analysis/refeval_cx333/ref/`); no arm was submitted.

**To decide it.** Build the backbone and accessory contigs (cheap), then run
the REFEVAL harness: 35 genomes, assembly-derived truth, only the reference
varied. Compare it against H37Rv, MTBC0 and the matched references, and
measure the per-sample cost of each arm. Cost the evaluation before
submitting; it is on the scale of the earlier REFEVAL run. Not started.

### Lab libraries reviewed for scaling (2026-10-03, read-only)

Two of Peter Culviner's libraries were read, not run, to see whether they
help with the items above. Neither is a drop-in fix.

**samarray** (`/n/netscratch/sfortune_lab/Lab/mchase/samarray`, v0.1.0, last
commit 2026-08-27). "SAM" means SAM/BAM/CRAM; it is not a job-array tool. It
stores a cohort's calls, alignments, ancestral states and events as TileDB
arrays.

| item | does it help? |
|---|---|
| 1. Job submission | Partly. Work runs on a few long-lived Dask workers on Slurm, so the job count stays small. But one failed task stops the run, it cannot resume, and it runs only its own Python functions, not shell passes. |
| 2. States array, pairwise validator | Partly. `VariantArray` stores calls sparsely, with zstd compression. Queries load dense into memory unless chunked. There is no pairwise-distance code. |
| 3. Event matrix | Partly. `EventArray` stores per-branch events sparsely; `AncestorArray` runs chunked multistate Fitch with bounded memory. There is no permutation or association code. |
| 4. Tree | Indirectly. It builds no trees, but `writeSNPAlignment` exports MAPLE format (single contig only), and MAPLE is built for very large datasets. |
| 5. Storage | No. `CRAMArray` copies existing CRAMs into its own store. |

Risks:

- no test suite (removed as stale);
- single-author research code, still changing;
- runs only in its author's separate conda environment;
- pickled metadata, which ties arrays to a Python version;
- some functions raise NotImplementedError.

It is used once in the working tree (`is6110/bin/is6110_clipend_scan.py`,
for its Dask depth computation). `bin/ancestral_alleles.py` and
`is6110/bin/is6110_junctions.py` chose not to adopt it.

**mtbvartools** (v0.1.0dev, 2024-2025, no tests) is a toolkit for download,
mapping, calling, VCF merge, PastML reconstruction and tree metrics.

- **Item 1, job submission.** `scripts/slurm_dask_splitter.py` runs one shell
  command per sample on a small pool of Dask workers on Slurm. A failed task
  cancels nothing else. But the scheduler must stay alive for the whole run,
  and there is no retry or record of what finished.
- **Item 3, event matrix: the weak point.** The event matrix is stored in its
  `CallBytestream`, a pair of `KeyedByteArray` files.
  - Every row is dense and zlib-compressed on its own.
  - The matrix is stored twice, once in each orientation.
  - The whole index is loaded into memory.
  - Reading a column decodes every row.
  - A file cannot be appended to once closed.
  - This is workable at 200k × 2k, not at about 1.3M × 80k.
- **Reconstruction and tree metrics.** `writeAncestorCalls` runs PastML in
  100-variant blocks, very slow at 40k. The tree metrics are roughly quadratic
  in tip count.
- It builds no trees and has nothing for storage.

**Reproducibility issue.** The association environment (`MTB_PY_VT`,
`mtb_isolates_cluster`) does not load the netscratch copy of mtbvartools. It
loads an editable install at `/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbvartools`
with uncommitted changes: `conversions.py` adds `writeNodeStates`,
`writeAncestralFastas`, `writeNodeVariantStates_stream` and
`make_informative_variant_index`, and there are new PastML scripts. The event
matrix depends on code in no commit.

**Conclusions:**

1. **Job submission.** Both offer the same Dask worker-pool idea, without the
   retry, resume and per-sample completion records `refbias_run.sh` relies on.
   The smaller change is in our runner: each array task processes K samples
   and writes a done-marker per sample.
2. **Cohort-wide data at 40k** (states, events, ancestral states). samarray's
   sparse TileDB arrays are the right kind of design. Adopting them is a
   substantial rewrite onto an untested dependency; they can also serve as a
   model for our own sparse format.
3. **The tree.** samarray's MAPLE export points to a realistic option for 10k
   to 40k tips.
4. **Storage.** Neither library helps.

## 0d. CX333 pangenome against Behruznia et al. 2024 and MTBC0 (2026-10-03)

Analysis only, in `analysis/pangenome_compare/`. Read `README.md` there for the
whole account.

**Source.** The CX333 graph's VCF against H37Rv
(`graphs/CX333.s10k.k23.K15/all_variants.decomposed.vcf.gz`), reduced to
structural events of 250 bp or more (`compare.py`). Three corrections were
needed:

1. **All snarl levels.** The classic RDs are nested records.
2. **Pooling.** The file is decomposed (one record per allele sequence), so
   records are pooled by snarl, position, kind and size.
3. **Parent and child.** A parent and a nested record of the same event count
   once.

The supplementary tables are in the working tree as `584580_file08` to
`file11`.

**Genome overlap** (`genome_overlap.py`):

- Their 339 against CX333's 333: **186 shared**, by BioSample.
- 139 of theirs are not in RefSeq at all; mostly the Marin 2024 Zenodo
  long-read assemblies.
- The rest of theirs are contig-level, excluded by CX333's QC, or M. orygis,
  which the panel's organism-name selection never listed.
- Their table has two BioSamples each listed under two lineages.

**Their RDs in CX333.** Their file 11 RDs fall in the right lineages:

- RD239, RD105, RD750/RD316 and RD312 cover all of lineages 1, 2, 3 and 5
  respectively;
- RD702/RD303 cover all of lineages 6 and 9;
- RD9 and RD10 cover lineages 5, 6 and 9 and the animal lineages;
- RD711 is in 5 of 8 lineage 5;
- ND1 is in lineage 7.

RD207 is in the CRISPR region, which the graph represents as a complex snarl.
The La3 RDs cannot be tested: CX333 has no La3.

**The accessory genomes** (`gene_detail.py`; `accessory_detail.md`):

- **Theirs:** 394 accessory genes (Panaroo, corrected), 298 accessory
  Pangraph blocks, 111 and 116 sub-lineage-specific genes, and no LSPs.
- **CX333, H37Rv genes** (absent when 50% or more is deleted): 176 accessory
  genes, 48 sub-lineage-specific absences (the classic RDs). Mobile elements
  (38%), ESX (24%) and PE/PPE (16%) are the most often accessory classes.
- **CX333, by lineage:** absent genes per genome are highest in lineage 6 and
  the animal lineages, matching their Figure 1. Their "L2 smaller than L1, L3
  and L4" is not reproduced.
- **CX333 insertions** (1,159 events against H37Rv):

| class | events |
|---|---|
| IS6110 | 509 |
| copies or expansions of H37Rv sequence | 461 |
| ancestral (in MTBC0) | 43 |
| novel | 124 (100 in a single genome) |

**Conclusion: the paper's main claim holds.** Of 77 lineage-specific
insertions, 53 are IS6110 and 17 are H37Rv-sequence copies. The 7 others are
all PE/PPE loci: lineage-specific alleles at PPE50/51, PPE38 and PPE57 to 59.

- No lineage gains new genes.
- A 1,442 bp element (4 copies in H37Rv) has an extra copy in all of lineage 1.
- About 3.1 kb near pknH, in neither H37Rv nor MTBC0, is carried by lineages
  1, 5 and 6 and the animal lineages. It is unidentified.

**Withdrawn.** An earlier pass reported the pknH insertion as specific to
lineage 1.1.3.1, and named two "candidate genuine insertions". Both were wrong.
The specificity test skipped genomes the barcode cannot place; it now
disqualifies them.

**MTBC0.** The imputed MRCA of Harrison, Kapur and Behr 2024 was already
rejected as a backbone (`REFERENCE_BIAS.md` section 3b-bis, REFEVAL). Measured
against all 333 assemblies:

- **MTBC0 holds only 153 bp that no CX333 genome has.** It is a DR plus spacer
  arrangement inside the CRISPR array.
- **MTBC0 has one, truncated IS6110** (1,167 of 1,355 bp).
- **CX333 holds 366,587 bp that MTBC0 lacks**, in 767 sequences:

| class | bp |
|---|---|
| novel | 232 kb |
| IS6110 fragments, where MTBC0's copy is truncated | 97 kb |
| H37Rv-sequence copies | 22 kb |
| H37Rv itself, mostly the DR region | 9 kb |
| partly ancestral | 6 kb |

  They are written out in `analysis/pangenome_compare/mtbc0_missing/`, all
  together and also split with and without IS6110, plus a TSV index. MTBC0 is
  essentially a subset of CX333.

**Open:**

- identify the pknH-region sequence and the 1,442 bp and 1,458 bp elements;
- verify, against the source assemblies, that the lineage-specific PPE alleles
  are not assembly artifacts.

The alignments `gene_detail.py` and `mtbc0_missing.py` read are in the
session's scratch directory, not in the tree.

## 0e. The insertion gap: assessment, not scheduled (2026-10-04)

Full notes: `analysis/insertions/INSERTION_GAP_ASSESSMENT.md`, assessing
`insertion_gap.txt` (repo root).

**Checked against the pipeline:**

- **The IS6110 arm already does the IS-element trick.** Native copies are
  excised, and the element is kept as its own contig, so junction reads clip and
  their SA tags point at it.
- **Unmapped reads are kept** in P2 BAMs and in the CRAM archives.
- **The read collection is a transient netscratch copy.** That is the input
  risk for any read-level method.
- **Nothing calls a new copy of sequence H37Rv already has, at its new site.**
  This appears correct and should be verified before building anything.
- **The text's class counts are panel against H37Rv.** The gap that matters is
  relative to each sample's matched reference; measure it first.

**Suggested order, when taken up:**

1. Benchmark recall by class. First the simulated framework
   (`analysis/is6110_simbench/`, widened to all classes); then the Marin 2024
   hybrid assemblies with real reads (139 are not in CX333).
2. Identify the 1,442 bp (all lineage 1) and 1,458 bp (lineage 7) elements in
   ISfinder.
3. Prototype an SA-pair scanner for templated insertions, on scale200 BAMs
   (the elements in step 2 are positive controls).
4. Genotype recovered insertions in every sample, with REF and ALT junctions,
   on the SV interval arm's pattern. Without it they fail the callability
   floor.
5. Decoy contigs for novel and mosaic sequence only; T3's read-stealing
   caution applies.
6. Pooled-carrier targeted assembly for recurrent novel sites.
7. External tools (GRIDSS, panISa and others) on the benchmark, for comparison
   only.

**Side items:**

- A classifier for IS6110-mediated deletions (18 of 29 RD mismatches sit near
  IS6110).
- DR-region IS6110 keying: noted only, user's decision; the flank-consensus fix
  was reverted.

Nothing reaches the chain before it is prototyped, benchmarked and approved
under the one-rerun rule.

## 0f. External long-read truth sets: Marin 2024 and Behruznia 2024 (2026-10-04)

Purpose: an independent check of assembly reliability, and later the
insertion-gap benchmark (section 0e). Nothing downloaded yet except Marin's
supplementary zip; download location and isolate set are pending the user's
decision.

**Where the data is:**

- **Marin et al., Bioinformatics 2025 (PMC12119186):** 151 Mtb clinical
  isolates.
  - Assemblies (hybrid and short-read) on Zenodo 10846276.
  - Run accessions in Supplementary File S2 (zip via Europe PMC; PMC itself
    serves a JavaScript page instead of the file).
  - Analysis repo: github.com/farhat-lab/mtb-pg-benchmarking-2024paper.

| set | data | isolates | long + short reads public |
|---|---|---:|---:|
| A | PacBio RS II + Illumina | 27 | 26 |
| B | PacBio Sequel II + Illumina (PacBio runs found at ENA by BioSample) | 21 | 20 |
| C | PacBio HiFi, "TRUST PB Set 1" | 8 | 0: PRJNA1250160 still unreleased at ENA and SRA (checked 2026-10-04) |
| D | ONT, Hall 2022 Madagascar | 78 | 62 |
| E | ONT + Illumina, Peker 2021 | 17 | 17 |
| **total** | | **151** | **125**: about 101 GB long-read + 50 GB short-read FASTQ; median 150x and 96x |

  R21770's only Illumina run is 2.7x. Most of these isolates are not in CX333.

- **Behruznia et al. 2024, the 200 non-Marin genomes:**
  - 24 have long and short reads under one BioSample; 4 duplicate Marin, so
    **20 are new** (about 57 GB). They include 7 La1, lineage 8, M. orygis and
    M. microti; 13 are CX333 genomes, which suits a leave-one-out test.
  - Mb3601's PacBio run has no FASTQ at ENA.
  - The 11 ITM genomes Behruznia sequenced show long-read runs only.

**Assembly reliability, from Marin's own statistics**
(`analysis/marin_assembly_qc/`: a per-isolate table, `README.md`, flag rules).

Median Pilon changes per assembly, by data type:

| data | median | max |
|---|---:|---:|
| HiFi | 0 | 1 |
| PacBio RS II (ChinerOms) | 20 | |
| Sequel II (TB Portals) | 54 | |
| RS II (Farhat Peru) | 190 | |
| ONT Hall | 474 | |
| ONT Peker | **3,294** | **11,817** |

Flags:

| flag | isolates |
|---|---:|
| homopolymer-indel signature | 104 |
| 1,000 or more corrections | 35 |
| no public short reads | 18 |
| length outlier (crude) | 12 |
| no flags (natural controls) | 30 |

The long-read F2 measures read error, not mixture, so it is not used. Pilon
fixes only what short reads see; heavily corrected assemblies are the most
likely to keep repeat, collapse and misjoin errors.

**Suggested check order:** Peker (17), flagged Hall isolates with short reads,
the length outliers, then the 30 clean PacBio isolates as controls.

The source metadata are in the session scratch. `build_table.py` documents how
to rebuild them.

**External assembly QC, run 2026-10-04** (`analysis/external_assemblies/`:
`QC_PLAN.md`, `QC_RESULTS.md`, `qc_summary.tsv`).

- **Assemblies:** 155. Marin's 151 hybrid assemblies (Zenodo) and Behruznia's 4
  long-read genomes outside CX333, run through the CX333 selection screens
  unchanged.
- **Read-free long-read checks**, against the 332 CX333 genomes as the baseline:
  indel and private-homopolymer-indel excess (the SNP-excess rule), IS6110 full
  and partial copies, the DR array, PGAP pseudogenes, each hybrid against
  Marin's own short-read assembly, and PPE38 locus typing.
- **Cost:** about 6 billing-hours.

**CX333 screens:**

- **153 pass, 2 review** (SNP outliers: mada_1-38, TB3237), **0 exclude.**
- None is scrambled, reference-guided or BCG; all are single, dnaA-started
  chromosomes.
- 42 are the same isolate as a CX333 genome.
- Lineage came from tb-profiler's barcode on each assembly's own calls (not
  tb-profiler itself, about 50 billing-hours saved). It agrees with Marin's for
  every isolate. Mixed lineage is not assessed.

**Foreign-DNA screen: method fault, fixed for this run only.**

- **The two docs disagree.** RUNBOOK.md's command uses the whole 1.48 Gb
  panel as the background (1,167 false "foreign" inserts). QC_PIPELINE.md
  section 1.4 requires a small one.
- **Even with a 69-genome background, 952 inserts were called foreign.** The
  script's minimap2 `-x asm10` frequency filter drops IS6110 minimizers.
- **A re-check with the filter off** (`bin/foreign_recheck.py`): all 952 are
  IS6110. **No foreign DNA.**
- **Open:** fix RUNBOOK.md, and check whether the same seed filter affected
  CX333's own foreign calls (QC_PIPELINE.md reports 0.9% of IS6110-sized
  inserts misclassified there).

**Long-read errors that CX333's screens do not catch.** All of these pass
CX333:

| assembly | evidence |
|---|---|
| CP010333 (M. microti), Behruznia, 2015 PacBio RS II (PRJNA270004) | 417 private homopolymer indels |
| CP010329, same set | 241; barcode lineage 4.9, consistent with Behruznia's "H37Rv" label, though GenBank says strain F1 |
| CP010337, same set | 117 |
| Peker ONT (4549-04, 696-05, 8129-04, 8651-04, 702-06, QC-7 and others), and Hall's R25048 and 18_0621851 | indel excess and/or disagreement with their own short-read assemblies (QC-7: 343 SNPs; R25048: 131) |

- The other PacBio sets (ChinerOms, TB Portals, Lee 2020) and the TRUST HiFi
  assemblies show no error signal.
- **Not errors:** truncated IS6110 in lineage 2.2.1 (CX333 has the same), and
  PPE38 class E (3.9% here against 2.7% in CX333).
- **For a truth set:** prefer TRUST HiFi (reads unreleased), then the clean
  PacBio sets, then unflagged Hall. Avoid, or reassemble from reads, the three
  PRJNA270004 genomes and the flagged ONT isolates.

**CX333's own foreign-DNA calls re-checked, and two artifacts found
(2026-10-04)** (`analysis/external_assemblies/qc/cx333_foreign_recheck/`:
`README.md`, `foreign_insertions.cx333.recheck.tsv`; script
`bin/cx333_foreign_recheck.py`; the working tree was read, not written).

CX333's screen called 54 inserts FOREIGN in 24 genomes. All are accounted for:

| explanation | inserts |
|---|---:|
| IS6110 copies (the asm10 seed-filter false positive) | 27 |
| M. canettii divergence (CIPT 140010059: 14; ET1291: 9) | 23 |
| attB vectors in the two genomes already excluded as engineered | 2 |
| **kilobase single-base runs: assembly artifacts** | **2** |

No engineered DNA was missed. **The two artifacts are in CX333 panel genomes.**
Both genomes are otherwise clean (SNP, indel and private homopolymer-indel
counts normal for their sublineage; no scrambling; no Ns).

| genome | strain | lineage | artifact | assembly metadata (NCBI) |
|---|---|---|---|---|
| GCF_045348265 | ITM-2003-01539 (one of Behruznia's 11 sequenced genomes; SRR34323114, 156x) | 4.6.1.1 | 1,243 bp of A at H37Rv 545,413 | Flye 2.9.2, ONT MinION: a genuine ONT assembly. The poly-A is most likely a low-quality or adapter-derived stretch Flye kept and mismatch polishing could not remove |
| GCF_050259585 | SY-1 (M. caprae) | La2 | 1,121 bp of G at H37Rv 1,594,703 | **long-read status doubtful:** "PacBio", but assembled with SOAPdenovo v2.04, a short-read assembler. Long poly-G is the two-colour Illumina (NextSeq/NovaSeq) no-signal artifact. Likely a short-read assembly with PacBio gap-filling, or wrong metadata. CX333's provenance screen took the label and passed it (5 or more insertions of 50 bp or more) |

- **Effect now: small.** Each artifact is a private insertion on one genome's
  path in the graph; no read maps to a kilobase single-base run.
- **One isolate is affected:** gwas1000 SAMEA112806781 uses GCF_045348265 as
  its matched reference and gets about 1.2 kb of NOCALL at a position H37Rv
  lacks. No scale200 isolate uses either genome.
- **Decision (2026-10-04): do not remove them now.** Removal means a graph
  rebuild (3 days, 48 cores) and a full rerun of every cohort, against the
  one-rerun rule, for a negligible effect.

**At the next panel rebuild:**

1. **Mask both runs with N,** or replace the assemblies if corrected versions
   exist.
2. **Add a run-length rule** to the assembly QC: mask or reject any single-base
   run over about 100 bp. No current CX333 screen tests for it; it would catch
   both.
3. **Do not trust assembler or technology metadata alone.** SY-1 shows a
   "long-read" label can sit on a short-read assembly. Add the read-free
   long-read checks used on the external assemblies (indel and
   private-homopolymer-indel excess) to the panel QC.
4. **Fix the foreign-DNA screen:**
   - RUNBOOK.md: use a small one-per-sublineage background (QC_PIPELINE.md
     section 1.4);
   - re-check FOREIGN calls with the minimap2 frequency filter off and against
     IS6110 (`bin/foreign_recheck.py`), so IS6110 copies stop reading as
     foreign;
   - exempt the outgroup.
5. **Optional, only alongside a rerun:** exclude the two genomes from
   matched-reference selection. That changes SAMEA112806781's reference, so
   it follows the one-rerun rule.

**Projection precompute: measured, and declined (2026-10-04).** Job 50407681
filled GCF_000023625 with all 244,005 store positions:

| | value |
|---|---|
| wall time | 1 h 44 min on 4 cores |
| CPU time | 5 h 26 min (87% efficient) |
| peak memory | 0.66 GB |
| cost | **about 6.9 billing-hours per reference** |
| store growth | 19 MB |

gwas1000's 1.05 per reference was for its 58k SV probes only; odgi cost
scales with positions. Filling the other 173 references would be about 1,200
billing-hours and 3.3 GB. **Decision (user): do not precompute.** The store
fills lazily the first time a cohort needs a reference and is reused after.
GCF_000023625 stays filled.

## 0g. Forward strategy: a lineage 1-4 panel, and lessons for the next build (2026-10-04)

A discussion note, **no change now**: `analysis/strategy/L1_4_GRAPH_STRATEGY.md`.
Related: `analysis/external_assemblies/PANEL_CANDIDATES.md`.

**Why consider it.** The genotyping need is lineages 1 to 4. In CX333, the 46
of 332 genomes outside L1-4 (lineages 5 to 9, animal lineages, M. canettii)
supply a large share of the variation:

| | share from non-L1-4 genomes only |
|---|---:|
| SNP records | 43% |
| indel and SV records | 30% |
| SV events of 250 bp or more | 37% |
| allele sequence | 29% |

M. canettii is the largest single contributor.

**Assessment.** A lineage 1-4 graph plus a separately held outgroup is the
better long-run design for L1-4 genotyping.

- **Gains:**
  - a simpler graph: fewer and shallower nested snarls, the complexity that
    defeated vg giraffe;
  - cheaper per-reference projection;
  - a cleaner accessory catalogue and RD set;
  - panel slots freed to rebalance within L1-4;
  - no loss for L1-4 matched references (no gwas1000 isolate was matched
    across lineages).
- **Costs:**
  - non-L1-4 isolates have no close reference; screen them out by barcode
    lineage before P2;
  - **the outgroup must move out of the graph:** keep one or two canettii
    genomes aligned separately to H37Rv, and give `add_outgroup.py` that source
    instead of the graph VCF;
  - MTBC-root and animal-lineage questions are no longer answerable from the
    graph.
- **Timing:** do it at the next panel rebuild, with the rebuild items in
  section 0f. Any panel change is a rebuild and a full rerun.

**Panel-building inputs gathered so far:**

- **90 clean candidates** from the QC'd Marin and Behruznia assemblies
  (`panel_candidates.tsv`).
  - 27 are about 300 or more differences from any CX333 genome.
  - The useful gaps they fill: M. orygis (none in CX333), lineage 1.1.2,
    lineage 2.2.2, lineage 4.8, lineage 4.1.3 and 4.6.
  - They add almost nothing to lineages 5 to 9.
- **Behruznia's 11 gap-filling genomes are already all in CX333.**
- **GCF_045348265** (one of them, ITM-2003-01539) carries the poly-A artifact
  (section 0f).

**Lessons to carry forward** (the note gives detail):

1. Assembly QC beyond CX333's screens: the long-read indel and
   private-homopolymer-indel checks, a single-base run-length rule, and
   sequencing technology verified from the sequence, not the metadata.
2. Foreign-DNA screen: a small one-per-sublineage background, the IS6110
   re-check, the outgroup exempt; fix RUNBOOK.md.
3. Sample for balance within L1-4:
   - cap lineage 2.2.1;
   - fill L1 and L3 sublineages and L2.1 and L2.2.2;
   - collapse clonal clusters (`snp_nonredundant.py`);
   - add GenBank-only complete genomes, not just RefSeq.
4. Prefer assemblies with public reads (HiFi where possible), so they can be
   rebuilt with your own assembly pipeline.
5. Record panel provenance at build time.
6. The decomposed VCF splits events across records and nests most SNPs: pool
   records and read all levels.
7. Projection is per reference and bounded by the panel. Keep the lazy store;
   test pairwise-alignment liftover against odgi.
8. Keep matched references. Core accuracy is limited by PE/PPE and by GATK's
   blindness at 500 bp and above, not by the reference.
9. Every new variant class needs absence genotyping (REF and ALT observed).
10. Everything that produces a result goes in the chain script, with the audit
    listing every expected product. The scan and burden omissions were all
    hand-run steps.
11. Positive controls in every association report.
12. Run documented settings, and reconcile the run documentation
    (RUNBOOK.md against QC_PIPELINE.md).

## 0h. GenBank scan for the sparse lineages (2026-10-04)

Analysis only: `analysis/genbank_scan/README.md`; one row per genome in
`genbank_typed.tsv`. Cost under 1 billing-hour (job 50558650).

**What was scanned:** complete MTBC genomes in GenBank that are not in CX333,
**173 in all** (68 released in 2026, after CX333). Each was typed by barcode
lineage, checked for indel and homopolymer-indel excess, and given its mash
distance to the nearest CX333 genome.

**Sparse lineages: barely helped.**

| lineage | CX333 | GenBank, not in CX333 |
|---|---:|---:|
| 5 | 8 | 3 |
| 6 | 5 | 5 (one error-rich) |
| 7 | 2 | 1 |
| 8, 9 | 1 each | 0 |
| La3 (M. orygis) | 0 | 2, clean |

- **Seven of the new lineage 5, 6 and 7 genomes are one 2026 set**
  (CT2018-00095 to -00101). Each is 31 to 94 differences from a CX333
  genome.
- **The one real gap filled is M. orygis.**
- **Conclusion:** public genomes cannot fix lineages 5 to 9; that would need
  sequencing.

**Lineages 1 to 4:**

- **82 of 173 genomes are flagged for indel or homopolymer-indel excess.**
  - **PRJNA994284 (Oman, 58 genomes):** treat the whole project as
    error-rich.
  - **PRJNA270004:** the 2015 PacBio project.
  - **PRJNA1254888:** Masan XDR.
- **43 clean genomes add new diversity** (300 or more differences from any
  CX333 genome), mostly lineage 2.2, 4.1, 4.3, 1, 2.1 and 3.

## 0i. Proposal: a new lineage 1-4 panel and an outgroup outside the graph (2026-10-04)

**A proposal, nothing built:** `analysis/strategy/L1_4_PANEL_PROPOSAL.md`.

1. **Timing: rebuild before scaling to 10,000+ samples.** A panel change
   forces a full rerun, which is cheapest now.
2. **The gap, measured on gwas1000's 768 lineage 1-4 isolates:**
   - median 251 SNPs to the matched reference; 328 isolates over 300;
   - worst: lineage 1.1.3 (613), 1.1.2 (528), 1.2.2 (416) and 3 (348);
   - lineage 2.2.1 has 87 panel genomes and a median of 129.
3. **Panel:** about 300 lineage 1-4 genomes, chosen from a QC'd pool of about
   450 (CX333 286, Marin 88, GenBank 72, plus new genomes):
   1. full QC, including the long-read checks;
   2. a single distance matrix;
   3. collapse clones within about 50 SNPs;
   4. a sublineage floor;
   5. greedy selection that minimises the number of isolates over 150 SNPs
      from their nearest panel genome;
   6. long-read sequencing for the gaps the pool cannot fill (1.1.3, 1.2.2,
      3, 3.1.2).
4. **Outgroup (corrected after check 2, section 0j):**
   - **lineages 5 and 6 as the nearest outgroups**, then lineage 8
     (RW-TB008), then M. canettii ET1291;
   - lineages 1-4 plus 7 are one clade, and lineages 5, 6, 9 and the animal
     lineages are its sister clade (an earlier version of this section said
     the opposite);
   - the kit is outside the graph:
     - SNP states by direct whole-genome alignment to H37Rv;
     - event states by running the kit genomes through P1 to P5 as simulated
       pseudo-isolates;
     - ancestral polarity takes the nearest called outgroup first.
5. **Cost:**
   - graph build up to about 1,500 CPU-hours (CX333's build log: 1,478);
   - QC and selection about 20 billing-hours;
   - then reruns: scale200, then gwas1000 (about 350+).
6. **Three cheap checks can run now in `analysis/`:**
   - the kit against the graph VCF for canettii;
   - lineage 8 rooting on scale200;
   - read sketches of gwas1000 isolates against the new candidates.
7. **Open decisions:** strictly lineages 1-4 or including lineage 7;
   exclude or flag non-lineage 1-4 isolates (228 of gwas1000's 996); panel
   size; sequencing; whether to run the checks.

## 0j. Three pre-rebuild checks: outgroup, rooting, panel coverage (2026-10-04)

Analysis only: `analysis/panel_checks/README.md`. Cost about 1 billing-hour
(jobs 50572387 and 50574023). Every genome's states in H37Rv coordinates come
from direct alignment (`states/`; 596 genomes).

**Correction to section 0i: the outgroup topology.**

- **What the trees show** (the panel tree and both cohort+CX333 trees agree):
  - lineages 1-4 plus 7 are one clade, with lineage 1 the first split inside;
  - its sister clade is lineages 5, 6, 9 and the animal lineages;
  - then lineage 8, then canettii.
- **So the nearest outgroups for a lineage 1-4 tree are lineage 5 and 6
  genomes.** The proposal and section 0i are corrected.

**Check 1: outgroup states from the graph VCF against direct alignment.**

- **The panel alignment is sound:** `cx333.snps.fasta`, the ancestral-allele
  input, agrees at 99.97% of sites; the median genome differs at 19.
- **Fault A (production): `assoc/bin/add_outgroup.py` lets the last duplicate
  record win.**
  - **The mechanism:** the decomposed VCF has one record per allele path,
    and the outgroup's ALT is in only one of them. At H37Rv 1845, ET1291's
    genotype is 1, 0, 0 across three `G>C` records, and the script writes
    REF.
  - **The scale:** the ET1291 column is wrong at 527 sites in scale200
    (1.3%) and 461 in gwas1000 (0.5%).
  - **Where it lands:** `data/trees/<cohort>.og.fasta`, used for tree
    rooting, and `write_event_matrix.py`'s outgroup fallback, used only where
    `AA` does not resolve. The cohort roots still look right.
- **The same fault, analysis only:** `analysis/combined_tree/build_alignment.py`
  has the same line, about 200 wrong cells per panel genome in the combined
  trees.
- **The kit works:** direct alignment reproduces the graph's genotypes (lineage
  8: 24 disagreements in 63,650 sites).

**Check 2: ancestral alleles.**

- **Sanity check:** recomputing the production rule reproduces
  `assets/ancestral.tsv` at all 72,615 resolved sites.
- **Fault B (production): the `AA` node also contains the second canettii.**
  `bin/ancestral_alleles.py` takes "the root's child that is not ET1291", and
  in the panel tree that node includes CIPT 140010059.
  - **Against the true MTBC ancestor:** 132 of 35,530 sites variable within
    lineages 1-4 differ. At 108 of them, H37Rv's allele is called derived
    where it is ancestral.
  - **Against the lineage 1-4,7 ancestor:** 169 differ (0.48%).
  - **Where it lands:** the build asset, the `AA` tag in every cohort's merged
    VCF, and the direction of events in `write_event_matrix.py`.
- **The REFEVAL backbone already excluded both canettii.** Production was
  never changed.

**Check 3: isolate-to-panel distance.** 992 lineage 1-4,7 isolates; the method
matches P1's choice for 88% of isolates (r = 0.94 with `snps_vs_matched`;
absolute values run lower, so read them as relative).

- **All 161 clean public candidates added to CX333:**
  - median 155 → 148; isolates over 300: 97 → 85;
  - only lineage 1.1.2 gains much;
  - lineage 1.1.3 does not move (12 of 16 still over 300).
- **Greedy demand-based selection, cross-validated:**
  - about 150 genomes match CX333's 288 on held-out isolates;
  - about 300 do slightly better;
  - selecting on 802 isolates works far better than on 190.
- **Implications:**
  - choose the panel after a large cohort's P1;
  - 150 to 300 genomes are enough;
  - the remaining gaps need sequencing.

**Status of faults A and B: open, not fixed.** Both affect the current
scale200 and gwas1000 outputs at about 0.4 to 1% of sites, through event
polarity. Options:

- (a) record them and fix them at the rebuild;
- (b) measure their effect on events and association results first
  (analysis only);
- (c) fix them now and do one planned rerun of the affected passes (P0
  ancestral, P5 finish, the outgroup step, the chain), scale200 first, with
  costs before any submission.

Your decision.

## 0k. Outgroup in or out of a lineage 1-4 graph (2026-10-05)

Discussion and a measurement, nothing built. The measurement uses the
direct-alignment states from section 0j (`analysis/panel_checks/states/`;
SNPs outside the repeat mask).

### Recovering the root when the outgroup is not in the graph

- **The outgroup's SNP-tree row comes from direct alignment.** The cohort
  alignment is in H37Rv coordinates, so the row is read from the outgroup
  assembly aligned to H37Rv ("the outgroup kit"):
  - align once per panel version (`genome_states.py`, seconds per genome);
  - store the base at every H37Rv position, N where the genome is
    unaligned, aligned more than once, deleted, or next to an insertion;
  - at each cohort column, write REF, ALT or N.
- **Only the input of `add_outgroup.py` changes;** its output and the chain
  stay the same.
- **Validated in check 1:** lineage 8 agrees with the graph at all but 24 of
  63,650 sites. The kit also avoids fault A and the "no record = REF" rule.
- **Limits:**
  - accessory (node_*) columns are N for the outgroup unless the insert
    sequences are mapped to it;
  - canettii gives more N.
- **Event states** (IS6110, SV, accessory) for polarity: simulate reads from
  the outgroup and run them through P1 to P5 as a pseudo-isolate.
- **Cross-check:** lineage 1 must be the first split inside the lineage 1-4,7
  clade.

### What an outgroup costs inside the graph

| added to the 289 lineage 1-4,7 panel genomes | new SNP sites (on top of their 31,427) | runs of 50 bp or more not uniquely aligned to H37Rv |
|---|---:|---:|
| **M. canettii ET1291** | **+10,770 (+34%)** | 134 kb |
| M. canettii CIPT 140010059 | +20,167 (+64%) | 164 kb |
| lineage 8 (RW-TB008) | +918 (+3%) | 90 kb (about the baseline for any MTBC genome) |
| lineage 6 (GCF_000253355) | +1,059 (+3%) | — |

**Assessment:**

- **One canettii in the graph would undo about a third of the lineage 1-4
  simplification at the SNP level, and more structurally.**
  - Its recombined, divergent regions form the nested snarls that defeated
    vg giraffe.
  - It would serve only rooting and polarity, since no lineage 1-4 isolate is
    ever matched to it.
- **A lineage 5 or 6 genome is the cheap in-graph outgroup:**
  - it is the nearest one (the sister clade);
  - it adds about 3%;
  - it is genotyped by the same code;
  - and it gives stray lineage 5/6 isolates a closer reference.

**Recommendation:**

1. Build the graph from lineages 1-4 and 7, plus one or two clean lineage 5/6
   genomes.
2. Keep lineage 8 and canettii outside the graph, in the kit, for questions
   about the root of the whole complex.
3. **Whatever is in the graph, take the ancestral node as the MRCA of the
   ingroup lineages,** not "the root's other child" (fault B, section 0j).

**Not measured:** graph complexity itself (snarl nesting). That would need test
builds: lineages 1-4 alone, plus lineage 5/6, plus canettii. They cost real
compute and should get a cost estimate before submission.

## 0l. Toward a 10K cohort on a new graph; test graphs running (2026-10-05)

**Plan:** `analysis/strategy/ROAD_TO_10K.md`. The test graphs are in
`analysis/graph_tests/`.

**State after the cluster maintenance:**

- netscratch survived;
- the working tree is intact (882 GB);
- the colleague's CRAM collection is still present, with **54,461 CRAMs**.

### Test graphs (submitted 2026-10-05, about 200 billing-hours approved)

**Design:**

- **Base:** 50 lineage 1-4,7 CX333 genomes, the first 50 of a greedy
  demand-based selection over the 992 lineage 1-4,7 isolates of gwas1000 and
  scale200 (`base50.txt`).
  - Composition: 11 lineage 1, 11 lineage 2, 6 lineage 3, 21 lineage 4,
    1 lineage 7.
  - Median distance to the nearest genome is 174, against 155 with all 289
    CX333 lineage 1-4,7 genomes.
- **Settings:** CX333's pggb settings (`-s 10000 -l 30000 -p 95 -k 23 -K 15`),
  pggb container v0.7.4.
- **Measurement:** `vg deconstruct -a` against H37Rv, plus `odgi stats`.
- **Inputs and outputs:** sequences come from the CX333 panel FASTA. Output
  goes only to `analysis/graph_tests/<arm>/metrics.tsv`.

| arm | genomes | job |
|---|---|---|
| A | lineages 1-4,7 | 50713714 |
| B | A + GCF_022870225 (lineage 5.1) + GCF_022870205 (lineage 6.3.1) | 50713715 |
| C | A + M. canettii ET1291 | 50713717 |
| D | B with `-x auto` (sparse mapping) | 50713719 |

**The first submission (50712874 to 50712881) failed in seconds.** The
script's relative list path was read after changing directory. Fixed and
resubmitted.

**Why sparse mapping is tested:** CX333's build was 97% wfmash all-against-all
alignment (1,438 of 1,478 CPU-hours, about 94 CPU-seconds per genome pair).
That cost grows with the square of the panel size.

### Route to 10K (proposed order)

1. **Test graphs:** composition and build cost.
2. **Settle the code before any 10K spend:**
   - faults A and B (section 0j);
   - the panel QC lessons (0f, 0g);
   - scheduler changes A, B and E (0c). E is needed at 10K anyway: one pass
     would hit the 10,100-job limit;
   - keep P1's per-isolate H37Rv VCFs (gwas1000's were not kept);
   - **decide on a cheaper P1 selector (D) first,** or the 10K P1 would be
     redone.
3. **Choose the 10K isolates** from the 54,461 CRAMs (criteria pending from
   the user).
4. **Run P1 on the 10K set first.** Its H37Rv alignment and calls do not
   depend on the panel. They give:
   - the demand data for panel selection;
   - the lineage 1-4 screen;
   - the direct-calling arm.

   Cost: about 2,800 billing-hours right-sized, about 6,000 as currently
   requested.
5. **Build the new panel by demand:**
   - 150 to 300 genomes;
   - one or two lineage 5/6 genomes in the graph if arm B is cheap;
   - lineage 8 and canettii in the outgroup kit;
   - gap sequencing as soon as P1 shows the gaps;
   - then the graph, the assets and P0.
6. **Validate and scale:** scale200, then gwas1000, then the 10K in batches of
   about 1,000. After P1: about 5,000 to 10,000 billing-hours. The cohort-level
   steps (0c, item I) need design before about 2,000 samples.

**Open decisions for the user:**

1. the 10K selection criteria, and lineage 1-4 only or not;
2. the cheaper P1 selector, adopt or not;
3. approval of the step 2 code work (repository only, no fairshare).

## 0m. Full pipeline audit: 97 findings (2026-10-05)

Seven read-only audits, one per area. The consolidated list is in
`analysis/audit/CONSOLIDATED.md`, with the per-area reports beside it.
**Totals: 10 HIGH, 39 MEDIUM, 48 LOW.** The main findings were re-checked in
code and, where marked, on production data.

**The current scale200 and gwas1000 association results are provisional,**
above all the IS6110, SV and branch-null results.

### A. Findings that change current results

- **P3IS-1 (HIGH, verified on data):** IS6110 carriers whose matched
  reference already holds a copy are written GT=0 on H37Rv-frame records.
  - gwas1000: 5,103 wrong cells against 2,424 correct ALT; scale200: 1,098
    against 469.
  - Polarity is inverted where H37Rv itself has a copy.
- **ASSOC-1 (HIGH, verified):** a read-based M. canettii isolate (`canettii`)
  is a tip of both association trees, with phenotype 0.
  - Its branch is 25.8% (scale200) and 11.4% (gwas1000) of tree length.
  - Removing it cuts q_branch passes from 138 to 31 and from 68 to 38.
  - **Deferred by the user (2026-10-05):** keep the isolate in for now and
    decide after the test-graph data, once the costs are clear. No fix
    touches it.
- **PGB-6 (verified):** a third artifact genome. GCF_039770655, the only
  lineage 9 genome and the reference for all lineage 9 isolates, has a
  377 bp poly-T. It produces false calls in all 16 lineage 9 isolates. A
  single-base-run rule at 100 bp flags exactly the three artifact genomes.
- **SV genotyping:**
  - **P4P5-4:** delly records genotyped 0/0 with FILTER PASS become ALT
    (6,604 rows in gwas1000);
  - **P4P5-3:** one clip cluster "confirms" an inherited deletion despite
    full depth (7,321 cells);
  - **P4P5-1:** any inexact odgi projection is written ABSENT (about
    109,000 gwas1000 SNP cells; the total ABSENT on scale200 SNP records is
    134,308, the real-deletion share to be settled);
  - **P4P5-2:** REF stated at reference-carried core indels;
  - **P4P5-5, P4P5-6, P4P5-8:** widened intervals, duplicate insertion
    records, overlapping inherited records.
- **Accessory and IS6110:**
  - **P3IS-2:** 674 of 802 presence loci cannot be seen by the read route
    but are written absent;
  - **P3IS-3:** node-frame IS6110 keys on repeated nodes merge unrelated
    insertions;
  - **P3IS-4:** DR-rescue offset error.
- **Event direction:**
  - faults A and B (section 0j);
  - **TP-1:** further REF-for-N and REF-for-ALT cases in `add_outgroup.py`;
  - **TP-2:** `X,*` SNP records (19% of scale200 SNPs) dropped from the tree
    alignment and the outgroup;
  - **TP-3:** the H37Rv tip is set to REF at node-frame records;
  - **GRAPHVCF-3:** `panel_polarity.py` (hand-run) has the last-duplicate
    fault;
  - **GRAPHVCF-5:** decompose leaves REF padded, so 376 SNPs lose their
    ancestral allele.
- **Other:**
  - **GRAPHVCF-4:** P1 selection ignores the ALT base (10 isolates
    affected);
  - **ASSOC-2 to ASSOC-5:** null and burden details;
  - **PGB-8:** 79 panel genomes never had the SNP-outlier screen.

### B. Hazards for the rebuild and rerun

- **Skip guards test file existence only:** P1, P2, P0 ancestral, the
  association chain and the accessory tables. **The rerun goes into new
  output folders.**
- **Paths hard-coded to CX333:** panel SNPs, node tables and panel allele
  frequencies, the build ID stamp, the IS6110 tie-break files.
- **The repeat mask changed** under the build's recorded checksum.
- **`pggb_build.sh` cannot start a clean build** (provenance is written
  before the emptiness check), and its defaults differ from CX333's.
  vcfwave settings also differ.
- **Parts of CX333's panel construction were never coded:**
  - the 484 → 333 step;
  - the production provenance rule;
  - the foreign screen as run.
- **Duplicate graph-VCF records come from vcfwave.** The collapse step
  removes them almost losslessly, so new-graph readers should use the
  collapsed VCF.

### C. Code location

- **The association chain, the tree builders, `panel_polarity.py` and all
  graph and panel build code run from the working tree only,** and were never
  reviewed.
- **The working tree also holds stale copies** of `cohort_assoc_tail.sh` and
  `audit_chain.py`. Run from there, they bring back the omissions fixed in
  0b.

### Order of work (the user approved steps 1-3 on 2026-10-05)

1. Bring all production code into the repository, with no change in
   behaviour. **Done 2026-10-05:**
   - 34 scripts copied from the working tree;
   - `runroot/assoc/bin` linked to the repository;
   - the repository copies, run from runroot, reproduce scale200's
     `snps.fasta`, `sites.tsv` and `og.fasta` byte for byte;
   - the working tree's stale copies stay as they are (read-only), and
     production runs only from runroot.
2. Fix section A, with a test for each finding. ASSOC-1 waits for the user's
   decision.
3. Fix section B.
4. A second review of the fixes.
5. Build the new panel and graph.
6. Rerun into new folders: scale200 first, checked against these findings,
   then gwas1000. Cost estimate before each.

**Test graph D (sparse mapping) finished in 82 minutes.** A, B and C are still
running.

## 0n. Test graphs, step 1 done, fixes in progress (2026-10-05)

### Test graphs: results

`analysis/graph_tests/RESULTS.md`. Base: 50 lineage 1-4,7 genomes, built with
CX333's pggb settings.

| change from the base (A) | nodes | VCF records | records at nesting level 4+ | build time |
|---|---:|---:|---:|---:|
| B: + one lineage 5, one lineage 6 | +11% | +12% | +15% | +5% |
| **C: + M. canettii ET1291** | **+60%** | **+63%** | **+101%** | **+26%** |
| D: B with sparse mapping (`-x auto`) | +10% | +12% | +15% | −20% against B |

- **Canettii in the graph is about five times the cost of two lineage 5/6
  genomes.** It adds 175 kb of new sequence and doubles deep nesting. **Keep
  it outside the graph** (the outgroup kit, section 0k).
- **A lineage 5/6 in-graph outgroup is affordable.**
- **Sparse mapping gives the same graph within 1%, 20% faster,** with the
  saving growing as n² with panel size. **It is adopted only after a
  record-by-record comparison of B and D's variants.**
- **Build model:** a 150-300 genome panel costs roughly 300-1,200 CPU-hours
  with full mapping. CX333's build was 1,478.
- **Cost of the four arms:** about 300 billing-hours, above the 200 estimated,
  because builds ran 82-123 minutes.

### Step 1 (code into the repository): done

Commit 19c039a. 34 working-tree-only scripts were copied byte for byte, and
`runroot/assoc/bin` now links to the repository. The repository copies, run
from runroot, reproduce scale200's `snps.fasta`, `sites.tsv` and `og.fasta`
exactly. The test suite passes (33 tests).

### Step 2 (section A fixes): in progress

Five groups, each in its own git worktree, with a regression test per finding
and before/after counts on real data copies. Brief:
`analysis/audit/FIX_BRIEF.md`.

| group | findings |
|---|---|
| is6110 | P3IS-1, -3 to -7 |
| polarity | faults A and B, GRAPHVCF-2/3, TP-1 to TP-3, TP-6 to TP-11 |
| p4p5 | P4P5-1 to -6 and -8, P3IS-2 |
| assoc | ASSOC-2 to -6 and LOWs (not ASSOC-1) |
| selection | GRAPHVCF-4 to -6, PGB-10/11, the P2 summary path, line endings |

The branches are merged and reviewed only after every group reports.

### The canettii isolate (ASSOC-1)

**Deferred by the user.** The test graphs concern a canettii *genome* in the
graph, not this cohort *sample*.

**Requested (2026-10-05): a full before/after list of which associations
change** without the isolate. The user stressed that:

- downstream analysis is almost all tree based;
- ancestral reconstructions are crucial;
- right answers outrank speed.

The comparison holds everything else equal: the same production code (main at
19c039a, whose association scripts are byte-identical to what production
ran), the same seeds and settings. The only change is the isolate, pruned from
the same tree. Both arms are rerun, so Monte Carlo noise is not mistaken for
an effect. It covers the scan, all three burdens and the ancestral
reconstruction at internal nodes.

## 0o. Audit fixes complete on branch `audit-fixes`; canettii comparison (2026-10-06)

### Status

- **Branch `audit-fixes` (ffd7bb9, pushed; NOT yet merged into main)** holds
  steps 1-3 of 0m:
  - step 1: code into the repository;
  - step 2: section A fixes, from 5 groups;
  - step 3: section B fixes, from 4 groups, plus 2 clean-up passes.
- **Size:** 95 files changed. 254 tests pass, none skipped: 33 original and
  221 new regression tests, each failing on the old code. `tests/run_tests.py`
  also runs `tests/test_audit_*.py`.
- **Open before main:**
  - the user's confirmation of the 44 decisions in
    `analysis/audit/DECISIONS.md`;
  - the D41 decision;
  - step 4, the independent second review (not yet approved).
- **Briefs and reports:** `analysis/audit/FIX_BRIEF.md`,
  `analysis/audit/CONSOLIDATED.md` and `analysis/audit/DECISIONS.md`.

### What the fixes change (measured on copies of production data)

**IS6110 (P3IS-1):**

- carriers at H37Rv-empty sites: scale200 469 ALT + 1,098 wrongly REF →
  **1,567 ALT, 0 REF**; gwas1000 2,424 + 5,103 → **7,528 / 0**;
- insertions on graph nodes their reference visits more than once are no
  longer merged, so they leave the cohort records (D8): scale200 456 of 772,
  gwas1000 1,945 of 3,585.

**P5 states:**

- ABSENT is written only where the reference lacks the position (about half
  of the old ABSENT SNP cells are confirmed deletions);
- R's own allele decides REF; cells where R carries a third allele are NOCALL
  (D3, about 600 per sample);
- clip-only SV confirmations are dropped (scale200 1,918, gwas1000 14,087
  cells).

**Polarity:**

- fault A: 386 / 323 outgroup cells fixed (about 1,400 per cohort in all,
  mostly REF → N);
- fault B: the ancestral node is the MRCA of the non-canettii genomes, and
  132 lineage 1-4 sites change;
- `X,*` SNPs enter the alignment: scale200 46,379 → 54,903 columns.

**Association:**

- the region-null leave-one-out works;
- SV evidence tiers are rebuilt and passed;
- deletions are credited to every gene they overlap: SV units 179 → 263,
  small deletions too;
- the burdens get the 80% floor;
- the survivor rule needs q_branch (D15). This changes no survivor today.
- **The DR positive controls hold:** rpoB, rpoC, embB, gyrA, pncA, ethA pass;
  katG narrowly misses the lineage null.

**Selection:** P1 matches on alleles. Final references change for 4 / 200
scale200 and 28 / 997 gwas1000 isolates, mostly closer, by about 6 SNPs.

**Graph VCF:**

- the collapsed file has 0 duplicate keys and 0 padded records;
- every production reader reads it;
- the duplicate bug had changed 26 of 802 accessory loci (for example
  ACC_2165937: 216 → 2,078 bp).

**Bugs found during the fixes (not in the audit):**

- **P4 / P4b odgi strand off-by-one at inverted path steps.** The placed REF
  matched H37Rv at 22 of 59 SNPs before the fix and 59 of 59 after.
  - Affected: scale200 234 records in 52 samples; gwas1000 about 1,660 in 276
    samples.
  - P4b: 1,285 scale200 SV keys move by 1-2 bp.
- **Node offsets counted along the walk direction.** One base got two keys,
  and P5 read 536 scale200 and 15,972 gwas1000 cells at the wrong base. Keys
  now use the node's forward offset.

### Build and rerun safety (now enforced)

- **Guards are keyed on build ID:** outputs from another build are refused,
  never silently reused or overwritten.
- **Assets are copied into the build and checksum-verified before every run.**
- **Every input comes from the build or an explicit argument;** no silent
  CX333 or pilot defaults remain.
- **New P0 steps:** nodes, catalogue, is6110_intervals, panel_polarity.
- **The outgroup is set in one place** (`MTB_OUTGROUP`, recorded in
  `build_info.tsv`).
- **`pggb_build.sh` builds clean** with CX333's settings by default, and
  records every setting plus checksums.
- **The panel construction is in code and reproduces CX333's 514 → 333**
  with all 181 reasons. It adds the single-base run rule (it flags exactly
  GCF_039770655, GCF_045348265 and GCF_050259585), the long-read error checks
  and the production provenance rule. See `docs/PANEL_BUILD.md` and
  `docs/PANEL_TREE.md`.

**The current production build 7713a8d71d8e fails the new verification** (its
assets are links outside the build, and the repeat mask changed). This is moot
under the plan for a new graph and build.

**Rerun prerequisites (new build):**

1. Regenerate the collapsed graph VCF with the fixed `vcf_collapse.sh`.
2. Run P0 with all the new steps; build the panel tree (`docs/PANEL_TREE.md`)
   before the ancestral step.
3. Run p1iv in full (`node_occ`).
4. Twoframe and svgt rerun, because interval IDs change.
5. Write everything into new output folders.

### Canettii isolate: full before/after (`analysis/canettii_effect/README.md`)

**Setup:** the same code (production), seed and inputs, with only the
`canettii` tip pruned. The rerun with the isolate reproduced production byte
for byte.

- **Survivors unchanged:** scale200 3 → the same 3; gwas1000 28 → the same
  28. All are known DR mutations, and all pass q_branch either way.
- **The isolate inflated only the branch null:** q_branch passes scale200
  137 → 29, gwas1000 62 → 38; SV units 39 → 8 and 39 → 3. Every change is
  significant → not.
- **Rows testable only through it:** 95 / 188. The embA promoter small burden
  becomes a survivor without it (scale200).
- **The ancestral reconstruction is unchanged within the MTBC,** except at
  the MTBC root: 242 / 326 variants there, mostly resolved → unknown (the
  isolate was breaking ties).
- **The decision stays with the user** (ASSOC-1).

### Decisions pending (44, `analysis/audit/DECISIONS.md`)

- **Those that change results most:** D3, D8, D9, D15, D16, D37.
- **D41 (recommended):** write node-frame alleles on the node's forward
  strand, so the same event from opposite walks gets one key, not two.

## 0p. Step 4: independent second review of `audit-fixes` (started 2026-10-06)

Approved by the user 2026-10-06. The brief is
`analysis/audit/review2/REVIEW_BRIEF.md`; reports go to
`analysis/audit/review2/<area>.md`.

**Rules for the reviewers:** read-only for code and data, no cluster jobs, no
writes to the working tree. They verify the fix groups' claims rather than
trust them.

| reviewer | priorities |
|---|---|
| genotyping | the strand corrections in P4 (R→H37Rv) and P5 (H37Rv→R) for every path-orientation combination; forward node offsets in every producer and reader; rerun P5 on 2-3 scale200 samples and audit changed cells by hand |
| is6110_accessory | P3IS-1 carrier semantics from writer to event direction on the tree; D9's effect on polarity and burdens; whether the D8 repeat-node exclusions are truly ambiguous; the D6 unmeasurable loci; node_locus with forward offsets |
| trees_ancestral_assoc | Fitch reconstruction, root, outgroup polarity and gains re-derived by hand; the fixed ancestral alleles against check 2's MTBC-ancestor column; the fixed association tail end to end on scale200 in scratch, each change traced to a fix; the DR controls; the statistics of the nulls and BH |
| build_graph_panel | a full fresh P0 on test graph arm B (no wrongful refusal, no silent CX333 input); stale outputs refused; `docs/PANEL_BUILD.md` and `docs/PANEL_TREE.md` followed literally; collapse union and trimming |
| integration | every HIGH/MEDIUM finding and every D1-D44 against the code; every producer→consumer hand-off; test skips by environment; mutation tests; cost at 10,000 samples |

**After the review:**

1. Verify and consolidate the findings.
2. Fix them, with tests.
3. The user confirms the decisions.
4. Merge `audit-fixes` into main.
5. Then the new panel and graph.

The working tree was synced to durable storage on 2026-10-06.

## 0q. Step 4 review results: 12 items before any rerun (2026-10-06)

Full list: `analysis/audit/review2/CONSOLIDATED.md`. The five reports are
beside it.

### Verified independently on real data

- **The P4 and P5 strand corrections:** right in every orientation (31-mer
  homology).
- **Forward node offsets:** consistent in every producer and reader.
- **Fitch reconstruction, gains and losses:** match a separate implementation
  on all 85,575 comparable scale200 records.
- **The fixed ancestral alleles:** equal check 2's MTBC ancestor at 72,986 /
  72,986 sites.
- **The fixed association tail on scale200:**
  - every change traces to a fix;
  - survivors unchanged (rpoB 761155, embB 4247429/4247730);
  - DR controls hold;
  - leave-one-out, BH and the nulls are correct.
- **IS6110 P3IS-1:** correct end to end.
- **A fresh build** works through P0; stale outputs are refused.
- **D1-D44:** implemented as described.

### To fix before any rerun

| # | Severity | Item |
|---|---|---|
| 1 | HIGH | Indel keys in the panel polarity table are not left-aligned like the cohort VCF's, so ancestral alleles are missed for 2,951 of 8,711 scale200 indels (103 inverted) |
| 2 | HIGH (process) | The chain builds a cohort-only tree (against the CX333-included rule); the guards refuse the combined tree; `build_alignment.py` still has fault A and TP-1 (TP-5 unaddressed) |
| 3 | MEDIUM | False ABSENT (`*`) calls from repeat-copy anchors in `p5_states.deleted_in_ref`; read downstream as unknown |
| 4 | MEDIUM | D6 "unmeasurable" judged at any identity; 3 real variable accessory loci wrongly dropped |
| 5 | MEDIUM | Variants conditioned on unmeasurable accessory loci silently leave the scan (about 249 of 308) |
| 6 | MEDIUM | Reciprocal overlap not applied when matching caller deletions to intervals; 65 scale200 deletions dropped |
| 7 | MEDIUM | Chain provenance omits code, tree, outgroup and node-locus |
| 8 | MEDIUM | The manifest does not require `ancestral`; `build_info.tsv` is not hashed |
| 9 | MEDIUM | The I/O contract refuses a fresh build |
| 10 | MEDIUM | `p3acc` has no dependency on P1 (pre-existing) |
| 11 | decision | D8 repeat-node exclusion removes about 19% of IS6110 sites, about half of them keepable |
| 12 | decision | The permutation floor of 1/20,000: a lone true hit in the small gene burden can reach at best q = 0.177 |

There are also about 20 LOW items in the reports.

### Progress on the 12 items (on `audit-fixes`)

- **Item 1, done (88b27a9). Panel polarity keys are left-aligned like the
  cohort VCF.**
  - `panel_polarity.py` normalises keys with `mtb_norm.normalise` against the
    build's H37Rv and pools collisions.
  - Indels that find a key: scale200 1,553 → 4,452, gwas1000 1,953 → 5,690.
    159 / 161 of the newly matched have ALT ancestral; none is lost or
    changed.
  - P0 step `panel_polarity` needs `refs`, and records the H37Rv and
    `mtb_norm` checksums.
- **Item 2, done (2554f71). The chain builds the cohort + panel tree** (the
  user chose to move it into the chain).
  - **New steps:**
    - 2, `assoc/bin/combined_alignment.py`: every panel genome's row via
      `add_outgroup.outgroup_allele`, with the outgroup as one of them, plus
      panel-only SNP columns only where no cohort record covers the position;
    - 3, IQ-TREE rooted on the outgroup;
    - 3b, `assoc/bin/prune_for_cohort.py`.
  - The event writer reads the outgroup from the combined alignment.
  - **scale200 (scratch):** 533 taxa × 90,812 columns in 31 s. Panel rows
    against direct alignment: 648 discordant cells (0.002%), against 68,935
    (0.25%) in the old combined alignment.
  - **Not yet run:** IQ-TREE on the combined alignment (a job; part of the
    rerun).
- **Item 3, done (31de8ce). No ABSENT where R carries the position.**
  - `p5_states.deleted_in_ref` first checks whether p's H37Rv context (12
    bases each side, any base at p) occurs in R within 2,000 bases plus
    odgi's distance of the target, either strand. If it does, the cell is
    NOCALL.
  - **On three scale200 references,** false ABSENT positions go 7/654 →
    0/624, 59/667 → 0/564 and 7/668 → 0/655. Every removed call had p's
    context in R; no real deletion was lost.
  - **Tried and rejected:** unique anchors, which lost about 100 real
    deletions per sample.
- **Definitions given to the user (2026-10-06):**
  - **ABSENT** is a positive claim, inferred from R's sequence, that the
    position is deleted: GT 2, `*`, and unknown in the tree and event
    matrix.
  - **NOCALL** makes no claim.
  - **Polarity:** only `AA` (Fitch on the panel tree) and the cohort-tree
    reconstruction are ancestral reconstructions. The panel polarity table
    (one outgroup's allele), presence-is-derived and ALT-is-derived are
    readings or assumptions.
- **Item 4, done (eaaa255). Accessory read-route blindness is judged at 95%
  identity.**
  - `merge_catalogues.py` writes `h37rv_cov95` (BLAST; minimap2's assembly
    presets miss the 30 sequences under 200 bp), and
    `locus_presence.read_route_blind` uses it: blind if ≥0.9 at ≥95%
    identity. A catalogue without the column keeps the old rule.
  - P0 step `catalogue` needs `refs`, and uses `MTB_BLASTN` (default
    `$MTB_QC_BIN/blastn`).
  - **CX333:** blind loci 674 → 604, including measurable ACC_2867346,
    ACC_2165937 and ACC_0334653. The FASTA and every other column are
    unchanged.
  - **Not changed:** contiguous H37Rv copies with mixed old PRESENT/ABSENT
    calls stay blind. Their pool coverage comes from homologous copies
    elsewhere, because the pool is aligned to the whole catalogue.
- **Item 5, done (4b5f33a). No level-2 conditioning on a never-measured
  locus.**
  - For a node-frame variant inside an accessory locus no sample was
    measured at (every cell UNMEASURABLE), `write_event_matrix.py` now leaves
    it unconditional on P5's own states. The locus goes in a new
    `variants.tsv` column, `acc_locus_unmeasured`.
  - A measured locus with no carrier still makes its variants inapplicable.
  - The scan, burden and audit follow, because they read `variants.tsv`.
  - **Production event matrices with item 4's rule:** 167 of 254 (scale200)
    and 862 of 1,251 (gwas1000) conditional variants are now kept,
    unconditional; 87 / 389 remain conditioned.
- **Item 6, done (33aff55). One same-deletion rule.**
  - `sv_intervals.same_deletion`: near in position and length, and at least
    half of the longer deletion's bases shared.
  - It is now used for the graph clustering, for adding caller deletions to
    the catalogue, and for the merge's "the catalogue supersedes this caller
    row". The caller's first deleted base (h37rv_pos + 1) is compared with
    the catalogue's start.
  - Caller deletions no longer wrongly dropped: 98 (scale200), 312 (gwas1000).
- **Item 7, done (2ac6da3). Chain provenance is a chain of checksums.**
  - Each `.prov` of the association chain holds the build, the VCF, the code
    and every input's checksum, upstream products included (`_args_for` in
    `cohort_assoc_tail.sh`). A change upstream or a code fix makes the
    downstream products stale, and they are refused, naming the inputs that
    differ.
  - `MTB_CHAIN_PRINT_PROV=<product>` prints the expected record.
  - **After any code change,** move the chain products aside (or use a new
    cohort name).
- **Item 8, done (b730b24). The P0 manifest is strict.**
  - It requires `ancestral`.
  - `manifest.done` is cleared at the start of each manifest run, and
    whenever any other step is marked done after it.
  - `p0_check verify` checks `build_info.tsv` against its manifest row.
- **Item 9, done (fb9997d). The I/O contract accepts a fresh build.**
  - P0 step `assets` makes and indexes `accessory_panel.fasta` (identical to
    production's on CX333). P3 requires it and never writes into the build.
  - The contract declares p1g's `isclean.bam` and p1iv's cohort key table as
    outputs.
- **Item 10, done (ca121eb).** `p3acc` waits for P1's summary job (it reads
  `refmap.tsv`). The I/O contract lists the refmap as an input of p3.
- **278 tests pass.** Items 1-10 are done.
- **User's decisions (2026-10-06), both option (b):**
  - **Item 11 (D8), IS6110 repeat-node rows:** exclude only the ambiguous
    ones. Keep a row when every carrier's flank placement agrees on one
    locus (about half of the 456 / 1,945 excluded rows).
  - **Item 12, the permutation floor:** adaptive. Rerun only the rows at the
    1/20,000 floor with more permutations (for example 1,000,000), so a lone
    true hit can be significant.
- **Item 11, done (ec2cda7). IS6110 repeat-node rows (option b).**
  - In `is6110/bin/is6110_write_vcf.py`, a repeated node keys its rows when
    every site on it is `on_path`, its H37Rv span is at most 1 kb
    (`REPEAT_SPAN`), and no sample appears on it twice.
  - On scale200, 37 of the 64 repeated nodes pass. That recovers 156 of the
    456 excluded sites, all at the same H37Rv position (span 0).
  - The other 27 nodes are still excluded as ambiguous.
- **Item 12, done (4524fac). Adaptive permutations (option b).**
  - Every row starts with 20,000 permutations. A null whose count of permuted
    values at least as extreme is 10 or fewer is rerun with 1,000,000 fresh
    permutations, before BH (`--refine-permutations`, `--refine-below`;
    0 disables).
  - This applies to the scan (branch, region, lineage and conditional nulls)
    and to the burdens (branch, region and lineage). A new `refined` column
    marks the nulls that were rerun.
  - **scale200 small-variant burden:**
    - units passing all three nulls go from 6 to 8: katG and the embA promoter
      now pass;
    - rpoB, rpoC, embB, gyrA, pncA and ethA pass more strongly (q_region
      0.022 to 0.0006);
    - the fabG1 promoter still fails the lineage null;
    - the refinement adds about 4 minutes.
- **All 12 items are done on `audit-fixes`; 283 tests pass.** The branch is
  not yet merged into main.
- **D41, done (f835497; the user chose option b).** Node-frame alleles are
  written on the node's forward strand.
  - **P4:** complements the alleles where R reads the node reverse
    complemented (odgi's walk flag combined with R being stored flipped in
    the panel); re-anchors indels on the forward left
    (`mtb_norm.node_forward_alleles`).
  - **P5:** the reversion test complements to match.
  - **scale200** (P4 rerun locally, 199 samples):
    - reverse-complement twin keys 3 to 0;
    - SNP keys whose REF is the node's own base: 10,932 to 22,008 of 22,008.
  - **Edge, mostly closed (746e57d):** 244 records (1.1%) are reverse-read
    indels or MNPs whose forward start falls on the neighbouring node (64%
    of CX333 nodes are 1 bp).
    - 96 are now keyed on that node, found from the build's node table.
    - 148 stay as R reads them, counted `off_node`:
      - 138 were counted as anchored on a node on H37Rv's path. That count
        used the node-FORWARD anchor (R's base after the span), not the
        H37Rv-strand one. The user chose option (b): key them as a forward
        reader does.
        - **(b) is done for the clean group (0100c0c):** R reads against
          H37Rv's strand, so the H37Rv-strand anchor really is on the path.
          This covers 22 records in scale200: 15 are keyed `h37rv:<t>` with
          H37Rv's anchor, and 7 in core are left to the direct arm. After
          left-normalisation, 9 of the 15 match keys that other samples
          wrote for the same event; the other 6 are events no other sample
          carries. p4_place.sh now projects the base after every indel
          (6,729 extra positions in scale200).
        - **The along-strand group, done (16715bb; the user took the
          recommendation, 2026-10-07):** R reads along H37Rv's strand and
          the node in reverse, so R's own anchor is off the path. These are
          keyed with on_strand's t-len(D)-1 only where the REF is H37Rv's
          own sequence (always true for an insertion). Deletions of
          off-path sequence keep their node key.
          - scale200, both groups together: 73 records keyed in H37Rv, 30
            in core left to the direct arm, 35 kept as node keys
            (`off_node_not_h37rv_sequence`), 45 still `off_node`
            (inversion junctions, nodes visited twice, anchors themselves
            off the path).
          - Most of the 58 along-strand records have no event of the same
            size within 15 bp in any other scale200 sample, so there is
            nothing to compare against. Where there is one, it usually
            agrees (CCCG>C at 2,802,268).
        - **Existing, not caused by (b):** 407 of 1,702 composed-arm
          h37rv-frame deletions (24%) have a REF that differs from H37Rv.
          These are deletions of R's own sequence that H37Rv lacks. 3 of
          (b)'s 15 records are of this kind.
      - 5 are at inversion junctions;
      - 5 are on nodes visited twice.
    - Homopolymer indels read in opposite directions can still sit at
      different anchors: node-frame keys are not left-aligned on the
      forward strand.
  - The explanation is in `analysis/audit/D41_EXPLAINED.md`. 308 tests pass.
  - **Homopolymer: fixed (16715bb).** Node-frame indels are left-aligned
    on the node's forward sequence (`mtb_norm.node_left_align`), spelled
    from R's own sequence via R's node table; the shift stops at the
    node's first base.
    - scale200: 59 records shifted. Checked against the GFA's node
      sequences:
      - node indel keys not left-aligned: 42 to 0;
      - events with two keys: 3 to 0.
    - 31 keys could shift further onto the previous node; that is not
      followed.
    - Records whose node R visits twice are left unchanged
      (`no_node_sequence`).
- **New item (the user's request, 2026-10-06): an ancestral allele (`AA`)
  for node-frame variants.**
  - `bin/ancestral_alleles.py` writes `AA` only for H37Rv-frame panel SNPs.
    Node-frame variants get none; their polarity rests only on the cohort
    tree's Fitch.
  - **Plan:**
    - reconstruct the state at node-frame sites from the 333 panel genomes'
      paths through each node, on the panel tree, as for H37Rv sites;
    - record each node's ancestral orientation as an annotation.
  - **When:** with the new panel and graph build, since node ids are
    graph-specific. Not started.
- **Waiting on the user before the merge:**
  - the 44 decisions in `analysis/audit/DECISIONS.md` (D41 is decided: b);
  - the canettii isolate;
  - R2-TREES-6: 80 variants derived at the MTBC root. Should they be pinned
    to the `AA` allele?
- **A running chat log** for the user is at `analysis/CHAT_LOG.md`
  (committed on main), updated each turn.

### Working mode (the user's preference, 2026-10-06)

- Fewer agents, one bug at a time, so the user can respond as bugs emerge.
- Items 1-12 are to be fixed sequentially in the main session, on
  `audit-fixes`, each with a test, starting with item 1.

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
- **Not started:** the scaling work in section 0c, the per-reference
  coordinate map from the 10k scaling plan, and a review of about 15 GB of "possible junk".
- **The 50 (sample, key) pairs** that occur twice in gwas1000's key table
  predate today's changes. They are worth a look, but they do not block.

## 8. Working rules for the next session

- Code changes go only in the repository. The working tree is data, written
  only by approved pipeline runs.
- Conda environments go under the home directory. `gh` is in `~/.conda/envs/gh`.
- Any output-affecting code change makes the affected passes stale for every
  cohort. Record it here, and settle it before the next rerun, not after.
