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

**Why it matters for scaling.** A single reference removes every
per-reference cost:

- **SV projection:** run per reference, it was 158 of gwas1000's 348
  billing-hours (section 0c);
- **P1 reference selection;**
- **the per-reference P3/P4 coordinate projection** and the per-reference
  IS6110 projection store;
- **one bwa index** instead of one per panel genome.

The trade is per-isolate false-positive suppression and some of the
interpretive gain, for a simpler and cheaper pipeline at 10k or more.

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
