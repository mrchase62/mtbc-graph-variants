# Can canettii root the tree without being in the graph?

2026-10-07. Analysis only. Nothing in the pipeline reads these files.

**Question (the user's):** canettii is the outgroup of the cohort trees. A
lineage 1-4 graph would leave it out, because one canettii adds about 60% to
the graph's nodes and variant records (test graphs, `../graph_tests/`). Can
the trees and the association still use canettii as the outgroup?

**Answer: yes, on every measure tested.** Reading canettii from its own
assembly aligned to H37Rv, instead of from the graph, gives:

- the same root;
- the same lineage structure;
- 99.98% of the same ancestral alleles;
- the same association survivors and drug-resistance controls.

## The test

The two arms differ only in where the two canettii genomes' rows come from:

- **G (graph):** read from today's fixed graph VCF (review 2's
  `graph_collapsed.vcf.gz`). This is canettii inside the graph.
- **K (kit):** each canettii genome's assembly aligned directly to H37Rv
  (`../panel_checks/states/`, `genome_states.py`). Any site the alignment
  cannot call cleanly is N.
  - Columns where only canettii varies are dropped, because without canettii
    in the graph they would not exist.
  - This is canettii outside the graph.

The read-based `canettii` isolate in scale200 is a cohort sample, not a
panel genome, and is the same in both arms.

| step | what | where |
|---|---|---|
| 1 | scale200 + CX333 SNP alignment, arm G | `G.fasta`, `G.sites.tsv` (`../combined_tree/build_alignment.py`) |
| 2 | arm K from arm G | `K.fasta`, `K.sites.tsv` (`make_kit_arm.py`) |
| 3 | one tree per arm, with the production tree's settings (`-m GTR+F+ASC+G4 -B 1000 -alrt 1000 -o GCF_035581225`) | `G.contree`, `K.contree` (`run_tree.sbatch`; jobs 51184235, 51184236) |
| 4 | the two trees compared | `compare_trees.out` (`compare_trees.py`) |
| 5 | panel ancestral alleles (AA) per arm, same panel tree and code | `ancG.tsv`, `ancK.tsv` |
| 6 | scale200 event matrix and scan per arm (review 2's `rr2` inputs, MTBC node pinned as in main) | session scratchpad |

**Compute:** two IQ-TREE runs of about 31 minutes on 24 cores, 15.7 CPU-hours
in all. That is above the roughly 11 CPU-hours estimated, because the runs
used more CPU time than the production tree did. Everything else ran
locally.

## Results

### 1. Canettii's rows: graph vs direct alignment

ET1291 (GCF_035581225), the outgroup, over the 94,383 columns of the
scale200 + CX333 alignment:

| | sites |
|---|---:|
| identical | 92,468 (98.0%) |
| called by the graph, N by direct alignment | 1,266 |
| N in the graph, called by direct alignment | 332 |
| opposite allele | 317 (0.34%) |

The second canettii (GCF_000253375) differs more, as expected for the more
divergent genome: 1,205 opposite alleles (1.3%). 5,340 columns vary only through
canettii and are dropped in arm K.

Over the 54,903 columns of the cohort-only tree used by the association,
ET1291's row is identical at 54,324. 35 sites have the opposite allele, 240
are N only by direct alignment and 304 are N only in the graph.

### 2. The trees

| | G (graph) | K (direct alignment) |
|---|---|---|
| alignment | 533 taxa × 94,383 columns | 533 × 89,043 |
| MTBC clade | 530 tips, support 100; outside it the two canettii genomes and the canettii isolate | **the same** |
| each lineage (1-7, 9 and the animal lineages) | monophyletic, support 100 | **the same** |
| lineages 1-4,7 | one clade, support 100 | **the same** |
| first split inside lineages 1-4,7 | lineage 1 vs the rest | **the same** |

**Robinson-Foulds distance: 12 of 1,062 splits (1.1%).**

- All 12 are inside lineage 4, among 2 to 14 closely related genomes.
- Their UFBoot support is 16 to 79 (median about 40). None is at 95 or
  above.
- These are poorly resolved near-clonal clades, far from the root. Where the
  outgroup's rows come from cannot place them. No replicate run was made to
  measure IQ-TREE's own run-to-run variation at these clades.

### 3. The ancestral alleles (AA)

Computed on the same panel tree with each arm's canettii rows:

| | sites |
|---|---:|
| in both arms | 49,179 |
| same AA | **49,168 (99.98%)** |
| opposite AA | 7 |
| resolved in only one arm | 4 (3 only G, 1 only K) |
| only in arm G (canettii-only columns, constant across the MTBC) | 24,077 |

The 24,077 sites only in arm G are constant across the MTBC, so they carry no
information about the MTBC ancestor.

### 4. The association (scale200, review 2's inputs)

Each arm's AA was re-joined into the merged VCF, and each arm's ET1291 row was
used as the outgroup row. The same cohort tree and the current code (MTBC node
pinned) were used for both.

| | G | K |
|---|---:|---:|
| scan rows | 2,729 | 2,732 |
| testable variants (2 or more gains) | 3,679 | 3,682 |
| q_branch < 0.05 | 141 | 142 |
| **survivors** | rpoB 761155, embB 4247429, embB 4247730 | **the same 3** |
| drug-resistance controls (rpoB, rpoC, katG, fabG1, embB, gyrA, rpsL, rrs, pncA and others) | | **every status unchanged** |

**18,773 records lose their AA in arm K.** 18,669 of them are carried by the
`canettii` isolate alone. They had an AA in arm G only because the graph's
canettii genomes made the panel vary there. A variant on one tip produces no
event inside the MTBC.

**582 variants change their event counts.** 389 of them are on a single tip.
The other 193 lost their AA and fell back to outgroup or default polarity.

## What this does not cover

- **Polarity for IS6110, SVs and accessory sequence.** These need canettii's
  state on sequence H37Rv lacks. The plan is to simulate reads from the
  canettii assembly and run them through P1-P5 as a pseudo-isolate. That
  route is untested.
- **Whether leaving canettii out changes the other genomes' calls in the
  graph.** The test graphs can answer this cheaply: arms A (lineages 1-4,7)
  and C (A plus canettii) share 50 genomes, and their VCFs can be compared
  genotype by genotype (`../graph_tests/`). Not done yet.
- **Inside the new graph itself.** This test swaps canettii's rows in an
  alignment built from the CX333 graph. It does not build a graph without
  canettii.
- **gwas1000.** Only scale200 was run.

## Files

`make_kit_arm.py`, `run_tree.sbatch`, `compare_trees.py` and this README are
the analysis. `G.*`, `K.*`, `panelG.*`, `panelK.*`, `ancG.tsv`, `ancK.tsv`,
`compare_trees.out` and `logs/` are its outputs.
