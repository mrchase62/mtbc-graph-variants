# Step 4 second review: consolidated (2026-10-06)

Five reports in this folder: `genotyping.md`, `is6110_accessory.md`,
`trees_ancestral_assoc.md`, `build_graph_panel.md`, `integration.md`.

**What holds up (verified independently):**

- **The strand corrections:** in P4 and P5, in every orientation (31-mer
  homology on real data).
- **The forward node offsets:** in every producer and reader (57,966 P4 records
  match the graph base).
- **The Fitch reconstruction and gains/losses:** match a hand-written
  independent implementation on all 85,575 comparable scale200 records.
- **The fixed ancestral alleles:** equal check 2's MTBC ancestor at 72,986 /
  72,986 sites.
- **The fixed association tail on scale200:**
  - every change traces to a fix;
  - survivors unchanged;
  - DR controls hold;
  - leave-one-out, BH and the nulls are correct.
- **IS6110 P3IS-1:** carriers are correct end to end.
- **Build:** a fresh build works through P0, and about 25 stale-output cases
  are refused.
- **Decisions:** all D1-D44 are implemented as described; no prohibitive cost
  at 10,000 samples.

## Findings to fix, by severity

| # | Severity | Finding | Where |
|---|---|---|---|
| 1 | **HIGH** | **R2-INT-1:** indel keys in the panel polarity table are trimmed but not left-aligned, while the cohort VCF is left-aligned, so the exact lookup misses. scale200: 2,951 of 8,711 indels miss, 103 are inverted. | `panel_polarity.py`, `write_event_matrix.py:508-521` |
| 2 | **HIGH (process)** | **R2-TREES-1 / TP-5:** the chain still builds a cohort-only tree (against the standing rule to include CX333's genomes); the guards refuse the combined tree; `build_alignment.py` (the combined-tree route) still has fault A and TP-1 for all 332 panel genomes. | `cohort_assoc_tail.sh` step 3, `analysis/combined_tree/build_alignment.py` |
| 3 | MEDIUM | **R2-GENO-1:** false ABSENT (`*`) calls. Deletion anchors can sit in another repeat copy; 59 of 667 on GCF_000193185, inherited by every sample on it. Read as unknown downstream, so no false event. | `p5_states.py:260-362` |
| 4 | MEDIUM | **R2-IS-1:** D6 "unmeasurable" accessory loci are judged at any identity; 60 of 674 are divergent enough to measure. 3 are clearly real and variable in gwas1000 (ACC_2867346, ACC_2165937, ACC_0334653). | `locus_presence.py`, the merge |
| 5 | MEDIUM | **R2-INT-2 / R2-IS-4:** variants (and some IS6110 node records) conditioned on an unmeasurable accessory locus get no carriers and silently leave the scan (about 249 of 308 conditional variants). | `write_event_matrix.py:826-834`, `assoc_scan.py:303-320`, `audit_chain.py` |
| 6 | MEDIUM | **R2-GENO-2:** P4P5-5 only partly fixed. Reciprocal overlap is not applied when matching caller deletions to intervals; 65 scale200 caller deletions are dropped as covered. | `sv_intervals.py:185`, `merge_cohort_vcf.catalogued` |
| 7 | MEDIUM | **R2-INT-4 / R2-TREES-2:** chain provenance records hold only the build, VCF and polarity hashes, not code, tree, outgroup or node-locus, so a code fix or new tree is silently skipped. | `cohort_assoc_tail.sh` |
| 8 | MEDIUM | **R2-BUILD-2 / R2-INT-3:** the manifest does not require `ancestral`; `build_info.tsv` is not hashed; a stale `manifest.done` can survive. | `p0_prepare.sh:697`, `p0_check.py` |
| 9 | MEDIUM | **R2-BUILD-1:** the I/O contract refuses a fresh build (`accessory_panel.fasta` not made by P0; two outputs undeclared). | `refbias/io_contract.tsv`, P0 |
| 10 | MEDIUM | **R2-INT-5 (pre-existing):** `p3acc` is submitted without a dependency on P1 but reads P1's refmap. | `refbias_run.sh:452` |
| 11 | MEDIUM (decision) | **R2-IS-2 / D8:** the repeat-node exclusion removes about 19% of IS6110 sites; about half are consistent keys that could be kept. | IS6110 writer |
| 12 | MEDIUM (decision) | **R2-TREES-3 / ASSOC-11:** the permutation floor of 1/20,000 means a lone true hit in the 3,548-gene small burden can reach at best q = 0.177. | `assoc_scan.py`, burden |

**LOW items** (about 20) are listed in the reports. Among them:

- D3's NOCALLs are often decidable (lost information, not errors);
- a blank offset is read as 0 for IS6110 on H37Rv-path nodes;
- shard-edge differences;
- `retier` uses the working-tree mask;
- the H37Rv tip at node-frame alignment columns;
- doc commands that do not run as written;
- 18 tests skip silently without the env;
- the D15 test does not catch the rule's removal;
- D41's description is off (node-frame alleles already follow the projected orientation).

**Still unaddressed from the audit:** TP-5 (item 2), PGB-11, ASSOC-11/12, and
PGB-6 (the artifact genomes are detected, not removed; they go at the new
panel).
