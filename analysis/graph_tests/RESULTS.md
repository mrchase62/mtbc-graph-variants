# Test graphs: what an outgroup and sparse mapping cost (2026-10-05)

**Setup:**

- **Base:** 50 lineage 1-4,7 CX333 genomes chosen by demand (`base50.txt`).
- **Settings:** CX333's pggb settings (`-s 10000 -l 30000 -p 95 -k 23 -K 15`).
- **Resources:** 32 cores and 48 GB each.
- **Measurement:** `vg deconstruct -a` against H37Rv. Per-arm detail is in
  `arm_*/metrics.tsv`.

| arm | genomes | build (min) | graph length | nodes | edges | VCF records | SNP alleles | other alleles | max nesting level |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A: lineages 1-4,7 | 50 | 98 | 4,691,189 | 64,649 | 89,214 | 20,557 | 16,883 | 6,791 | 6 |
| B: A + lineage 5.1 + lineage 6.3.1 | 52 | 102 | 4,697,999 | 71,754 | 98,603 | 23,054 | 19,097 | 7,318 | 6 |
| C: A + M. canettii ET1291 | 51 | 123 | 4,866,582 | 103,212 | 140,803 | 33,489 | 28,159 | 9,344 | 6 |
| D: B with `-x auto` | 52 | 82 | 4,707,710 | 71,337 | 98,139 | 23,014 | 19,008 | 7,406 | 5 |

## Against the baseline A

| | B: +2 lineage 5/6 | C: +1 canettii | D: B, sparse |
|---|---:|---:|---:|
| nodes | +11% | **+60%** | +10% |
| VCF records | +12% | **+63%** | +12% |
| SNP alleles | +13% | **+67%** | +13% |
| graph length | +0.1% | **+3.7% (175 kb of new sequence)** | +0.4% |
| records at nesting level 4 or more | +15% | **+101%** | +15% |
| build time | +5% | **+26%** | −16% (−20% against B) |

## Reading

1. **One canettii costs about five times what two lineage 5/6 genomes do.**
   - It adds about 60% to the graph's nodes, variant records and SNPs on a
     50-genome base, plus 175 kb of new sequence.
   - It doubles the deeply nested records, the complexity that hurts
     genotyping.

   This confirms keeping canettii out of the graph, in the outgroup kit
   (HANDOFF 0k).
2. **A lineage 5/6 in-graph outgroup is affordable:** about 12% on 50 genomes.
   On a 150-300 genome panel its share would be smaller (about 3% against all
   289 CX333 lineage 1-4,7 genomes, HANDOFF 0k).
3. **Sparse mapping (`-x auto`):**
   - **The same graph, within 1%:** nodes 0.6% fewer, records 0.2% fewer,
     SNP alleles 0.5% fewer, one fewer nesting level.
   - **20% faster at 52 genomes.** Its saving should grow with panel size,
     because full all-against-all mapping grows with the square of the genome
     count.
   - **Before adopting it, check that the variants agree:** compare the B and
     D VCFs record by record (cheap, local).
4. **Build cost model:**
   - at 50 genomes and 32 cores, about 100 minutes (about 53 CPU-hours of
     wall time × cores);
   - CX333's 333-genome build was 1,478 CPU-hours.

   A 150-300 genome lineage 1-4 panel with full mapping would cost roughly
   300-1,200 CPU-hours (pairs scale as n²), less with sparse mapping. This is
   to be refined with a mid-size test if needed.

**Cost of the four arms:** about 7 hours of wall time at 32 cores and 48 GB,
roughly 300 billing-hours at FASRC weights. That is higher than the
200-billing-hour estimate, because builds ran 82-123 minutes instead of about
75.

# odgi lookup timing (2026-10-12, job 52474931)

User: "do we have a good estimate on how much time a reduced complexity
graph saves runtime of the pipeline?" then "run the odgi timing on the test
graphs". Script: `odgi_timing.sbatch`; table: `timing/timing.tsv` (120 runs,
4 h 37 min on 2 cores, about 10 billing-hours).

**Method.** The same `odgi position` queries on every graph, one thread, one
run at a time:
- three references present in every graph (L1 GCF_009730215, L2
  GCF_000193185, L4 GCF_013010385);
- r2h: 20,000 random reference positions to H37Rv (the P4/P4b form);
- h2r: 20,000 random H37Rv positions to the reference (the P3 and SV
  projection form);
- load: 10 positions, i.e. graph loading alone;
- two repeats; CPU = user + system; median over references and repeats.

| graph | genomes | nodes | CPU-s per 20,000 lookups (r2h / h2r) | ms per lookup | vs A | load (CPU-s) | memory |
|---|---:|---:|---:|---:|---:|---:|---:|
| A: L1-4,7 | 50 | 64,649 | 85 / 92 | 4.4 | 1.00 | 0.1 | 47 MB |
| B: A + L5.1 + L6.3.1 | 52 | 71,754 | 97 / 94 | 4.8 | 1.08 | 0.1 | 52 MB |
| C: A + canettii | 51 | 103,212 | 131 / 134 | 6.6 | 1.50 | 0.1 | 67 MB |
| D: B, sparse mapping | 52 | 71,337 | 89 / 92 | 4.5 | 1.02 | 0.1 | 51 MB |
| CX333 | 333 | 275,621 | 858 / 862 | 43 | 9.7 | 0.5 | 640 MB |

Every run returned all 20,000 positions. Repeat-to-repeat spread was up to
about 25% on the first repeat (B, CX333), so differences under about 10% are
not meaningful.

## Reading
1. **Lookup cost follows graph complexity.** Two L5/6 genomes add about 8%
   (D, the same graph built with sparse mapping, about 2%); one canettii adds
   50%. CX333 is about 10 times slower per lookup than the 50-genome graph.
   Loading is negligible on every graph (under 1 s, under 1 GB).
2. **How much of the pipeline is odgi.** gwas1000 (997 samples, 2,130
   billing-hours in all) spent about 290 billing-hours in the stages that
   query the graph: SV projection 158 (odgi is its whole cost), P3 66 (with
   bwa against the panel), P4b 49, P4 14. That is about 13% of the total, of
   which odgi is the larger part.
   - scale200_fix lookups per sample: P4 about 670 (two passes), P3 about
     710, P4b about 1,700 (before our caller was wired in). About 3,000 in
     all: about 2 CPU-minutes per sample on CX333, about 13 s on graph A.
3. **What a ~300-genome L1-4,7 panel would save is not measured here.** The
   test graphs have 50 genomes; the proposed panel has about 300. Its graph
   would be close to CX333 in genome count but without the 44 divergent
   genomes, which add 30-43% of CX333's variation. If lookup cost falls in
   proportion, odgi stages get about 25-40% cheaper, which is **about 3-5% of
   the per-sample cost** (about 300-500 billing-hours at 10k, at current
   requests). This revises the earlier 5-10% guess down.
4. **Panel size matters more than composition.** At equal size, dropping
   canettii and the animal lineages saves a third of the lookup cost; going
   from 333 to 50 genomes saves 90%. The build cost estimate (300-1,200
   CPU-hours against 1,478) is unchanged.
5. **The compute saving is small either way.** The case for the reduced
   graph rests on genotyping quality (fewer nested records) and build cost,
   not on per-sample runtime.

# B vs D: does sparse mapping change the variants? (2026-10-13)

User: "yes, run the B vs D check". B and D are the same 52 genomes and
settings; D was built with pggb `-x auto` (each genome aligned to a subset
of the others). Both deconstruct VCFs were run through the production
decomposition (`decompose_arm.sbatch`, jobs 52595664/5, about 10 min each,
under 1 billing-hour), then compared genome by genome (`compare_arms.py`,
`compare_B_D.out`/`.tsv`) and the differences classified
(`explain_B_D.py`, `explain_B_D.out`).

Records after collapsing: B 24,111, D 24,176.

| class | ALT calls in both | only B | only D |
|---|---:|---:|---:|
| SNP | 75,078 | 322 | 303 |
| MNP | 2,015 | 14 | 26 |
| indel < 50 bp | 10,214 | 513 | 518 |
| SV >= 50 bp | 4,944 | 239 | 278 |

Summed over 51 genomes: 92,251 calls agree, 2,213 (2.3%) differ.

**What the differing calls are** (each genome's one-arm-only calls):
- **in the repeat mask:** 1,468 (66%), mostly indels and SVs;
- **the same variant written differently:** the other arm has a call of
  the same class and the same length change within 50 bp in that genome
  (560, 25%). For SNPs this test is loose: any SNP within 50 bp counts;
- **unexplained:** 185 (0.2% of all calls; 65 only in B, 120 only in D).
  71 are in one genome, GCF_008761675, and 61 of those form one dense
  D-only cluster at H37Rv 4,224,126-4,224,342 (a local alignment
  difference in D). The rest are spread thinly.
- SNP bases after splitting MNPs: 80,938 the same, 2 different.

**Reading.** Outside the masked repeats and representation shifts, the two
builds agree at all but about 0.2% of calls, and these are clustered in
local regions rather than systematic. There is no baseline for how much two
full builds of the same input differ, so part of this may be ordinary pggb
run-to-run variation. Sparse mapping gives the same variants for practical
purposes, at a 20% faster build on 52 genomes (more on about 300). A
rebuild can use `-x auto`; the 4.22 Mb cluster can be looked at in the
rebuilt graph's own QC.
