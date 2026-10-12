# Which loci need genomes outside lineages 1-4 and 7 for their ancestry? (2026-10-12)

User: "I am still on the fence about removing lineages from the graph that
may provide key ancestry for specific loci." Then: "check which loci need
non-L1-4 genomes for ancestry". Report only.
- Script: `ancestry_need.py` (local, about 10 s per run).
- Outputs:
  - `out/records.tsv`: every record that varies inside lineages 1-4,7,
    with the best L5/6 pair;
  - `out/genes.tsv`;
  - `out/l56_sweep.tsv`, `out/l56_sweep.sorted.tsv`: all 81 L5/6 choices;
  - `out/report_best_pair.txt`, `out/report_testgraph_pair.txt`.

## Method
- **Inputs:** the build's panel tree (the one production's AA uses) and the
  graph VCF (every CX333 genome's genotype at every SNP, indel and SV
  record against H37Rv).
- **Ingroup I:** the MRCA of the 289 lineage 1-4 and 7 genomes. Nothing
  else lies under it.
- **Records:** those that vary inside I, 47,280 in all (37,031 SNPs,
  6,270 indels, 3,256 SVs of 50 bp or more, 723 MNPs). Only these can give
  events inside lineages 1-4.
- **Step 1.** Fitch parsimony over I's own genomes. Where I's set is one
  allele, the ingroup decides the ancestral state alone.
- **Step 2.** Otherwise, I's state is taken from Fitch over the whole tree
  (the downward pass) in three panels:
  - full CX333 (333 genomes);
  - the proposed panel: I + one or two L5/6 genomes + the kit (L8
    GCF_012923765, canettii ET1291). The kit genomes' graph genotypes stand
    in for kit rows;
  - the kit only, with no L5/6 genome.
- **"Differs" and "unresolved" are measured against the full panel's
  parsimony.** That is the current answer, not a truth set.

## Results

| | records |
|---|---:|
| vary inside lineages 1-4,7 | 47,280 |
| ancestral state decided by lineages 1-4,7 alone | 46,464 (98.3%) |
| need an outside genome | 816 |
| ... proposed panel (best L5/6 pair) agrees with full CX333 | 799 |
| ... differs | 8 |
| ... unresolved | 7 |
| ... tied even with full CX333 | 2 |

1. **For 98% of what varies in lineages 1-4,7, no outside genome matters.**
2. **Of the 816 that need an outside genome, the proposed panel gives the
   full panel's answer for 799.**
   - **All 15 that differ or stay unresolved are in sequence the pipeline
     already masks:**
     - Rv0094c/Rv0095c (paralog);
     - PE_PGRS10, 28, 53 and 57 (PE/PPE);
     - Rv2258c (tandem).
   - The 2 tied even with the full panel (PE_PGRS38, whiB6-Rv3863) would be
     tied under any panel.
3. **No drug-resistance gene is affected.** 1,114 records in the
   resistance genes vary inside I. For 1,107 of them I decides alone, and
   for the other 7 the proposed panel agrees.
4. **The L5/6 genomes matter, and so does which ones.**
   - With the kit alone (L8 and canettii, no L5/6), 44 records differ or
     are unresolved.
   - Across all 81 choices of one or two L5/6 genomes, the count is 15-43.
     The best is 15: any element-carrying L5.1 genome (GCF_022870265,
     GCF_045348065 or GCF_022870245) with L6.1.1 GCF_000253355
     (*M. africanum* GM041182).
   - Single genomes do worse (up to 43) than pairs.
5. **The tatC-helY element (H37Rv 2,352,065) depends on the choice.**
   - **Present (by sequence, flank check):** canettii ET1291, L8, L9,
     *M. microti* (5), GCF_039555375, 6 of 8 L5, and 215 of 289
     L1-4,7.
   - **Absent:** all 5 L6, all 21 *bovis*/*caprae* (La1, La2), and L5.1
     GCF_022870225 and L5.3 GCF_021559895.
   - Full CX333 parsimony calls the element present in the L1-4,7 ancestor.
     So does the kit alone.
   - **Test graph B's pair (GCF_022870225 L5.1 + GCF_022870205 L6.3.1)
     lacks it in both genomes and flips it to absent.** That would invert
     its polarity in every association test, as the tatC-helY polarity
     fault did before panel_polarity (31 deletions read as 11 gains and 31
     losses).
   - The element comes out present in 63 of the 81 choices. It is absent
     exactly when no chosen L5/6 genome carries it.

## What this means for the reduced graph
- **For lineage 1-4 events, the ancestry the removed genomes provide can be
  kept:**
  - choose the L5/6 pair with a check like this one;
  - and keep L8 and canettii in the parsimony as kit rows.
- **Not only nearest-first.** Production's `ancestral_alleles.py` breaks a
  tie with the nearest single outgroup that decides it. With an L5/6 genome
  lacking the element as the nearest, that rule would call it absent. The
  reconstruction for a reduced panel should use parsimony over all the
  outgroups (as here), or the L5/6 choice must be checked locus by locus.
- **The kit's SV and indel states need OUT-1** (canettii, and here also
  L8, as simulated pseudo-isolates). Here their graph genotypes stood in.
- **Not covered:**
  - questions about the MTBC root itself (TbD1, the RvDs) beyond the kit's
    states;
  - SVs absent from every panel genome;
  - whether the removed genomes improve the graph's representation of L1-4
    loci (no gwas1000 isolate was matched across lineages).
