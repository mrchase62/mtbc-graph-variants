# ABSENT calls checked against each sample's own reads (scale200_fix, 2026-10-08)

**Why:** in ctpV, SAMEA2297133 had an intact gene (72-107x, no clip or split-read
pileups) but was ABSENT at four records, inherited from its matched reference's
deletion (`analysis/ctpV_check.md`). This counts how often that happens.

**How:** `count_false_absent.py`, run by `run_scale200_fix.sbatch`. For every
ABSENT cell of `refbias/scale200_fix/p5/merged.vcf.gz` that has an H37Rv span,
it takes the sample's MAPQ>=20 depth from its P1 H37Rv BAM and compares it with
the sample's genome-wide median:
- contradicted: >= 0.5x the median and no zero base;
- supported: < 0.1x;
- partial: anything in between.

Results are in `out/scale200_fix/{summary,records}.tsv`. This is a report only.

## Result

| | ABSENT cells |
|---|---:|
| all | 2,284,798 |
| not testable (node-frame, IS6110, accessory: no H37Rv span) | 2,199,333 |
| testable | 85,465 |
| supported by the reads | 70,576 (83%) |
| **contradicted (region covered at normal depth)** | **10,247 (12%)** |
| partial | 4,642 (5%) |

- **Spread:** 199 of 200 samples have at least one contradicted cell (median 26,
  maximum 277); 81 have 50 or more.
- **Sites:** 3,197 distinct sites. By record class: small 9,967, SV 280. By
  region: core 6,245, PE/PPE 2,327, masked 1,395.
- **Concentrated:** the top 20 2-kb windows hold 55% of them.
  - phiRv1 prophage (Rv1573-Rv1588c, 1.778-1.790 Mb): about 3,000 cells in 37
    samples, mostly lineage 4 matched to GCF_000153685 and others;
  - plcA/plcB (2.630 Mb);
  - PPE57/58 (3.842 Mb);
  - Rv3766-Rv3770 (4.212 Mb);
  - wag22 (1.990 Mb);
  - an insertion record at 3,552,705 (180 samples).
- **More with a poorer match:**
  - Spearman rho 0.44 between contradicted cells and the SNP distance to the
    matched reference.
  - In the 7 samples whose reference comes from a deeper sublineage than their
    own lineage call, the median is 72, against 24 for the rest.

## Caveat: 12% is an upper bound

H37Rv depth shows the sequence is somewhere in the sample, not that it is at
this locus. Reads from a copy at another place still map uniquely to H37Rv's
single copy. phiRv1 and the IS6110 hotspots (plcA/plcB, PPE57/58) are exactly
the regions where sequence moves, so some of these cells may be right to say
"not here". Clips and split reads at the reference's deletion junctions decide
it. In ctpV they did: 0 junction reads in SAMEA2297133, against 30-65 per
breakpoint in the carriers.

Cost: one 4-core job, a few minutes.

## Junctions and mappability (2026-10-08)

`junction_check.py`, run by `run_junction_scale200_fix.sbatch` in 2 minutes.
For each ABSENT cell it finds the matched reference R's deletion over the cell
in the build's collapsed graph VCF. It then counts, in the sample's H37Rv BAM
(MAPQ >= 20):
- reads that run across each end of the deletion;
- reads clipped at each end;
- reads that join the two ends (split reads, or a CIGAR deletion).

Verdicts:
- present_here: reads run across both ends and fewer than 3 join them;
- deleted_here: 3 or more reads join the ends;
- mixed: both;
- unresolved: neither;
- R_no_genotype: R has no genotype at the site (its path does not cross it in
  place);
- no_R_deletion: nothing in R explains the call.

**Masking.** A cell counts as masked when its span, or either end of R's
deletion within 150 bp (a read length), overlaps the build's
`repeat_mask.bed`: 594 kb of PE/PPE, paralog, tandem and IS-element sequence,
which includes all PE/PPE. In masked sequence, depth, clips and split reads
are all unreliable, so these cells cannot be called either way and are kept
apart.

Outputs: `out/scale200_fix/junction.tsv` and `junction.masked.tsv`.

**All testable ABSENT cells:** 64,713 of 85,465 (76%) are masked.

**The 10,247 contradicted cells:**

| junction verdict | core | masked |
|---|---:|---:|
| present_here (ABSENT is wrong) | **2,279** | 4,760 |
| deleted_here (sequence is elsewhere; ABSENT is right at this locus) | 71 | 375 |
| mixed | 1 | 24 |
| unresolved | 109 | 561 |
| R_no_genotype | 151 | 880 |
| no_R_deletion | 175 | 861 |
| **total** | **2,786** | **7,461** |

- **Masked:** 73% of the contradicted cells. They should become "not callable"
  rather than ABSENT or REF.
- **Confident errors in core sequence:** 2,279 cells in 110 samples, from 94
  distinct R deletions.
  - median 5 per affected sample, maximum 212;
  - these include the four ctpV cells of SAMEA2297133;
  - the biggest clusters are 4.212 Mb (Rv3766-70, 5 samples), 2.786 Mb (7.5 kb
    R deletion, 4 samples), 0.660 Mb (8.4 kb, 1 sample), 1.536 Mb (5.2 kb, 3
    samples) and 1.332 Mb (2.8 kb, 3 samples).
- **Correctly ABSENT despite the depth:** 71 core cells. The locus is deleted
  and the reads come from a copy elsewhere.
- **Not resolved in core:** 435 cells (R_no_genotype, no_R_deletion,
  unresolved). They need local reassembly.

The first run failed on one sample: SAMEA7526648 is matched to H37Rv itself,
which is not a column of the graph VCF. That case now returns no R deletions;
the sample has no testable ABSENT cells.

## Local reassembly (2026-10-08)

**Scripts:**
- `make_events.py`: groups cells into regions per sample.
- `local_assembly.py`: for each region, collects reads within 1.5 kb plus
  their mates, assembles them with SPAdes `--isolate`, aligns the contigs to
  H37Rv with minimap2 `asm5`, and requires 300 bp of flank on both sides.
- `annotate_complex.py`: blastn of the contig stretches that do not align to
  H37Rv, against H37Rv, the matched reference and the accessory catalogue.
- Run by `run_local_assembly.sbatch <set>`.
- Tools from the lab `autocycler` env: SPAdes 4.3.0, minimap2 2.31.

**Calibration** (128 regions, job 51392271, 10 min, about 1.4 billing-hours):
- the 5 ctpV regions with known answers;
- every core region the junction check called deleted_here (45);
- a seeded random 49 of the core present_here regions;
- 29 core zero-depth regions with junction reads.

The first classification counted deletions of up to 20 bp as present. That
called 22 deletions of 1-22 bp present, a bug in the rule, not a
disagreement. Present now allows at most 10% of the region deleted, counted
over the region +-10 bp when it is under 50 bp. Calibration was then
reclassified from the saved alignments (`results.v2.tsv`), without
reassembling:

| junction verdict | assembly: deleted_here | present_here | unresolved |
|---|---:|---:|---:|
| contradicted / deleted_here (45) | 43 | 0 | 2 |
| contradicted / present_here (49) | 0 | 48 | 1 |
| supported / deleted_here (29) | 29 | 0 | 0 |

- **Agreement:** all 120 decided regions agree with the junction check; none
  contradict it.
- **ctpV:** SAMEA2297133 is present_here at both regions; SAMEA5542103 is
  deleted_here.
- **The two 235-bp-deletion carriers come out unresolved**, because their next
  event starts 55 bp after the deletion, which is less than the 300-bp flank.
  `annotate_complex.py` resolves them. The contig reads H37Rv to 1,078,811,
  then about 340 bp of other sequence, then H37Rv again from 1,080,104 at
  99.9% identity. That is the 1,290-bp deletion with the graph's "337-bp
  insertion" in its place.
- **That insertion is not accessory.** It matches nothing in the accessory
  catalogue, and is an 85% diverged copy of ctpV 1,079,759-1,080,109
  containing a 97.7% match to 1,079,652-739. The H37Rv-frame SNP/indel cluster
  at 1,079,879-927 comes from reads of this copy.

The 109 unresolved core regions were then submitted (job 51395857).
