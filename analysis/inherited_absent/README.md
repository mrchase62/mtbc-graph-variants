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
