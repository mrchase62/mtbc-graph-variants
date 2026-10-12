# selfcal54: fold-back cut-off, rate vs rate plus insert size (2026-10-11)

The test is pre-specified in `../CALIBRATION_PLAN.md` ("Added 2026-10-11").
- Script: `calib.py`, run from runroot.
- Outputs: `out/calib_points.tsv` and `out/calib_report.txt`.

## Data
- **Points: 59.** All are QC'd reads mapped to their own genome, so every
  call is false.
  - 39 selfcal54 isolates that passed QC, from 19 studies.
  - 11 clean Phase B controls.
  - 9 chimeric controls (`../removal`, all reads).
- **QC losses (15 of 54):**
  - 14 below 60x. Five PRJEB49093 runs sit at 59.5x, and two lose more than
    80% of reads in fastp.
  - 1 at 80% MTBC.
- **Outcome y:** false calls in callable sequence, our typed calls plus
  dysgu at AF >= 0.2.
- **Clean controls:** y has mean 0.73 and SD 0.90, so the threshold
  T = mean + 2 SD = 2.54.

## Pre-specified result: insert size adds nothing
| rule | cut-offs | points kept (of 59) | selfcal kept (of 39) |
|---|---|---|---|
| 1. rate <= c | c = 5.80 per 1,000 (v2) | 40 | 29 |
| 2. rate <= c, or rate <= c2 with median insert >= i | c = c2 = 5.80, so the insert term is unused | 40 | 29 |

- **Leave-one-study-out.** Rule 2 kept 2 more held-out isolates than rule 1,
  both from marin2024: N0145 (y 292) and N0155 (y 39). The insert term lets
  bad libraries through. This is because two of the worst chimeric
  libraries have long inserts (N0145 median 501 bp, N0155 389 bp).
- **Correlations.** Median insert does not correlate with y (Spearman
  -0.02 for all points, -0.12 for selfcal).
- **Negative binomial fit.** It did not converge (dispersion went to
  infinity), so there is no model estimate.

**Conclusion:** rule 2 is rejected. Insert size should not enter the
exclusion rule.

## The pre-specified criterion is too lenient
The criterion "kept isolates average <= T" lets individually noisy
isolates through when many clean ones dilute the mean. At c = 5.80, all 3
isolates in the 3-6 band exceed T, and 2 have y >= 10.

False calls by rate band:

| v2 rate per 1,000 | n | median y | max y | y > T | y >= 10 |
|---|---|---|---|---|---|
| < 1 | 31 | 0 | 5 | 3 | 0 |
| 1-3 | 6 | 3.5 | 14 | 4 | 1 |
| 3-6 | 3 | 20 | 21 | 3 | 2 |
| 6-13 | 8 | 11 | 43 | 8 | 5 |
| >= 13 | 11 | 107 | 967 | 11 | 6 |

- **Below 1:** no isolate reaches 10 false calls, across 10 studies and 7
  instruments.
- **1-3 (the gap the plan wanted filled):** false calls start rising here.
  - The one isolate at 10 or more is GCF_035581225 (L6, 62x, rate 2.39):
    12 dysgu calls and 2 of ours.
  - It is the lowest-depth passing isolate and has not been checked for
    genome/read identity (plan limit 3). It should be checked before it
    moves the cut-off.
- **Choosing a cut-off.** A per-isolate criterion (for example, few kept
  isolates with y >= 10) points to a cut-off between 1 and 3 per 1,000
  (v2). The earlier proposal of 2.5 (v2) lies inside that range. Changing
  the criterion is a departure from the pre-specification, so the user
  decides.

## Not pre-specified: junction recurrence
Recurrence is the number of junctions with >= 5 fold-back reads.
- **Correlation with y:** about equal to the rate.
  - All points: 0.73 for recurrence, 0.75 for rate.
  - Selfcal only: 0.57 for recurrence, 0.71 for rate.
- **Where recurrence explains outliers the rate gets wrong:**
  - PRJDB10923 is MiSeq 2x300 TruSeq with about 400 bp inserts. Its rate
    is 11-28 per 1,000, but it has 0-3 recurrent junctions and y of only
    4-8.
  - RW_TB008 is similar: rate 10.9, y 3.
  - Fold-backs in these libraries are spread out rather than stacked at
    hotspots.
- **Rule 3 (recurrence <= h, same criterion):** h = 11. It keeps 44 of 59
  points, and 44 vs 39 held out under leave-one-study-out.
- **Rule 3's extra keeps:**
  - Most are the PRJDB10923 isolates (y 4-8).
  - In one leave-one-study-out fold it also kept GCF_026727615 (y 21).
- **Depth dependence:** recurrence scales with depth, and some passing
  isolates are at 62-77x.

**Conclusion:** recurrence is a promising second metric. It is not yet
better than the rate as a single rule.

## Where the mid-range false calls come from
- **The y 11-21 isolates at rates 3-12** are mostly MiSeq Nextera (Tn5)
  libraries with very short inserts:
  - PRJNA902023, PRJNA930491, PRJNA1089215, PRJDB11640;
  - median insert 130-165 bp; nearly all pairs overlap.
  - This is the same library type as the user's within-host samples.
- **Our caller is the main source** in these isolates (11-19 calls vs dysgu
  0-7). In the chimeric Phase B controls, dysgu dominated.
- **The "share of pairs shorter than 2x read length" predictor is
  confounded.** It reaches 1.0 for every 2x300 run whatever the insert,
  so its selfcal correlation (0.71) mixes that effect with short Tn5
  inserts. It should not be used as is.

## Not done here
- Part A attribution and part B AF bins on selfcal54 (`../fb_effect.py`).
- An identity check for GCF_035581225.
- Duplicate rate. P2 marks no duplicates; see `../removal/RESULTS.md`, where
  the duplicate rate turned out to be a depth signal, not a fold-back one.
