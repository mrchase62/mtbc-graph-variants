# Read-removal experiment: fold-back reads cause most SV false calls (2026-10-11)

Part C of `../FOLDBACK_SV_EFFECT_PLAN.md`.
- **Isolates:** 9 chimeric Phase B controls (N0004, N0054, N0145, N0155,
  N1176, N1202, N1272, N1274, RW_TB008), fold-back v2 6.5-66.5 per 1,000.
- **Reads:** QC'd as the cohort standard and mapped to the isolate's own
  genome, so every SV call is false.
- **Four read sets per isolate:**
  - all reads;
  - 50% of fold-back pairs removed;
  - 100% of fold-back pairs removed;
  - the same number of random non-fold-back pairs removed (the depth control).

  Fold-back pairs are 1.3-10.7% of pairs.
- **Check:** after removal the v2 rate is 0.0000-0.0004 (fb100), half
  (fb50) and unchanged (rand).
- **Analysis:** `../fb_effect.py`, writing `out/effect.calls.tsv` and
  `out/effect.summary.tsv`. Calls are counted in callable sequence only.
- **Cost:** QC 9 tasks, about 1 h; P2 45 tasks; post about 1 h; about 20
  billing-hours in all.

## False calls, 9 isolates summed

| caller | AF | all | fb50 | fb100 | rand | fb100 vs rand |
|---|---|---|---|---|---|---|
| ours, typed | all | 781 | 207 | 19 | 919 | -98% |
| ours, BND | all | 795 | 226 | 74 | 990 | -93% |
| ours | < 0.1 | 894 | 249 | 48 | 1038 | -95% |
| ours | 0.1-0.2 | 616 | 168 | 33 | 786 | -96% |
| ours | 0.2-0.5 | 55 | 10 | 7 | 71 | -90% |
| dysgu PASS | all | 3999 | 3107 | 455 | 4523 | -90% |
| dysgu | < 0.1 | 735 | 856 | 306 | 694 | -56% |
| dysgu | 0.1-0.2 | 1561 | 1335 | 122 | 1696 | -93% |
| dysgu | 0.2-0.5 | 1279 | 821 | 26 | 1556 | -98% |
| dysgu | >= 0.5 | 424 | 95 | 1 | 577 | -100% |
| delly PASS | all | 16343 | 10528 | 6373 | 16237 | -61% |
| delly | < 0.1 | 12560 | 8343 | 5551 | 11965 | -54% |

## What it shows
1. **Causal: fold-back reads cause nearly all of our caller's and dysgu's
   false calls.** Removing them, with insert size, library, depth and
   genome held fixed, removes 95% of ours and 90% of dysgu's. For dysgu at
   AF >= 0.2 it removes 98-100%. Removing the same number of random pairs
   does not reduce false calls.
2. **The effect is steeper than the dose.** Halving the fold-back reads cuts
   our false calls by 72% and dysgu AF >= 0.2 by 46-78%. A false call needs
   several chimeric reads at the same junction, so fewer chimeras cross the
   support threshold less often. A cut-off should therefore act sharply over
   a modest range of rates.
3. **Allele fraction does not rescue chimeric libraries.**
   - Ours: 57% of these false calls are below AF 0.1 and 39% at 0.1-0.2.
     A 10% floor keeps 43% of them; a 20% floor keeps 4%.
   - dysgu: 1,703 of its 3,999 are at AF >= 0.2, and 424 are >= 0.5
     (mostly N1176, N0054, N1274). AF >= 0.2, our current dysgu rule, is
     not enough in chimeric libraries.
   - This agrees with the user's within-host experience: Tn5 short-insert
     libraries were too noisy below 10% AF.
4. **What remains after full removal.** Per isolate this is about 2 typed
   calls from ours, about 3 dysgu calls at AF >= 0.2 and about 34 below it.
   Clean controls show 5 and 3 in total across 11 isolates. The low-AF dysgu
   residue may be chimeras the scan misses (synthetic sensitivity 0.74-0.88)
   or other library noise. It is not separated here.
5. **delly is noisy apart from fold-backs.** 6,373 calls remain after full
   removal, mostly below AF 0.1. This supports dropping delly from SV calling
   (production P2/P4b/P5 still use it).
6. **Junction attribution (part A) is uninformative in heavily chimeric
   libraries.** With junctions from all reads, 54-81% of callable sequence
   lies within 20 bp of some junction, so chance level is high (92-100%
   observed vs 54-75% chance). The removal experiment is the decisive test
   here. Attribution remains useful for low-rate isolates (selfcal54).
7. **Unexplained:** random removal raises false calls by 13-21%. The likely
   cause is fraction-based thresholds (e.g. `--bnd-min-frac 0.2`): removing
   other pairs raises the fold-back share at each junction. Not checked.

## Per isolate (all / fb50 / fb100 / rand)

| isolate | fold-back v2 | ours typed | dysgu AF>=0.2 | delly |
|---|---|---|---|---|
| N0004 | 7.7 | 22/7/3/31 | 21/6/2/21 | 613/435/347/614 |
| N0054 | 29.9 | 51/11/4/66 | 277/114/2/331 | 2504/1521/817/2479 |
| N0145 | 13.0 | 95/20/0/95 | 197/84/4/203 | 1243/912/817/1257 |
| N0155 | 6.5 | 9/3/1/9 | 30/6/0/33 | 421/275/178/421 |
| N1176 | 66.5 | 463/132/8/568 | 504/438/14/763 | 6246/3801/1949/6210 |
| N1202 | 24.0 | 23/5/0/31 | 178/50/0/221 | 1862/1240/701/1831 |
| N1272 | 20.2 | 35/8/2/35 | 72/15/0/84 | 731/327/120/722 |
| N1274 | 32.3 | 81/20/1/82 | 423/203/5/475 | 2620/1980/1429/2602 |
| RW_TB008 | 10.9 | 2/1/0/2 | 1/0/0/2 | 103/37/15/101 |

The rate alone does not predict the count. RW_TB008 (10.9, MiSeq, long
inserts) has almost no false calls, while N0145 (13.0) has 95 typed calls.
So junction recurrence and insert size probably matter as well as the rate.
selfcal54 will test this across isolates.
