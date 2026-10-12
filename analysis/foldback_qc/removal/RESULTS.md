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

## Duplicate marking (2026-10-11, `DEDUP_PLAN.md`, job 52282456)
Question: P2 does not mark PCR duplicates. Would marking them remove
fold-back false calls? Method: samtools markdup -S on each `_all` BAM.
- `dedup`: duplicates removed.
- `randdup`: the same number of random pairs removed.
- `dupmark`: duplicates flagged only.

dysgu ran on all three; our caller on dedup and randdup only. Outputs are in
`out/dedup/` (`self_noise.tsv`, `effect.*.tsv`, `dups/`).

### Duplicates are not concentrated in fold-back reads
| isolate | pairs duplicated | fold-back pairs duplicated | other pairs duplicated | fold-back /1,000 before → after |
|---|---|---|---|---|
| N0004 | 42.8% | 32.3% | 42.9% | 7.7 → 9.2 |
| N0145 | 24.7% | 16.1% | 25.0% | 13.0 → 14.6 |
| N0155 | 19.5% | 13.0% | 19.6% | 6.5 → 7.1 |
| N0054 | 6.8% | 6.2% | 6.8% | 29.9 → 30.2 |
| N1176 | 6.1% | 6.7% | 6.0% | 66.5 → 66.2 |
| N1202 | 6.0% | 5.0% | 6.0% | 24.0 → 24.3 |
| N1272 | 7.0% | 5.8% | 7.0% | 20.2 → 20.5 |
| N1274 | 5.2% | 4.8% | 5.3% | 32.3 → 32.5 |
| RW_TB008 | 1.2% | 0.4% | 1.2% | 10.9 → 11.0 |

Duplicate shares were measured after downsampling to 100x, as production
does.

### False calls, 9 isolates summed, callable sequence
| caller | AF | all | dedup | randdup | fb100 |
|---|---|---|---|---|---|
| ours, typed | all | 781 | 623 | 748 | 19 |
| ours, BND | all | 795 | 590 | 809 | 74 |
| dysgu PASS | all | 3,999 | 4,862 | 4,451 | 455 |
| dysgu | < 0.1 | 735 | 737 | 740 | 306 |
| dysgu | 0.1-0.2 | 1,561 | 1,939 | 1,757 | 122 |
| dysgu | 0.2-0.5 | 1,279 | 1,671 | 1,498 | 26 |
| dysgu | >= 0.5 | 424 | 515 | 456 | 1 |

dysgu on `dupmark` matched `dedup` to within 1 call per isolate, so dysgu
already skips flagged reads.

### What it shows
1. **Duplicate marking does not fix fold-back false calls.**
   - dysgu: duplicate removal raised false calls by 22% overall and by 28%
     at AF >= 0.2 (1,703 to 2,186). Random removal of the same number of
     pairs raised them by 11% and 15%. The rise is largest in the
     high-duplicate libraries: N0004 21 to 124 at AF >= 0.2, N0145 197 to
     355, N0155 30 to 63.
   - Ours: duplicate removal cut typed calls by 20% (random removal: 4%).
     Most of the cut is in N1176 (463 to 345) and N1274 (81 to 50), both
     libraries with low duplication. The remaining 623 are still 33 times
     the 19 left after full fold-back removal.
2. **Why.** Fold-back pairs are duplicated no more often than other pairs,
   and in high-duplicate libraries less often (N0004 32% vs 43%). Removing
   duplicates therefore removes ordinary reads at least as fast as chimeric
   ones, so the chimeric share at each site goes up (N0004 7.7 to 9.2 per
   1,000). It also lowers depth (N0004 loses 43% of pairs), and dysgu
   gains false calls with lower depth even under random removal.
3. **This corrects the explainer page.** The page had said PCR copies
   inflate the hotspots. The 42% distinct molecules at N0004 hotspots is
   simply that library's overall duplication. Recurrence at hotspots comes
   from sequence (inverted repeats), not PCR. The page was corrected
   (version 4).
4. **For production:**
   - No case for adding duplicate marking to P2 as a fold-back fix.
   - Duplicate marking for other reasons (HaplotypeCaller, depth) is a
     separate question that this test does not answer.
   - A high duplicate rate leaves fewer distinct molecules than nominal
     depth, which is a sample-quality signal in itself (N0004: 100x
     nominal, about 57x distinct). It is worth adding to the selfcal54 QC
     table.
