# scale200: old run vs fixed rerun (report only, 2026-10-08)

- **Old:** `refbias/scale200` and `assoc/scale200`. Code of 2026-10-01/02,
  build `7713a8d71d8e`.
- **New:** `refbias/scale200_fix` and `assoc/scale200_fix`. All audit fixes,
  build `7713a8d71d8e-fix1`.
- **The same** 200 isolates, CX333 graph and phenotype file (57 RRDR carriers of
  202 tips).
- **Script:** `compare_old_new.py`, tables in `compare/`. It reads outputs only.
- **Not in this run:** ABS-1, the inherited-ABSENT check (fix list, HANDOFF).

## Summary

- **Every DR positive control is still found.** The fixes add katG and rpsL:
  - variant-level survivors go from 3 to 6 (adding katG S315T, rpsL K43R and
    embB 4,247,431);
  - small-variant burden survivors go from 6 to 8 (adding katG and rpsL).
- **q_branch passes fall from 137 to 92.**
  - Most of the drop is SV intervals (100 to 65), which never survive the
    region or lineage null.
  - In core sequence (36 to 25), 8 of the lost passes were IS6110 junction
    fragments double-counted as small insertions, and the rest moved because
    their calls changed or sat at the edge.
- **Calls are 98% unchanged cell for cell** on the 77,312 records in both runs.
  The main change is about 245,000 more NOCALL cells, 85% of them in PE/PPE
  and masked sequence.
- **Accessory presence** records fall from 802 to 204, by design. Loci whose
  sequence H37Rv carries are now UNMEASURABLE rather than ABSENT.
- **Ancestral alleles:** the number of records where the reference carries the
  derived allele (AA_INVERTED) falls from 6,451 to 1,079, the ancestral fix
  (Fault B). AA coverage is about the same (35,182 and 35,292).

## 1. Matched reference (P1)

- 30 of 200 samples get a different reference.
- In all 30 the new one is closer: median SNP distance 203 -> 148
  (`compare/refmap_changes.tsv`).

## 2. Records in the merged VCF

| class | old | new | shared IDs | old only | new only |
|---|---:|---:|---:|---:|---:|
| small | 79,429 | 79,007 | 72,759 | 6,670 | 6,248 |
| sv | 6,693 | 6,223 | 3,838 | 2,855 | 2,385 |
| is6110 | 778 | 752 | 663 | 115 | 89 |
| accessory_presence | 802 | 204 | 52 | 750 | 152 |

- **IDs change without the event changing** in three classes:
  - indels move with left-alignment (D22);
  - SV interval IDs follow the rebuilt interval catalogue;
  - accessory locus IDs carry their position.
- **Old-only and new-only counts are therefore upper bounds on real change.**
  Matching by coordinates (same kind, within 20 bp, length within 10%) pairs
  162 of the 762 old-only scan rows with a new ID.

## 3. Calls

| class | state | old | new | change |
|---|---|---:|---:|---:|
| small | REF | 12,711,323 | 12,541,702 | -169,621 |
| small | ALT | 362,844 | 372,688 | +9,844 |
| small | ABSENT | 2,415,969 | 2,262,319 | -153,650 |
| small | NOCALL | 395,664 | 624,691 | +229,027 |
| sv | ALT | 29,696 | 26,235 | -3,461 |
| sv | NOCALL | 1,039,529 | 957,112 | -82,417 |
| is6110 | ALT | 597 | 2,043 | +1,446 |
| is6110 | ABSENT | 21,056 | 20,126 | -930 |
| accessory_presence | REF | 147,697 | 29,597 | -118,100 |

**On records shared by ID** (15.5 M cells):

| old state | new state | share of cells |
|---|---|---:|
| same state | | 98.2% |
| ABSENT | NOCALL | 0.86% |
| REF | NOCALL | 0.72% |
| REF | ALT | 0.07% |
| ABSENT | REF | 0.06% |
| NOCALL | REF | 0.05% |
| ABSENT | ALT | 0.03% |

New NOCALLs on shared small records, by region:

| region | old REF | old ABSENT |
|---|---:|---:|
| PE/PPE | 78,156 | 50,100 |
| masked | 19,312 | 54,619 |
| core | 11,005 | 28,542 |

The core ABSENT -> NOCALL cells are sites that the old run called absent from
the reference alone. That is the direction ABS-1 continues.

## 4. Ancestral alleles

- With AA: 35,182 -> 35,292 small records.
- AA_INVERTED: 6,451 -> 1,079.
- On 72,759 shared small records, AA differs on 5,912: the base changes on
  5,723, it is gained on 169, and lost on 20. This is the ingroup-MRCA fix
  (Fault B).

## 5. Association

| | old | new |
|---|---:|---:|
| variants tested | 2,602 | 2,690 |
| q_branch < 0.05 | 137 | 92 |
| of which SV intervals | 100 | 65 |
| of which core | 36 | 25 |
| survivors (all three nulls) | 3 | 6 |
| smallest p_branch | 5e-5 | 1e-6 |
| small-variant burden survivors | 6 | 8 |
| SV burden survivors | 0 | 0 |
| IS6110 burden survivors | 0 | 0 |

The smallest p_branch is lower in the new run because it uses more
permutations.

**Survivors:**
- **old:** rpoB S450L (761,155), embB M306 (4,247,429), embB G406 (4,247,730);
- **new:** those, plus katG S315T (2,155,168: 28 gains, 59 carriers), rpsL
  K43R (781,687) and embB 4,247,431.
- rpoB S450L has 23 independent gains, against 16.

**Burden survivors:**
- **old:** embB, ethA, gyrA, pncA, rpoB, rpoC;
- **new:** those, plus katG and rpsL.

**Old core q_branch passes that no longer pass (14 of 36):**
- **8 are IS6110 junction fragments** recorded as small insertions (sequence
  GAACCGCCCCGG...), at 1,657,016, 1,986,638/41 and 3,797,827/8. They are
  not small-variant records in the new run, so they are no longer counted
  twice.
- **2 changed with the genotyping fixes:**
  - 2,372,436 (a 42-bp insertion) goes from 32 carriers in 20 gains to 116
    carriers in 3 gains;
  - 2,266,550 G>T goes from 118 carriers to 85.
- **1 is a GC insertion at 2,266,613** with no matching new ID.
- **3 are at the edge:** 1,481,337, 2,122,395 and 2,123,146, at q 0.055-0.09.

## Cost

| part | billing-hours |
|---|---:|
| chain | 577 |
| tree | 53 |
| tests | 1 |
| **total** | **about 631** (approved 380) |

The overruns:
- P5 states (362), because of the empty projection store of a new build id;
- the 24-core tree.
