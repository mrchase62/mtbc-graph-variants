# Accessory genome in scale200_fix (2026-10-08, report only)

"Accessory" here means sequence H37Rv does not have.
- **Data:** build `7713a8d71d8e-fix1`'s accessory catalogue
  (`assets/accessory_catalogue.tsv`), the per-isolate presence tables
  (`accessory/scale200_fix/*.presence.tsv`, 200 isolates) and the merged VCF.
- **Durable copy:** `results/cohorts/scale200/2026-10-08_audit_fixes/`.

## 1. The catalogue: 806 loci, insertions of 500 bp or more found among the 332 panel genomes

| | copy_number | mosaic | novel | all |
|---|---:|---:|---:|---:|
| polymorphic (varies among panel carriers) | 129 | 74 | 47 | 250 |
| reference_gap (a gap in H37Rv relative to the panel) | 550 | 4 | 2 | 556 |
| **all** | 679 | 78 | 49 | 806 |

What the categories mean:
- **copy_number:** the inserted sequence is also elsewhere in H37Rv, for
  example an IS6110 copy or a duplication.
  - 481 of the 806 loci are IS6110-sized (1,340-1,370 bp). They are new
    IS6110 insertion sites, which the IS6110 arm types.
- **novel:** the sequence is not in H37Rv.
- **mosaic:** part of the sequence is in H37Rv and part is not.

So the accessory sequence proper is the 49 novel and 78 mosaic loci.

## 2. Presence calls in the 200 isolates

Each locus gets one state per isolate:

| state | locus x isolate cells |
|---|---:|
| PRESENT | 1,735 |
| ABSENT | 29,597 |
| UNCERTAIN | 9,468 |
| UNMEASURABLE | 120,400 |

- **UNMEASURABLE** means the locus's sequence is also in H37Rv, so reads
  cannot tell presence at this locus from the other copy. It applies to
  602 of the 806 loci, nearly all of them copy_number (audit fix P3IS-2; the
  old run called these ABSENT).
- **Measurable loci:** 204, the `accessory_presence` records in the merged VCF.
  - 86 vary in the cohort, present in some isolates and absent in others;
  - 118 are absent from every measured isolate;
  - none is present in all.
- **The 86 variable loci:**
  - by type: 31 novel, 32 mosaic, 23 copy_number;
  - most are rare: 52 in 5% of measured isolates or fewer, and 18 in 20-50%;
  - 54 have all their carriers in one lineage.

**Most common variable loci**
(carriers / measured isolates of that lineage in the cohort; H37Rv gene at the
insertion point):

| locus | type | length | present / absent | carriers by lineage | at |
|---|---|---:|---|---|---|
| ACC_0742634 | novel | 533 | 176 / 23 | L4 80/104, L2 46/46, L1 19/19, L3 17/17 | Rv0647c |
| ACC_2268721 | mosaic | 6,359 | 174 / 23 | L4 82/104, L2 43/46, L1 19/19, L3 17/17 | Rv2024c |
| ACC_1594703 | novel | 1,121 | 151 / 7 | L4 73, L2 36, L1 15, L3 14 | uvrC |
| ACC_3846777 / _3846796 | mosaic | 3,780 / 9,768 | 144 / 55 | L4 51/104, L2 46/46, L1 18, L3 17 | |
| ACC_3527837 | novel | 655 | 115 / 68 | L2 41, L4 29, L3 17, L1 17 | PPE53 |
| ACC_2209604 / _2219417 | novel | 4,403 / 4,496 | 85 / 109 | L2 44/46, L1 19/19, L3 17/17, L7 4/4, no L4 | mce3A |
| ACC_1761789 | mosaic | 3,511 | 33 / 167 | L1 19/19, L7 4/4, L6 4/4, L5 3/3, none in L2/3/4 | mmpL6 |

**Positive control: TbD1.** ACC_1761789, at mmpS6-mmpL6, is the TbD1 region.
- It is present in every L1, L5, L6 and L7 isolate and absent from every L2,
  L3 and L4 isolate.
- This is the known split between "ancient" and "modern" lineages: modern
  lineages, H37Rv among them, lost TbD1.

**ACC_2209604 disagrees with the panel** (noted in `analysis/ctpV_check.md`).
- The reads call it in all L1, L2, L3 and L7 isolates but in no L4 isolate.
- The catalogue lists only 5 panel carriers.
- The read pattern (everything but L4, which includes H37Rv) is what a
  deletion in L4's ancestor would give. So the panel carrier list is the
  suspect, but this is not resolved.

## 3. Variants inside accessory sequence (merged VCF, small variants)

| region | SNPs | indels |
|---|---:|---:|
| off_path_accessory (inside catalogue insertions) | 326 | 136 |
| off_path_near (graph sequence next to H37Rv's path, not in H37Rv) | 10,832 | 203 |

- The 462 off_path_accessory variants are mostly private: median 1 carrier,
  366 singletons, 8 with 5 or more carriers.
- They are scored only within carriers of the insertion ("level 2,
  conditional").
- off_path_near records are alternative haplotypes beside H37Rv sequence
  rather than accessory sequence proper.

## 4. Association (RRDR phenotype)

- **Not tested under the nulls:**
  - 53 accessory_presence records have 2 or more independent gains on the
    cohort tree, and 8 pass the callability floor;
  - the accessory null pool has 53 events, below the 200 needed to draw a
    null, so they are not tested under the region and lineage nulls
    ("no null");
  - 0 pass the branch null.
- **off_path_accessory:** 1 variant is testable, and does not pass.
- No accessory locus survives.

## 5. Change from the old run

- **Records:** 802 -> 204 accessory_presence records, and ALT cells 2,447 ->
  1,735. This is by design: loci whose sequence is also in H37Rv are now
  UNMEASURABLE instead of ABSENT.
- **Locus IDs changed with left-alignment (D22),** so old and new loci are
  matched by sequence and carriers, not by ID.
