# Does duplicate marking remove fold-back false calls? (2026-10-11, user request)

**Why:** P2 does not mark PCR duplicates. At fold-back hotspots in these
libraries only 42-97% of the supporting reads are distinct molecules, so one
chimeric molecule copied by PCR can count as several reads of support. The
explainer page (`../explainer/`) suggested duplicate marking might cut false
calls. This tests it before any production change.

**Isolates:** the 9 chimeric controls of the read-removal experiment
(`isolates.txt`). Start from their P2 BAMs of all reads (`<isolate>_all`):
QC'd, downsampled to 100x, and mapped to their own genome, so every SV call is
false.

## Steps, per isolate (`03_dedup.sbatch`, 9-task array)
1. **Mark duplicates** with samtools markdup (established; same tool family
   as P2):
   - `collate`, then `fixmate -m`, `sort` and `markdup -S`;
   - a pair is a duplicate when both mates start at the same positions and
     orientations as another pair;
   - `-S` also marks the supplementary (split) pieces of a duplicate read,
     which carry the fold-back junction.

   Output: `_dupmark` (duplicates flagged, all reads kept).
2. **Remove the flagged reads.** Output: `_dedup`.
3. **Depth control (`_randdup`).** From `_all`, remove as many random pairs
   (seed 7) as step 2 removed. This separates the effect of losing those
   particular reads from the effect of lower depth. Part C showed that
   random removal can raise false calls by 13-21%.
4. **dysgu** (the P2 command) on `_dedup` and `_randdup`. Also on `_dupmark`,
   to check whether dysgu ignores flagged reads by itself. If it does,
   production could flag duplicates without removing them.
5. **Our caller v6 plus local assembly** on `_dedup` and `_randdup`. Our
   caller does not check the duplicate flag today, so `_dupmark` is not run
   with it.
6. **Fold-back scan v2** on `_dedup`, giving the rate and junction sites after
   duplicate removal.
7. **Duplicate rates.** For each isolate:
   - the duplicate share among fold-back pairs (names from `names/`);
   - the duplicate share among all other pairs.

   If fold-back pairs are duplicated more often, PCR amplification is
   inflating the hotspots.

## Then
- Run `self_noise.py --all-own` on the new sets. HaplotypeCaller is not
  rerun, so its columns are copied from `_all` and are not part of this
  test.
- Report false calls in callable sequence by caller and allele-fraction bin
  (`../fb_effect.py`). Compare the sets `all`, `dedup`, `randdup` and
  `fb100` (part C), and `dupmark` for dysgu only.

## Limits
- The reads were downsampled to 100x before mapping. This keeps each PCR
  copy with the same probability, so some duplicates are already gone.
  Production downsamples the same way, so this is the production-relevant
  number.
- Short inserts mean two independent molecules can share both end positions
  by chance, so some true molecules are counted as duplicates. The random
  control does not fully cover this.
- Reads are only 9 chimeric libraries. Effects on clean libraries
  (recall) are not tested here.

**Cost:** about 0.5 billing-hours per isolate (sorting and marking 0.1;
dysgu x3 0.15; our caller plus assembly x2 0.2), about 5 in total.
