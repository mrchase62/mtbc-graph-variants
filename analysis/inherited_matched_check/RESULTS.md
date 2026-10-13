# Can the matched-genome alignment confirm inherited SVs? (2026-10-13)

User: "Can't read metrics from the matched genome confirm or deny the SV in
the sample?" then "yes, run the matched-genome check". Analysis only.

- Scripts:
  - `check.py`: places each inherited SV on the matched genome R and reads
    the sample's P2 BAM there;
  - `score_truth.py`: scores against truth7's finished assemblies;
  - `compare_h37rv.py`: compares with production's H37Rv-frame genotypes.
- Runs:
  - truth7, locally;
  - scale200_fix, job 52598492 (6 min, 1 core).
- Outputs: `out/`.

## Method
1. **Events.**
   - scale200_fix: P4b's inherited rows (15,180; 7,229 distinct
     reference-event pairs).
   - truth7: every SV of 50 bp or more that each isolate's R carries in the
     current graph VCF (507). truth7's own P4b rows came from an older graph
     VCF.
2. **Placing on R.** R's local haplotype (300 bp of H37Rv each side + R's
   allele from the graph VCF) is aligned to R's own FASTA with minimap2
   (asm10). The event must sit inside the alignment with 100 bp to spare,
   which also checks the placement. 95% of events are placed.
3. **Reads in R's frame** (the P2 BAM, MAPQ >= 20):
   - **inherited DEL** (R has a junction where H37Rv has extra sequence):
     - "carries": >= 5 reads span the junction cleanly, < 15% disrupted;
     - "lacks": >= 3 reads clipped at, or with a >= 10 bp indel at, the
       junction, >= 30% disrupted.
   - **inherited INS** (R has extra sequence Y):
     - "carries": depth across Y >= 0.5 x flank depth, and both of Y's
       junctions spanned cleanly;
     - "lacks": depth <= 0.1 x flank at any mapping quality, plus >= 2
       clipped reads at Y's ends;
     - "repeat": low MAPQ>=20 depth but reads present at MAPQ 0.
4. **Tandem flag.** At least 50% of the changed sequence's 20-mers also
   occur within 300 bp either side (a VNTR unit or duplication).

**The INS rules were revised after the first truth7 score:**
- *lacks* needs clipped reads: GC-rich PE_PGRS inserts the isolates carry
  showed near-zero depth;
- *carries* needs clean junctions: tandem repeats.

So the INS numbers below are tuned on truth7, not blind. The DEL rules were
not changed.

## Truth: the seven isolates with finished assemblies (`out/truth7_score.txt`)
An isolate truly carries an event when its assembly (paftools against H37Rv)
has the same type within 200 bp and within 50% in size.

| | today (all reported, no check) | matched-genome "carries" |
|---|---|---|
| DEL, all | 192 right / 28 wrong (0.87) | 181 / 191 (0.95) |
| DEL, outside the repeat mask | 27 / 1 (0.96) | 27 / 27 (1.00) |
| INS, outside the mask, not tandem | 36 / 7 (0.84) | 24 / 26 (0.92) |
| INS, outside the mask, tandem | 12 / 20 (0.38) | 12 / 27 (0.44) |
| INS, inside the mask | 68 / 89 (0.43) | 36 / 66 (0.55) |

**"lacks":**
- DEL: right 10 of 13 times, catching 10 of the 28 truly absent;
- INS: never called on a truly absent event. Those went to unclear or
  repeat.

## scale200_fix: against production's H37Rv-frame genotype (`out/scale200fix_compare.txt`)
**Inherited deletions outside the repeat mask** (in cells, columns = the
state in the cohort VCF):

| matched-genome call | ALT | NOCALL | REF |
|---|---:|---:|---:|
| carries | 581 | 40 | 4 |
| lacks | 0 | 8 | 78 |

The two frames agree almost completely (659 of 663 where both decide), and
the matched genome decides 48 of the 56 that are NOCALL today.

**Inside the mask they conflict:**
- 903 "carries" against H37Rv-frame REF;
- 948 "carries" against NOCALL;
- 25 "lacks" against ALT.

truth7 has no H37Rv-frame states, so these conflicts are not adjudicated
here. truth7 puts matched-genome "carries" inside the mask at 154 of 164
right.

**Inherited insertions** (all reported present today):

| class | events | carries | lacks | unclear / repeat | low depth / unplaced |
|---|---:|---:|---:|---:|---:|
| outside the mask, not tandem | 1,827 | 759 | 28 | 994 | 46 |
| tandem | 1,074 | 786 | 2 | 272 | 14 |
| inside the mask | 5,812 | 2,644 | 32 | 2,350 | 786 |

## Reading
1. **Deletions: yes.**
   - The matched-genome reads confirm or reject inherited deletions
     reliably. Outside the mask they agree with the H37Rv-frame check at all
     but 4 cells, and with the assemblies at every cell.
   - They can resolve most of the deletions left NOCALL today.
   - Inside the mask the two frames disagree in about 1,900 cells, and
     neither is validated there.
2. **Insertions: only outside repeats and tandem repeats.**
   - There, "carries" is right 92% of the time (24/26), against 84% with no
     check. But more than half of these events end up unclear or repeat
     rather than decided.
   - In tandem repeats (VNTR loci such as MIRU40) short reads aligned to R
     look the same whether or not the isolate has R's extra unit. Neither
     depth nor junctions can tell.
   - Inside the mask the check is weak.
   - No truly absent insertion was ever called "lacks", so the check can
     confirm insertions but not yet reject them.
3. **Small truth set.** Seven isolates, all close to their R. Few events
   are truly absent, so "lacks" is barely tested.

## Possible next steps (not done)
- Add the matched-genome junction test to P5 SV genotyping for inherited
  deletions, as a second confirmation beside the H37Rv frame.
- Insertions: report as confirmed only non-tandem insertions outside the
  mask that the check calls "carries"; mark the rest unconfirmed.
- Tandem repeats need read pairs or repeat-length genotyping (e.g. a VNTR
  caller), not this test.
