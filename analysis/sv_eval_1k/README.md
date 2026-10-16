# SV / insertion evaluation of the new 1K cohort VCF (2026-10-16)

Analysis only. OLD = `runroot/refbias/gwas1000/p5/merged.vcf.gz` (graph 7713a8d71d8e,
before ABS-1/ABS-2/AA-1/SPLIT-1). NEW = `runroot/refbias/new1k_p200/p5/merged.vcf.gz`
(P200 graph 3893a3762ed1, with the fixes). 525 isolates are in both.

## What was done, in standard terms

1. **Record extraction** (`extract.py`). H37Rv-frame records of four classes:
   - **DEL**: `SVTYPE=DEL` of CLASS sv (any COMPONENT), CLASS is6110 DEL (NEW's IS6110 copy losses) and
     CLASS small deletions >= 50 bp. OLD keeps >= 50 bp HaplotypeCaller indels as CLASS small, NEW writes
     them as CLASS sv `SVSOURCE=small_allele`, so both are included on both sides.
   - **INS**: `SVTYPE=INS` of CLASS sv, plus CLASS small insertions >= 50 bp.
   - **IS**: CLASS is6110 `<INS:ME:IS6110>`. Node-frame IS6110 records (519 OLD / 470 NEW) cannot be
     matched across graphs and are left out.
   - **ACC**: CLASS accessory_presence. `<REPL>` (51 NEW records) are left out.
2. **Event matching** (`compare.py`). Candidate pairs by geometry: DEL reciprocal overlap >= 50% of the
   deleted bases; INS start within 50 bp and length within 20%; IS same site within 20 bp; ACC same
   ACCLOCUS id, else start within 50 bp and length within 20% (only 16 of the 99 NEW ids exist in OLD).
   Pairs are then made one-to-one, greedily, by the Jaccard index of ALT carriers over the shared
   isolates (geometry breaks ties), so that where a locus has several allele records, each OLD record
   pairs with the NEW record that best represents it. "secondary" = a record whose candidates were all
   taken by a better pair (split records, extra alleles). Cell transitions are counted on paired records.
3. **Read check** (`sample_cells.py`, `readcheck.py`, `context.py`). Fixed-seed random sample
   (seed 20261016) of disagreeing cells on paired records, ALT vs REF/ABSENT, 30 per class and direction
   (13 for OLD-only IS6110, which is every such cell after a cap of 2 cells per record pair). ACC has no
   ALT<->REF flips at all, so its 20 cells are ALT<->NOCALL. Each cell is checked in the isolate's own P1
   H37Rv BAM (MAPQ >= 20, no duplicates), at the coordinates of the build that called ALT:
   - DEL: depth inside vs flanks, reads with the deletion in the CIGAR, reads soft-clipped at either
     breakpoint, reads running unbroken across a breakpoint (REF evidence); a depth-only ALT is not
     accepted where MAPQ>=0 depth is high (repeat).
   - INS: reads with the insertion in the CIGAR or clipped at the site (when the record carries the
     sequence, the clipped bases must occur in it), unbroken spanning reads; tandem duplications are
     also checked at the far junction and by the depth of the duplicated copy.
   - IS6110: junction reads (clipped at the site, clipped bases in IS6110, or SA tag into an H37Rv copy)
     vs unbroken spanning reads.
   - ACC: reads that are unmapped, clipped, MAPQ < 20 or NM >= 4 re-aligned (bwa mem) to H37Rv + the
     locus representative; breadth over the locus positions unique against H37Rv.
   A REF verdict is downgraded to unclear when >= 10 reads are clipped in one 100 bp bin within 1 kb of
   the site (a shifted or complex allele the site check cannot see). Extra checks: 60 ALT cells of
   unpaired records (`sample_unpaired.py`) and 30 OLD ALT -> NEW NOCALL deletion cells (`sample_nocall.py`).
4. **Tree consistency** (`homoplasy.py`, `homoplasy_summary.py`). Fitch parsimony change count per record,
   leaves ALT=1, REF=0, all else unknown; scored only when a record has >= 2 ALT and >= 2 REF leaves
   (NOCALL-heavy records otherwise score 0). Trees: the existing cohort IQ-TREE tree
   `analysis/combined_tree/itol/gwas1000/tree.nwk` pruned to the 525 (OLD vs NEW), and a new FastTree
   (`build_tree.py`, `tree.sbatch`: GTR+gamma on the NEW VCF's biallelic SNPs outside the repeat mask,
   all 1,000) for NEW on all isolates and as a check of the 525 numbers. SPLIT-1 is simulated with the
   event writer's own `split_sv_groups()` (assoc/bin/write_event_matrix.py; same type, starts within
   50 bp, lengths within 5%, not a copy-number step) on DEL+INS; merged genotype = ALT if any member ALT,
   else REF if any REF. **Not simulated**: the SPLIT-1 extension that maps truncated small-allele
   insertions to their graph allele.

## 1. OLD vs NEW on the 525 shared isolates

Record counts differ a lot because OLD wrote one record per caller call (10,510 `called` DEL and
6,606 `called` INS records, 99.9% NOCALL) and 802 accessory loci; NEW writes 824 / 1,942 caller
records and 99 accessory loci.

Records with >= 1 ALT among the 525:

| class | OLD paired | OLD secondary | OLD only | NEW paired | NEW secondary | NEW only |
|---|---|---|---|---|---|---|
| DEL | 1,529 | 333 | 4,002 | 1,482 | 320 | 317 |
| INS | 1,364 | 74 | 2,758 | 1,352 | 375 | 1,262 |
| IS6110 | 814 | 4 | 8 | 978 | 0 | 8 |
| ACC | 27 | 4 | 58 | 27 | 4 | 13 |

The OLD-only side is mostly singleton caller calls: 3,573 of the 4,002 OLD-only DEL records are
`called` records (4,365 ALT cells; 3,607 of all unpaired OLD DEL records are singletons), and 2,452 of
the OLD-only INS are `called` (3,855 cells). The NEW-only INS are mostly `small_allele` records
(923 records, 33,095 ALT cells), largely extra alleles at loci OLD also has.

Cell transitions on paired records (rows OLD, columns NEW):

**DEL** (1,956 pairs)

| OLD \ NEW | ALT | REF | ABSENT | NOCALL |
|---|---|---|---|---|
| ALT | 26,657 | 3,593 | 0 | 12,106 |
| REF | 3,099 | 520,419 | 0 | 146,850 |
| ABSENT | 90 | 791 | 200 | 3,859 |
| NOCALL | 6,871 | 79,660 | 72 | 222,629 |

**INS** (1,507 pairs)

| OLD \ NEW | ALT | REF | ABSENT | NOCALL |
|---|---|---|---|---|
| ALT | 14,305 | 1,967 | 0 | 8,783 |
| REF | 6,665 | 118,332 | 167 | 50,916 |
| NOCALL | 6,783 | 56,255 | 53 | 526,683 |

**IS6110** (1,110 pairs)

| OLD \ NEW | ALT | REF | ABSENT | NOCALL |
|---|---|---|---|---|
| ALT | 1,280 | 10 | 4 | 25 |
| REF | 2,956 | 502,277 | 23,360 | 8,961 |
| ABSENT | 10 | 15,047 | 18,797 | 1,612 |
| NOCALL | 97 | 1,973 | 1,370 | 4,971 |

**ACC** (57 pairs)

| OLD \ NEW | ALT | REF | NOCALL |
|---|---|---|---|
| ALT | 3,364 | 0 | 170 |
| REF | 0 | 12,753 | 552 |
| NOCALL | 238 | 81 | 12,767 |

NOCALL share over all records of a class, shared isolates (`out/state_counts_shared.tsv`):
DEL 65% OLD -> 40% NEW; INS 92% -> 75%; IS6110 1.4% -> 2.3%; ACC 6.8% -> **43%**.

What ABS-1/ABS-2 did (NEW EV and component of the changed cells, `out/new_change_ev.tsv`):
- **OLD ABSENT on deletions** (4,940 cells) is now ABSENT in only 200: 3,859 NOCALL (1,910
  `absent_masked`, 1,081 on `called` records with no EV, 225 `masked`, 186 `r_other_allele`, 176
  `inherited_unrecorded_neither_frame`), 791 REF (689 with no EV, 51 `absent_checked_present_ref`) and
  90 ALT. NEW ABSENT on DEL is 360 cells in total on the 525 (`junction_joined` / `absent_junction_joined`).
- **NEW ALT gained from OLD REF** (3,099): 950 small_allele, 741 interval_graph with no EV, and 1,321
  inherited cells confirmed by ABS-2 (`inherited_unrecorded_matched_only` 462, `inherited_reference_two_frames`
  300, `inherited_unrecorded_h37rv_only` 263, `inherited_reference_matched_only` 167, ...).
  NOCALL -> ALT (6,871) is mostly interval_graph with no EV (3,906).
- **OLD ALT lost to REF** (3,593): small_allele 1,484, interval_graph 950, interval_caller 468, IS6110 copy
  losses written as small alleles 434. Lost to NOCALL (12,106): small_allele `no_coverage` 2,773,
  interval_graph `ambiguous_depth` 1,837, `called` 1,561, small_allele `masked` / `r_path_unresolved` /
  `r_other_allele` ~2,370, `uncovered__and_so_are_the_flanks` 667.
- On INS the changes are almost all small_allele records: REF -> ALT 6,535 (no EV), NOCALL -> ALT
  6,192, ALT -> NOCALL 5,038, and REF -> NOCALL 36,132 (26,897 `r_other_allele`: the isolate's matched
  reference carries another allele at the site). IS6110 and ACC cells carry no EV at all.

## 2. Read check of disagreeing cells

Who the reads favour (paired records; "ALT side" = the build that wrote ALT; `out/readcheck_results_ctx.tsv`):

| class | ALT side | cells | reads ALT (ALT side right) | of which the other build is ALT at a nearby record | reads REF (other side right) | unclear |
|---|---|---|---|---|---|---|
| DEL | NEW | 30 | 5 | 2 | 6 | 19 |
| DEL | OLD | 30 | 11 | 7 | 3 | 16 |
| INS | NEW | 30 | 16 | 6 | 1 | 13 |
| INS | OLD | 30 | 15 | 13 | 3 | 12 |
| IS6110 | NEW | 30 | 29 | 0 | 0 | 1 |
| IS6110 | OLD | 13 | 8 | 0 | 4 | 1 |
| ACC (vs NOCALL) | OLD | 10 | 10 | - | 0 | 0 |
| ACC (vs NOCALL) | NEW | 10 | 10 | - | 0 | 0 |

"Nearby" = the build that wrote REF/ABSENT gives the same isolate ALT at another record of the same
class within 100 bp: the disagreement is then which allele record carries the call, not whether the
event is there.

- **IS6110: NEW is right.** 29/30 NEW gains are real insertions (16-96 junction reads, no spanning
  reads); OLD wrote REF at those sites (e.g. 1694664, 475214, 3797827, 850087, 2166312, which have
  119-144 carriers among the 525 in NEW and 0 in OLD). Of the 13 OLD-only IS6110 ALT cells, the reads
  support OLD in 8: three 5-23 bp from an H37Rv copy end, one at 888,763 with junction reads in 25% of
  reads (mixed), and four where NEW writes ABSENT but the reads show an IS6110 junction and no reads
  crossing the site (an IS6110-flanked deletion edge). NEW is right in 4 (no junction reads, 47-91
  spanning reads).
- **INS: NEW gains are mostly real** (16 ALT, 1 REF); 8 of the 13 unclear are the 2.1-3.5 kb graph
  alleles at 3,528,851-3,528,974 (PPE region), where reads cross the anchor unbroken but have 40-96
  clipped reads 300-700 bp away (a complex allele, not a simple insertion at the anchor). OLD-only INS
  ALT cells are mostly real too (15 vs 3), but 13 of those 15 have NEW ALT at a neighbouring allele.
- **DEL: mostly undecidable from short reads** (35/60 unclear: repeat / MAPQ-0 interiors, flank depth
  < 5, or partial depth). Decidable NEW gains split 5 right / 6 wrong. OLD-only ALT: OLD right 11 vs 3,
  but 7 of the 11 have NEW ALT on a nearby deletion record; 4 are deletions NEW now misses outright.
- **ACC: the ALT side is right in all 20**: the locus sequence is present (breadth 0.88-1.00 over its
  H37Rv-unique positions). So the 170 OLD ALT -> NEW NOCALL cells are lost calls, and the 238
  NOCALL -> ALT cells are recovered ones.
- Unpaired records (ALT with no counterpart, >= 2 carriers; `readcheck_results_unpaired_ctx.tsv`):
  OLD-only DEL 5 ALT / 0 REF / 10 unclear, NEW-only DEL 4 / 3 / 8, OLD-only INS 8 / 2 / 5,
  NEW-only INS 9 / 0 / 6.
- OLD ALT -> NEW NOCALL deletions (`readcheck_results_nocall.tsv`): `no_coverage` small alleles 5 ALT /
  2 REF / 13 unclear; interval `ambiguous_depth` 2 / 1 / 7. NOCALL is mostly defensible (repeat or
  uncovered regions), but about a quarter are clear deletions with covered flanks.

## 3. Tree consistency (Fitch changes)

525 shared isolates, cohort IQ-TREE tree (`out/homoplasy_summary_gwas1000tree.tsv`); records with
>= 2 ALT and >= 2 REF leaves:

| set | class | records | median changes | >= 2 changes | >= 5 changes |
|---|---|---|---|---|---|
| OLD | DEL | 1,468 | 6 | 81% | 56% |
| NEW | DEL | 1,213 | 4 | 74% | 48% |
| OLD | INS | 300 | 1 | 46% | 14% |
| NEW | INS | 1,113 | 1 | 49% | 18% |
| OLD | IS6110 | 209 | 2 | 58% | 5% |
| NEW | IS6110 | 322 | 2 | 53% | 5% |
| OLD | ACC | 67 | 5 | 99% | 52% |
| NEW | ACC | 41 | 7 | 98% | 54% |

Outside the repeat mask (`out/homoplasy_by_mask_gwas1000tree.tsv`): DEL OLD 210 records, 51% >= 2 ->
NEW 152, 58%; INS OLD 46 records, 89% >= 2 -> NEW 184, 46%. Inside the mask DEL is 86% (OLD) / 76% (NEW)
>= 2 and 1,061 of NEW's 1,213 scored DEL records are masked. The DEL and ACC homoplasy is dominated by
masked / PE-PPE / IS6110-adjacent loci in both builds.

SPLIT-1 simulation (event = SPLIT-1 group; `out/split1_groups_gwas1000tree.tsv`):

| set | type | merged groups (records) | scored events | >= 2 changes, records -> events | total changes in merged groups: members summed -> merged |
|---|---|---|---|---|---|
| OLD | DEL | 362 (751) | 1,329 | 81% -> 86% | 4,088 -> 3,800 |
| NEW | DEL | 506 (1,143) | 963 | 74% -> 83% | 11,022 -> 11,149 |
| OLD | INS | 208 (595) | 242 | 46% -> 70% | 721 -> 954 |
| NEW | INS | 395 (1,339) | 656 | 49% -> 56% | 2,635 -> 2,957 |

Merging does not lower the change count: merged loci are where the member records sit in different
clades, so one event over their union needs as many changes as the members did separately (and
"REF if any member REF" turns some unknown leaves into REF). SPLIT-1 will reduce the number of tests
(NEW: 1,143 DEL and 1,339 INS records into 506 and 395 events) but it will not make the SV events less
homoplastic; the share of events with >= 2 changes goes up.

Cross-check on the new FastTree (`out/homoplasy_summary.tsv`, `out/homoplasy_by_mask.tsv`), and NEW on
all 1,000 isolates (records with >= 2 ALT and >= 2 REF; share with >= 2 changes):

| class | OLD 525 | NEW 525 | NEW 1,000 | NEW 1,000 outside mask |
|---|---|---|---|---|
| DEL records | 80% (n 1,468) | 71% (1,213) | 66% (1,577) | 41% (252) |
| DEL SPLIT-1 events | 85% (1,329) | 79% (963) | 75% (1,209) | 40% (244) |
| INS records | 51% (300) | 49% (1,113) | 51% (1,357) | 47% (228) |
| INS SPLIT-1 events | 74% (242) | 55% (656) | 56% (803) | 42% (165) |
| IS6110 | 55% (209) | 42% (322) | 44% (607) | 43% (475) |
| ACC | 99% (67) | 98% (41) | 98% (46) | 92% (13) |

The two trees agree on the direction of the OLD -> NEW changes except INS records, which are flat on both (46 -> 49% IQ-TREE, 51 -> 49% FastTree). Outside the mask the FastTree gives
lower DEL homoplasy than the IQ-TREE tree for both builds (OLD 40% vs 51%, NEW 43% vs 58% >= 2), so the
unmasked DEL comparison is within tree noise; the INS improvement outside the mask (OLD 87-89% -> NEW
45-46%) holds on both.

## 4. Problems worth fixing (things the fixes did not cover or introduced)

1. **IS6110 ABSENT is unchecked and site-wide.** IS6110 cells carry no EV, so ABS-1 did not reach them.
   46 NEW IS6110 sites are ABSENT in >= 90% of isolates (45,024 cells; OLD had 36 such sites), e.g.
   4,230,303: ALT 3, ABSENT 997, while SAMD00089643 has 116x MAPQ>=20 depth across the site and SNPs
   are called beside it. On the 525, 23,360 REF -> ABSENT and 15,047 ABSENT -> REF cells, whole sites
   switching. GT=2 at these sites removes ~all isolates from any IS6110 test there. Pre-existing, but
   larger in NEW.
2. **Accessory presence: NOCALL rose from 7% to 43%** (shared isolates; EV empty), and only 99 loci are
   emitted (OLD 802; the build catalogue has 564). Mosaic loci are ~60% NOCALL. All 20 checked
   ALT<->NOCALL cells have the locus present, so NOCALL is discarding real presence calls.
3. **Long "small_allele" insertions are graph alleles, not HaplotypeCaller calls.** 486 NEW
   small_allele INS records are >= 500 bp (21,459 ALT cells on 1,000), up to 3,501 bp, and several
   carry the same position with near-identical sequence (five 2,143 bp alleles at 3,528,974). At the
   3.53 Mb PPE locus the reads do not show an insertion at the anchor but a complex rearrangement
   300-700 bp away. These calls are not positional insertions as written; SPLIT-1 will merge them
   into one event, which is the right unit, but the anchor position should not be read as the site.
4. **Deletions written NOCALL `no_coverage` include real deletions.** 2,773 OLD ALT cells on small_allele
   deletions are now NOCALL with EV `no_coverage`; in 20 checked, 5 are clear deletions (covered
   flanks, no reads inside), 2 are REF, 13 undecidable. No coverage inside a deletion with covered
   flanks is the deletion itself.
5. IS6110-flanked deletions: 4 of 13 OLD-only IS6110 ALT cells are NEW ABSENT with a clear IS6110
   junction in the reads; the ABSENT (region deleted) and the IS6110 presence at the edge are both true,
   but the IS6110 call is lost.

## Files

- `extract.py`, `compare.py` -> `out/{old,new}.pkl`, `out/pairs.tsv`, `out/events.tsv`, `out/transitions.tsv`,
  `out/new_change_ev.tsv`, `out/disagree_cells.tsv.gz`
- `sample_cells.py`, `sample_unpaired.py`, `sample_nocall.py`, `readcheck.py`, `readcheck.sbatch`,
  `context.py` -> `out/readcheck_*.tsv`, `out/is6110_absent_sites.tsv`, `out/state_counts_shared.tsv`
- `build_tree.py`, `tree.sbatch` -> `out/snps.fasta`, `out/new1k.fasttree.nwk`
- `homoplasy.py`, `homoplasy_summary.py` -> `out/homoplasy_*`, `out/split1_groups*.tsv`
- `tables.py`: markdown tables from `out/`
