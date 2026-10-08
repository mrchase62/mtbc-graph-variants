# Novel structural events: development plan (proposal, 2026-10-08)

The user: absorbing SVs from best-match genomes helps, but external SV
callers are not working well and novel events must be captured, so develop
this soon, within the current pipeline testing.

Status: proposal only. Nothing is run. Compute and any download need the
user's approval.

## What "novel" means here

P1 picks each sample's matched reference R from the panel, and P2 aligns
the sample's reads to R. An event the panel already carries is mostly
absorbed: it is in R, or in the graph. A **novel event** is one that differs
between the sample and R. So it is found in the **P2 alignment to R**, not
the H37Rv one, and placed afterwards on H37Rv and the graph by P4b, as
caller SVs already are.

Today the only novel-event callers are delly and dysgu, run in P2 against R.
In the ctpV check, delly called nothing and dysgu called fragments, and
neither can call an insertion longer than a read (`analysis/ctpV_check.md`).

## The method: a breakpoint-evidence caller

This is the same evidence as the ABS-1 diagnosis, read genome-wide instead of
at known deletion ends:

| signal in the P2 BAM (reads vs R) | event |
|---|---|
| clip pile at A, SA partners at B (same strand), low depth between | deletion A-B |
| clip piles on both sides of one position, no joining reads, clipped sequence and unmapped mates not in R | insertion |
| SA partners on the opposite strand | inversion |
| SA partners far away or on another contig | translocation, or a duplicated or relocated segment |
| read pairs too far apart or too close, or wrongly oriented | support for all of the above |

Steps for each sample:
1. **Collect** the signals above from the P2 BAM, MAPQ >= 20 for anchors.
2. **Cluster** them into candidate breakpoints.
3. **Assemble** each candidate locally: the audit assembly tool
   (`analysis/inherited_absent/local_assembly.py`; SPAdes, then minimap2),
   pointed at R instead of H37Rv. This gives the exact breakpoints, and the
   inserted sequence.
4. **Classify** the inserted sequence with blastn against IS6110, the
   accessory catalogue, R and H37Rv: known element, known accessory,
   duplication of R's own sequence, or novel.
5. **Output** a VCF in R's frame, for P4b to place.

Across the cohort: recurrent events get a REF/ALT junction test in every
sample, the way SV intervals are genotyped now.

One module, `bin/breakpoint_evidence.py`, would serve both ABS-1 and this
caller. Building it for ABS-1 is the first step of this work.

## Testing within the current pipeline testing

### Phase A: simulated reads from real genomes outside the panel (no download)

- **Genomes:** the external assemblies already on disk
  (`analysis/external_assemblies/`: 151 Marin hybrid assemblies and 4
  Behruznia long-read genomes, all outside CX333), plus the 20 simbench
  genomes. Start with about 40 across lineages, using the Marin isolates
  with no QC flags first (`analysis/marin_assembly_qc/`).
  - Panel genomes cannot be used: P1 would pick the genome's own assembly as R.
- **Truth:** each assembly against its matched reference R
  (minimap2 `-x asm5` + paftools call; syri for rearrangements), 50 bp or
  more. These are real, natural novel events relative to the best panel
  match, which is the quantity that matters.
- **Reads:** wgsim, as in `analysis/is6110_simbench/`: 150 bp pairs, fragment
  450 +/- 50, 0.2% base error, about 80x.
- **Run:** P1 and P2 through the normal runner, as a separate test cohort
  (the simbench pattern), then the new caller. dysgu and delly come free
  from P2.
- **Score:** recall and precision by class (DEL/INS/INV/other), size
  (50-500, 500-read pair span, longer), repeat context (masked or not) and
  inserted-sequence type. The new caller against dysgu and delly, alone and
  in union.
- **Cost:** about 0.75 billing-hour per genome (P1 + P2, from scale200_fix),
  so about 30 billing-hours for 40 genomes. The caller adds little.

### Phase A2: planted events

Insert or delete random and foreign sequence of chosen sizes at chosen
positions in an assembly (unique and repeat contexts), and simulate reads.
This gives recall curves by size where natural events are few. The cost is
small.

### Phase B: real reads (needs download approval)

Simulated reads have no real error profile, GC bias or library artefacts.
Marin sets A, B and E have public Illumina and long reads: 63 isolates, with
their assemblies already local. About 25 GB of short reads is needed for the
caller test, but not the long reads. Same scoring as Phase A.

### On scale200_fix, without a truth set

Run the caller on the 200 existing P2 BAMs; no P1/P2 rerun is needed.
- **Tree consistency:** a true event should fit the tree, in one clade or a
  few. Scatter across unrelated tips suggests false calls.
- **Agreement** with dysgu and delly calls.
- **Positive controls:** the ctpV events, against each sample's R.

## How it joins the chain

- **New step P2n after P2:** reads the P2 BAM and writes the novel-event VCF
  in R's frame.
- **P4b** places it, P5 genotypes recurrent events, and the merged VCF
  carries them with SVSOURCE=breakpoint.
- **Gate:** it enters the chain only after Phases A and B show it beats
  dysgu/delly at equal precision, and the user approves.
- **Timing:** it reads only existing P2 BAMs. If it is ready when the fix list
  is complete, it can join the same P5-onward rerun as ABS-1; if not, it can
  be added later as an extra pass, without rerunning P1/P2.

## Order

1. `bin/breakpoint_evidence.py`, shared with ABS-1 (fix-list work).
2. Phase A truth sets: assembly vs R for about 40 external genomes. Small
   compute.
3. Prototype caller and Phase A scoring.
4. Phase A2 planted events.
5. scale200_fix tree-consistency check.
6. Phase B, if the user approves the download.
7. Proposal to add P2n, for approval.
