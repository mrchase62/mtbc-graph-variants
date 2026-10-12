# Homopolymer tracts in the graph's panel genomes (2026-10-12)

User: "we need to redo the graph because of homopolymeric tracts in some of
the genomes". This is a report only; nothing was changed. Script: `scan.py`
(run from runroot, about 3 minutes on one core). Report: `out/report.txt`.
Tables:
- `out/runs.tsv`: every pure run of 15 bp or more, placed on H37Rv;
- `out/tracts.tsv`: near-pure tracts;
- `out/loci.tsv`: runs grouped by H37Rv locus;
- `out/per_genome.tsv`.

Genomes: the 333 references of build 7713a8d71d8e-fix1, which are the
CX333 graph's panel. Coordinates are those of the deposited assemblies
(build `refs`).

## What was measured
1. **Pure runs:** every run of a single base of 15 bp or more.
2. **Near-pure tracts:** stretches in which one base makes up at least 90%
   of every 50 bp window. These catch tracts with a few interrupting bases.
3. **Where each pure run sits on H37Rv,** found from its 40 bp flanks, and
   how long H37Rv's own run is there. A locus where many genomes carry runs
   of varying length is a natural variable tract. A run found in one genome
   only, and far longer than any other, is an assembly artifact.
4. **Which cohort samples use an affected genome** as their matched
   reference (every cohort in `refbias/cohorts.tsv`).

## Result: three genomes, as found by the 2026-10-04 QC

| genome | lineage / strain | tract (deposited coordinates) | longest pure run | likely source |
|---|---|---|---|---|
| GCF_050259585 | La2, M. caprae SY-1 | 1,585,213-1,586,349: 1,137 bp, 99% G | 1,015 bp G | two-colour Illumina no-signal poly-G; "PacBio" label, SOAPdenovo assembly |
| GCF_045348265 | L4.6.1.1, ITM-2003-01539 | 546,692-547,661: 970 bp, 98.5% A; then 547,680-547,904: 225 bp, 94% T | 484 bp A, 464 bp A, 150 bp T | ONT (Flye) |
| GCF_039770655 | L9, ITM-2021-00504 | 2,070,519-2,071,038: 520 bp, 97% T, inside (CCATT)n | 377 bp T | Nottingham Trent set (assembly method not checked) |

**The separation is clean.**
- Across the other 330 genomes the longest pure run is 48 bp. 284 genomes
  have none of 15 bp or more.
- No run falls between 48 and 107 bp. A run-length rule of 100 bp (the
  default in `bin/assembly_qc_stats.py`) catches exactly these three.
- A pure-run rule alone under-masks GCF_045348265. Its artifact is two
  tracts of about 1.2 kb in all, with A- and T-rich stretches between the
  pure runs. Masking should use the tract coordinates.

## Natural variable G/C tracts (not artifacts)
59 of the 80 runs of 15 bp or more were placed on H37Rv, at 16 loci. They
are G or C runs that recur across lineages, with lengths varying between
genomes:
- **2,338,189:** C runs of 15-42 bp in 23 genomes; H37Rv has 2 C there.
- **1,623,086:** G runs of 18-46 bp in 9 genomes; H37Rv has 3.
- **976,897:** G runs of 15-22 bp in 6 genomes.
- **1,700,416 and 854,253:** G and C runs of 15-27 bp in 4 and 3 genomes.

These are variable mononucleotide tracts, whether the length differences are
real or long-read length errors. They are not the kilobase artifacts and
were not flagged.
- **GCF_039555365 (La2) has several private G runs up to 48 bp,** e.g. at
  3,840,780, where H37Rv has none. These are within the natural range seen
  at other loci, so not flagged; noted as a genome to watch.

## Which cohorts are affected (matched references)

| cohort | GCF_039770655 (L9) | GCF_045348265 | others |
|---|---|---|---|
| pilot, scale100, pilot_rerun, scale200, scale200_fix | 2 each | | |
| gwas1000 | 14 | 1 | GCF_965121955 9 (its "tract" is a 57 bp A-rich stretch, 88% A, longest run 18: natural, not an artifact) |

GCF_050259585 is no cohort sample's reference.

## Decision needed for the rebuild (user)
Per genome:
- **mask the tract with N,**
- **replace the assembly** (e.g. reassemble from its reads), or
- **exclude the genome.**

Things that bear on the choice:
- **GCF_039770655 is CX333's only lineage 9 genome.** Excluding it removes
  L9 from the graph, and 14 gwas1000 isolates would lose their matched
  reference.
- **GCF_050259585 (SY-1)** has a close La2 neighbour in the panel,
  GCF_039555365.
- **GCF_045348265** has a close L4.6.1.1 neighbour, GCF_045348245.

The 2026-10-04 note proposed masking with N and adding the run-length rule
to the assembly QC. The rule exists (`assembly_qc_stats.py`, flag
`single_base_run`); `build_panel.py` refuses to write a panel until each
flagged genome has a recorded decision. Masking must be done on the
rebuild's input FASTA, by sequence, since panel coordinates may be rotated.
