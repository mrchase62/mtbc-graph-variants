# Panel selection by demand (new1k), 2026-10-13

Analysis only; nothing in the pipeline changed. Cost: ~2 CPU-hours
(QC 52755519, low-depth 52755541, selection 52757589).

1. `qc_candidates.py` -> `candidates.tsv`: 596 candidates (CX333 333, Marin
   90, GenBank 173). Excluded: single-base run >= 100 bp (4: GCF_039770655,
   GCF_045348265, GCF_050259585, GCA_061390545), indel / private-homopolymer
   excess (31 CX333), GenBank flags and PRJNA994284, non-L1-4,7 lineages.
   Pool 426 (CX333 266, Marin 88, GenBank 72). Forced: H37Rv,
   GCF_022870265 (L5.1.1), GCF_000253355 (L6.1.1), 3 lineage 7 genomes.
2. `lowdepth.sbatch`: depth < 5 intervals from P1's kept H37Rv BAMs.
3. `select_panel.py` -> `out/`: SNP distance (outside the repeat mask) from
   each new1k isolate's own P1 calls to each genome's direct-alignment
   states; greedy k-median; cross-validated on lineage-stratified halves.

## Results (`out/summary.txt`)

- Validation: nearest current CX333 genome = P1's choice for 776/1000;
  r = 0.98 with P1's snp_distance (values run higher; read as relative).
- Gain stops at 174 genomes on all 1,000 isolates: beyond rank 174 no
  genome is the nearest for any isolate, and the order is arbitrary ties.
- Held-out (cross-validated) median: CX333 now 207.5 / 221.5; 150 genomes
  199.5 / 229; 250 genomes 189.5 / 217; no gain past ~250-300.
- Gaps that no candidate closes (median > 300 at any size): lineage 1.1.3
  (558), 1.2.2 (307), 1.1.1 (303), unresolved "lineage3" (286, 41 of 122
  over 300), 4.1, 4.3.1, 4.6.3-4.6.5. These need sequencing.
- Lineage 4.6.1 gets worse (182 -> 260) because GCF_045348265 (poly-A run)
  is excluded.

## Selection on the whole passing pool (`select_pool.py` -> `out_pool/`)

User: "yes, start it" (pool = standard screen, L1-4, fold-back <= 2.5;
33,364 isolates). Isolate calls from the collection's Mutect2 pickles; no
new alignment. Lineages weighted equally. Job 52760614, 8 min at 8 cores.

- Validation on new1k: same nearest genome as the P1-based distance for
  99.2%; r = 1.000. Missing per-base depth does not matter.
- Held-out median: CX333 now 150 / 152; 150 genomes 147 / 149; 175 146 /
  148; 200 145 / 147; flat after 200. Over 300: 13.1% now, 11.9% at 175+.
- Gaps (most isolates > 300 at any size): L1.1.3 (88%), L4.1 / 4.1.4 /
  4.3 / 4.3.1 / 4.4 / 4.6.x (90-100%), L4.5 (74%), L1.2.2 (53%),
  L1.1.2 (44%). These need sequencing.
