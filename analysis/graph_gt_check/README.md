# Graph genotype check, P200 build 3893a3762ed1 (2026-10-14)

User: "yes, run the graph genotype check". Analysis only. Each panel
genome's graph VCF genotype vs its own direct alignment to H37Rv (SNPs from
analysis/panel_checks/states; indels/MNPs/SVs from outgoup kit route A,
`kit.sbatch`, matched with panel_ancestral.kit_event_states). `score.py` ->
`out/summary.txt`, `out/records.tsv`, `out/sv_records.tsv`. ~2 CPU-hours.

- Discordance (cells called by both), outside / inside the repeat mask:
  SNP 0.014% / 0.24%; MNP 0.02% / 0.36%; indel 0.11% / 0.30%; SV 0.18% / 0.48%.
- One outlier genome: QC-1 (Marin set E, Peker 2021 ONT, heavily polished)
  2.0% (867 cells, nearly all SNPs); every other genome <= 0.11%, median
  0.067%. Not investigated yet.
- Representation: 422 loci with >= 2 SV records of one type within 50 bp,
  411 with disjoint carriers (one event possibly written as several records).
- SV homoplasy barely moves when discordant cells are removed: >= 2 changes
  739 -> 680 (all), 175 -> 165 (outside the mask); 599 of the 3,082 SV
  records with no discordant cell still change >= 2 times. So genotype
  errors do not explain the homoplasy.
