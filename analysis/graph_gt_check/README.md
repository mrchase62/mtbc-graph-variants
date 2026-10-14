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

## QC-1 and the split loci (2026-10-14)

- **QC-1 was our bug, not the graph's.** `analysis/panel_checks/genomes.tsv`
  pointed QC-1 at QC-10.fna.gz and mada_1-1 at mada_1-10.fna.gz (prefix
  match). Fixed; their states regenerated (old ones in
  `analysis/panel_checks/states_old/`). Re-scored: QC-1 now ordinary; SNP
  discordance outside the mask 0.004% (was 0.014%).
- **Panel impact:** QC-1 was selected (rank 117) on QC-10's sequence.
  Re-run with corrected states (`analysis/panel_select/out_pool_fixed/`): the
  top 200 would swap one genome (QC-6 out, GCF_977011625 in); median
  distance 147 and 11.95% over 300 are identical. No rebuild needed.
- **Split loci** (`split_loci.py` -> `out/split_loci.tsv`; 419 loci with
  disjoint carriers):

  | class | outside mask | inside mask |
  |---|---:|---:|
  | copy-number series (VNTR repeats, the 1,674 bp Rv1318c-Rv1320c duplication) | 10 | 72 |
  | IS6110 at one site, copies differing | 58 | 28 |
  | one event, alleles within 5% length (internal SNPs/indels) | 20 | 40 |
  | other (mostly one large insertion with a few length outliers) | 12 | 179 |

  Examples outside the mask: a ~5 kb insertion at 2,268,721 written as 57
  records (5,000 or 6,358 bp); 4.5 kb at 2,219,417 as 23 records. Each
  record carries a scattered subset of genomes, so per-record homoplasy is
  inflated. A large insertion whose copies differ by a base becomes several
  records; the graph cannot merge them.
