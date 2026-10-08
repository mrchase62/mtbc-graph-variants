# ctpV (Rv0969) deletions in the CX333 graph (2026-10-08)

Question: does the CX333 graph contain the ctpV deletion reported by
"Genome graphs reveal the importance of structural variation in
*Mycobacterium tuberculosis* evolution and drug resistance" (Nat Commun 2025,
doi 10.1038/s41467-025-65779-9, PMC12663572)?

What the paper reports (text; exact SV IDs are only in its Fig. 5b image):

- several different deletions of ctpV, found only in L1.2.1 and absent from its
  deepest branch;
- the deletions get longer further from that branch;
- every one overlaps the start codon, so the gene is likely lost;
- the isolate used for RNA-seq carries a confirmed 297-bp deletion.

## Answer: yes

The graph has a 297-bp deletion of the same size at the gene's start. The
whole pattern the paper describes is also there.

ctpV is H37Rv 1,078,743-1,081,055 (+ strand), so its start codon is
1,078,743-745. Rv0968 (1,078,391-687) and csoR (1,077,975-8,334) lie just
upstream, in the same operon. The positions below are deleted H37Rv bases.
They come from the fixed collapsed VCF of build 7713a8d71d8e-fix1
(`refbias/build/7713a8d71d8e-fix1/assets/graph_collapsed.vcf.gz`, records
`>126446>126583_*`). Lineages are from `analysis/pangenome_compare/cx333_panel_lineages.tsv`.

| genomes | sublineage | events in the region | covers ctpV start |
|---|---|---|---|
| GCF_014900005, GCF_022870445, GCF_030323805, GCF_040209215, GCF_040209395 | 1.2.1.2 | **297-bp DEL 1,078,519-1,078,815** | yes |
| GCF_040209205 | 1.2.1.2.1 | 1,536-bp DEL 1,078,519-1,080,054; 274-bp INS at 1,080,093 | yes |
| GCF_009730215, GCF_014899965, GCF_040208985, GCF_040209325 | 1.2.1.2.1 | 235-bp DEL 1,078,522-1,078,756; 337-bp INS at 1,078,811; 1,290-bp DEL 1,078,814-1,080,103 | yes |
| GCF_040208995 | 1.2.1.2.1 | 235-bp DEL 1,078,522-1,078,756; 1,296-bp DEL 1,078,811-1,080,106 | yes |
| GCF_033124885, GCF_033125685 | 1.2.1.2.1 | 235-bp DEL 1,078,522-1,078,756; 872-bp DEL 1,078,871-1,079,742; 87-bp INS at 1,079,743 | yes |
| GCF_014899985 | 1.2.1.2.1 | 239-bp DEL 1,078,518-1,078,756; 872-bp DEL 1,078,871-1,079,742; 87-bp INS at 1,079,743 | yes |
| GCF_965124535 | lineage 6 | 82-bp DEL 1,079,801-1,079,882 (inside the gene, start kept) | no |

**How this matches the paper:**

- All 14 panel genomes from L1.2.1.2 and below carry a deletion over the start
  codon. They carry five different deletions.
- None of the other 34 lineage-1 panel genomes has one, including the two from
  L1.2.1.1, the deepest branch of L1.2.1.
- The deletions in the deeper sublineage, L1.2.1.2.1, remove more of the gene
  (up to about 1.5 kb), so the deletions grow away from the branch split, as
  the paper says.

**Earlier catalogue:** the SV catalogue from the old collapse
(`analysis/pangenome_compare/cx333_features.tsv`, 2026-10-03) has the 297, 1,536,
1,296, 1,290 and 872-bp deletions and the 337 and 274-bp insertions, with the
same carriers. It has no 235/239-bp deletions. The fixed collapse (GRAPHVCF-5
trim-and-union) splits those haplotypes into the upstream 235/239-bp deletion
plus the downstream event. Under either collapse, every L1.2.1.2 genome lacks
the start codon.

**Not checked:** whether the scale200 or gwas1000 short-read calls genotype
these deletions. An earlier cross-check, `analysis/rd_crossref/novel_stable_gwas1000.tsv`,
lists the 297, 1,587 and 872-bp deletions as graph-only (not in the published
RD list), all in lineage 1.
