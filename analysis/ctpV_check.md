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

## Calls in scale200_fix (2026-10-08)

The cohort has 19 lineage-1 isolates. Four of them are in L1.2.1:
SAMEA112800746, SAMEA5542103 and SAMN07766100 (lineage caller 1.2.1.2.1),
and SAMEA2297133 (lineage caller only 1.2.1). All four were mapped in P2 to
an L1.2.1.2.1 panel genome that lacks the ctpV start: GCF_040208995 for three
of them and GCF_040209325 for SAMEA5542103.

**The cohort catalogue clusters the panel deletions.** It has no separate
297-bp interval. The 297, 239 and 235-bp deletions form one interval,
`svi:DEL:1078518:239` (13 panel carriers), the "start-codon" interval. The
1,536, 1,296 and 1,290-bp deletions form another, `svi:DEL:1078519:1536`
(6 panel carriers).

**Read depth on H37Rv:** mean depth and zero-depth positions per segment, from
the P1 H37Rv BAM, MAPQ>=20.

| sample | start-codon segment 1,078,519-756 | 1,078,814-1,079,742 | merged VCF: 239 interval | 1536 interval | 872 / 82 intervals |
|---|---|---|---|---|---|
| SAMEA5542103 | 0x (236/238 zero) | 4x (847/929 zero) | ALT | ALT | ABSENT |
| SAMEA112800746 | 1x (228/238 zero) | 4x (845/929 zero) | ALT | NOCALL | ABSENT |
| SAMN07766100 | 2x (14/238 zero) | 4x (426/929 zero) | ALT | NOCALL | ABSENT |
| SAMEA2297133 | 72x (0 zero) | 107x (0 zero) | REF | REF | **ABSENT (wrong)** |
| SAMEA1119809 (L1.2.2.2, control) | 56x | 67x | REF | REF | REF |

**Results:**

- **The three L1.2.1.2.1 isolates are correctly called ALT for the deletion
  over the ctpV start codon.** Their reads cover neither the start codon nor
  most of the gene, which fits the 235-bp plus 1,290-bp haplotype of their
  matched references.
- **SAMEA2297133 has an intact ctpV.** It has full depth with no gaps across
  the whole gene. It is probably a basal L1.2.1 isolate, the branch the paper
  excludes. Its start-codon and 1536 intervals are correctly REF, from the
  two-frame check ("inherited_reference_contradicted").
- **Its 872 and 82-bp intervals are wrongly ABSENT.** They lie inside the
  1,296-bp deletion of its matched reference, so no probe projects and
  `p5_sv_genotype.py` calls ABSENT. Unlike the enclosing interval, these are
  not checked against H37Rv-frame depth.
- **Its small variants have the same problem in part.** P4b also places the
  reference's 1,296-bp deletion in this sample as INHERITED. Small-variant
  states are mostly NOCALL inside it (27 NOCALL, 1 ALT). Two small records in
  the start-codon segment, 1,078,644 and 1,078,715, are ABSENT although the
  reads cover them at 72x.
- **The other 15 lineage-1 isolates (L1.1.x, L1.2.2.x) are REF on every ctpV
  interval.**

**Not looked at:** whether other samples get the same wrong ABSENT where their
matched reference has a deletion the sample does not. In the comparison, this
belongs with the NOCALL/ABSENT counts.
