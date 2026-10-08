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

## Read-level evidence in the four L1.2.1 samples (2026-10-08)

### HaplotypeCaller on H37Rv (P1 VCF, 1,078,500-1,081,100)

- **SAMEA2297133:** two SNPs, both outside the start-codon segment:
  1,079,927 C>A (DP 117) and 1,080,192 G>A (DP 126). No deletion. Nothing at
  1,078,644 or 1,078,715, so REF at about 72x. This agrees with an intact gene.
- **SAMEA5542103:** only 1,080,192 G>A, past the deletion; no reads inside it.
- **SAMEA112800746 and SAMN07766100:**
  - 14 SNPs and 1-12 bp indels at 1,079,879-1,079,927, at depth 4 rising to
    50, inside the "deleted" segment;
  - these look like reads from rearranged ctpV pieces (below) forced onto H37Rv;
  - HaplotypeCaller's local reassembly spans about one read length and cannot
    represent the event.

### Soft clips (>= 5 bp, piles of >= 3 reads) and split reads (SA) on H37Rv

| breakpoint | SAMEA5542103 | SAMEA112800746 | SAMN07766100 | SAMEA2297133 | SAMEA1119809 (control) |
|---|---|---|---|---|---|
| ~1,078,519-522 | 65 | 51 | 42 | **0** | 4 |
| 1,078,757 | - | 12 | 24 | **0** | 0 |
| 1,078,811 | - | 12 | 24 | **0** | 0 |
| 1,079,652 / 1,079,733 | 42 / 38 | 42 / 45 | 29 / 27 | **0** | 0 |
| 1,079,759 / 1,079,827 / 1,079,921 | - | 39 / 39 / 33 | 33 / 39 / 37 | **0** | 0 |
| ~1,080,103-104 | 56 | 49 / 47 | 35 / 31 | **0** | 0 |

- **The carriers have sharp two-sided clusters with split-read partners.**
  - In SAMEA112800746 and SAMN07766100, the 235-bp deletion joins 1,078,522 to
    1,078,757: clips at one breakpoint have their SA at the other.
  - SAMEA5542103 joins 1,078,519 straight to about 1,080,104, a single larger
    deletion like the panel's 1,536 bp, not its reference's two-piece version.
  - In all three, the "1,290-bp deletion" is not clean. Pieces at about
    1,079,652-733, 1,079,759-827 and 1,079,921-1,080,103 remain, joined in a
    different order. These are probably what the graph records as the 337-bp
    insertion.
- **SAMEA2297133 has no cluster at any breakpoint.** Its only pile is about 27
  one-sided clips spread over 1,079,080-1,079,163, with no SA partners. The
  control shows the same at 1,079,080, so it is background.
- **Together:** SAMEA2297133's reads (even depth, no junctions, HaplotypeCaller
  REF) say ctpV is intact. Its four ABSENT cells come only from the reference.

### delly and dysgu (P2, against each sample's matched reference)

P2 calls SVs against the matched reference R, not H37Rv. On R, ctpV is at about
2,814,000-2,819,600 in GCF_040208995 (reverse strand; csoR at 2,818,290-649)
and about 2,470,000-2,475,500 in GCF_040209325.

| sample (R) | delly | dysgu |
|---|---|---|
| SAMEA2297133 (GCF_040208995) | INV 2,816,189-3,869,989 and INV 2,817,770-811, both LowQual, GT 0/0 | DEL 42 bp, INS 111 bp, INV 33 bp at 2,818,015; all lowProb |
| SAMEA112800746 (GCF_040208995) | none | INS 270 bp at 2,818,015, PASS, 1/1 |
| SAMN07766100 (GCF_040208995) | none | INS 173 bp at 2,818,015, PASS, 0/1 |
| SAMEA5542103 (GCF_040209325) | INV 2,474,533-559, LowQual, GT 0/0 | DEL 179 bp at 2,471,469-648, PASS, 1/1; INV 86 bp lowProb |

- **delly finds none of it.** Its calls are low-quality inversions with
  genotype 0/0.
- **dysgu finds fragments of the carriers' differences from R.**
  - The PASS insertions at 2,818,015, the junction of R's two deletions, are
    most likely a partial view of the 337-bp insertion the two carriers have
    and GCF_040208995 lacks.
  - SAMEA5542103's PASS 179-bp deletion fits it losing more of the gene than
    its reference does.
- **Neither sees SAMEA2297133's real difference from R.** Relative to
  GCF_040208995, SAMEA2297133 carries about 1.5 kb more sequence: the 235-bp
  and 1,296-bp segments its reference lacks. On R that is an insertion longer
  than a read. Short-read callers rarely call insertions that long, and delly
  hardly at all; dysgu reports only lowProb noise near the junction.
- **By design, a deletion the sample shares with R is invisible to both.**
  Showing it is the job of the H37Rv frame and the graph.

### Cohort-wide

The same false-ABSENT pattern across all 200 scale200_fix samples is counted in
`analysis/inherited_absent/README.md`. Of the testable ABSENT cells, 12% (an
upper bound) are covered at normal depth; they concentrate in phiRv1,
plcA/plcB, PPE57/58, Rv3766-70 and wag22.
