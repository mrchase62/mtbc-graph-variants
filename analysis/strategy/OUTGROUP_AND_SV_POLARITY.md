# How the tree and the SVs use the outgroup, now and in an L1-4 graph

2026-10-12. A refresher, written for the rebuild decision (scope: lineages
1-4 and 7 plus one or two lineage 5/6 genomes; canettii and lineage 8
outside the graph). Nothing here is new design: it gathers HANDOFF 0g, 0i,
0j, 0k, OUT-1 and `analysis/outgroup_test/README.md`, and says where the
code stands.

## 1. The tree is built from SNPs only

`assoc/bin/cohort_assoc_tail.sh`, steps 1-3b:
1. **Cohort SNP alignment** from the merged cohort VCF.
2. **Combined alignment:** the cohort plus every panel genome's SNP row
   (`combined_alignment.py`, rows by `add_outgroup.py`'s allele rule from
   the panel VCF). The outgroup, canettii ET1291 (GCF_035581225), is one of
   these rows. The panel genomes fix the deep branches the cohort alone
   cannot place.
3. **IQ-TREE,** rooted on the outgroup (`-o GCF_035581225`).
4. **Pruned** to the cohort, H37Rv and the outgroup.

SVs, IS6110 and accessory events play no part in building the tree.

### Without canettii in the graph
Canettii stays the root; only the source of its row changes.
- **The outgroup kit:** canettii's assembly aligned directly to H37Rv
  (`genome_states.py`, seconds per genome). It gives the base at every H37Rv
  position, and N where the alignment is unclear.
- **`add_outgroup.py` reads the kit instead of the graph VCF.** Its output
  and the rest of the chain stay the same.
- **Tested 2026-10-07:**
  - same root and lineage structure; 12 of 1,062 splits differ, all inside
    lineage 4 at low support;
  - 99.98% the same ancestral alleles;
  - the same scale200 association survivors and drug-resistance controls.
- **Lineage 8 goes in the kit too,** for questions about the root of the
  whole complex.
- **The one or two lineage 5/6 genomes are inside the graph.** They are the
  nearest outgroup to the L1-4,7 clade, so their rows come from the panel
  VCF like any panel genome.

## 2. SVs enter as columns of the cohort VCF and are scored on the tree

1. **Calling and placing, per isolate:**
   - P2 calls dysgu and, from 2026-10-12, our caller;
   - P4b places both breakpoints on H37Rv through the graph and adds the
     inherited half (the matched reference's own differences, from the
     graph VCF);
   - P5 clusters events across isolates (`sv_matrix.tsv`);
   - P5 SV genotyping (`p5_sv_genotype.py`) checks reference and alternative
     by read depth at probe positions, so "no call" and "reference" are
     kept apart.
2. **The merged cohort VCF** carries SVs as records of class `sv`, beside
   SNPs, small indels, `is6110` and `accessory_presence`.
3. **Ancestral and derived states per tree node** (`write_event_matrix.py`):
   - Fitch parsimony on the pruned SNP tree reconstructs every column's
     state at every node;
   - each state change is a gain or a loss on a branch.
4. **The convergence test** (`assoc_scan.py` / phyoverlap2) asks whether
   isolates below independent gains of a variant are enriched for the
   phenotype. So it matters which allele is the derived one: getting it
   backwards turns repeated gains into repeated losses, and the variant
   becomes untestable.

## 3. How an SV's direction (polarity) is decided today
In `write_event_matrix.py`, in this order:
1. **`AA`** (Fitch over the 333-genome panel tree). This covers SNPs only,
   since it comes from a SNP alignment. SVs never have it.
2. **`panel_polarity.tsv`** (`bin/panel_polarity.py`): the outgroup's own
   genotype in the graph VCF, for every variant class, at sites that are
   graph VCF records.
   - An outgroup ALT counts as ancestral only if the panel also carries it
     at 5% or more.
   - This is what fixed the tatC-helY 68 bp insertion. Canettii and 229
     panel genomes carry it and H37Rv does not, so it is really 31
     deletions, not 11 gains and 31 losses.
3. **Otherwise ALT is taken as derived.** For an insertion that is a
   definition: presence is derived. It is wrong wherever H37Rv itself
   carries the derived state.
4. **`--polarise-by-root`, a second pass,** can flip variants whose
   reconstruction puts the derived allele at the root.

The outgroup leaf in the event matrix has data only for SNP columns. For SV
columns its state is unknown, so an SV's root state comes from step 2 or 3
and the cohort's own parsimony. The MTBC-node pin (R2-TREES-6) applies only
to variants with `AA`.

**What step 2 covers.** It covers SVs that are records of the graph VCF:
differences between panel genomes. That includes the inherited half and any
caller SV that matches a graph record. A new SV found only in cohort reads
(not in any panel genome) has no panel record. It falls to ALT-as-derived,
in both the current graph and a new one.

## 4. What breaks when canettii leaves the graph, and the plan

**Step 2 loses its source.** With canettii out of the graph VCF, the
panel-polarity table has no outgroup genotype to read for SVs, indels,
IS6110 or accessory loci. (SNPs are unaffected: the kit gives their rows.)

**The plan has two sources, nearest outgroup first** (HANDOFF 0i point 4
and 0k):
1. **The lineage 5/6 genome(s) in the graph.** They are genotyped in the
   graph VCF at every panel variant, SVs included, so `panel_polarity.py`
   can read the outgroup state from them instead of from canettii. They are
   the nearer outgroup, so they are the better first choice.
   - **Code change needed:** `panel_polarity.py` reads one outgroup name
     from the build's `outgroup` key. It would take an ordered list:
     L5/6 first, then canettii's event states.
   - **The ancestral node also moves:** take it as the MRCA of the L1-4,7
     ingroup, not "the root's other child" (fault B, 0j).
2. **OUT-1, the canettii pseudo-isolate** (HANDOFF, added 2026-10-08, user):
   - simulate 150 bp read pairs at about 80x from the two canettii
     assemblies (wgsim, as `analysis/is6110_simbench/`);
   - run them through P1-P5 in their own registry row, like any isolate;
   - their SV, IS6110 and accessory calls are canettii's event states, and
     serve as the second polarity source and as the outgroup leaf's data in
     the event matrix;
   - checks: agree with scale200's read-based canettii isolate where both
     are called, and with the in-graph canettii's polarity on a build that
     still has it (the current CX333 build);
   - status: planned, not built. It can be built and tested on the current
     CX333 graph before the rebuild.

**Not covered by either source:** SVs absent from every panel genome and
from canettii. They stay ALT-as-derived, as today. For a new insertion this
is right by definition. For a deletion of sequence that canettii and L5/6
also lack, the pseudo-isolate would supply the state.

## 5. Open points to settle before the rebuild
- **Which lineage 5/6 genome(s) go in the graph.** Decided 2026-10-13
  (user): GCF_022870265 (L5.1.1, carries the tatC-helY element) +
  GCF_000253355 (L6.1.1), the best pair in `analysis/ancestry_need/`.
  Earlier text: the proposal says one or
  two clean ones. L5 has 8 in CX333 and L6 has 5.
- **Build OUT-1 on CX333 now,** or with the new graph.
- **Extend `panel_polarity.py` to an ordered outgroup list,** and move the
  ancestral node to the L1-4,7 MRCA.
- **Non-target isolates in cohorts** (L5/6 beyond the in-graph genomes, L8,
  L9, animal lineages): exclude, or flag, at P1 by barcode lineage
  (0i open decision).
