# D41: which strand node-frame alleles are written on (2026-10-06)

## What a node-frame variant is

- Most variants are keyed in H37Rv coordinates.
- Some variants sit on graph sequence that H37Rv's path does not pass
  through. Examples are lineage-specific insertions and regions deleted in
  H37Rv.
- These have no H37Rv position, so they are keyed to the graph node instead,
  as `node:<id>:<offset>` with an allele such as `C>T`. This is the "node
  frame".

## The problem

- **A graph node has a forward strand.** A reference genome's path can pass
  through the node forwards (+) or backwards (-).
- **Offsets are already fixed** (clean-up pass, fb199f4). Offsets are
  measured on the node's forward strand, so one base always gets one
  `node:<id>:<offset>`.
- **Alleles are not fixed.** They are still written in whichever direction
  the sample's matched reference walks the node.
- **So one event can get two keys.** One reference walks the node +, another
  walks it -, and the same mutation is written two ways:
  - `node:200110:0 C>T`, from references walking +;
  - `node:200110:0 G>A`, the reverse complement, from references walking -.
- These do not match, so the cohort VCF gets two records for one biological
  event.

## How often it happens

- **scale200:** 11 bases have two keys this way, out of about 11.5K
  node-frame small records. It is rare.
- **gwas1000:** not counted separately. There are more node-frame records,
  and more references in a bigger graph, so the count is expected to grow
  with the cohort and the graph.

## Why it matters for trees and ancestral reconstruction

- **The event matrix sees two characters for one mutation.** Depending on
  how P5 calls samples at each key, the effect is one of two things:
  - **the carriers are split** between the two columns. Each column looks
    like a smaller, possibly homoplastic event, the reconstruction places it
    on the wrong branches, and the power of an association test is diluted;
  - **or the event is counted twice,** and the two rows are perfectly
    correlated. That inflates the number of tests and double-weights the
    event.
- The current numbers do not say which of the two happens. Fixing D41
  removes both, so it does not matter for the decision.
- **A lineage-specific event is exactly the kind that lives in node frame,**
  and references from different lineages are the ones likely to walk a node
  in opposite directions. The error therefore falls on events that matter
  for lineage-level reconstruction.

## The choice

- **(a) As implemented: offsets fixed only.**
  - Nothing more to do.
  - The duplicate keys stay. They are few now, but grow with scale.
- **(b) Recommended: write node-frame alleles on the node's forward strand.**
  The change has three parts:
  1. **P4 (`bin/p4_place.py`).** When the reference walks the node `-`,
     reverse-complement REF and ALT. `odgi position -v` already reports the
     strand.
  2. **Indels.** Re-anchor them after the complement. The VCF padding base
     must be the base to the left on the forward strand, so the anchor
     moves.
  3. **P5 (`bin/p5_states.py`).** Make the same flip when genotyping a
     sample at a node key through a reference that walks it `-`, so P5
     reads the allele on the same strand the key is written in.
- **Cost of (b):**
  - one code change, with tests: a node walked +, a node walked -, a SNP and
    an indel each way, and a check that the two walks give one key;
  - no extra compute, because the rerun is happening anyway.
- **Risk of (b):** getting the indel re-anchoring wrong. The tests above are
  aimed at that, and real data can be checked by confirming the 11 duplicate
  bases collapse to one key each.

## Recommendation

**(b), done before the merge, as one item in the same sequential style as
items 1-12.** All data will be regenerated anyway. Trees and ancestral
reconstructions are the priority, and this defect splits or doubles events
along lineage lines.

## Follow-up: use an ancestral sequence to choose the strand?

The user asked: could an ancestral composite sequence built from the graph
choose the strand, instead of one rule for every node?

### The strand is only a label

- **The key just has to be the same for every sample.**
  - Complementing every allele of an event, in every sample, changes nothing
    in the tree or the reconstruction. Fitch gives the same result whether a
    site is written C/T or G/A.
  - The only thing that breaks trees is inconsistency: one event written
    both ways.
- **The node's forward strand is consistent by construction.** It is a fixed
  property of the graph, so every reference and sample gets the same key.
- **An ancestral orientation would also be consistent.** But it adds nothing
  to the tree, and it costs a new method (below).

### Two different questions

- **(1) Which strand the key is written on.** This is naming, and D41 is
  about it. Any fixed rule works.
- **(2) Which allele is ancestral** (polarity). This is biology, and it is
  what an ancestral sequence is really for.
- The two are independent. A node read in either direction still needs its
  ancestral allele determined separately.

### Where the suggestion points at a real gap

- **`bin/ancestral_alleles.py` makes `AA` only at panel SNP sites in H37Rv
  coordinates.** Node-frame variants get no `AA`.
- Their polarity comes only from the cohort tree's Fitch, with the outgroup
  rooting (and D14: the H37Rv tip is unknown there).
- **Extending the panel-level ancestral reconstruction to node-frame
  variants would be a real improvement for ancestral reconstruction.**
  - It would use the 333 panel genomes' paths through each node, on the panel
    tree, to give an `AA` at node-frame sites as at H37Rv sites.
  - That is question (2). It can be done on top of either strand rule.
- **If it is built, it would also be natural to report each node's ancestral
  orientation** as an annotation, for example "the MTBC ancestor reads this
  node -". That keeps the key stable and lets the output be shown in
  ancestral orientation.

### Why not key on an ancestral composite sequence directly

- **No ancestral genome path exists yet.** Building one is a new method:
  - choose the ancestral state of every presence or absence and every
    inversion;
  - stitch it into one path;
  - handle segments the ancestor lacked, such as later insertions. Those
    would have no ancestral orientation, so a fallback rule would be needed
    anyway.
- **Inside inversions,** the ancestor reads the node one way and the inverted
  lineages the other. Either label is fine, as long as it is fixed.

## Revised recommendation

1. **D41 (b):** the node's forward strand as the key. It is consistent,
   cheap, and needs no new method.
2. **Add as a separate, later item: an `AA` for node-frame variants** from
   the panel paths and tree, with the ancestral orientation as an
   annotation.
   - It belongs with the new panel and graph build, since node ids are
     graph-specific.
