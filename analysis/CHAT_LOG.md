# Chat log

A running record of this working session, so nothing is lost when the terminal
scrolls. Newest entries are at the bottom. Detailed write-ups are in the files
named here, and the formal record is `HANDOFF.md` (sections 0h-0q).

---

## 2026-10-04: GenBank scan for the sparse lineages

**Question:** do complete genomes in GenBank fill the lineages CX333 samples
thinly?

**Answer (`analysis/genbank_scan/README.md`, HANDOFF 0h): barely.**

- **What the scan found:** 173 complete genomes not in CX333.
- **Sparse lineages:**
  - lineage 5: 3; lineage 6: 5; lineage 7: 1; lineages 8 and 9: none;
  - mostly close relatives of CX333 genomes;
  - the one real gap filled is M. orygis (2 clean genomes).
- **Error-rich submissions:** 82 of 173 are flagged for indel or homopolymer
  errors.
- **Inside lineages 1-4:** 43 clean genomes add new diversity.

## 2026-10-04: proposal for a lineage 1-4 panel and outgroup

**Request:** propose a new lineage 1-4 panel and how to get an outgroup for
trees and ancestral reconstruction.

**Answer:** `analysis/strategy/L1_4_PANEL_PROPOSAL.md`, HANDOFF 0i.

- About 300 genomes, chosen by isolate demand.
- The outgroup is held outside the graph.

## 2026-10-04: the three pre-rebuild checks (HANDOFF 0j)

**Correction found by the checks:** lineages 1-4 plus 7 form one clade, and
lineages 5, 6, 9 and the animal lineages are its sister. So **lineage 5/6
genomes are the nearest outgroups**, not lineage 8. The proposal was corrected.

**Two production faults found:**

- **Fault A:** `add_outgroup.py` let the last duplicate record win.
- **Fault B:** the ancestral-allele node included the second canettii.

**Panel coverage:** the public candidates barely close the reference gap; a
demand-chosen panel of about 150 genomes matches CX333's coverage.

## 2026-10-05: outgroup in or out of the graph (HANDOFF 0k)

**Recovering the root without the outgroup in the graph:** align the outgroup
genome to H37Rv directly; this was validated in check 1.

**Cost of an outgroup inside the graph:**

| outgroup added to lineages 1-4,7 | extra SNP sites |
|---|---:|
| one canettii | +34% |
| lineage 5/6 | about 3% |

## 2026-10-05: test graphs (HANDOFF 0n, `analysis/graph_tests/RESULTS.md`)

Four 50-genome builds were compared.

| arm | effect |
|---|---|
| + canettii | +60% nodes and records, deep nesting doubled |
| + lineage 5/6 | +12% |
| sparse mapping (`-x auto`) | the same graph within 1%, 20% faster |

**Recommendation:** lineages 1-4,7 plus a lineage 5/6 outgroup in the graph;
lineage 8 and canettii outside it.

## 2026-10-05: road to 10K (HANDOFF 0l, `analysis/strategy/ROAD_TO_10K.md`)

Run P1 on the 10K isolates first. Its H37Rv alignment does not depend on the
panel, and it gives the demand data for choosing the panel.

## 2026-10-05: full audit, 97 findings (HANDOFF 0m, `analysis/audit/CONSOLIDATED.md`)

**Scale:** seven read-only audits found 10 HIGH, 39 MEDIUM and 48 LOW
findings.

**The current association results are provisional.** The main problems:

- IS6110 carriers were written as non-carriers;
- a read-based canettii isolate is in the cohort trees;
- a third artifact genome (lineage 9, a 377 bp poly-T);
- several SV genotyping errors;
- outgroup and ancestral-node faults.

**User's questions answered:**

- **Where did these come from?** Faults A and B came from checks the user
  approved on 2026-10-04. B came from our own earlier review fix (2026-10-01).
  The rest are lessons for the panel, cost items or decisions.
- **Why so many agents?** It was my choice, for scale and independence. The
  user prefers fewer agents and slower work (now a standing preference).

## 2026-10-05/06: fixes on branch `audit-fixes` (HANDOFF 0o)

**Steps 1-3 are done:**

1. code brought into the repository;
2. section A fixes;
3. section B fixes and two clean-up passes.

**Bugs found along the way:**

- the P4 strand off-by-one at inverted path steps;
- node offsets read at the wrong base.

**44 decisions** are listed in `analysis/audit/DECISIONS.md`, awaiting the
user (D41 especially).

## 2026-10-05: canettii isolate before/after (`analysis/canettii_effect/README.md`)

**Survivors are unchanged with or without the isolate:** 3 in scale200 and
28 in gwas1000, all known DR mutations.

**What the isolate does change:**

- it inflated only the branch null (scale200 137 → 29 passes);
- the ancestral states move only at the MTBC root (242 / 326 variants, mostly
  resolved → unknown).

**Status:** deferred by the user until the costs are clear.

## 2026-10-06: step 4, the second review (HANDOFF 0p/0q, `analysis/audit/review2/CONSOLIDATED.md`)

**Verified:**

- the strand fixes;
- the node offsets;
- Fitch reconstruction, against an independent implementation on 85,575
  records;
- the ancestral alleles, against check 2;
- the association tail on scale200.

**12 items remain before any rerun.** They are being fixed one at a time, as
the user asked, with no new agents.

## 2026-10-06: items 1-6 fixed (all on `audit-fixes`, pushed)

| item | what | commit | effect |
|---|---|---|---|
| 1 | panel polarity keys left-aligned like the cohort VCF | 88b27a9 | indels finding their ancestral allele: scale200 1,553 → 4,452, gwas1000 1,953 → 5,690 |
| 2 | the chain builds the cohort + panel tree (user: "move it into the chain") | 2554f71 | panel rows wrong in 0.002% of cells, against 0.25% in the old combined alignment |
| 3 | no ABSENT where the reference carries the position | 31de8ce | false ABSENT 7/59/7 → 0/0/0 on three references; no real deletion lost |
| 4 | accessory "unmeasurable" judged at 95% identity | eaaa255 | blind loci 674 → 604; three real variable loci recovered |
| 5 | no conditioning on a never-measured accessory locus | 4b5f33a | 167 (scale200) and 862 (gwas1000) variants back in the scan |
| 6 | one same-deletion rule (position, length and half the bases shared) for clustering and both caller matches | 33aff55 | caller deletions no longer wrongly dropped: 98 (scale200), 312 (gwas1000) |

**270 tests pass.**

**User's questions answered on 2026-10-06:**

- **ABSENT vs NOCALL** (in P5 states):

  | state | meaning | how it is written | downstream |
  |---|---|---|---|
  | ABSENT | a positive claim, inferred from the matched reference's sequence, that the H37Rv position is deleted in the sample | GT 2, `*` | unknown in the tree and event matrix |
  | NOCALL | no claim | missing | unknown |

  A false ABSENT asserts a deletion that is not there.
- **Polarised calls:** polarity means knowing which allele is ancestral.
  - **True ancestral reconstructions:** only the `AA` tag (Fitch parsimony on
    the panel tree, at the MTBC ancestor) and the cohort-tree reconstruction
    in the event writer.
  - **Readings or assumptions, not reconstructions:** the panel polarity
    table (one outgroup genome's allele), "presence is derived" for
    insertions, and "ALT is derived" for everything else.
  - **Open question (R2-TREES-6):** 80 variants reconstruct as derived at the
    MTBC root; should they be pinned to `AA`?

**Next:** item 7. The chain's provenance records should include the code,
tree, outgroup and node-locus table, not just the VCF and build.

**Still waiting on the user:**

- the 44 decisions (D41 especially);
- item 11, the IS6110 repeat-node exclusion;
- item 12, the 20,000-permutation floor;
- the canettii isolate.

## 2026-10-06: item 7 fixed (2ac6da3)

**Chain provenance is now a chain of checksums.** Each association-chain
product's `.prov` records:

- the build;
- the merged VCF;
- the code that made it;
- the checksum of every input, upstream products included.

It used to hold only the build and the VCF (plus the polarity table for the
event matrix). A rebuilt tree, new presence tables or a code fix could leave
old products looking current.

**Behaviour now:** a change anywhere upstream makes the chain refuse the stale
product, naming the input that differs.

- **Consequence:** after any code change, the chain's existing products must
  be moved aside (or the cohort run under a new name).
- **Print-only mode:** `MTB_CHAIN_PRINT_PROV=<product>` prints the record a
  product must carry.

**271 tests pass.** Next: item 8, the build manifest must require the
`ancestral` step and hash `build_info.tsv`.

## 2026-10-06: item 8 fixed (b730b24)

**The P0 manifest is now strict.**

- **`ancestral` is required:** a build without its ancestral-allele table can
  no longer be marked complete.
- **No stale completeness marker:**
  - `manifest.done` is cleared at the start of every manifest run;
  - it is also cleared whenever any other build step is marked done
    afterwards (that step changed what the manifest vouches for).
- **`build_info.tsv` is verified:** `p0_check verify` now checks the build's
  identity record against the manifest, as well as `assets/` and
  `annotation/`.

**Tests:** 4 new ones, all failing on the old code. **275 tests pass.**

**Next:** item 9. The I/O contract refuses a fresh build:
`accessory_panel.fasta` is not made by P0, and two outputs are undeclared.

## 2026-10-06: item 9 fixed

**The file-dependency check (`refbias/io_contract.tsv`) no longer refuses a
fresh build.**

- **The accessory panel is a build asset now.** `accessory_panel.fasta` (the
  128 contigs P3 aligns to) is made and indexed by P0's assets step. P3
  requires it and never writes into the build.
  - Before, P3 built it on its first task, changing the build after its
    manifest.
  - On CX333 the P0-made FASTA is identical to production's.
- **Two outputs are now declared:** p1g's IS-clean BAM (which p5svgt reads)
  and p1iv's cohort key table (which p1is reads).

**Tests:** 2 new ones, failing on the old code. **277 tests pass.**

**Next:** item 10. `p3acc` is submitted with no dependency on P1, but it
reads P1's reference map.

## 2026-10-06: item 10 fixed

**The accessory-presence array (p3acc) now waits for P1's summary job.**

- **Why it must:** it reads P1's `refmap.tsv` for each task's sample and for
  the sample's matched reference.
- **What went wrong before:** it was submitted with no dependency. On a fresh
  cohort its tasks would fail and the final VCF pass would never run.
- **Contract:** the refmap is now listed as an input of p3 in the I/O
  contract.

**Tests:** 1 new test (a dry run shows the dependency), failing on the old
code.

**Items 1-10 are all done.** Remaining are two decisions for the user:

- **Item 11:** the IS6110 repeat-node exclusion. It removes about 19% of
  IS6110 sites, about half of them keepable.
- **Item 12:** the 20,000-permutation floor. A lone true hit in the
  small-variant gene burden can reach at best q = 0.177.

## 2026-10-06: user chose option (b) for items 11 and 12; item 11 done (ec2cda7)

**Item 11: IS6110 sites on repeated graph nodes.**

- **New rule:** a repeated node is now keyed when every site on it projects
  onto the H37Rv path within 1 kb of the others and no sample has two sites
  there. Otherwise the whole node stays excluded.
- **scale200:**
  - 37 of 64 repeated nodes are keyed, recovering 156 of the 456 excluded
    sites;
  - every recovered key has all carriers at the same H37Rv position (span 0
    bp), with no sample twice, and up to 29 carriers per key.

**Tests:** 1 new test, failing on the old code. **279 tests pass.**

**Next:** item 12. Adaptive permutations: rows at the 1/20,000 floor are
rerun with more permutations.

## 2026-10-06: item 12 done (4524fac); all 12 review items fixed

**Adaptive permutations in the scan and the gene burdens.**

- **How it works:**
  - every row is still tested with 20,000 permutations;
  - any null whose count is 10 or fewer is re-tested with 1,000,000 fresh
    permutations before BH;
  - a new `refined` column records which nulls were refined.
- **scale200 small-variant burden** (production event matrix; 225 s, 1.4
  GB):
  - gene units passing all three nulls go from 6 to 8;
  - **katG and the embA promoter now pass.** They were the two drug-resistance
    controls that had narrowly missed.
  - The original six (rpoB, rpoC, embB, gyrA, pncA, ethA) pass more strongly
    (q_region 0.022 → 0.0006).
  - The fabG1 promoter still fails the lineage null.

**Tests:** 4 new ones. **283 tests pass.**

**Status:**

- All 12 items from the second review are fixed on `audit-fixes` (pushed).
- Still waiting on the user:
  - the 44 decisions in `analysis/audit/DECISIONS.md` (D41 especially);
  - the canettii isolate;
  - approval to merge `audit-fixes` into main.

## 2026-10-06: items 11 and 12 recorded in HANDOFF (main 3ca5776, pushed)

- **Section 0q now records:**
  - items 11 and 12, with their scale200 results;
  - that all 12 review items are fixed on `audit-fixes` (283 tests).
- **The "Start here" block is updated** to match.
- **Next, waiting on the user:**
  - the 44 decisions (D41 especially);
  - the canettii isolate;
  - R2-TREES-6;
  - then the merge of `audit-fixes` into main.

## 2026-10-06: D41 explained

- **Full write-up:** `analysis/audit/D41_EXPLAINED.md`.
- **The problem:** node-frame variants (those on sequence not on the H37Rv
  path, keyed `node:<id>:<offset>`) are written in the direction of the
  reference's walk. When two references walk a node in opposite directions,
  one event gets two keys, for example C>T and G>A.
- **How often:** 11 bases in scale200.
- **Recommended:** option (b). Write the alleles on the node's forward strand
  (P4 complement, indel re-anchor, matching flip in P5), done as one item
  before the merge.
- **Awaiting the user's choice.**

## 2026-10-06: should an ancestral composite pick the strand? (D41 follow-up)

**The user's question:** could an ancestral sequence from the graph choose
each node's strand, rather than one rule for all?

**Answer** (details in `analysis/audit/D41_EXPLAINED.md`, follow-up
section):

- **The strand is only a label.** Any rule that is the same for all samples
  gives the same trees and reconstruction. The node's forward strand is
  already fixed per graph.
- **What an ancestral sequence really answers is polarity,** a separate
  question. Here there is a real gap: `ancestral_alleles.py` writes `AA`
  only for H37Rv-frame panel SNPs, and node-frame variants get none.
- **Proposal:**
  - keep D41 (b);
  - add a later item: a node-frame `AA` from the panel paths and tree, with
    the ancestral orientation as an annotation, built with the new graph.

## 2026-10-06: D41 (b) done (audit-fixes f835497); AA item added (main)

**The user's request:** "yes, go with (b) and add the AA item".

**D41: node-frame alleles are now written on the node's forward strand.**

- **P4:**
  - complements the alleles when the reference reads the node reverse
    complemented (odgi's walk flag combined with the 22 genomes stored
    flipped in the panel);
  - re-anchors indels on the forward left.
- **P5:** the reversion test complements to match.
- **scale200 check** (P4 rerun locally from production inputs, 199 samples,
  about 7 minutes, output in the scratchpad only):
  - node-frame records 22,709 before and after, so nothing is lost;
  - duplicate reverse-complement key pairs 3 to 0;
  - **SNP keys whose REF is the node's own base: 50% (10,932) to 100%
    (22,008 of 22,008).** Before, alleles were on H37Rv's strand.
  - **Open edge:** 244 records (1.1%) are reverse-read indels or MNPs whose
    forward start falls on the neighbouring node. 64% of the graph's nodes
    are 1 bp. They are left in the reference's orientation and counted.
- **Tests:** 10 new; 293 pass.

**AA item:** added to HANDOFF 0q (main e49c3e4, pushed). It covers ancestral
alleles for node-frame variants from the panel paths and tree, with the
ancestral orientation as an annotation, and is built with the new graph.

**Question to the user:** close the 244-record edge now? P4 can find the
neighbouring node from the build's node table (no new odgi runs) and key the
event there, as a forward-reading reference already does.

## 2026-10-06: the 244 node-boundary records (audit-fixes 746e57d; main a86d315)

**The user's request:** "yes, close the 244 now".

**What changed:** P4 finds the neighbouring node from the build's node table
(no new odgi pass) and keys the event there.

**scale200 check** (P4 rerun locally, using a full node table built from the
CX333 graph in 40 s; the old CX333 table only listed a subset of nodes):

- **96 of the 244 now get a single key.** Multi-base keys whose REF matches
  the node's own sequence go from 457 to 553 of 701.
- **3 keys now merge** across references.
- **3 false merges are now split.** For example, +C read on one strand and
  +G read on the other had both been written as A>AC.
- **148 remain:**
  - **138 are anchored on H37Rv's own path.** A forward reader keys these in
    H37Rv coordinates (or the direct arm owns them in core), while a reverse
    reader keys them to a node. This is a question of which coordinate frame
    the event belongs in, not of node naming. **It needs the user's
    decision.**
  - 5 are at inversion junctions;
  - 5 are on nodes visited twice.
- **Also noted:** node-frame indels are not left-aligned on the forward
  strand, so a homopolymer indel read in opposite directions can still get
  two keys.

**Tests:** 3 new; 296 pass.

---

## 2026-10-06: D41 option (b) done; homopolymer measured (audit-fixes 0100c0c; main 296e42d)

**What you asked:** go with (b) and measure the homopolymer issue.

**(b) is committed (0100c0c on `audit-fixes`, 299 tests pass).** Validated by
rerunning P4 on all 199 scale200 samples:

- **22 records changed:** 15 now get `h37rv:<t>` keys and 7 in core are left
  to the direct arm.
- **9 of the 15 match** (after left-normalisation, as P5 does) keys that
  other samples wrote for the same event. The other 6 are events no other
  sample carries.

**Correction:** the earlier "138" counted the wrong anchor. Only 22 have
their H37Rv-strand anchor on the path.

**Still open, about 116 records: your call.** Here R reads along H37Rv's
strand, so H37Rv also reads the node in reverse, and R's own anchor is off
the path.

- Keying these "as a forward reader does" means using `on_strand`'s
  t-len(D)-1 position.
- That matches other samples' keys for short events (TGG>T at 976,895).
- For long off-path deletions it lands about 100 bp away from where others
  key them.
- For now they stay as node keys.

**Existing, not caused by (b):** 24% of composed-arm H37Rv-frame deletions
(407 of 1,702) have a REF that differs from H37Rv, because R's own sequence
differs there.

**Homopolymer:**

- 410 node-frame indel keys, 331 checkable;
- 42 are not left-aligned on the node (59 sample records);
- only 3 events actually get two keys, each from two different references;
- 31 more would shift past the node start.

A fix is optional: left-align node-frame indels against the node sequence.

---

## 2026-10-07: along-strand group keyed; node indels left-aligned (audit-fixes 16715bb; main e79c93e)

**What you asked:** go with the recommendation, and left-align the node
indels.

**Committed:** 16715bb on `audit-fixes`, 308 tests pass. Checked by rerunning
P4 on all 199 scale200 samples.

**The along-strand group:**

- These events are keyed in H37Rv coordinates only where the key's REF is
  H37Rv's own sequence. That is always true for an insertion. For a
  deletion it means the deleted bases are H37Rv's.
- Long deletions of off-path sequence stay as node keys (35 records).
- Totals for D41 (b), both groups:
  - 73 records keyed in H37Rv;
  - 30 in core left to the direct arm;
  - 45 still keyed as R reads them (inversion junctions, nodes visited
    twice, anchors themselves off the path).
- Most of these events occur in no other scale200 sample, so there is
  little to cross-check against. Where there is a match, it agrees.

**Left-alignment:**

- Node indels are now shifted left along the node's own sequence, stopping
  at the node's first base.
- 59 records moved.
- Checked independently against the graph's node sequences:
  - keys not left-aligned: 42 to 0;
  - events with two keys: 3 to 0.

---

## 2026-10-07: is the chat log still being updated?

Yes. Every substantive turn gets an entry here. The two headings above were
brought into the same format as the earlier ones, with commit references.
It is now committed on main (the user asked, 2026-10-07).
