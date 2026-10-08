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

---

## 2026-10-07: chat log and audit folder committed to main

- Chat log: b6bf039.
- `analysis/audit/` (19 Markdown files, including DECISIONS.md and
  D41_EXPLAINED.md): committed next.
- The rest of `analysis/` stays uncommitted.

---

## 2026-10-07: going through DECISIONS.md, group 1 (D1-D7)

**You confirmed all seven as implemented.** Each one writes NOCALL rather
than guessing:

- D1: no clean junction;
- D2: no projection;
- D3: the reference carries a third allele;
- D4: a reference-only deletion of 50 bp or more;
- D5: depth plus a clip cluster;
- D6: no record for accessory loci H37Rv carries;
- D7: duplicate insertions are dropped.

They are recorded in a new answers table at the top of DECISIONS.md. Next is
group 2, IS6110 (D8-D10).

---

## 2026-10-07: group 2, IS6110 (D8-D10)

**You confirmed all three as implemented.**

- **D8:** insertions on repeated nodes are not keyed.
- **D9:** at sites where H37Rv has a copy, non-carriers are NOCALL.
- **D10:** conflicting duplicate rows stop with an error.

**Tracked for the new graph build:** keying D8's insertions by the nearest
single-copy node. It is listed in HANDOFF 0q beside the AA item.

Next is group 3, ancestral states (D11-D14).

---

## 2026-10-07: group 3, ancestral states (D11-D14)

**You confirmed all four as implemented.**

- **D11:** an uncalled outgroup is N.
- **D12:** `*` and sibling-allele samples count as missing.
- **D13:** panel allele frequency comes from genotypes.
- **D14:** the H37Rv tip stays REF at node-frame IS6110 records. H37Rv
  lacks the site, so it carries no element there.

Next is group 4, association tests (D15-D19).

---

## 2026-10-07: group 4, association tests (D15-D19)

**Confirmed:**

- **D15:** the all-three-nulls survivor rule;
- **D16:** any overlap credits a gene for a deletion;
- **D17:** burdens use the 80% callability floor;
- **D19:** each cohort is tiered from its own genotypes.

**Switched, D18:** a point variant in overlapping genes is credited to every
gene, not just the earlier-starting one. That matches D40 and SV deletions.

The code change is queued with any other switched decisions, to be made one
at a time after the review.

Next is group 5, reference selection and the graph VCF (D20-D24).

---

## 2026-10-07: group 5, reference selection and the graph VCF (D20-D24)

**Confirmed, D23:** `0` wins over `.` in duplicate records.

**Switched, to implement before the new-panel rerun:**

- **D20 + D24:** pick the matched reference by mismatches per site that both
  sides called. Isolate coverage comes from the P1 BAMs.
- **D21:** H37Rv becomes a candidate reference, tested through P2-P5.
- **D22:** left-align the collapsed graph VCF, then re-collapse the 41 events
  that split.

**Queue so far:** D18, D20/D24, D21, D22.

Next is group 6, the graph build (D25-D27).

---

## 2026-10-07: group 6, graph build (D25-D27)

**You confirmed all three as implemented.**

- **D25:** sparse mapping is off by default. Revisit only if the new
  panel's cost estimate is too high.
- **D26:** keep the decomposed graph VCF.
- **D27:** `--rank-by` is required.

Next is group 7, panel construction (D28-D33).

---

## 2026-10-07: group 7, panel construction (D28-D33)

**Confirmed:** D28, D30, D31 and D33.

**D29 confirmed, with an addition:** a UniVec vector check. You approved the
download, which I'll do when implementing it.

**Switched, D32:** the foreign-screen background is chosen by quality, using
the same ranking rule as D27.

**Queue:** D18, D20/D24, D21, D22, the UniVec check, D32.

Next is group 8, build safety (D34-D38).

---

## 2026-10-07: group 8, build safety (D34-D38)

**Confirmed:** D34-D37 (panel_af and add_outgroup from the build's collapsed
VCF, the off-path node table, deposited references).

**Switched, D38:** refuse all pre-guard outputs. The rerun regenerates
everything.

**Queue:** D18, D20/D24, D21, D22, the UniVec check, D32, D38.

Next is group 9, leftovers (D39-D40).

---

## 2026-10-07: group 9, leftovers (D39-D40)

**Switched, D39:** complex and MNP records are credited to every gene their
REF span touches. The open item about indels credited by their anchor base is
folded in: indels follow the span rule too. It will be done together with
D18.

**Confirmed, D40:** small deletions are credited to every overlapping gene.
This is now the same rule as D18.

**Queue:** D18 with D39 and the indel item, D20/D24, D21, D22, the UniVec
check, D32, D38.

Next is the last group, D42-D44.

---

## 2026-10-07: group 10 (D42-D44); the decisions review is finished

**Confirmed:**

- **D42:** drop and count records with no node length (0 in scale200);
- **D43:** catalogue inputs are used only when given;
- **D44:** canettii is the default outgroup, and it is fatal if they are not
  tree leaves.

**Queued, the D43 fix:** a path that is given but missing becomes fatal.

**Closed:** the "found during the fixes" section. Five items are already
fixed and the burden-indel item is folded into D39.

**All 44 decisions are answered.** Queue of 8 fixes, each one at a time with
a test:

1. D18 + D39 + indels (the span rule);
2. D20 + D24;
3. D21;
4. D22;
5. UniVec;
6. D32;
7. D38;
8. D43.

Still with you: ASSOC-1 (the canettii isolate) and R2-TREES-6.

---

## 2026-10-07: fix 1 done, D18 + D39 + burden indels (audit-fixes 27618e3)

**One rule for every variant: the bases it changes decide its genes.** It is
used by the gene burdens and by the scan's genic/intergenic stratum, from one
shared helper (`assoc/bin/variant_span.py`).

- **SNPs in overlapping genes** count for every gene, not just the
  earlier-starting one.
- **MNPs and complex records** count for every gene their changed bases
  touch.
- **Insertions** (small, SV and IS6110) count for a gene only if they fall
  inside it. One just past a gene's end goes to the neighbouring promoter.

**Effect on the current events** (callability floor not applied in this
count):

- scale200: 271 records change unit, 260 of them gaining an overlapping
  gene;
- gwas1000: 584, including 9 IS6110 records;
- SVs: none.

**Tests:** 6 new and 1 updated (the old test encoded the anchor rule); 314
pass. Two mutations of the new rule are both caught.

**Not changed:** P6's annotation still names one gene per position. It is
descriptive output, and changing its column format would affect everything
downstream of it.

Next is fix 2, D20 + D24 (reference selection).

## Fix 2 done; how the fixes change the results (2026-10-07)

**User:** "will you compare the previous runs with the fixed runs to show
how the fixes change the results? I haven't seen any results from these
fixes yet."

**Answer:** each fix has so far been checked against the production
outputs at the level it changes. The end results (P5 matrix, trees,
association survivors, DR positive controls) can only be compared after the
post-fixes scale200 rerun, which still needs a cost estimate and approval.
Proposed: (1) now, a before/after document collecting every fix's measured
effect; (2) after the rerun, an old-vs-new comparison at every stage.

**Fix 2 (D20 + D24), commit a30bcd4:** the new selector was run on every
isolate of both cohorts (outputs in the scratchpad) and compared with the
production choices:

| cohort | isolates | rank-1 changes | new pick was old rank 2 | old choice trails by ≤2 SNPs | ≤5 | >10 |
|---|---|---|---|---|---|---|
| scale200 | 200 | 24 (12%) | 18 | 13 | 21 | 1 |
| gwas1000 | 997 | 112 (11%) | 87 | 59 | 91 | 5 |

Panel sites covered per isolate: median about 74,100 of 75,587 (minimum
about 71,600). Most changes are near-ties between closely related
references, but a few move by more than 10 SNPs. Tests: 324 pass.

**Before/after document:** `analysis/audit/FIX_EFFECTS.md` collects the
measured stage-level effect of D41, fix 1 and fix 2. Part 2 (end results)
will be added after the scale200 rerun. Next: fix 3 (D21).

## 2026-10-07: fixes 3 and 4 done (audit-fixes 1c6b365, cf46f41)

**User:** "yes, write the document, then start fix 3"; "how many isolates
pick H37Rv now?"; "commit fix 3 when gwas1000 finishes, then start fix 4".

**Fix 3 (D21), H37Rv as a candidate reference:**

- scale200: 1 isolate (SAMEA7526648) now picks H37Rv, at 29 vs 36 SNPs.
- gwas1000: none; no other choice changes.
- **Found on the way:** a path projected onto itself through the graph is
  not the identity in repeats. With R = H37Rv, 729 of 51,139 P4 records
  moved, by up to 1 kb; in the P5 direction, 16 of 52,007. In the current
  IS6110 arm (any R), 22 of 2,707 same-reference lookups landed on another
  copy.
- Fixed once in `frame_convert.py`, plus the two IS6110 scripts that read
  odgi directly.

**Fix 4 (D22), graph VCF left-aligned before the collapse:**

- 8,304 of 93,214 keys move; 49 records fold into 41 keys; no carrier
  lost.
- **Also found:** the production graph VCF predates the earlier GRAPHVCF-5
  collapse fix (1,461 padded keys, 61 missing carrier cells). The new build
  regenerates it.

**Tests:** 344 pass. Both fixes are recorded in
`analysis/audit/FIX_EFFECTS.md`. Next: fix 5 (UniVec).

## 2026-10-07: fix 5 done, UniVec vector check (audit-fixes 6dec4b2)

**User:** "yes, start on fix 5".

- **Download:** UniVec_Core build 10.0 (3,155 sequences), fetched into the
  repository's `data/univec/` (git-ignored) by the new
  `bin/fetch_univec.sh`, which records the build and checksum.
- **The check:** every foreign-screen insert, native or not, is searched
  with VecScreen's blastn settings. A strong match gives `VECTOR`.
- **Real data:** of 2,703 inserts (CX333 449, external assemblies 2,254),
  exactly the two known constructs match: pJEB in GCF_044324775 and the attB
  vector in GCF_021535155. Nothing else matches, not even moderately.
- **Also fixed:** a failed minimap2 call was read as "no inserts" and marked
  done.

**Tests:** 349 pass. Next: fix 6 (D32, the background chosen by quality).

## 2026-10-07: fix 6 done, the foreign screen's background (audit-fixes cba2135)

**User:** "yes, start on fix 6".

- Each sublineage's background genome is now the best by `--rank-by`, the
  same file and rule `snp_nonredundant.py` uses (D27). Before, it was the
  first accession in sorted order.
- A competing genome with no score stops the run.
- In CX333, 44 of the 65 slots have a choice. The score file itself comes
  with the new panel's clone-collapse step, so the actual picks are decided
  then.

**Tests:** 355 pass. Next: fix 7 (D38, refuse outputs from before the
guards).

## 2026-10-07: fix 7 done, no adoption of old outputs (audit-fixes 6afb220)

**User:** "yes, start on fix 7".

- **Four reuse paths removed:** P1 (H37Rv VCF stamp), P2 (VCF stamp plus
  GATK reference), the accessory presence tables (newer than the catalogue)
  and the chain runner (a folder adopted on 20 VCF stamps; the presence
  folder on no evidence).
- **Now:** a sample is done only by its marker; anything else is made
  again. A folder without a build record is refused.
- **Real data:** all 8 production cohort folders lack a build record, so
  the scale200 rerun will need a new registry row with new output folders.

**Tests:** 360 pass (4 old tests rewritten, 6 added). Next: fix 8 (D43).

## 2026-10-07: fix 8 done; all 8 queued fixes complete (audit-fixes 1b359b6)

**User:** "yes, start on fix 8".

- **D43:** a catalogue input that is named but missing now stops the run,
  in `merge_catalogues.py` and in P0 (whose input record wrote a missing
  path as "absent", the same as none).
- **Real data:** with CX333's real insgt tables the catalogue is
  byte-identical. With a typo, the old code exited 0 and dropped the insgt
  cluster from 214 of 802 loci.

**Tests:** 364 pass.

**All 8 queued fixes are done.** Each one's measured effect is in
`analysis/audit/FIX_EFFECTS.md`.

**Still waiting on the user:**

- the canettii isolate (ASSOC-1);
- R2-TREES-6;
- the merge of `audit-fixes` into main.

**Then:** the new panel and graph build, then the scale200 rerun into new
output folders, with a cost estimate first.

## 2026-10-07: MTBC node pinned; ASSOC-1 deferred; audit-fixes merged (main 3f0260a)

**User:** "pin the MTBC node, keep ASSOC-1 deferred, then merge".

**R2-TREES-6, the MTBC node pinned to the panel AA** (audit-fixes a4c01e1):

- The cohort tree's MTBC node (the common ancestor of every leaf but the
  outgroup and the canettii genomes) is pinned to AA wherever AA resolves.
  `mtbc_pinned.tsv` lists every change.
- **Checked on review 2's scale200 run:**
  - 95 variants change: the 80 review 2 found, plus 15 that parsimony
    left tied. No other variant changes.
  - Testable variants among them go from 28 to 87.
  - The scan's survivors are identical (rpoB 761155, embB 4247429, embB
    4247730); no drug-resistance control changes status.
- **Along the way:**
  - My first validation command pointed at the wrong polarity table and
    stopped with a FATAL; nothing was written.
  - The rerun was blocked for an unguarded `rm -rf`; it was rewritten to
    use fresh folders.
  - A bookkeeping bug (a NumPy view) made the pin report 0 changes; it was
    fixed and is covered by a test.

**ASSOC-1** stays deferred.

**Merge:** `audit-fixes` was merged into main (3f0260a). There were no
conflicts, and 369 tests pass on main.

**Next:** the new lineage 1-4 panel and graph build, then the scale200
rerun into new output folders, with a cost estimate first.

## 2026-10-07: can canettii stay out of the graph? The outgroup test

**User:** "I need to be convinced that we can do this without canettii in the
graph as this will be the outgroup in the tree." Then "yes, run the test" and
"write up the results when the trees finish".

**Full write-up:** `analysis/outgroup_test/README.md`.

**The test:** two arms that differ only in where the two canettii genomes'
rows come from, the graph (G) or each assembly aligned straight to H37Rv
(K). Columns where only canettii varies are dropped in K.

**Results:**

- **Canettii's rows:** ET1291 is identical at 92,468 of 94,383 columns.
  Only 317 (0.34%) have the opposite allele; the rest differ only in which
  route leaves a site uncalled.
- **Trees** (scale200 + CX333, production IQ-TREE settings):
  - the same root (MTBC clade at 100% support with the canettii outside);
  - every lineage monophyletic;
  - lineages 1-4,7 one clade, with lineage 1 the first split.
  - Robinson-Foulds distance 12 of 1,062 splits (1.1%). All 12 are inside
    lineage 4, with support 16-79.
- **Ancestral alleles:** 49,168 of 49,179 identical (99.98%).
- **Association (scale200):**
  - the same 3 survivors (rpoB 761155, embB 4247429, embB 4247730);
  - every drug-resistance control unchanged;
  - testable variants 3,679 → 3,682.

**Not covered yet:**

- IS6110, SV and accessory polarity, which needs a simulated canettii
  pseudo-isolate;
- whether a graph without canettii calls the other genomes differently
  (test graph arms A vs C, a cheap local comparison);
- gwas1000.

**Compute:** 15.7 CPU-hours, above the roughly 11 estimated.

## 2026-10-07 (evening): fixes first; graph work on hold

- **User:** run the test-graph comparison; are there scale200 results with
  CX333 before and after the fixes? "We don't need to do this with a new
  graph now."
- **Found:** netscratch was purged at 14:24 (files older than about 90
  days): the pggb and vg containers, the H37Rv FASTAs, known_RDs beds and a
  few data tables. All are on the mirror. Restoring them is the user's call.
  The CRAMs (dated 2026-07-28) face the same purge around 2026-10-26.
- **Answer:** no full fixed-code scale200 run exists yet. Rerun estimate
  about 500 billing-hours (about 350 right-sized; about 150-200 reusing P1).
- **User:** "How can I evaluate the fixes? ... I have not seen any output
  yet." Published a page with what is measured:
  https://claude.ai/artifact/S5zrD6Qzr5viBctJ3hGeRH. The association scan
  with old vs fixed association code: the same 3 survivors and every DR
  control in place; each earlier-step fix measured on its own.
- **Test graphs A vs C:** canettii in the graph changes about 2-4% of the
  other genomes' SNP calls, mostly near where canettii varies; 21 SNP
  positions change base.
- **User:** "Please hold off on the new graph comparisons until I get some
  data to evaluate the fixes. write this to handoff and we will get back to
  it." Written to HANDOFF "Start here (2026-10-07, evening)".
- **User:** "yes, restore the purged files from the mirror." Done: 13 files
  copied back, each identical to the mirror, with fresh dates; every
  configured path resolves.
- **User:** "How are the jobs sized now? what does sized mean?" Explained:
  the cores and memory each job reserves, which is what fairshare charges.
  P1 and P2 reserved 8 cores / 16 GB and used about 2 cores / 1.4 GB.
- **User:** resize P1 and P2 to 4 cores, 8 GB; "use slightly more detail in
  your descriptions". Checking the tools found three needed companions
  (bwa's thread-dependent batches, P1 not passing its core count, GATK's
  8 GB Java heap). **User:** "make all three and run the test."
- **Test:** scale200's largest sample at 4 cores / 8 GB gives alignments and
  calls identical to production's September P1 output. Peak memory 2.2 GB,
  20% longer, about 40% cheaper. Rerun estimate now about 380 billing-hours.
  One test assumed the purged H37Rv file was absent; made independent of it.
  369 tests pass.

## 2026-10-08

- **User:** "yes, run the scale200 rerun at 380 billing-hours."
- **Setup** (`analysis/rerun_scale200/README.md`): a new build,
  `7713a8d71d8e-fix1`, on the same CX333 graph, and a new registry row,
  `scale200_fix`, with new output folders. The old run is untouched for the
  comparison. The registry is now a repository file that `runroot` links to,
  as `io_contract.tsv` already was. `vcf_split_classes.sh` gained `--outdir`,
  so the build's files do not replace those in the graph folder.
- **Graph VCF with the fixed collapse:** 79,474 SNP/MNP, 9,027 indel and
  4,652 SV records; no duplicate keys; indels left-aligned.
- **Accessory candidates:** the same 1,643 insertions, but D22 moves most of
  them 1-14 bp left, so accessory locus ids change between old and new.
