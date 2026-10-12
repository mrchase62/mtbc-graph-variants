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
- **Build `7713a8d71d8e-fix1` complete:** references and frame table
  identical to the old build's; IS6110 crossmaps reproduce (0 of 1,000 files
  differ); accessory panel 806 loci (old 802); panel tree rebuilt from the
  fixed VCF; 73,256 ancestral sites.
- **Chain submitted** (p1 to p5vcf, jobs 51285817-51285924). The association
  tail and the comparison follow when it finishes.
- **Chain finished 11:10** (started 06:30), every task COMPLETED, no
  failures or requeues. Measured cost of the chain by sacct billing weights:
  **577 billing-hours**, over the approved 380. P5 states alone took 362
  (median 49 min per sample, against 15 for both P5 submissions of the old
  scale200): the new build id starts with an empty projection store, so 143
  of the 200 tasks projected their reference's positions through odgi before
  storing them. The estimate did not include this. P1 61, P2 89, P1g+P1i 39,
  P3 13, P4b 9, P4 3, rest under 1 each. Association tail waiting for the
  user's cost decision.
- **ctpV check** (user asked whether the CX333 graph has the L1.2.1 ctpV
  deletion of Nat Commun 2025, s41467-025-65779-9): yes. A 297-bp DEL at
  1,078,519-815 over the ctpV start codon in all 5 L1.2.1.2 panel genomes
  (the paper's confirmed deletion is 297 bp); the 9 L1.2.1.2.1 genomes carry
  longer deletions over the start (235/239 bp plus 872-1,296 bp downstream,
  or one of 1,536 bp); no other lineage-1 genome. `analysis/ctpV_check.md`.
- **ctpV in scale200_fix:** the three L1.2.1.2.1 isolates are correctly ALT
  for the deletion over the start codon (depth 0-2x there). SAMEA2297133
  (basal L1.2.1, matched to GCF_040208995) has an intact ctpV at 72-107x but
  is ABSENT on two intervals nested in its reference's 1,296-bp deletion, and
  on two small records. Inherited absence is not checked against H37Rv-frame
  depth for nested intervals and small states. Reported, not fixed; details
  in `analysis/ctpV_check.md`.
- **Inherited ABSENT checked cohort-wide** (`analysis/inherited_absent/`).
  Of the 85,465 ABSENT cells in scale200_fix that have an H37Rv span, 83%
  are supported by the sample's own depth, 12% (10,247; 3,197 sites; 199 of
  200 samples) are covered at normal depth, and 5% are partial. They cluster
  in phiRv1, plcA/B, PPE57/58, Rv3766-70 and wag22, and rise with distance to
  the matched reference. 12% is an upper bound: relocated sequence also gives
  depth, so the junction clips and split reads decide. 2.2M node-frame
  ABSENT cells are not testable this way.
- **ctpV read evidence** added to `analysis/ctpV_check.md`: HaplotypeCaller,
  clip/split-read clusters and delly/dysgu on the matched reference. delly
  calls nothing over ctpV in the four L1.2.1 samples. dysgu finds PASS
  fragments of the carriers' differences from their references. Neither sees
  SAMEA2297133's ~1.5 kb of extra sequence relative to GCF_040208995: an
  insertion longer than a read.
- **Junction check of contradicted ABSENT cells** (`analysis/inherited_absent/`).
  73% (7,461 of 10,247) are in repeat-masked sequence (PE/PPE, paralog,
  tandem, IS) and cannot be called either way; the user asked that these be
  distinguished. In core sequence: 2,279 cells (110 samples, 94 R deletions)
  are confidently wrong (bridging reads at both ends, no junction reads); 71
  are correctly ABSENT (junction reads, depth from a copy elsewhere); 435 are
  unresolved. 76% of all testable ABSENT cells are masked.
- **Local assembly of the 109 undecided core regions:** 77 present (161
  ABSENT calls wrong), 3 deleted, 6 complex, 23 unresolved. The leftover
  sequence is IS6110 copies and canettii divergence; none is accessory loss.
  The user lost track of the method because of loose terms ("realign", "local
  reassembly", "accessory"), and I had dropped the proposed read-realignment
  step without saying so. Now explained in standard terms: SPAdes local de
  Bruijn assembly, then contigs aligned with minimap2. Saved as a working rule.
- **ABS-1 added to the fix list** (HANDOFF). The user chose NOCALL with a
  reason for masked sites. Open: whether the assembly step joins the
  standard chain.
- ABS-1: the user decided local assembly stays an audit run on request.
- **ABS-1 diagnosis step specified:**
  - evidence: depth, crossing reads and joining reads at R's deletion ends;
  - calls: ABSENT confirmed / no reads, own REF/ALT call when contradicted,
    NOCALL masked / mixed / unresolved / R has no genotype.
- **Novel-insertion extension** (genome-wide two-sided clip scan, then blastn
  sorting): the user asked to keep it separate from ABS-1, on the
  insertion-gap plan (HANDOFF 0e). Not scheduled.
- **Novel events raised in priority** by the user ("external SV callers are
  not working very well ... we need to be able to capture novel events").
  The plan is `analysis/novel_events/PLAN.md`: a breakpoint-evidence caller on
  the P2 BAMs vs the matched reference, sharing a module with ABS-1. Phase A
  uses simulated reads from about 40 external assemblies already on disk,
  about 30 billing-hours. Waiting for approval.
- **Novel events Phase A started** (user approved 30 billing-hours).
  - 40 external assemblies, simulated with wgsim at 80x, run through P1/P2 as
    test cohort `novelA40` (registry row; outputs `refbias/novelA40/run`).
  - Truth set: 993 events of 50 bp or more, assembly vs matched reference
    (minimap2 asm5 + paftools).
  - Prototype `breakpoint_caller.py` and `score.py` written.
  - Spend 17.5 billing-hours with P2 running; expected about 35-40, over the
    approved 30. Reported to the user.
- **Phase A first results** (`analysis/novel_events/README.md`, 24.5
  billing-hours).
  - On simulated reads: dysgu PASS typed recall 0.51 / precision 0.65; delly
    0.26 / 0.49; prototype v2 0.33 / 0.62.
  - Prototype + dysgu finds 0.65 of breakpoints.
  - Misses: 40% repeats (MAPQ), about a third tandem copy-number changes that
    leave no clipped reads (need a depth scan).
  - Next options: local assembly of candidates, a depth scan, Phase B real
    reads.
- **scale200_fix association done.**
  - Variant survivors 3 -> 6 (adds katG S315T, rpsL K43R, embB 4247431);
    burden survivors 6 -> 8 (adds katG, rpsL); q_branch passes 137 -> 92.
  - All DR controls found.
  - The tree cost 53 billing-hours; the rerun's total is about 631 against
    380 approved.
- **Old vs new comparison** (`analysis/rerun_scale200/COMPARISON.md`):
  - 98% of shared cells unchanged; NOCALL +245k, mostly PE/PPE and masked;
  - accessory presence 802 -> 204 (unmeasurable loci, by design);
  - 30 references changed, all closer;
  - AA_INVERTED 6,451 -> 1,079;
  - survivors 3 -> 6, burden 6 -> 8, all DR controls found;
  - q_branch 137 -> 92, mostly SV intervals and IS6110 fragments no longer
    double-counted.
- **Novel events steps 1-2 done** (Phase A total 30.0 billing-hours).
  - prototype + assembly + depth + dysgu PASS: typed recall 0.77, breakpoints
    0.81, precision 0.58 (dysgu alone 0.51 / 0.65).
  - The depth scan carries the tandem copy-number classes; assembly raises
    precision.
  - Phase B compute awaits approval; the download is running.
- **Phase B approved at 45 billing-hours.**
  - Cohort `marinB63` (63 Marin isolates, real Illumina, subsampled to about
    120x; 28 runs were deeper, up to 1,150x).
  - Read preparation 51425836 is queued after download 51420772.
  - Then P1/P2, truth, callers, scoring.
- The user asked where the scale200 outputs are; all are on netscratch
  (listed in the reply), with no durable copy.
- **scale200 key outputs copied to the mirror**, `results/cohorts/scale200/{2026-10-01_original,2026-10-08_audit_fixes}` (178 MB, verified). Descriptive names, at the user's request.
- **NAME-1 added to the fix list:** rename the `refbias/` working folders (`results/` or `out/`), deferred by the user until the pipeline is stable.
- **AA explained to the user:** Fitch parsimony on the CX333 panel tree,
  reporting the MTBC ancestor (MRCA of non-canettii), with canettii breaking
  ties nearest first.
  - Measured: canettii decides 1,134 of 73,256 sites (1.5%); its main role is
    rooting.
  - The 7 Oct kit test gave 99.98% the same AA without canettii in the graph.
- **OUT-1 added to the fix list:** a canettii pseudo-isolate (simulated reads
  through P1-P5) for IS6110/SV/accessory polarity.
- **AA-1 added to the fix list:** AA_CHANGES (Fitch count), AA_FLAG HOMOPLASTIC/OUTGROUP, AA_SOURCE for assumed polarity. Measured: 1,741 panel sites (2.4%) homoplastic, 196 of 1,084 alt-ancestral.
- **Phase B P1/P2 submitted** (p1 51453199 -> p2sum 51453203), with 62
  isolates. MT_0080 (an in-panel control) was dropped after ENA reset its
  1.8 GB read-2 transfer three times. Groups: novel, in_panel (false-call
  controls), setE_polished (Peker, less reliable truth).
- **Accessory summary for scale200_fix** (`analysis/accessory_scale200_fix/SUMMARY.md`): 806 catalogue loci (481 IS6110-sized); 204 measurable, 86 variable; TbD1 (ACC_1761789) behaves as expected; 462 small variants inside accessory sequence, mostly singletons; no accessory association survivor.
- **Phase B caller job 51469186 hit its 3 h limit** after finishing 1 of 62
  isolates (prototype done on 5). On real reads (about 120x) the prototype
  finds about 1,100 clip clusters per isolate, against about 10-20 on the
  simulated reads. Most are unpaired single-side clusters (BND) where reads are
  clipped at scattered positions in GC-rich sequence. These are probably
  untrimmed low-quality ends or library chimeras, not breakpoints. Every one
  becomes a local-assembly window, so assembly took about 3 h per isolate.
- **Phase B cost overrun:** prep 1.2 + P1 28.7 + P2 45.1 + callers 24.1 =
  **99 billing-hours against 45 approved** (plus 18.5 for the download).
  P1/P2 alone cost 74, about twice the 38 estimate, because the real reads
  are deeper than the novelA40 simulations. No further compute until the user
  approves.
- **Real-read changes approved (steps 1-4) and tested** (job 51555392, 5 isolates):
  - Caller v3 trims low-quality clip bases (stops at the first base below
    Q20) and keeps a cluster only if half its reads clip within 1 bp of the
    main position and their clipped sequences agree (60%). Every cluster is
    logged to `<out>.clusters.tsv`.
  - Assembly runs only for typed candidates and for one-sided clusters with
    support of at least 20% of median depth, with at most 40 windows per isolate.
  - The depth scan is GC-corrected by GC-percentage bin.
  - Regression check: v3 caller on the 40 Phase A simulations.
- **5-isolate test results** (jobs 51555392 and 51556221, 1.8 billing-hours; Phase B total about 101 against 45):
  - Real reads have binned quality scores (scattered Q14 bases), so first-low-base
    trimming cut good clips. Changed to a sliding-window trim (v3b): stop where
    the mean quality of 5 bases first falls below Q20.
  - Clip clusters on mar_QC_6: 1,297 (v2) -> 92 (v3b). Typed calls fell from
    108 to 24, precision rose from 0.17 to 0.58, and breakpoint recall went from
    0.51 to 0.47.
  - Phase A simulations are unchanged (typed recall 0.333 -> 0.330, breakpoint
    recall 0.556).
  - Novel isolates (TB1236, TB1612), prototype + assembly + dysgu PASS: typed
    recall 0.39 / 0.50, precision 0.28 / 0.48.
  - In-panel controls (should have no events), typed calls:
    - prototype + assembly: 0 and 2;
    - dysgu PASS: 4 and 11;
    - delly PASS: 2 and 9;
    - depth scan: 161 and 170.
  - The depth scan stays unusable on real reads even with GC correction:
    precision 0.02-0.10 and about 20 strong runs in each control. Decision on
    it is left to the user.
- **Phase B full run approved (57 isolates, about 9 billing-hours, depth scan
  option 1: run it, report it separately)**: job 51557679,
  `phaseB/06_full_v3b.sbatch`, outputs in `phaseB/out/v3b/` (the 5 test
  isolates copied in).
- **User decision: call novel events only in mappable sequence and flag the
  rest as uncallable.** Written (not yet run):
  - `analysis/novel_events/callable_mask.py`: a mask in each matched
    reference's own coordinates. It unions two sources:
    - R's own repeats, measured with build_repeat_mask.py's tests (50-mer
      paralogy and 9-mer tandem recurrence);
    - the H37Rv repeat mask lifted to R (minimap2 asm5 + paftools liftover).
    The gene extension and name list are dropped for R (no annotation) and
    come in only through the lift.
  - `score.py --refmap --mask-dir`: true events and calls within 50 bp of the
    mask are uncallable. Recall is measured over callable events, uncallable
    calls are flagged, and both counts are reported.
- **Uncallable masks built** (job 51560887, 3 min, about 0.6 billing-hours
  against 5 approved): 51 matched references (Phase A and B), each 9.3-10.6%
  uncallable (median 10.0%). In each, about 150 kb is measured on the
  reference itself and about 410 kb is the lifted H37Rv mask (594 kb in
  H37Rv). Outputs in `analysis/novel_events/out/masks/`.
- **GRIDSS** (user prefers established tools): the user approved the install
  and a 2-isolate timing test. Installing v2.13.2 from bioconda into
  `~/.conda/envs/gridss`, with the package cache in `~/.conda/pkgs`. The
  timing job is `analysis/novel_events/gridss/01_time.sbatch` (mar_QC_6,
  mar_TB1612; P2 BAMs; the R mask as the exclude list).
- **Phase B full run done** (job 51557679, 64 min, 9.6 billing-hours).
  Rescored on callable sequence; review written in
  `analysis/novel_events/REVIEW_2026-10-09.md`.
  - 84-87% of true events (simulated and real) have breakpoints in
    uncallable sequence.
  - In callable sequence, ours + dysgu typed recall: 0.74 (novel), 0.77
    (setE), 0.90 (simulated). Precision on real reads is low (0.11-0.34),
    driven by dysgu.
  - In-panel controls: median 1 false call per isolate for our caller, but
    about 7 controls (mostly the N series) have hundreds.
  - The depth scan is unusable.
  - Phase B total about 111 billing-hours against 62 approved in all.
- **GRIDSS timing** (job 51564761, 201 s, about 0.45 billing-hours): about
  100 s per isolate with 4 threads (mar_QC_6 108 s, mar_TB1612 91 s) on the
  un-QC'd P2 BAMs. Full 102 isolates: about 23 billing-hours at this setting.
  - PASS records: mar_QC_6 189 (177 of them single breakends), mar_TB1612 21
    (13 single).
  - Not scored yet: GRIDSS reports breakend pairs, which need converting to
    DEL/INS/etc before score.py can use them.
- **User raised read preprocessing:** reads should be QC'd and downsampled to
  100x. Currently there is no QC at all (no adapter or quality trimming, in
  Phase B or in the pipeline), and Phase B is subsampled to about 120x
  (47-120x). Proposed: fastp (adapters, 3' sliding-window Q20 trim, min
  length 50), then rasusa to 100x; a 5-isolate test (about 4 billing-hours)
  before a full rerun (about 80). Awaiting the user. Open question: were the
  production CRAMs trimmed upstream?
- **Read QC written to match the user's standard** (their hybrid-assembly
  pipeline settings, plus the user's statement: 100x, minimum 60x):
  - `phaseB/08_qc.sbatch`: fastp with the user's options plus `--cut_right`
    (4-base window, Q20) and `--trim_poly_g` (user approved these); TB-Profiler
    6.7.0 from the container; rasusa 2.2.2 to 100x, seed 7; unaligned BAM.
  - `09_kraken.sbatch`: kraken2, lab standard database (59 GB, memory-mapped,
    one job for all samples).
  - `qc_summary.py`: the user's pipeline rules, ported (MTBC >= 85% of
    classified, any other taxon <= 5%, mixed = two top-level lineages each
    >= 10%, final depth >= 60x).
  - `gridss/02_array.sbatch`: 2 cores, 6 GB, 4 GB heap, 15 min per isolate
    (from seff: 50% CPU, 3.4 GB peak).
  - Corrected cost for the 5-isolate test: about 17 billing-hours, not 4
    (kraken2's 72 GB about 6; QC about 5; P1+P2 about 6). Awaiting approval.
- **Read QC test (user: without kraken2):** QC array 51567729, about 2.4
  billing-hours. All 5 PASS: TB-Profiler single-lineage; fastp removed
  6-19% of bases (3' trim); final depth 64-100x (`refbias/marinQC5/qc_table.tsv`).
  Cohort `marinQC5` registered; P1/P2 chain 51568650 -> 51568653. Then callers
  51568688 and GRIDSS 51568689 on the QC'd reads, and GRIDSS 51568691 (tuned
  resources) on the same 5 without QC.
- **score.py reads GRIDSS** (`--gridss-dir`): breakend pairs typed as
  DEL/INS/dup/INV/REARR by standard breakend notation; single breakends untyped.
  First look, 2 isolates without QC, callable sequence: GRIDSS PASS
  breakpoint recall 0.54, typed recall 0.15 (13 events), mostly single
  breakends.
- **Read QC test results** (5 isolates, callable sequence, 17 true events):
  - P1 picked the same references after QC.
  - QC cost about 6.5 billing-hours in all (QC 2.4, P1 1.1, P2 1.7, callers
    0.5, GRIDSS before and after QC 0.8), against 11 approved.

  | | before QC | after QC |
  |---|---|---|
  | ours + dysgu PASS, typed recall / precision | 0.71 / 0.38 | 0.65 / 0.57 |
  | dysgu PASS calls | 45 | 21 |
  | ours + assembly, typed recall / precision | 0.53 / 0.67 | 0.47 / 0.64 |
  | GRIDSS single breakends (PASS) | mar_QC_6 177, TB1236 19, TB1612 13 | 10, 2, 5 |

  - Controls: 0-2 calls per caller after QC, against 0-7 before.
  - GRIDSS PASS typed recall is 0.18 either way: it finds breakpoints (0.53
    before QC, 0.59 after) but types few, and makes 2-3 typed calls in each
    control.
  - With tuned resources GRIDSS ran in 55-141 s at billing 3, about 0.1
    billing-hours per isolate.
- **Tool decisions (user, 2026-10-09):**
  - GRIDSS dropped ("doesn't seem so promising"); scripts kept in
    `novel_events/gridss/`.
  - breseq not pursued: too slow, and it misses events without saying why.
  - Lower bwa seed length (-k 15) not tested.
  - fastp settings reported: the user's (`--detect_adapter_for_pe --correction
    --qualified_quality_phred 20 --length_required 36`) plus `--cut_right`
    (4-base window, mean Q20) and `--trim_poly_g`.
- **Depth scan dropped** (user, 2026-10-09): removed from `phaseB/10_test_qc.sbatch`,
  the template for the QC rerun; README notes it.
- **Noisy in-panel controls** (`novel_events/controls/`, user request; no compute yet):
  - `false_calls.py` -> `false_calls.tsv`: typed false calls per control in
    callable sequence. 10 noisy: all 9 N-series (ERR27046xx-7xx, one ENA
    study, HiSeq 2500) and RW_TB008 (MiSeq), with dysgu PASS 135-505 and ours
    9-424. The other 15 have 0-7.
  - Our caller's false calls there are almost all inversions of about 60-220
    bp, typed from clipped reads alone ("not assembled").
  - Reads at one (N1274, 69,083): split reads whose other part maps on the
    opposite strand about 100 bp away (fold-back chimeras), some repeated as
    identical copies (duplicates; P2 does not mark them).
  - Quick count on 300 kb (login node): inverted local chimeras per 1,000
    reads are 4.8-59 in 4 noisy controls and 0.0-0.1 in 4 clean ones.
    The duplicate share does not separate them.
  - Written, not run: `01_artifacts.sbatch` (genome-wide split-read counts,
    GATK CollectSequencingArtifactMetrics for 8-oxoG) and `kraken_list.txt`
    (10 noisy + 5 clean) for `phaseB/09_kraken.sbatch`.
- **Control artifact results** (jobs 51622822 artifacts, 51623580 Picard;
  `controls/summary.tsv` via `summary.py`):
  - Fold-back chimeras per 1,000 reads, whole genome: noisy 4.8-56,
    clean 0.006-0.31.
  - Picard CollectAlignmentSummaryMetrics PCT_CHIMERAS (PAIR; mates on
    different contigs, insert over 100 kb, or not FR): noisy 1.4-10.3%,
    clean 0.23-1.05%. Same ranking, narrower gap (RW_TB008 1.37% vs
    M0003941_3 1.05%).
  - Oxidative damage (GATK pre-adapter G>T, Phred): noisy 48-59, clean
    39-47. Not the cause; no library looks damaged.
  - kraken2 (51622823) still running.
- **User: Peter's fastp settings for the production CRAMs likely differ** from
  the hybrid-assembly settings used in Phase B QC. CRAM @PG lines show
  fastp_trimming then bwa mem 0.7.19 to H37Rv, without the fastp options.
- **fastp:** user doubts production used --cut_right/--trim_poly_g; include them
  only if they help. `08_qc.sbatch` now takes an optional third argument
  `extra` (default: user's base settings). Comparison approved: cohort
  `marinQC5b` (same 5 isolates, base settings), QC 51635721 (final depth
  70.6-100x vs 64-100x with extras), P1/P2 51637341-51637354, callers
  51637363 -> `phaseB/out/qc5b`. `make_qc_cohort.py` builds the runner tables
  (reproduces marinQC5's).
- **kraken2 51622823 stalled** memory-mapping the 59 GB database from boslfs02
  (23 s CPU in 36 min); cancelled, about 13 billing-hours lost.
  `09_kraken.sbatch` now copies the database to /dev/shm once; resubmitted as
  51636061.
- **Fold-back chimeras in production CRAMs** (`controls/03_production.sbatch`,
  job 51631401, 500 random of 54,461, seed 7, about 1.5 billing-hours;
  `controls/production_chimera.tsv`): per 1,000 reads median 0.31, p75 3.6,
  p90 9.6, max 87; 40% above 1. Bimodal: one mode near 0.1 (like the clean
  controls) and one near 3-5 (like the noisy controls), trough about 0.3-1.
  SAMEA median 0.17 (95 of 282 above 1), SAMN median 0.92 (103 of 214).
  So a 1-per-1,000 fail cutoff would drop about 40% of production; this
  looks like a library-prep difference, not rare bad samples.
- `split_reads.awk` holds the fold-back count, shared by 01 and 03.
- **kraken2 on controls** (51636061, 14 min: database copy 2 min, then about
  25 s per isolate; about 5 billing-hours, about 18 with the stalled job):
  all 15 pass the user's rules. MTBC 98.5-99.7% of classified reads; largest
  other taxon at most 1.8% (Homo, in clean control 01_R1134). Noisy controls:
  Staphylococcus at most 1.0% (N1176), the rest under 0.6%. Contamination
  does not explain the noisy controls. `controls/summary.tsv` updated.
- **fastp comparison** (marinQC5 with --cut_right/--trim_poly_g vs marinQC5b
  base settings; same 5 isolates, same references picked; 17 callable true
  events; 6.7 billing-hours):

  | | with extras | base only |
  |---|---|---|
  | ours + dysgu PASS, typed recall / precision | 0.65 / 0.57 (35 calls) | 0.71 / 0.56 (39) |
  | ours + assembly | 0.47 / 0.64 (14) | 0.53 / 0.67 (15) |
  | dysgu PASS | 0.59 / 0.52 (21) | 0.65 / 0.50 (24) |
  | dysgu, all calls | 48 | 402 |

  Differences are 1 event and 3-4 calls: no measurable gain from the extras.
  Without 3' trimming dysgu makes many more non-PASS calls, but its PASS
  filter removes them. Decision: base settings (closer to production).
- **Correction: the 500-CRAM chimera sample ignored the standing inclusion
  screen** (`bin/select_isolates.py`: paired-end, meandepth >= 60, not mixed
  at any frequency, complete run, error_rate <= 0.01). 422 of 500 pass (52 low
  depth, 28 high error rate, 6 mixed). Among the 422: median 0.23 per 1,000,
  35% above 1, 83 above 4.8. By read length (avglen_1): 110-160 bp median
  0.18 (27% above 1); 160-260 bp and 260-400 bp all 29 above 1 (medians
  9.3, 11.0). The results table has no contamination (kraken2) column.
- **Fold-back study approved (synthetic study + ENA query).** `analysis/foldback_sim/`.
  - ENA (`ena_library.py` -> `ena_library.tsv`, portal API, metadata only):
    498/500 production samples, 63/63 Phase B runs. Rate by instrument
    (median per 1,000; share above 1): NextSeq 500 0.15, 12%; HiSeq 2500
    1.89, 64%; MiSeq 7.1, 85%; HiSeq 2000 0.06, 7%; NextSeq 2000 0.28, 32%.
    Protocol text is mostly empty. Strong study effect: in the 21 studies
    with 5 or more isolates, 89% of a study's isolates fall on the same side
    of 1 per 1,000; studies run from 0% to 100% above 1. Nextera XT studies
    differ (PRJEB9680 median 0.16; PRJEB6273 8.99), so the kit name alone
    does not decide it. N-series = PRJEB27802 (Swiss TPH, HiSeq 2500).
  - Synthetic design: Phase A genomes R27252 (L1), QC_6 (L4), R37765 (L2)
    with their novelA40 references; wgsim 1.47M pairs (100x), 150 bp, frag
    450+/-50, 0.2% error, seed per genome; `foldback_sim.py` replaces 0,
    0.5, 2, 6, 15% of pairs (inverted-repeat model, k=10 IR 20-300 bp apart,
    about 7,000 sites per genome) and 2, 15% (random model). Cohort
    `foldsim` (P2 only, reference pinned; run folder stamped with the build
    because the hand-made refmap predates the runner's record). Login test:
    at 5% replaced, 15.8 junction-crossing reads per 1,000 and 10.9 detected
    by the fold-back count (about 70%).
  - Jobs: reads 51676142 (done), P2 51676410/51676412, callers+score
    51676462 (`02_call_score.sbatch`, `analyze.py`).
  - Polishing (`make_draft.py`, `03_polish.sbatch`: 300 one-base homopolymer
    indels + 30 SNPs; Polypolish as the user's 05_polish, Pilon one round as
    Marin; scored with minimap2 + paftools call against truth) written and
    tested (draft vs truth = 330 differences), not submitted: cost check.
- **User clarified the design; simulation cancelled.** The user meant: SV calls
  from the real Illumina reads used for polishing vs from synthetic reads
  simulated from the same isolate's complete genome. My chimera-injection
  simulation was a different design; cancelled at the user's request
  (51676410/12/62). It cost 12.2 billing-hours (reads + 18 of 21 P2 tasks),
  not the 7 I had stated.
- **Real vs synthetic, 25 in-panel controls** (`novel_events/realsyn/`):
  `01_reads.sbatch` (51679094): wgsim from each control's Marin assembly,
  matched to its real run from samtools stats on the marinB63 P2 BAM (pairs,
  read length 100-151, insert mean/sd), 0.2% error. Cohort `marinSyn25`
  (P2 only, reference pinned to marinB63's P1 choice; run folder stamped).
  P2 51679401/02, then `02_call_compare.sbatch` 51679403 (our caller +
  assembly on synthetic BAMs, `compare.py`: real-only and synthetic-only
  typed calls in callable sequence per caller; small variants real vs
  synthetic).
- **marinSyn25 P2: 11 of 25 tasks failed, both P2 robustness bugs (fix list):**
  - DYSGU-1 (8 tasks): dysgu's own `--clean` work-folder removal crashed on
    an NFS lock file (`.nfs...`) after calling finished; p2_call.sh then
    deletes the output and fails the sample.
  - DYSGU-2 (3 tasks: M0016737_0, 01_R1430, TB3368 synthetic): dysgu found
    no events and wrote a VCF whose header puts ##contig after #CHROM with
    no sample column; stamp_build_id.sh could not parse it.
  For this analysis, without touching production code: the 3 count as 0
  dysgu calls (their BAM, HaplotypeCaller and delly outputs are complete);
  for the 8, `02_call_compare.sbatch` reruns dysgu alone with P2's options on
  the P2 BAM into `realsyn/out/dysgu/`. P2 cost 10.3 billing-hours with the
  failures; comparison job resubmitted (no dependency).
- **Real vs synthetic results** (51681813; `realsyn/out/compare.tsv`; total
  about 11 billing-hours). Typed SV calls in callable sequence, real /
  synthetic:
  - Synthetic reads of every control: 0-2 calls per caller.
  - 10 chimeric controls (fold-back 4.8-56 per 1,000): ours 2-424,
    dysgu 135-505, delly 206-6,586 with real reads; essentially all real-only.
  - 15 clean controls (0.006-0.31): ours 0-2, dysgu 0-7, delly 0-3.
  - Small variants: real-only SNVs 0-13 and indels 0-21 in both groups;
    chimeras do not visibly add SNP/indel calls. Side finding: in the
    N-series, RW_TB008, N0004 and TB3251 the synthetic reads (from Marin's
    assembly) call 21-43 indels and share 13-42 calls with the real reads
    against the panel genome R: Marin's assembly and the panel GCF assembly
    of the same isolate differ at those sites. Not investigated.
- **DYSGU-1 and DYSGU-2 fixed in `bin/p2_call.sh`** (user request):
  - DYSGU-1: dysgu's working folder is now node-local (mktemp under
    $TMPDIR or /tmp) and removed by P2 with a guarded rm, not by dysgu's
    --clean. Test (51684684): dysgu with the new command on mar_TB3251's P2
    BAM gives records identical to its existing VCF (364 of 364).
  - DYSGU-2: when dysgu writes no records, P2 appends the sample name to
    the #CHROM line (dysgu writes FORMAT with no sample column). Checked on
    syn_mar_TB3368's empty VCF: bcftools rejects the original and parses the
    fixed one; stamp_build_id.sh stamps it.
  - Affects only new P2 runs; existing outputs are unchanged.
- **Fold-back on cohorts in use** (user request): `analysis/foldback_qc/`,
  `cohorts_in_use.tsv` = union of pilot, scale100, l49, l7, scale200,
  gwas1000 crams tables (1,200 isolates; 16 already in the 500 random
  sample). `01_measure.sbatch` job 51684612 (8 cores; about 3.3
  billing-hours expected, above the 2 I quoted).
- **Fold-back on cohorts in use** (51684612, 24 min, 3.6 billing-hours; 1,200
  of 1,200 measured; `foldback_qc/cohorts_in_use_rates.tsv`,
  `cohorts_in_use.summary.txt`). The job's output overwrote its input list
  (named the same); results moved to `_rates.tsv`, list restored from git,
  script fixed to write `<name>_rates.tsv`. Share above 1 per 1,000: pilot
  26%, scale100 39%, l49 15%, l7 0%, scale200 30%, gwas1000 35%. gwas1000 by
  lineage above 1: L1 41%, L2 39%, L3 32%, L4 34%, L5 54%, L6 19%, L7 0%, L9 0%;
  RRDR 0 35%, RRDR 1 34%. By read length: under 110 bp 41%, 110-160 bp 27%,
  160-400 bp 72%.
- **Fold-back vs Peter's Mutect2 columns** (no compute; 1,684 isolates with a
  measured rate: cohorts in use + 500 random). Within lineage, median per
  isolate, rate <=0.31 vs >4.8 per 1,000: m2_low_alt L4 48 vs 784, L2 74 vs
  591, L3 78 vs 1,252 (3-16x); m2_fail about 3x; m2_pass +8% (L1) to +57%
  (L4); m2_mix_call +2% to +32%; mean depth similar. The lineage-based mixed
  flag (conflict_lineages / n_mixed / mix_freq) is 0-0.7% in every group.
  Association only (studies differ in other ways too); column definitions
  are Peter's and not checked. HaplotypeCaller in P2 is haploid, which is
  why the real-vs-synthetic test showed no SNP effect.
- **Fold-back on Phase B** (`foldback_qc/02_phaseb.sbatch`, 51704033, 53 s,
  about 0.13 billing-hours; `foldback_qc/phaseB_rates.tsv`, all 62):
  - novel: 19 of 20 at 0.01-0.04; **N0153 14.4** (PRJEB31443, HiSeq 2500).
    N0153 alone carries 457 of the novel group's 669 false dysgu PASS calls
    and 19 of our 38 false calls in callable sequence. Without it: dysgu
    precision about 0.20 (was 0.07), ours about 0.67 (was 0.50); its 0
    callable true events are lost. Clean novel isolates still have 15-21
    false dysgu calls each.
  - setE_polished (one study, PRJNA720906, NextSeq 500): all 0.25-1.21,
    between the clean (<=0.31) and chimeric (>=4.8) controls. 8651_04 (1.21)
    and QC_6 (1.06) are not outliers in false calls, so a cutoff at 1 would
    split this study for no visible benefit. Data now support a cutoff
    somewhere in 1.2-4.8; cohorts in use above 2: 30%, above 3: 24%.
  - in_panel: as before (10 above 4.8).
- **Fold-back cutoff set at 2 per 1,000 (user, 2026-10-09).** Phase B rerun
  list `foldback_qc/phaseB_keep.txt`: 51 isolates (19 novel, 17 setE, 15
  clean controls); excluded N0153 and the 10 chimeric controls.
- **Firm rerun cost** from measured per-isolate costs (marinQC5b chain and
  kraken2 51636061), QC scaled by raw size (rerun mean 123x vs test 113x):
  QC 0.53 x 1.09 x 51 = 29.5; kraken2 2 min copy + 51 x 25 s at billing 22
  = about 9; P1 0.26 x 51 = 13.3; P2 0.40 x 51 = 20.5; callers + assembly
  0.13 x 51 = 6.8; scoring under 0.5. Total about 79 billing-hours.
- **Phase B rerun submitted (user approved, about 79 billing-hours).** Cohort
  `marinQC51` (registry row added; `runroot/refbias/marinQC51` -> netscratch).
  QC array 51710211 (51 tasks, base fastp settings), kraken2 51710219. Next,
  by hand when both finish: `make_qc_cohort.py` (PASS only, kraken2 rules
  included), runner P1/P2, then `phaseB/11_call_score.sbatch` (callers +
  assembly, score per group). Known limit of score.py: isolates with no truth
  rows are not counted in precision, so control false calls are counted
  separately (as `controls/false_calls.py`).
- **Rerun QC done** (51710211 QC 61.5 billing-hours, estimate 29.5;
  kraken2 51710219 9.1, as estimated). QC tasks averaged 12.1 min (test
  5.3) at about 1.3 of 4 cores busy and up to 6 GB memory; probably
  filesystem contention with 51 tasks at once (not confirmed). Spent so far
  70.7 of the approved 79.
  - Result (`refbias/marinQC51/qc_report.tsv`, `qc_report.py`; page
    https://claude.ai/artifact/3fhVaqUwwf8222jQKHWPpg): 47 PASS, 4 FAIL on
    depth after QC (M0011368_9 45x, M0016395_7 49x, TB3368 58x: clean
    controls; TB3386 55x: novel, 4 callable events), 11 excluded for
    fold-back. TB-Profiler: none mixed. kraken2: MTBC 96.2-99.7% of
    classified; largest other taxon at most 3.5% (mostly Homo).
  - Paused before P1/P2 for cost approval: remaining about 37
    (0.79 x 47), projected total about 108.
- **QC slowdown diagnosed from task logs** (51710211 vs test 51635721):
  TB-Profiler's freebayes step took a median 310 s (16-634 s) against 13 s in
  the test; its alignment median 84 s against 59; fastp 48 s against 18.
  Freebayes time barely tracks depth (r = 0.25); the fastest tasks were the
  first few at low depth. Most likely cause: 51 tasks reading the
  TB-Profiler database (118 MB, boslfs02) by random access at once, as with
  kraken2. Not confirmed. Options: copy the database to node-local /tmp per
  task; throttle arrays (runner `--throttle N`).
- **TB-Profiler database copied to node-local /tmp per QC task** (user approved;
  `08_qc.sbatch`, `cp -rL` into the task folder, `--db_dir` points there).
  Untested until the next QC run.
- **Rerun continued (user approved, throttle 15):** cohort tables 47 PASS;
  P1 51741763/76, P2 51741784/86 (`--throttle 15`, arrays 1-47%15), callers
  and scoring 51741788.
- **Phase B rerun on clean data: results** (`phaseB/rerun_compare.py` ->
  `out/qc51/compare.tsv`). Cost: P1 14.6, P2 22.5, callers 6.1 = 43.2
  (estimate 37); rerun total 113.9 (QC 61.5, kraken2 9.1). 47 passed QC:
  18 novel, 17 setE, 12 clean controls (I earlier said 17/17/13: wrong).
  QC_10's P1 reference changed on clean reads, so it is left out of both
  sides (truth is in the old reference's coordinates): 46 compared.
  Same isolates, raw reads (before) vs clean (after), callable sequence,
  typed recall / precision (typed calls):
  - novel (70 events): ours+dysgu 0.77/0.28 (318) -> 0.77/0.48 (189);
    dysgu 0.69/0.20 (262) -> 0.71/0.39 (134); ours 0.50/0.66 (56) ->
    0.50/0.71 (55); delly 0.10 recall either way.
  - setE (78 events): ours+dysgu 0.81/0.33 (302) -> 0.74/0.39 (244);
    dysgu 0.65/0.27 -> 0.65/0.37; ours 0.51/0.46 -> 0.50/0.44.
  - 12 clean controls, typed calls (all false): ours 7 -> 6 (0-1 each),
    dysgu 37 -> 31 (1-7), delly 26 -> 21 (1-3).
  - Against the original full Phase B (raw, all isolates incl. N0153):
    novel ours+dysgu precision 0.11 -> 0.48 at recall 0.74 -> 0.77.
- **delly dropped from novel-event scoring** (user, 2026-10-10): removed from
  `score.py` (delly_all, delly_pass, dysgu+delly_pass) and
  `phaseB/rerun_compare.py`; compare.tsv regenerated (other rows unchanged).
  Production P2/P4b/P5 still run and use delly: not touched, decision open.
- **dysgu's remaining false calls on clean reads** (`phaseB/dysgu_false.py`
  -> `out/qc51/dysgu_calls.tsv`, `dysgu_filters.tsv`; existing outputs, no
  compute). Same 46 isolates; labels by score.py's rules. dysgu PASS false
  calls: novel 82, setE 95, controls 31.
  - dysgu's AF (share of reads at the site with direct evidence; true events
    in these clonal isolates have median 0.36, not ~1) separates them: true
    calls median 0.38 (novel) / 0.34 (setE), false 0.05 / 0.18 / 0.07
    (controls). SU separates less; PROB and MAPQP poorly.
  - Filter sweep, ours+asm + dysgu (AF >= x), recall / precision:
    none: novel 0.771/0.481, setE 0.744/0.393, controls 31 false;
    AF>=0.1: 0.771/0.689, 0.718/0.447, 9; AF>=0.2: 0.771/0.735,
    0.718/0.481, 3; AF>=0.3: 0.743/0.745, 0.615/0.464, 2.
    setE loses 2 events at 0.1 or 0.2: INS 87 bp at 28469 in QC_1 and QC_6
    (AF 0.04/0.06, our caller misses it).
  - What survives AF>=0.2 (59 false): 31 within 200 bp of a true event and
    5 within 1 kb (same-size insertions placed elsewhere in the same repeat,
    or dysgu INS at the end of a true deletion: scoring tolerance, not
    caller error); 15 are one 67 bp deletion of a palindromic sequence
    (GGGCTGGCGAGCAGACGCAAAATCCCCCGCACG...) at the homologous locus near
    2.35 Mb on several references, in 12 setE and 3 clean controls whose own
    genome is the reference (so artifact; AF 0.2-0.4, split reads, our caller
    calls it too); 8 other.
  - Caveat: cutoff chosen on the same isolates it is scored on.
- **Applied (user, 2026-10-10): dysgu AF >= 0.2, palindrome-locus mask,
  scoring fix.**
  - `score.py`: `dysgu_af20` (PASS and FORMAT AF >= 0.2) and
    `proto_asm+dysgu_af20` (the combination from now on; the unfiltered
    rows kept for reference). Matching: `matches()` adds truvari bench's
    default rule (starts within 500 bp, sizes within 70%; truvari's sequence
    similarity test not applied, no allele sequences); `at_event()` (for
    precision and breakpoint recall) also accepts a call within 50 bp of the
    true end. Call sizes from SVLEN / the tables' size column.
    `controls/false_calls.py` and `phaseB/dysgu_false.py` use the same rules;
    `realsyn/compare.py` (call vs call) keeps `near()`.
  - `artifact_mask.py`: masks the 67 bp palindrome site in every matched
    reference by its two flanking 30-mers (between flanks: 50 references;
    GCF_001870145 has one flank only, +-100 bp masked). Written as
    `out/masks/<ref>.artifact.bed`; `load_masks` unions it; `07_masks.sbatch`
    runs it. Element present in 41 of 51 references: the site varies between
    strains. Cost: TB3237's true 68 bp insertion (its reference lacks the
    element) is now uncallable; it had been found.
  - Regenerated: `out/qc51/compare.tsv`, `score_*`, `dysgu_*`. Earlier
    stages' score files keep the old rules.
  - Same 46 isolates, clean reads, ours+asm + dysgu, recall / precision
    (false calls); one change at a time:
    old rules:        novel 0.771/0.481 (98), setE 0.744/0.393 (148), controls 37
    + AF >= 0.2:      0.771/0.735 (31), 0.718/0.481 (94), 9
    + mask:           0.768/0.730 (31), 0.718/0.561 (68), 1
    + scoring fix:    0.797/0.826 (20), 0.744/0.677 (50), 1
    (scoring fix alone: 0.800/0.540, 0.756/0.484). Novel truth now 69
    callable events. By caller after all three: novel ours 0.536/0.778
    (54 calls), dysgu_af20 0.667/0.869 (61); setE ours 0.526/0.600 (80),
    dysgu_af20 0.603/0.760 (75). Controls: ours 1, dysgu_af20 0.
    Raw-read (before) side, same rules: novel 0.797/0.758, setE 0.795/0.667.
- **Artifact mask made opt-in** (user: "don't commit to this mask yet"):
  `score.py --artifact-mask`, `load_masks(..., artifacts=False)`; default
  scoring no longer applies it. Outputs regenerated without it. Clean-read
  rerun, no mask: ours+dysgu_af20 novel 0.800/0.829 (70 events), setE
  0.744/0.580; controls ours 6, dysgu_af20 3.
- **Our caller's setE false calls** (`phaseB/ours_false.py` ->
  `out/qc51/ours_false_setE.tsv`, `identity.tsv`; no compute). 46 false of
  94 typed calls, by cause:
  - palindrome site 14 (the 67 bp deletion; mask decision open);
  - truth mismatch 14, all QC_7: identity check (assembly SNPs vs read SNPs
    vs R, callable sequence) gives 276 discordant for QC_7 (218 reads-only,
    58 assembly-only) against 0-8 for 40 of 46 isolates. QC_7's reads are
    not the strain its assembly is from; several of these calls are found by
    dysgu at AF 0.3-0.4 too. Also flagged: TB1612 79 (68 assembly-only),
    TB3054 42 (41 reads-only), both novel;
  - short anchor 9 (+1 in novel, TB3162; 0 of 91 true calls): reads of
    sequence absent from R (QC_1's sits in its own assembly at ~1.98 Mb)
    aligned by a 19-base match at MAPQ 40, clipped on both sides; one read
    set makes both clip clusters, so our caller calls an INS. Recur at the
    same R positions in QC_1/QC_6/QC_7. Candidate fix: ignore clips from
    reads with < 30 aligned bases (bwa's default minimum score, -T 30);
  - small indel 6: real 24-29 bp indels in the assembly, below the truth's
    50 bp; our caller reports them as INS of unknown size (not assembled);
  - other 3: two tandem-repeat INS 72 bp where truth has 108 bp at the same
    site, one unexplained (3003_06, 1171103).
  - Projected setE ours / ours+dysgu_af20 (calls removed, not rerun):
    now 0.526/0.511, 0.744/0.580; - short anchors 0.526/0.565, 0.744/0.610;
    - QC_7 0.485/0.650, 0.727/0.664; + palindrome mask 0.485/0.830,
    0.727/0.809.
- **What the palindrome mask loses** (user asked; exploration, no code kept).
  Element state at the site (flank search) in all 63 Phase B assemblies vs
  their matched reference: 61 agree (QC_6 differs by 1 bp inside it). Two
  real events, both in the truth set: TB3237 INS 68 (L4.4.1.1; its
  reference lacks the element) and TB3386 DEL 68 (L4.3.3; failed read QC,
  so only in the raw-read run). Both found by both callers (ours exact; dysgu
  INS size 184/236, DEL 69). Element absent from 5 assemblies in L4.3.3,
  L4.5, L4.8, L6.2.1, L6.3.1, and from 10 of 51 references: scattered across
  lineages. In the 46 compared isolates (callable calls) the mask removes ours
  19 false (setE 14, controls 5) + 1 true; dysgu_af20 15 false (setE 12,
  controls 3) + 1 true. Real and artifact overlap on dysgu AF (real 0.20-0.45,
  artifact 0.2-0.4). Depth inside the element / flanks (MAPQ >= 20): TB3386
  0.00; every isolate carrying the element 0.08-0.66 (setE 0.08-0.21,
  controls 0.35-0.66): the element is poorly covered everywhere, never
  absent. The artifact only ever removes the element (all site calls in
  isolates carrying it are DEL). Alternative to the mask: at the site, keep
  INS; keep DEL only when depth inside the element is ~0. One real DEL to
  test it on.
- **The "palindrome site" is the user's tatC-helY hairpin element**
  (MtbPangenome/loci/hely_tatc; H37Rv 2,352,065; 68 bp, CTCGCC direct
  repeats, precise excision, ~30 independent deletions on the CX333 tree,
  230/332 carriers). The low-AF deletion calls are the same precise 68 bp
  excision product. Signal by batch, element carriers in marinQC51 (depth in
  element / flanks; dysgu DEL AF, max per isolate): setE NextSeq 500
  (SRR1419*, n=17) 0.14 (0.08-0.21), AF 0.23 (0.01-0.43); novel HiSeq 2500
  (SRR6356/6807/7516, n=15) 0.22-0.30, AF 0.07-0.11; controls HiSeq 2500
  (SRR1723*, n=10) 0.50 (0.35-0.66), AF 0.14 (0-0.37). Present in every
  batch, strength batch-dependent. Library protocol not in ENA metadata.
  Artifact (polymerase skipping the hairpin between the direct repeats) vs
  real excision in culture not resolved; decisive test: PCR-free long reads.
- **Repeat family members, short-read signal** (`phaseB/repeat_family.py`
  -> `out/qc51/repeat_family.{members,isolates}.tsv`; 46 isolates; 41 H37Rv
  copies from loci/hely_tatc/repeat_family_loci.tsv + tatC-helY + fadE22,
  located in each matched reference by H37Rv flanks / own sequence).
  - Coverage dropout is family-wide: depth inside / flanks median 0.25-0.42
    for 33 of 43 members, ~1.0 for the weakest folders (MFE > -0.5/nt).
  - Deletion-junction reads (clipped 15-mer continues exactly 30-150 bp
    further on): tatC-helY 885 reads in 41/42 carriers, length 68, per 100x
    depth setE 39, novel 18, controls 20 (medians). Others: Rv2522c 56 reads
    in 26 isolates (73 bp), nrdB 26 in 17 (75), fadE22 16 in 12 (68),
    3348486 12 in 8; medians 0-1.7 per 100x. 755253: 114 reads in 3
    isolates (QC_1, QC_6, QC_7) = a real 79 bp deletion in the truth set.
    22 members: no junction reads. dysgu calls a deletion only at tatC-helY
    (39 carriers) and at the real 755253 event (3).
- **setE long reads exist** (ENA PRJNA720906, "Short- and long-read WGS of
  M. tuberculosis"): Oxford Nanopore runs (MinION/GridION) on the same
  BioSample for all 17 setE isolates, 27-485x, 13.0 Gb / 13.8 GB fastq.
  Library protocol not in ENA metadata. Not downloaded (needs approval).
- **Long-read test submitted** (user approved download + run, 2026-10-10):
  `longread/runs.tsv` (17 setE Nanopore runs, ENA PRJNA720906),
  `01_download.sbatch` -> 51826243 (array 0-16%6, md5-checked, to
  netscratch data/external_reads/marin_ont), `02_align.sbatch` -> 51826244
  (minimap2 -ax map-ont to the marinQC51 matched reference, aftercorr),
  analysis `lr_element.py` (reads anchored >= 200 bp both sides: member
  present / absent (member-sized deletion or < 20% aligned) / other).
  Estimate ~5 billing-hours.
- **Long-read test, first pass (2026-10-10).** Download 51826243: 6 of 17
  completed; 11 failed md5 although full size, with no interruption logged
  and gzip CRC errors (silent corruption in transfer). `01_download.sbatch`
  now deletes a file failing md5 and fetches it again from scratch (up to 3
  times, curl `--fail`); corrupt partials deleted; the 11 resubmitted as
  51844221, alignment 51844222 (aftercorr). First pass cost: ~13 CPU-h
  download (1 CPU, 1 GB) + ~0.3 CPU-h alignment (4 CPU).
- **Interim result, 6 isolates (QC_1, QC_3, QC_6, QC_7, QC_9, QC_10;
  `lr_element.py` -> `longread/out/lr_element.tsv`):**
  - tatC-helY: present in 226/226 spanning Nanopore reads (14-57 per
    isolate), 0 absent, 0 partial. Short reads at the same isolates: 27-45
    deletion-junction reads per 100x, dysgu DEL AF 0.20-0.32. If the
    excised form were even 5% of molecules, ~11 absent reads would be
    expected (P(0) ~ 1e-5). The short-read signal is therefore made during
    library prep / sequencing, not present in the DNA.
  - Positive control: 755253 (real 79 bp deletion in the truth set) is
    absent in 60/60 reads (QC_1) and 29/29 (QC_6); present in the others.
  - QC_7 is mixed in the long reads: 755253 absent in 47/77, Rv2273 absent
    in 15/33 (0 in every other isolate). Short-read dysgu AF at 755253 is
    0.36 in QC_7 vs 0.67-0.73 in QC_1/QC_6, so the Illumina reads look
    mixed too. Fits QC_7's identity-check failure (276 discordant SNPs).
  - Other family members: no absent reads except single reads at 6 members.
- **Fold-back rate vs tatC-helY artifact** (`phaseB/foldback_vs_element.py`
  -> `out/qc51/foldback_vs_element.tsv`; 42 carriers, Spearman). Pooled:
  junction/100x rho +0.57 (p 9e-5), dysgu AF +0.43, depth loss +0.51. Within
  each batch: none (rho -0.35 to +0.28, all p > 0.18). setE has both the
  highest fold-back rates (median 0.54/1k) and the strongest artifact
  (junction 40/100x); HiSeq batches are low on both. So both depend on the
  library/sequencing batch, but a sample's fold-back rate does not predict
  its artifact level, and samples well under the 2/1k cutoff still show it:
  the fold-back QC filter will not remove the tatC-helY artifact. Other
  family members' summed junction signal: no relation (rho -0.22).
- **Long-read test, 12 of 17 isolates.** Rerun 51844221: 6 more complete
  (several needed 2-3 full downloads; md5 failures over https recur), 5
  largest (1.4-2.3 GB) timed out at 2 h; resubmitted over ftp with 8 h
  (download 51869430, align 51869436). tatC-helY: 4 absent of 843 spanning
  reads (0.5%; single reads in mar_702_06, mar_706_05, mar_8651_04, each
  101-187 reads deep), 0 in the other 9. Short reads at the same isolates:
  26-54 junction reads per 100x, dysgu AF 0.16-0.43. Other family members
  with no short-read signal show the same background of single absent reads
  (icl1 5/894, PE_PGRS52 3/831, Rv2522c 3/883): Nanopore noise. Conclusion
  unchanged: the short-read excision signal is a short-read artifact.
  Positive control and QC_7 mixture as before.
- **Long-read test closed at 12 isolates (user, 2026-10-10).** 51869430 /
  51869436 cancelled; mar_3003_06 and mar_696_05 fastq downloaded (md5 ok)
  but not aligned; 3 incomplete partials deleted. Total cost ~45 CPU-h, nearly
  all 1-CPU/1 GB download time (estimate was ~5 billing-hours; the overrun is
  repeated md5 failures over https). Final: tatC-helY present in 839/843
  spanning Nanopore reads (4 absent = Nanopore background), vs 15-40%
  apparent loss in short reads: the short-read low-fraction excision is a
  short-read library/sequencing artifact. Not explained by per-sample
  fold-back rate (batch-level only).
- **SV calling: three changes approved (user, 2026-10-10), "check TB1612
  and TB3054 first".**
  - Identity follow-up (`phaseB/identity_detail.py` ->
    `out/qc51/identity_detail.tsv`; read alt fraction at each SNP where reads
    and assembly disagree). QC_7: 272/276 discordant intermediate (0.2-0.8),
    and 70/159 shared: mixed. TB1612: 78/79 intermediate, genome-wide, 9/32
    shared intermediate: mixed. TB3054: 38/41 in one 1.3 kb block (2.973-
    2.975 Mb) at ~0.5 with double depth, all 447 shared SNPs clean: a
    duplicated, diverged region, not mixture; kept. TB3334 (clean) has the
    same pattern on 50 bp. Excluded: QC_7, TB1612 (`phaseB/exclude.tsv`).
  - tatC-helY genotype rule replaces the mask: `element_genotype.py`
    (depth inside element / 200 bp flanks; absent < 0.05) ->
    `phaseB/out/element_gt.{marinQC51,marinB63}.tsv`. TB3386 0.00 (absent);
    carriers 0.094-0.94 (present); 5 references lack the element.
    `score.py --element-gt` drops DEL calls at the site in carriers.
  - Short-anchor fix: `breakpoint_caller.py` v4, a read's clips count only
    with >= 30 aligned bases. Pre-fix outputs saved in `phaseB/out/qc51_v3b`.
    Rerun 51902242 (11_call_score.sbatch).
  - Exclusions + genotype rule on pre-fix outputs (44 isolates): combo
    (ours+asm + dysgu AF>=0.2) novel 0.806 / 0.825, setE 0.727 / 0.774;
    controls: ours 1 false call, dysgu_af20 0.
- **Short-anchor fix result** (51902242, 29 min, ~3 CPU-h). Calls changed vs
  `qc51_v3b`: 28 false calls removed (all 9 known short-anchor setE calls,
  plus 19 others in setE/novel/controls, e.g. TB3054 x5, TB3162 x5 incl. its
  known one), 0 true events lost. Two true calls changed form only:
  mar_4549_04's two overlapping DEL calls merged into one; TB3091's
  redundant INS at the far end of a true 140 bp REPL became a BND (the DEL
  call still matches the event). 1 false call gained (TB3162 DEL extended).
  Final, 44 isolates (QC_7, TB1612 excluded; QC_10 reference changed),
  tatC-helY genotype rule, recall / precision:
    novel: ours 0.537 / 0.788; ours+dysgu AF>=0.2 0.806 / 0.832
    setE:  ours 0.485 / 0.830; ours+dysgu AF>=0.2 0.727 / 0.809
    12 clean controls: ours 1 false call, dysgu AF>=0.2 0.
  Before these three changes (46 isolates, no rule): setE combo 0.744 /
  0.580, ours alone 0.526 / 0.511. Remaining setE false calls of ours
  (ours_false.py, before exclusion/rule): 14 tatC-helY (now dropped by the
  rule), 14 QC_7 (excluded), 6 real 24-29 bp indels, 3 other.
- **Our caller's 3 "other" setE false calls (after the fixes).**
  - mar_702_06 and mar_9050_05, INS at 801,398 (tandem_dup, 72 bp; dysgu
    AF>=0.2 also calls it, 54 bp at 801,399): the truth's 108 bp insertion
    of a repeated unit at 801,326. Nanopore reads of 702_06 (162 spanning)
    put net insertion at 100-120 bp, so the assembly size is right and both
    short-read callers underestimate it in the tandem repeat (local assembly
    failed there). Scored false and missed: 72 bp from the true start and
    72/108 = 0.67 < pctsize 0.7. Detection is right, size is wrong.
  - mar_3003_06, INS at 1,171,103 (dysgu also calls it, AF < 0.2): clipped
    sequences are the two IS6110 ends (IRL TGAACCGCCCCGGCATGTCCGGAGAC on one
    side, IRR on the other), clip positions 5 bp apart (target-site
    duplication). 36 junction reads vs 59 intact spanning reads: IS6110
    inserted in ~38% of the Illumina culture, absent from the (Nanopore)
    assembly. A real subclonal transposition the truth set lacks.
  - So all 9 remaining setE false calls of ours are real sequence
    differences: 6 indels of 24-29 bp, 1 event sized wrong, 1 subclonal
    IS6110 insertion.
- **3003_06 IS6110 at 1,171,103 in Nanopore reads** (user approved;
  02_align.sbatch 51934978, 3.6 min, ~0.25 CPU-h; `longread/lr_site.py`).
  Reads covering the site: 17 with a 1,329-1,433 bp insertion there, 23
  clipped there by >= 1 kb (rest of the read aligned elsewhere), 136 running
  straight through. All 40 inserted/clipped sequences align to the
  reference's own IS6110 copy end to end (minimap2, identity 0.79-0.94,
  Nanopore error level). Insertion in 11% of fully spanning reads, 23%
  counting the clipped ones; short reads 38%. Subclonal in both DNA
  preparations, so the assembly's consensus correctly lacks it, and the
  truth set (majority sequence) cannot contain it. A real event scored false
  by construction; no scoring change.
- **Do the 24-29 bp indels belong to HaplotypeCaller?** Yes: P2's GATK
  HaplotypeCaller (ploidy 1) calls all 6 at the assembly's exact size, PASS,
  QUAL 2,300-4,800, every read supporting (AD 0,52 to 0,107). Our caller
  reports them as INS of unknown size because bwa clips reads at a 24-29 bp
  indel. Across the 44 isolates, SV calls (ours, dysgu AF>=0.2) within 60 bp
  of a HC 10-49 bp indel: 11 false, 1 true. The true one (8129_04, 1,360 bp
  IS6110 insertion at 2,633,361) has a weak HC 18 bp call: QUAL 205, DP 17
  against ~100x, AD 6,11, since reads cannot align across 1.3 kb. Every
  false one: DP at full depth (78-107) and alt fraction 1.0. Proposed rule
  (not applied): drop an SV call when a HC indel of 10-49 bp lies within
  60 bp with DP >= 50% of the isolate's median depth and alt fraction >= 0.9
  (HC's small allele explains all reads at full depth; a large insertion
  cannot). Removes 11 false calls, 0 true, here.
- **HaplotypeCaller hand-off applied (user, 2026-10-10)**: `score.py
  --hc-dir` (rerun_compare passes P2). Drop an SV call when a PASS HC indel
  of 10-49 bp lies within 60 bp with DP >= 50% of the isolate's median SNP
  DP and alt >= 90% of AD. 44 isolates: recall unchanged everywhere; true
  calls unchanged (setE 89, novel 94). ours+dysgu AF>=0.2: novel 0.806 /
  0.879 (was 0.832), setE 0.727 / 0.848 (was 0.809); ours alone setE
  0.485 / 0.929.
- **Harmonization (user asked: should IS6110, SV and HC calls be
  harmonized?).** Production keeps three separate outputs with no
  cross-check: HC small variants -> P5 small-variant matrix; delly + dysgu
  -> P4b (merged across SV callers only, --min-sv 50) -> p5_sv_matrix.py, a
  separate table; IS6110 arm -> is6110_p5_merge.py, symbolic <INS> keys in
  P5's key space. Nothing in P4b / p5_sv_matrix refers to IS6110 or to HC
  indels, so one event can be represented more than once (an IS6110
  insertion as an IS6110 key and a dysgu INS; a 30-49 bp indel as a HC
  allele and a dysgu call, dysgu's minimum is 30 bp). Not measured yet.
- **Double representation in scale200_fix's merged VCF**
  (`analysis/harmonize/overlap.py` -> `analysis/harmonize/scale200_fix.tsv`;
  report only, no production change). 200 samples; 3,170 small indels >= 10
  bp, 752 IS6110 records, 2,029 SV INS, 4,194 SV DEL records. Counts are
  sample-level pairs (same sample carries both), by proximity and size:
  - small indel vs SV (within 60 bp, sizes within 70%): 8,629 pairs, 787
    distinct SV records (13% of SV records). By SV origin: inherited/mixed
    5,352, graph interval 3,009, called 268. So it is mostly the INHERITED
    reference-vs-H37Rv differences, written twice: as exact alleles by P4's
    composed arm (CLASS small, often >= 50 bp, PE/PPE and masked regions)
    and as SV inherited/interval records. Called SV vs HaplotypeCaller is
    only 3%.
  - small-class indels >= 50 bp: 18,988 sample calls (1,961 records), 10,932
    with an SV record within 60 bp. Not HaplotypeCaller: composed alleles.
    Some are whole IS6110 copies (e.g. h37rv:3552709).
  - IS6110: 368 pairs where a called SV INS (291 distinct records, 14% of SV
    INS; median size 234 bp, i.e. truncated IS6110 calls) sits within 200 bp
    of an IS6110 insertion the same sample carries. 1,544 sample calls of
    ~1.36 kb SV DEL graph intervals (18 records) at IS6110 sites: element
    copy losses written as SV DEL; whether the IS6110 arm also codes them
    not checked.
  Production code not changed (user). Parked in HANDOFF for the production
  update.
- **Missed true events by cause** (`phaseB/missed.py` ->
  `out/qc51/missed.tsv`; 44 isolates, all current rules; report only).
  novel 13 of 67 missed, setE 13+5=18 of 66.
  - **Truth errors, setE (8, confirmed):** 4 REPL (478-1,406 bp) and 4 INS
    (QC_4 x2, QC_6 x2, 276-670 bp). Short reads: normal depth, no clipped
    reads, no HaplotypeCaller SNP or indel inside. The isolate's assembly
    does differ from R there (only one 300 bp flank aligns), but Nanopore
    reads of the same isolates match R: mismatch inside = flanks (3-4%), no
    >= 30 bp indel, clipped long reads no more frequent than in windows
    10-20 kb away and at scattered positions. The setE assemblies are wrong
    at these 8 places.
  - **Same signature, novel (8 REPL, unconfirmed):** no long reads here.
  - **Tandem repeat expansions (12: setE 10, novel 2):** 2,397,748 (6 setE +
    TB3396; R has 2 copies of a 58 bp unit, isolates +2 units = +116 bp),
    801,326 x2, 4,348,983 x2, N1177 1,840,217. Our caller keeps a left-clip
    cluster at the array start A and a right-clip cluster at its end B
    (68-110 bp apart) but does not pair them: INS needs the two within 20
    bp, tandem_dup needs split reads, and the clipped bases continue one
    unit back inside the array, not at A. Result BND or an INS at B with
    the wrong size; dysgu calls ~half size at B. Depth across the array
    does not give the size reliably (estimates 19-166 vs truth 53-116).
  - **N1177 IS6110 at 1,312:** both clip clusters at 1,312/1,314, but split
    reads go to another IS6110 copy and our caller formed a 1.3 Mb DEL
    (1,315-1,298,112) instead of an INS; dysgu lowProb.
  - TB3334 2.6 kb INS (four breakpoints, dysgu ~330-365 bp pieces); TB3169
    truth INV (small inverted segment) called INS by both.
  So on events the reads support (truth errors removed), setE misses 10 of
  58, all tandem expansions; novel misses 5 of 59 (8 suspect REPL removed).
- **setE defined**: Marin et al. assembly set E (`MtbSetE`), 17 isolates
  from one study (ENA PRJNA720906, "Short- and long-read WGS of M.
  tuberculosis"), Illumina NextSeq 500 + Nanopore on the same BioSample;
  Phase B group `setE_polished` (assembly = Marin's genome for the isolate).
  A setE event = a structural difference >= 50 bp between that assembly and
  the isolate's matched reference R (truth.py: indels in alignments, or
  DEL/INS/REPL/INV/REARR between alignments), scored when callable (not
  within 50 bp of R's uncallable mask): 66 in the 15 scored isolates.
- **Caller v5 + truth flags (user approved all three, 2026-10-10).**
  - `breakpoint_caller.py` v5: (1) tandem expansion: left-clip cluster at A,
    right-clip cluster at B (A+20 < B <= A+500), not otherwise paired, most
    clipped 15-mers on both sides found in R between A and B -> INS at A,
    size unknown (needs --ref; 11_call_score.sbatch passes R); (2) split-read
    DEL > 5 kb kept only with depth inside < 30% of median, otherwise the
    cluster is left free (N1177's IS6110 had become a 1.3 Mb DEL); (3) found
    in testing: the split-read tandem_dup INS is now placed at the
    duplicated segment's start (left-aligned, as VCF and the truth set do),
    not its end, which had put 801,326 and 1,840,217 72 / 68 bp off.
    Test on 702_06 and N1177: IS6110 at 1,312 now INS; 2,397,739 called;
    three other tandem calls in 702_06 to be scored.
  - `phaseB/truth_flags.tsv`: 8 setE events refuted (removed from truth), 8
    novel REPL unconfirmed (kept; recall_confirmed leaves them out).
    `score.py --truth-flags`, rerun_compare passes it.
  - v4 outputs saved in `phaseB/out/qc51_v4`; rerun 51952533.
- **v5 result** (51952533, 37 min, ~3.7 CPU-h). 44 isolates, all rules,
  refuted truth removed. recall / precision (lenient = score.py's rule, a
  call within 50 bp of a true event counts true whatever its type):
    novel: ours 0.687 / 0.929 (v4 0.537 / 0.854); ours+dysgu AF>=0.2 0.836 /
           0.913 (v4 0.806 / 0.879); excl. 8 unconfirmed REPL: recall 0.949
    setE (58 events): ours 0.776 / 0.981 (v4 0.485 / 0.929, 66 events);
           ours+dysgu AF>=0.2 1.000 / 0.880 (v4 0.727 / 0.848)
    controls: ours 1 false call, dysgu_af20 0 (unchanged).
  Strict precision (type compatible and matches()): ours+dysgu_af20 novel
  0.766 -> 0.800, setE 0.705 -> 0.744. **All precision figures reported
  before today were lenient**, ~0.13 above strict.
  Call diff v5 vs v4 (44 isolates, typed): gained 97 true / 12 false, lost
  52 true / 23 false (many are the same event re-placed). Tandem expansion
  calls: 2397739 x6, 1897936 x4, 422355 x3, 4121239 x3, 2991082 x2,
  4348983, 3236299 x2, 3593748 true; 3943165 x2, 79189 false. Left-aligned
  tandem_dup: 801326 x2, 1840217, 28469 x2, 2076726 x2, ... now true; 5
  false (TB3054 x2, TB3237, TB3334 x2). Sizes stay unknown or under-
  estimated (matched by position).
  **New problem**: where the long-DEL guard freed an IS6110 insertion, the
  local assembly step (assemble_candidates.py) now writes a giant event
  from the contig: DEL 449 kb-2.3 Mb or INV 171 kb-3 Mb at IS6110 insertion
  sites in 15 isolates (contig's IS6110 aligns to another copy). Most are
  in uncallable sequence, so scores barely move (2 callable >100 kb calls,
  as in v4), but N1177's IS6110 at 1,312 (caller: INS, correct) was
  replaced by assembly with DEL 1,315-1,298,112. Needs the same depth guard
  in assemble_candidates.py.
  Remaining misses: setE 0 of 58; novel 11 of 67 = 8 unconfirmed REPL,
  N1177 IS6110 (assembly), TB3334 2.6 kb INS, TB3169 INV typed INS.
- **Assembly guard + strict precision (user approved, 2026-10-10).**
  `assemble_candidates.py`: a contig event (DEL/INV/REARR/REPL) spanning >
  5 kb is kept only if it is a DEL with depth inside < 30% of the median;
  if nothing is left the window keeps the caller's calls ("long contig
  event dropped"). `score.py`: precision is now strict (compatible type and
  matches()); precision_lenient alongside; rerun_compare prints both.
  Strict on v5 outputs: ours+dysgu_af20 novel 0.836 / 0.800, setE 1.000 /
  0.744; ours alone novel 0.687 / 0.821, setE 0.776 / 0.870. v5 outputs
  saved in `phaseB/out/qc51_v5`; rerun 51961508.
- **Assembly guard result** (51961508, 35 min, ~3.5 CPU-h). 24 windows
  dropped a long contig event and kept the caller's calls; calls > 100 kb
  67 -> 33. N1177's IS6110 at 1,312 recovered (INS). Strict precision /
  recall, 44 isolates:
    novel: ours 0.761 / 0.836 (v5 0.687 / 0.821); ours+dysgu_af20 0.851 /
           0.808 (v5 0.836 / 0.800); excl. unconfirmed 0.966
    setE:  ours 1.000 / 0.896 (v5 0.776 / 0.870); ours+dysgu_af20 1.000 /
           0.769 (v5 1.000 / 0.744)
    controls: ours 1 false call, dysgu_af20 0.
  Misses: setE 0; novel 10 = 8 unconfirmed REPL, TB3334 2.6 kb INS, TB3169
  INV typed INS. In setE, dysgu now adds only false calls to ours (ours alone
  1.000 / 0.896 vs combo 1.000 / 0.769).
  Still > 100 kb: 23 INV + 5 REARR from the caller's split-read rules
  (windows > 20 kb are not assembled; e.g. INV 1,980,173-2,363,063 in the 6
  setE isolates on GCF_001275565, both ends at IS6110 copies) and 5
  assembled INS of ~430 kb (contig tandem-dup reading). Same cause: reads
  in a repeat element linking two copies. Not yet guarded.

## 2026-10-11

- **383 kb INV 1,980,173-2,363,063 (GCF_001275565) checked with Nanopore**
  (user: check it first). `novel_events/longread/lr_inversion.py`; output
  `longread/out/lr_inversion_1980173.txt`. 4 of the 8 calling isolates have
  long reads (3003_06, 696_05, 702_06, 8651_04). Not an inversion:
    - end A: 0 reads in R's order across 1,980,173 (+-2 kb); 53-121 reads
      per isolate carry a 2.6 kb deletion at ~1,978,700 (R's IS6110 copy at
      1,978,793-1,980,148 plus ~1.2 kb of flank is absent);
    - end B: 46-63 reads per isolate carry a ~1.3-1.4 kb insertion at
      ~2,363,100 (a new IS6110 copy);
    - the "inverted" reads (28-72 per isolate) are reads through the new
      copy at B whose IS6110 part aligns, reversed, to R's copy at A
      (segments 1,978,792-1,980,173) - a repeat-placement effect, not an
      inversion junction.
  Matches the Marin assemblies exactly (truth: DEL 1,978,689-1,981,303 2,615
  bp; INS 2,363,130 1,358 bp; in all 6 setE isolates on this R; assemblies
  colinear with R 1.64-3.93 Mb). The short-read caller turned these into INV
  + BND 1,981,304 + DEL 2,361,348-2,362,692 (1,345 bp; R's copy at B is
  still present in long reads and assemblies) + DEL 2,362,820-2,363,129 +
  BND 2,363,132. None of this is scored: the uncallable mask covers both
  IS6110 copies, so the true events and the wrong calls are both set aside.
  So setE recall 1.000 says nothing about IS6110-copy loss/gain in repeat
  windows; these calls would go to production as they are.
- **Uncallable events now scored (report only; user: go ahead with both).**
  `score.py`: rows "uncallable:<all|repeat_copy|other>:<stratum>" for the
  truth events and calls the callable rows set aside; repeat_copy = touching
  a repeat-element copy (mask paralog stretches merged across <= 50 bp, >= 1
  kb: IS6110 and other multi-copy elements; breakpoint_caller.repeat_copies).
  `rerun_compare.py` adds groups <grp>:uncallable / :repeat_copy /
  :other_masked. Callable rows unchanged. On the assembly-guard outputs
  (proto_asm; recall / precision):
    setE: callable 58 events 1.000 / 0.896; uncallable 497 events 0.229 /
          0.663 (repeat_copy 148: 0.426 / 0.692; other 349: 0.146 / 0.630)
    novel: callable 67 0.761 / 0.836; uncallable 361 0.319 / 0.592
          (repeat_copy 113: 0.416 / 0.676; other 248: 0.274 / 0.545)
  So the callable scores cover ~10% (setE) and ~16% (novel) of true events.
- **Caller v6: repeat-copy join rule** (`breakpoint_caller.py --mask-bed`).
  Before it, all 52 of the caller's INV/REARR on Phase B were false (truth
  has one real INV, TB3169 1.2 Mb, never called as INV: its partner end is
  in a repeat at MAPQ < 20); 11 had an end within 100 bp of a repeat copy.
  Rule: no INV/REARR when either end is within 100 bp of a copy; cluster left
  free (INS pairing or BND "repeat_copy_join"). Not applied to DEL, against
  the approved wording: at copies 35 clip DELs are real copy losses and 16
  false, and depth inside cannot separate them (0.6-1.6 for real losses:
  reads from other copies fill a lost copy). The other 41 false INV/REARR
  (ends in tandem repeats, short paralogs, or unique sequence; often only
  2-9 of 10-60 clipped reads with the partner) are not addressed.
  Snapshot of the previous outputs: phaseB/out/qc51_v5g. Rerun 52093842.
- **Caller v6 result** (52093842, 40 min, 1.1 CPU-h). INV/REARR from the
  caller 52 -> 39; final calls > 100 kb 28 -> 18. The 383 kb INV is gone in
  all 8 isolates (BND "repeat_copy_join" at 1,980,173; 3 isolates now get
  DEL 1,980,174-1,981,303, which matches the true 2.6 kb copy loss by its
  end). Side effect: the INV had merged everything from 1.98-2.36 Mb into
  one window > 20 kb, never assembled; those candidates are now assembled.
  Scores (proto_asm, recall / precision):
    setE callable 1.000 / 0.892 (was 0.896: 2 duplicate true DELs at
         3,118,5xx in QC_1/QC_6 merged into one each; no new false calls)
    setE repeat_copy 0.459 / 0.773 (was 0.426 / 0.692); other masked 0.149
         / 0.619
    novel callable 0.761 / 0.836 (same); repeat_copy 0.425 / 0.701 (was
         0.416 / 0.676)
    controls: 1 false call (same).
  New: QC_1's true 1,359 bp IS6110 INS at 3,189,869 now assembled as INS
  of 76,741 bp (same contig reading as the ~430 kb INS: a contig gap read
  as a tandem duplication). Assembled INS > 20 kb: 6 (QC_10 429 kb, QC_1 77
  kb, TB2780 269 kb, TB3169 963 kb and 473 kb, TB3334 1.7 Mb). Not fixed.
  Previous outputs: phaseB/out/qc51_v5g.
- **Contig-length guard result** (52099486, 42 min, 1.1 CPU-h). Assembled
  INS > 20 kb: 6 -> 0; QC_1's IS6110 INS back at 3,189,869 (caller's INS,
  size unknown, "long contig event dropped"). Scores unchanged from v6
  (setE callable 1.000 / 0.892; novel 0.761 / 0.836; controls 1). Calls >
  100 kb left: 18, all the caller's INV (13) / REARR (5). Snapshot of v6
  outputs: phaseB/out/qc51_v6.
- User (2026-10-11): rules getting complicated, not robust; goal is reliable
  cross-sample calls for convergence testing and GWAS, staying within
  short-read limits. Proposal accepted: robust tier (DEL/INS in callable
  sequence, both breakpoints), everything else flagged and kept out; IS6110
  from its own caller; per-sample present/absent/no-call in a common frame.
  User: check heritability on the tree (outbreak / close relatives should
  share an event). Method to be written up before running.
- **Noise from reads mapped to the isolate's own genome**
  (`phaseB/self_noise.py`, out/qc51/self_noise.tsv; 11 in_panel isolates
  whose matched R is their own panel genome; raw output, no filters):
  callable / masked, totals over 11 isolates
    HC PASS SNPs           30 / 25   (20 + 16 in TB3251, 6 + 6 in
                                      M0010874_7; 0-2 in the other 9)
    HC indels 1-9 bp        7 / 46
    HC indels 10-49 bp      0 / 0
    our SV, typed           5 / 6    (callable: 4 tatC-helY 62-68 bp DEL,
                                      removed by the genotype rule; 1 CIGAR
                                      DEL 52 bp, 02_R1708; masked: all 6 in
                                      M0010874_7)
    our SV, BND            10 / 33   (27 + 6 in M0010874_7)
    dysgu PASS             29 / 23
    dysgu AF >= 0.2         3 / 8
  So in callable sequence, after the existing filters: ~0.1 false SV per
  isolate (ours), 0-2 SNPs in 9 of 11 isolates. M0010874_7 (P1 distance 103)
  and TB3251 look like reads and panel genome not being the same DNA.
- **Tree check method drafted** (user asked): analysis/tree_check/METHOD.md. scale200_fix (200 isolates, tree + p5 sv_matrix 5,845 SV keys in H37Rv frame); close pairs <= 12 SNPs (+13-50, 51-200); per-pair SV disagreement vs SNP distance, noise = intercept; split by same/different matched reference, inherited vs called, callable/masked, type/size/source/tier, IS6110; parsimony gains + lone-isolate calls; lineage RDs as positive control. Awaiting user decisions (threshold, cohort, inherited).
- **User (2026-10-11): before scale-up, QC samples with a more rigorous fold-back scan and harmonize HC / SV / IS6110 calls.** Plan written: analysis/qc_harmonize/PLAN.md (not run). ArtifactsFinder (asked earlier, not answered then): reference-side blacklist generator (exact inverted repeats >= 8 bp duplex, spacer >= 5; palindromes with 1 mismatch; hg19/hg38), no per-sample rate; complements split_reads.awk. Fold-back v2: all SA entries, MAPQ >= 20, plus reverse-complement soft clips >= 15 bp within 500 bp (current scan finds ~70% of synthetic junctions), recurrent reference sites separated; validate on foldsim + Phase B. Harmonization at merge stage (merge_cohort_vcf.py; ~1 billing-hour rerun vs 362 for P5 states): HC owns < 50 bp, SV block owns >= 50 bp (one record, exact allele in INFO), IS6110 arm owns insertions (any-size SV INS within 100 bp dropped), copy losses relabelled is6110 + flagged, SOURCES= tag. Awaiting decisions.
- **User approved (2026-10-11): fold-back v2 validation, harmonization rules
  1-5 in production code with a merge-only rerun of scale200_fix, and IS6110
  copy losses relabelled + flagged (LOWCONF).**
- **Fold-back scan v2** (`foldback_qc/foldback_v2.py`; 03_v2_validate.sbatch,
  52142746 after a refmap-column fix to 52142513; 80 s, 8.4 CPU-min; 25% of
  reads by name hash). Per read: split = any SA entry MAPQ >= 20, same
  contig, opposite strand, within 1 kb; clip = quality-trimmed soft clip >=
  15 bp whose 15 bases next to the alignment end, reverse-complemented,
  occur within 500 bp. Reads: primary, MAPQ >= 20, not duplicates. A 10 bp
  clip minimum was tested and rejected: +5-10 points sensitivity but the
  clean-sample background nearly doubled (0.12 -> 0.20 per 1,000).
  `validate_v2.py` -> `v2/validation.tsv`:
  - synthetic (3 genomes x ir 0.5-15%, random 2-15%): sensitivity v2
    0.74-0.80 (ir) / 0.85-0.88 (random) vs old 0.69-0.74 / 0.74-0.75;
    background at 0: v2 0.0. Reads whose junction is < 15 bp from an end
    cannot be seen by either method.
  - Phase B raw reads: v2 = ~1.27 x old, Spearman 0.991 with the old count;
    clean <= 1.55 (setE max, 8651_04), chimeric >= 6.45 (N0155): same
    isolates flagged. Old cutoff 2 corresponds to v2 ~2.5.
  - Spearman with false SV calls in callable sequence (62 isolates): v2
    0.49 (dysgu PASS) / 0.64 (ours), old 0.50 / 0.65: no gain.
    clip_other (clips matching nothing nearby) 0.53 / 0.17, but it also
    rises with distance to the reference (02_R0894, matched to another
    isolate's genome: 5.9), so it is not a clean library metric.
  - Recurrent reference sites (same reference, >= 3 isolates, >= 3
    carriers): <= 1% of junctions except TB3091 16%, TB3396 40% of tiny
    counts. Reference-driven fold-backs are negligible here.
  - Within the clean range, false calls do not follow the fold-back rate
    (setE 4-24 dysgu false calls at 0.48-1.55); QC_7 (mixed, excluded) has
    the most. Fold-back is one noise source, mixing another.
  Conclusion: v2 is somewhat more sensitive and better defined, but ranks
  isolates the same; the earlier QC decisions stand.
- **Harmonization in production** (`bin/harmonize_records.py`, called by
  `merge_cohort_vcf.py` before sorting; `--is6110-gff` passed by
  `p5_finish.sh` from `<build>/assets/is6110_elements.gff`, 16 H37Rv copies;
  `--no-harmonize` to skip; header: SOURCES, SMALLALLELE, IS6110LOSS,
  LOWCONF, ST=DUP). Tests: `tests/test_harmonize.py` (16), plus
  run_tests.py (369) and the audit suites, all pass. Merge-only rerun of
  scale200_fix (runner `--only p5vcf`; 52143216, then 52143721 and 52144022
  after two fixes found by the before/after check; ~45 s each; vcf_gate
  passed). Previous VCF kept as
  `scale200_fix/p5/merged.preharmonize_20261011.vcf.gz`.
  - Fix 1: rule 2's fill of no-call cells overwrote DUP cells set by rules
    1/3 (2 cells went ALT -> REF); now DUP cells are left alone.
  - Fix 2: rules 3 and 4 did not see small alleles that rule 2 relabelled
    CLASS=sv, so two H37Rv copy losses (3,551,230, 169 carriers; 2,635,577)
    and 16 IS6110-sized exact insertions were missed. Rule 3 now also takes
    relabelled insertions of 1,300-1,420 bp (exact sequence, so not
    truncated; shorter ones near IS6110 are other events); rule 4 takes
    relabelled deletions.
  Result (`analysis/harmonize/compare_scale200_fix.txt`):
    records small 79,007 -> 77,104 (400 merged into SV records, 1,503
    relabelled CLASS=sv); sv 6,223 -> 7,374 (319 caller records dropped,
    all calls absorbed); is6110 copy losses 0 -> 23 (LOWCONF).
    Cells: 30,605 SV no-calls now REF, 854 ABSENT, 476 ALT (from the merged
    small allele); 135 SV ALT -> DUP. 944 SV INS ALT cells where the IS6110
    arm says REF/ABSENT, kept and counted.
    Drug-resistance windows (gyrA, rpoB, rpoC, rpsL, rrs, fabG1-inhA, katG,
    pncA, eis, embB, gid, ethA): identical records and genotypes.
  overlap.py after (`scale200_fix_harmonized.tsv`) vs before:
    small_sv pairs 8,629 -> 64 (SV records 787 -> 15); small_ge50 with an SV
    record 10,932 -> 0; IS6110-near SV INS records 291 -> 103 (69 under 200
    bp, 29 other sizes: flank variants, not the element; 5 IS-sized where
    the arm disagrees or > 100 bp away); IS6110-sized DEL records 18 -> 14,
    pairs 1,544 -> 203: the rest are not within 50 bp of an H37Rv copy's
    ends (e.g. svi:DEL:2635522:1408 starts 54 bp before the copy at
    2,635,577, 179 carriers) or lie at copies H37Rv lacks.
  Not rerun: P6 annotation and the association tail of scale200_fix;
  scale200 and gwas1000 merged VCFs predate this code.
- **Open decision**: copy-loss tolerance 50 bp (approved) misses
  svi:DEL:2635522:1408 by 4 bp; 100 bp would take it.
- **Copy-loss tolerance widened to 100 bp (user).** On scale200_fix it adds
  2 deletions (svi:DEL:2635522:1408, 179 carriers; h37rv:3551190 small
  allele, 2 carriers), both removing the whole copy plus 40-80 bp of flank;
  no deletion within 300 bp of the 16 H37Rv copies removes only part of one.
  Merge rerun 52146472 (49 s; a watcher was first started on a mistyped id,
  52146893, and replaced): 25 IS6110LOSS records, DR windows unchanged,
  vcf_gate passed. Note: most of these "losses" have 159-199 of 200
  carriers, i.e. H37Rv-specific copies the isolates never had.
- **Fold-back calibration, step 1 (user: yes)**: `foldback_qc/panel_illumina.py`
  -> `panel_illumina.tsv` (NCBI assembly summary for BioSamples, ENA portal
  API metadata, no reads). Of 333 panel genomes, 69 have Illumina WGS runs
  of the same BioSample (68 paired); 13 of them are already Phase B in-panel
  controls; 2 are among the 54,461 production CRAMs. 56 new candidates,
  77 Gbases in all; instruments MiSeq 22, HiSeq 2500 8, HiSeq 4000 8,
  NextSeq 2000 8, NextSeq 550 3, others 7; largest studies PRJDB10923 (9),
  PRJNA1480832 (8), PRJNA1113402 (6), PRJEB49093 (5).
  Earlier answer (user asked about scanning Peter's ~54K): the 500 random
  production samples already give the distribution (<=0.3 50%, 0.3-1 10%,
  1-2 8%, 2-5 11%, >5 22%, old count); what is missing is the false-call
  rate between 0.4 and 6.4 per 1,000, where no control exists.
- **Fold-back calibration step 2 (user: straight to full reads, all 54;
  2026-10-11).** User's guidance (memory `sample-quality-from-wild-data`):
  learn sample filters from the distribution in real public data, not take
  all data because it exists; SVs need good reads (read length, depth,
  maybe insert size); a reusable junk-sample filter helps other projects.
  Parked for later: Mutect2 / non-haploid calling for within-host variants,
  and the effect of library fragmentation (sonication vs enzymatic) on low-
  frequency variants. Plan and limits: `foldback_qc/CALIBRATION_PLAN.md`
  (polishing bias, genome-vs-reads DNA mismatch, fold-back confounded with
  library era, unknown count in the 1-6 gap).
  - `selfcal/make_runs.py` -> `selfcal/genomes.tsv`: 54 genomes, one run each
    (largest paired Illumina WGS run of the BioSample, runs not pooled),
    35.6 GB FASTQ, with ENA library fields. Count corrected from 56 (one
    single-end set; Phase B controls counted by genome).
  - Cohort `selfcal54` (registry row; `runroot/refbias/selfcal54` ->
    netscratch), P2 only, reference pinned to the own genome
    (`selfcal/make_refmap.py`, run after QC).
  - `phaseB/08_qc.sbatch`, `09_kraken.sbatch`: QC_TABLE / QC_READS env
    overrides (defaults unchanged). `phaseB/self_noise.py`: `--all-own
    --refmap --p2 --asm-dir --out` (default output unchanged, checked).
  - Jobs: download 52206224 (54, 8 at a time), masks 52206464, QC 52206708
    (afterany download, 20 at a time), kraken2 52206709. Then by hand:
    make_qc_cohort.py, make_refmap.py, runner `--only p2`, then
    `selfcal/03_post.sbatch` (fold-back v2 at 25%, samtools stats, caller v6
    + assembly, self_noise, `collect.py` -> `selfcal/out/selfcal_table.tsv`;
    fold-back and stats also for the 11 Phase B own-genome controls on QC'd
    BAMs). Estimate about 100-110 billing-hours.
- **Library type from adapter read-through (user: all three sets).**
  `foldback_qc/adapters/adapter_count.py` + `01_count.sbatch` (jobs
  52211341 phaseB 32 s, 52211344 production 97 s, 52211345 selfcal after the
  download). Results in `adapters/RESULTS.md`:
  - Phase B: all 11 chimeric controls are Tn5 libraries (10 with very short
    inserts); clean are TruSeq set B, unclear set A, and Tn5 set E (0.33-1.55).
  - Production (500): fold-back > 5 in 26% of Tn5 and 23% of TruSeq. Within
    TruSeq, fold-back tracks short inserts (rho 0.61 adapter, 0.75 trimmed-
    length fraction); within Tn5, weakly (0.07, 0.29).
  - Adapter columns added to `selfcal/collect.py`.
- **Fold-back effect on SV calls: part C started, D queued (user, 2026-10-11:
  "start C now and queue D after"; also wants the major pipeline issues
  wrapped up).** Plan `foldback_qc/FOLDBACK_SV_EFFECT_PLAN.md`. User context
  (memory `within-host-tn5-short-insert`): their within-host Tn5 short-insert
  samples could not be called below 10% AF.
  - Part C isolates: 9 chimeric Phase B controls whose marinB63 reference is
    their own genome (N0004, N0054, N0145, N0155, N1176, N1202, N1272,
    N1274, RW_TB008; plan said 8). Cohort `fbremove` (registry row,
    `runroot/refbias/fbremove` -> netscratch). QC 52219006, kraken2 52219010.
  - `foldback_v2.py --names` (fold-back read names; test on N0004 5%:
    1,808 names, rate 7.55).
  - `removal/make_sets.py` (stage 1 `_all`, stage 2 adds `_fb50`, `_fb100`,
    `_rand`), `01_sets.sbatch` (names, read sets by samtools view -N ^),
    `02_post.sbatch` (fold-back, stats, caller v6 + assembly, self_noise).
    Order by hand: make_qc_cohort -> make_sets 1 -> runner p2 -> 01_sets ->
    runner p2 -> 02_post.
  - Part D after selfcal54: reads to the nearest other panel genome, truth
    from assembly differences.
- **Part C result (read removal), `foldback_qc/removal/RESULTS.md`.** 9
  chimeric isolates, 4 read sets each (jobs: P2 52222277 + 52227868, sets
  52222307, post 52227921). Removing fold-back pairs (1.3-10.7% of pairs)
  removes 95% of our false calls and 90% of dysgu's (98-100% at AF >= 0.2);
  random removal of the same number removes none (+13-21%, unexplained,
  likely fraction thresholds). Halving fold-backs cuts ours by 72%
  (non-linear). Chimeric false calls are mostly AF < 0.2 for ours, but
  dysgu has 1,703 at AF >= 0.2. delly keeps 6,373 after removal (its own
  noise). Junction attribution at 20 bp is uninformative in these libraries
  (chance 54-75%). RW_TB008 (10.9, few false calls) vs N0145 (13.0, 95)
  shows the rate alone is not the whole story.
- **delly dropped from production (user, 2026-10-11).** Basis: read-removal
  experiment (6,373 delly false calls left with every fold-back read
  removed). Changes:
  - `p2_call.sh`: no delly; dysgu now required (missing dysgu fails the
    sample); only the dysgu VCF stamped.
  - `p4b_place_sv.sh` / `.py`: dysgu only; `--delly` removed, so an old
    P2 folder's delly VCF is not read.
  - `p2_summary.py`: dysgu column, a sample without dysgu is incomplete.
  - `merge_cohort_vcf.py` help text.
  - `config`: MTB_DELLY_ENV removed.
  - P5 QUAL banding unchanged (per caller, from the calibration table;
    dysgu already calibrated).
  Tests updated (generic VCFs passed as --dysgu; p2_summary fixture writes
  a dysgu VCF; new test: --delly refused). All suites pass (run_tests 370,
  harmonize 16, 13 audit suites). First real use: selfcal54 P2.
  Existing cohorts' SV records still include delly calls until rerun
  (wrap-up item 5). dysgu AF >= 0.2 is not yet a production rule (item 3).
  User will ask the CRAM owner before any copy (purge about 2026-10-26).
- **Fold-back mechanism explainer (user asked for figures).** Page
  https://claude.ai/artifact/6mTMff14Ts66rLgU8maQdu (private; source
  `foldback_qc/explainer/foldback_explainer.html`). New measurements on the
  9 removal isolates (all-reads sets):
  - Calls removed by fold-back removal are almost all small INV (ours 753,
    dysgu 3,451) and our unpaired BND (721), mostly 50-500 bp.
  - Read geometry: the two pieces of a split fold-back read overlap
    (hairpin) in 59-71% of reads in most libraries; in N1176 they are
    usually 100-600 bp apart. Spaced fraction did not predict calls beyond
    the rate (n = 9).
  - Junction recurrence predicts false calls better than the rate:
    junctions with >= 5 reads vs dysgu AF >= 0.2, Spearman 0.95 (rate 0.88);
    ours 0.92 (0.77). RW_TB008: 478 such junctions, 1 call; N0145: 1,976, 197.
  - Recurrence comes from PCR copies (P2 marks no duplicates; distinct
    molecules at hotspots 42% N0004, 72% N0145, 82% N1176, 97% RW_TB008)
    and sequence hotspots: inverted repeats (10-mer stem within 120 bp) at
    14-39% of hotspots vs 7-8% of random sites.
  - Open: duplicate marking in P2 would remove part of the support (not
    tested); junction recurrence as a QC metric, to test on selfcal54.

## 2026-10-11 — duplicate-marking test on the 9 chimeric libraries
- User: "test duplicate marking on the 9 libraries".
- Plan: `analysis/foldback_qc/removal/DEDUP_PLAN.md`. On each `_all` P2 BAM:
  - samtools markdup -S (flags supplementary pieces too);
  - `_dedup` removes the flagged reads; `_randdup` removes the same number of
    random pairs (the depth control); `_dupmark` keeps them flagged, to test
    whether dysgu ignores flagged reads;
  - dysgu on all three; our caller v6 + assembly on dedup and randdup;
  - duplicate share among fold-back pairs vs the other pairs.
  - HaplotypeCaller not rerun.
- Job 52282456 (array 0-8, `03_dedup.sbatch`), about 5 billing-hours.
- Result: duplicate marking does not fix fold-back false calls. dysgu false
  calls rose by 22% with duplicates removed (random removal of the same
  number of pairs: 11%); our typed calls fell 20% (random: 4%), still 33x
  the fold-back-removed level.
  - Fold-back pairs are duplicated no more often than other pairs (N0004
    32% vs 43%), so dedup raises the chimeric share.
  - dysgu already ignores flagged duplicates (dupmark = dedup).
  - Corrected the explainer page (v4): hotspots are sequence-driven, not
    PCR copies.
  - Details in `foldback_qc/removal/RESULTS.md`.

## 2026-10-11 — selfcal54: QC done, P2 launched
- Downloads complete (54). Masks: all 54 present. The masks job reported
  FAILED only because 5 genomes already had masks from 2026-10-09, so its
  `mv` of a missing .tmp failed; callable_mask.py is unchanged since then.
- QC: 1 TB-Profiler crash (a temp file vanished, ps_GCF_042982005); rerun
  52293507 passed at 77x.
- 39 of 54 pass QC (the plan guessed about 44). Failures:
  - depth < 60x: 14, of which 5 are PRJEB49093 runs at exactly 59.5x (the
    submitters' own 60x) and 2 lose more than 80% of reads in fastp
    (GCF_037101195 155x to 28x, GCF_040209395 184x to 25x);
  - MTBC < 85%: 1 (GCF_000023625, KZN MDR 2009, 80.2%).
- P2 (dysgu only, first production run without delly): array 52294010,
  summary 52294011; post (`selfcal/03_post.sbatch`) 52294042, after the
  summary.
- selfcal54 calibration (`foldback_qc/selfcal/RESULTS.md`, `calib.py`):
  - 59 points (39 selfcal, 11 clean, 9 chimeric); clean y (ours typed +
    dysgu AF >= 0.2) mean 0.73, T = 2.54.
  - Pre-specified: rule 1 rate <= 5.80 (v2) keeps 29/39 selfcal; rule 2
    (plus median insert) adds nothing, and in leave-one-study-out it let
    N0145 / N0155 through. Insert size rejected (Spearman -0.02).
  - The mean-based criterion is too lenient: all 3 isolates at 3-6 exceed T.
    By band, < 1: max y 5; 1-3: median 3.5, one 14 (GCF_035581225, 62x,
    identity unchecked); 3-6: 20-21. A per-isolate criterion points to a
    cut-off of 1-3 (v2); the user decides.
  - Exploratory: recurrence (junctions >= 5 reads) rho 0.73 vs rate 0.75;
    explains PRJDB10923 2x300 (rate 11-28, y 4-8). Not better as a single
    rule yet.
  - Mid-range false calls: MiSeq Tn5 short-insert libraries (median insert
    130-165 bp), mostly from our caller.
  - Negative binomial fit did not converge.
