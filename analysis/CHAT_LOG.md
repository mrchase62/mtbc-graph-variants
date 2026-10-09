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
