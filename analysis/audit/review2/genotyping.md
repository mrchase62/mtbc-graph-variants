# Review 2: genotyping (P4, P4b, P5, SV genotyping, merge)

Reviewer scope: `bin/p4_place.{py,sh}`, `p4b_place_sv.{py,sh}`, `p5_states.py`,
`p5_states_io.py`, `p5_sv_genotype.py`, `p5_sv_matrix.py`, `sv_intervals.py`,
`proj_store.py`, `mtb_norm.py`, `merge_cohort_vcf.py` (SV, accessory and dedup
parts), `p5_merge.sh`, `p5_finish.sh`, `p5_svgt.sh`, `p5_sanity.py`,
`p5_matrix.py`, `sv2frame/bin/*`, on `audit-fixes` (ffd7bb9).

Method: code read on a `git archive` export of `audit-fixes` and of `main`.
Real-data checks used production scale200 inputs (read-only), with all outputs
in `/tmp/claude-12043/review2_genotyping/`:

- The odgi outputs already on disk: `work/scale200_p4/*.{h37rv,node}.pos` for
  200 samples, and the build projection store `build/7713a8d71d8e/proj/*`.
  Both carry odgi's own strand flag in column 5.
- Node tables regenerated with the branch's `node_path_membership.py` from the
  CX333 GFA: 95,828 nodes and 1,379,748 pairs, matching the fix group's figures.
- Old and new `p5_states.py` and `p4_place.py` run on three samples:
  - SAMEA111556136, reference GCF_014900175 (forward);
  - SAMN08612917, reference GCF_001870145 (stored flipped in the panel);
  - SAMEA1403838, reference GCF_000193185 (13,508 inverted-step keys).

  The old runs reproduce production exactly: 0 cells differ, and the placed
  tables are byte-identical.
- Changed cells checked against R's sequence and against the sample's own P2
  reads (`work/scale200_p2/*.bam`).

Scripts: `strandcheck.py`, `nodecheck.py`, `mm*.py`, `absent_check.py`,
`trace.py`, `readcheck2.py`, `svdup.py`.

---

## R2-GENO-1 (MEDIUM): `deleted_in_ref` calls ABSENT where the H37Rv base plainly exists in R

**Where:** `bin/p5_states.py:260-362` (`deleted_in_ref`), reached from `:753`.

**What is wrong.** The ABSENT test anchors on 12-mers found with
`h37.find/rfind` anywhere within ±50 kb of p. It does not require the 12-mer
to be unique or near its expected offset. It then accepts an empty R middle
(`rm` of length 0) as a deletion: with n = 0, `best = (0, 0)`, and
`c <= poff < m` holds for any p inside `hm`. So in GC-rich repeats
(PE_PGRS / PPE), two adjacent R k-mers that match two different repeat copies
in H37Rv "prove" a deletion of everything between those copies. One traced
case anchors 15 kb from p.

Traces (`trace.py`, GCF_000193185):

```
839092 (t 3574578, dist 153): anchor (357, 839160, after) ... other anchor (354, 839082)
                               rm 0  hm 75  poff 8   -> True
1095330 (t 3318739, dist 7):  anchor (218, 1080787, before) ... other anchor (220, 1095336)
                               rm 0  hm 14547 poff 14540 -> True
```

In both, the H37Rv 25-mer centred on p occurs exactly once in R, a few bases
from t. R carries p, and the sample's reads carry the H37Rv REF haplotype
(73 reads REF, 0 ALT at 1095330; 36 reads REF, 0 ALT at 964699).

**Count.** New ABSENT H37Rv-frame positions whose H37Rv 25-mer (centred on p)
occurs exactly once in R's whole genome and lies within ±(dist + 2 kb) of t,
so p certainly exists in R (`absent_check.py`):

| reference | new ABSENT positions | p exists uniquely in R | ambiguous (25-mer in several R copies) |
|---|---:|---:|---:|
| GCF_014900175 | 654 | 7 | 27 |
| GCF_001870145 (flipped) | 668 | 7 | 9 |
| GCF_000193185 | 667 | 59 (58 of its 154 `-`-strand positions) | 49 |

The false calls concentrate in the region R walks inverted. The count is per
reference, so every sample matched to GCF_000193185 inherits about 59 false
ABSENT cells.

**Effect.**
- The merged VCF writes GT=2 `*`, "Position deleted in this sample", where the
  base exists.
- The event matrix reads GT=2 as unknown, so for association this is lost
  information.
- An alignment or ancestral reconstruction that treats `*` as a gap, or a
  deletion, as a state would read a false deletion.
- The fix intended ABSENT to mean "measured deletion"; about 1% of new ABSENT
  positions for a typical reference, and 9% for GCF_000193185, are not.

**Fix.**
- Before ABSENT, test whether H37Rv's 25-mer around p (tolerating a mismatch
  at the centre) occurs in R within the window. If it does, the result is
  NOCALL (or genotyped), never ABSENT.
- In `deleted_in_ref`:
  - require each anchor k-mer to be unique in the H37Rv search window (or use
    k ≥ 20);
  - require the walked runs to be at least about 20 bp on each side;
  - reject an anchor more than `dist + reach` from p.
- Add a regression test on a tandem-repeat contraction.

## R2-GENO-2 (MEDIUM): the reciprocal-overlap rule (P4P5-5) is not applied when a caller deletion is matched to a graph interval

**Where:**
- `bin/sv_intervals.py:185`: caller deletions are dropped as "already covered" on `same_event` alone;
- `bin/merge_cohort_vcf.py:821-828` (`catalogued`): the same test decides which caller deletions the catalogue supersedes.

**What is wrong.** The fix added `overlap_frac >= 0.5` to graph-graph
clustering. Its own comment says that position tolerance alone "lets two
58 bp deletions 150 bp apart merge without sharing a base". The graph-caller
match still uses tolerance alone.

**Count.** Run on scale200 (`sv/new_iv.tsv`, `svdup.py`): of 693 caller
deletions dropped as covered, 65 share less than 50% with every matching
graph interval, and 25 share no base. They carry 164 ALT cells.

Such a deletion is neither catalogued nor written UNCATALOGUED, because the
merge skips it as `catalogued()`. Its carriers survive only if depth at a
different interval happens to call them.

**Fix.** Use the same predicate in both places:
`same_event(...) and overlap_frac(...) >= 0.5`, with `merge_cohort_vcf.catalogued`
importing it from `sv_intervals`.

## R2-GENO-3 (LOW, information loss, recorded as D3): REF→NOCALL at indel keys often has read support for REF

**Where:** `bin/p5_states.py:365-417` (`r_allele`), `:624`.

`r_allele` requires an exact match over the first d+1 bases. In a tandem array,
d runs to the end of the array, and the left 10-bp flank allows at most 1
mismatch. Any R SNP inside that window therefore gives OTHER, then NOCALL.

Per sample, 277 / 643 / 711 indel cells go REF→NOCALL. In a random 40 per
sample:

- reads in the sample's own BAM carry the discriminating H37Rv-REF haplotype
  (at least 3 reads, at most 1 for ALT) in 7 / 16 / 14 cells;
- ALT is supported in 0 / 2 / 1;
- the rest have no discriminating reads.

So about 25-35% of these new NOCALLs are decidable as REF. ABSENT→NOCALL is
similar: indels have about 20-25% REF and about 15% ALT read support; SNPs
10-30% REF. This is conservative, as D1 and D3 chose, not wrong; it is noted
so that the size of the loss is known.

**Possible improvement.** Compare R against the haplotypes allowing
mismatches outside the discriminating base(s), as `homologous` does.

## R2-GENO-4 (LOW): `p5_sv_genotype.load_projection` ignores odgi's inverted-step flag

**Where:** `bin/p5_sv_genotype.py:164-180`.

P4, P4b, IS6110 and p5_states all shift the target at column 5 `-`. The SV
genotyper does not, so its depth probes and flanks land 1 bp off at inverted
steps. For depth this is immaterial; it is noted for consistency.

## R2-GENO-5 (LOW, clarity): D41's description of node-frame alleles is inaccurate, and node keys change orientation in the new code

**Where:** `bin/p4_place.py:364`.

`on_strand` is applied to off-path (node-frame) records too, so their alleles
are in the orientation of the nearest H37Rv projection, not "each reference's
own orientation" as DECISIONS D41 states. Because that strand now includes
odgi's column 5, node-frame SNPs on locally inverted R steps change from R's
orientation to H37Rv's. Examples from SAMEA1403838:

- node:195962:0 T>C becomes A>G;
- node 198028, node 179906, node 149189 and node 149592 change the same way.

I checked them against the sequence: the new alleles and the shifted
`h37rv_pos` are the correct H37Rv homolog. The behaviour is better than
before, but:

- D41 should be corrected;
- node keys from the old and new code are not comparable, which a rerun
  handles anyway.

## R2-GENO-6 (LOW): sharded merge can drop or keep a caller INS differently from the unsharded merge

**Where:** `bin/merge_cohort_vcf.py`, `acc_emitted` and `dup_of`.

`acc_emitted` holds only this shard's level-1 loci (`cat` is filtered by
`in_shard`). A caller INS within 10 bp of a locus that sits in the next shard
is therefore kept, though the unsharded run drops it. The assembled VCF is
claimed byte-identical to the unsharded one. This is rare: it needs a locus
within 10 bp of a shard edge.

**Fix.** Build `acc_starts` from the whole catalogue, as `catalogued()`
already does.

## R2-GENO-7 (LOW, tests): the strand corrections are tested only for the forward-stored (+,-) case

`tests/test_audit_p4p5.py:200` and `tests/test_audit_leftovers.py:58-166`
exercise column 4 `+` with column 5 `-`. Two cases have no test:

- column 4 `-` (a reference stored flipped);
- (-,-), where p5_states uses t+1 and p4 uses t-1.

A sign error there would pass the suite. Add one flipped-reference fixture in
each direction.

There is also no test of `deleted_in_ref` on repeats (see R2-GENO-1).

---

## Checked and found sound

### (a) Strand corrections, verified on real data in both directions

Method: 31-mer homology between source and corrected target, flanks with at
most 1 mismatch each, tabulated against what the code chooses.

**R→H (p4_place `parse_pos_file`, p4b identical)**, 200 scale200 samples,
dist 0:

| col4 | col5 | reference | code's choice homologous | unique alternative | no homolog anywhere |
|---|---|---|---:|---:|---:|
| + | + | forward | 24,441 | 13 | 5,316 |
| + | - | forward | 165 | 0 | 62 |
| - | + | flipped | 2,511 | 0 | 509 |
| - | - | flipped | 7 / 7 | 0 | 0 |

With a ±5 bp search, every (+,-) "none" case has no homologous shift at all.
These are divergent PE/PGRS loci, not a wrong offset.

**H→R (p5_states)**, whole projection store, 160 references:

| col4 | col5 | code's choice homologous | wrong |
|---|---|---:|---:|
| + | + | 214,759 | 10 |
| + | - | 525,083 | 19 |
| - | + | 14,338 | 1 |
| - | - | 227 | 0 |

The 193 "none" cases in (-,-) have no local homology at any shift up to ±40
(paralog landings). Where homology exists, t+1 on the same strand is right
every time.

**Consistency check.** GCF_000193185 has 13,508 inverted-step keys. Under the
new allele check (`r_allele` at the corrected target), no SNP cell moved
REF→NOCALL. Had the correction been wrong, these would have shown as OTHER.

### (b) Forward node offsets

- Every P4 node record of 200 samples (57,966) was checked against the GFA
  node sequence at `forward_offset`, complemented by (walk strand XOR panel
  flip). 100% agree, for both strands, flipped and forward references, 1-bp
  and longer nodes.
- The producers all use `mtb_norm.forward_offset`:
  - p4_place (`:356`);
  - is6110_project_sites (`:304`);
  - p5_keys and p5_states `by_key` (pass P4's offset through);
  - is6110_p5_merge.
- The readers that turn a key back into a path position are:
  - p5_states node branch, `start + forward_offset(off, strand, L)`, with
    `start` the walk-entry offset from `node_path_membership.py`;
  - merge_cohort_vcf, which uses the offset only for the IS6110 POS.
- write_event_matrix, vcf_to_alignment and node_locus_from_p4 use only the
  prefix or the node id.
- L comes from `<build>/assets/node_positions.tsv`; p4_place.sh passes it, and
  p5_merge.sh requires it. A '-' walk with no length is dropped and counted.
  The table excludes H37Rv nodes, which is complete because node keys exist
  only where dist ≠ 0.

### (c) The new P5 rules

Old vs new on 3 samples; the old run matches production exactly. Transitions
at H37Rv-frame cells:

- ABSENT→NOCALL: 912 / 807 / 1,160;
- REF→NOCALL (indel): 277 / 643 / 711;
- REF→ALT (indel): 52 / 54 / 55;
- ABSENT→ALT: 19 / 18 / 15;
- ABSENT→REF: 8 / 1 / 9;
- no transitions to ABSENT.

Reads check (sampled 40 per category; `truth/reads_summary2.txt`):

**REF→ALT (P4P5-2).** The discriminating ALT haplotype is in reads in 23, 15
and 18 cells. The REF haplotype only is in 2, 4 and 1 cells (about 6%,
possible errors). The rest are ambiguous repeats or have no reads. The rule is
mostly right.

**ABSENT→REF.** All 14 checked cells have REF read support, or are genuine
reversions. For example, at 2715342 the sample's gVCF calls T→C back to
H37Rv's base. The rule is right.

**ABSENT→ALT.** 21 cells have ALT read support and 1 has REF. The rule is
mostly right.

**Retained ABSENT.** Reads mostly have no haplotype, consistent with a
deletion, apart from the R2-GENO-1 cases.

### Other items checked

- **P4P5-8** (inherited records overlapping a called record) is fixed with
  correct interval logic. The new p4 on 3 samples drops only 2-5 inherited
  records per sample beyond the old rule.
- **P4P5-4** (0/0 delly records) is fixed in `parse_sv_vcf`.
- **P4P5-3** (clip cluster): ALT only on an H37Rv-frame ALT; REF plus clips
  gives NOCALL. This matches D5.
- **P4P5-14** (depth inside a GATK-called deletion): the `gvcf_depth` span
  arithmetic is correct.
- **P4P5-11**: `qual_caller` is carried through to `p5_sv_matrix`.
- **P4P5-6** (caller INS de-duplication against IS6110 and level 1) is in
  place; see R2-GENO-6.
- **P4P5-7** (build assets):
  - graph, node tables, collapsed graph VCF and accessory catalogue are taken
    from the build;
  - `--build-id`, `--graph-vcf` and `--node-*` are required and passed by
    `p5_finish.sh`, `p5_merge.sh`, `p4_place.sh`, `p4b_place_sv.sh` and
    `p5_svgt.sh`;
  - no other chain caller of these scripts was found missing an argument.
- **Write-then-rename** in p4 and p4b.
- **p5_sanity** writes NA, not 0.
- **p5_matrix** panel AF per trimmed allele.
- **Tests:** `test_audit_p4p5`, `test_audit_cleanup` and `test_audit_leftovers`
  pass (50 tests). The P5 fixture tests behaviour (states), not just that the
  code runs.
- **Performance:** new p5_states takes 63-109 s per sample at 0.78 GB on 79k
  keys. `homologous`, `deleted_in_ref` and `r_allele` depend only on
  (reference, key), so at 10k samples they could be cached per reference.
  This is not a blocker.

## Audit findings not addressed (this area)

- **P4P5-9 (LOW):** the composed REF is still an equal-length H37Rv slice.
- **P4P5-10 (LOW):** node-frame small variants are still at POS=1.
- **P4P5-12 (LOW):** alleles are not oriented to the node's forward strand
  (D41, pending the user's decision).
- **P4P5-17 items:**
  - p5_merge.sh still passes the panel-frame anchors to `proj_store`; this is
    latent and harmless while H37Rv's panel offset is 0;
  - node-contig sort order is unchanged.
- **P4P5-5 is only partly addressed** (R2-GENO-2).

## Open questions

- R2-GENO-1's per-reference rate was measured on 3 references. A full count
  across all 160 references is a cheap rerun of `absent_check.py` per
  reference, but was not done. The rate is plausibly higher for references
  with large inverted segments. **UNVERIFIED** beyond the three.
- Whether `*`/ABSENT is ever read as a character state in the tree or
  ancestral code path decides whether R2-GENO-1 is HIGH there; that is the
  tree reviewer's area.
- The read-support classification uses exact 12-bp H37Rv flanks. "reads_none"
  therefore includes sites where the sample differs from H37Rv in the flank,
  and the percentages above are lower bounds on decidability.
