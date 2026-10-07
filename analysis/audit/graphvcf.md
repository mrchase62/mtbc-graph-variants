# Audit: readers of graph-derived VCFs and tables (area `graphvcf`)

Auditor: Claude (subagent), 2026-10-05. Read-only; no jobs submitted. Scratch:
`/tmp/claude-12043/audit_graphvcf/` (the scripts named below are there and
re-run in a few minutes each).

## How the graph VCFs relate (established, not assumed)

All under `graphs/CX333.s10k.k23.K15/`, one contig `GCF_000195955#1#NC_000962.3`,
332 samples.

| file | records | made by | notes |
|---|---|---|---|
| `all_variants.decomposed.vcf.gz` | 164,650 (164,976 allele records; 322 multi-allelic) | `bin/vcf_decompose.sh` steps 1-4 (working tree only) | 93,214 distinct (pos,ref,alt); **6,186 keys duplicated** over 77,948 records (4,205 of them SNPs; 3,895 LV=1, 2,291 LV=0). In **6,151** of them a sample has GT 1 in one duplicate and 0/. in another (**172,004 sample cells**). Taking the last duplicate would give 154,487 wrong carrier cells; the first, 154,435. |
| `all_variants.collapsed.vcf.gz` | 93,214 | `vcf_decompose.sh:151-153` (`norm -m +any \| norm -m -any \| fill-tags`) | 0 duplicate keys. See GRAPHVCF-5 for what the round-trip does wrong. |
| `all_variants.nolab.vcf.gz` | 91,304 | `bin/exclude_samples.sh` on collapsed | header shows 0 samples were dropped (332 before and after); only `view -c 1` removed 1,910 records |
| `snps.vcf.gz` = `refbias/build/7713a8d71d8e/assets/panel_snps.vcf.gz` (symlink) | 77,808 | `bin/vcf_split_classes.sh --in all_variants.nolab.vcf.gz` (header provenance) | 75,268 plain 1-bp SNPs, 319 SNPs with padded REF, 2,221 MNPs |
| `indels.vcf.gz`, `svs.vcf.gz`, `small_variants.vcf.gz` | 8,638 / 4,850 / 86,454 | same split | no production reader found |
| `*.inv64.*` | 164,639 / 91,528 | `bin/vcfwave_rerun.sh` | no production reader found (DR_ELEMENTS analysis only) |

Build assets: `panel_snps.vcf.gz` -> `snps.vcf.gz` (above);
`accessory_loci.tsv` -> `refbias/panel/panel_manifest.tsv` (from
`bin/build_accessory_panel.py`, whose candidates come from
`bin/t2_extract_candidates.py` on the collapsed/nolab VCF: 1,643 alleles, no
identical (pos, sequence) duplicates); `ancestral.tsv` from
`bin/ancestral_alleles.py` over `data/trees/cx333.snps.fasta` + `cx333.sites.tsv`
(72,986 sites, 0 duplicate keys, i.e. built from the collapsed lineage, not from
decomposed); `graph_frame_offsets.tsv` from `graphframe/bin/graph_frame_offsets.py`
(reads FASTAs, not the VCF).

## Production readers of graph-derived VCFs/tables (which copy runs)

| reader | input | copy that runs | verdict |
|---|---|---|---|
| `assoc/bin/add_outgroup.py` | decomposed | working tree (no repo copy), chain step 2 | GRAPHVCF-1, -2 |
| `bin/panel_polarity.py` -> `refbias/assets/panel_polarity.tsv` -> `assoc/bin/write_event_matrix.py:467-471` | decomposed | working tree; table hand-built 2026-09-25, read by chain step 4 by default | GRAPHVCF-3 |
| `bin/t8_select_reference.py` (P1) | `assets/panel_snps.vcf.gz` | repository (`runroot/bin` -> repo `bin`); identical md5 in working tree | GRAPHVCF-4 |
| `bin/vcf_decompose.sh` / `bin/vcf_collapse.sh` (producer of every non-decomposed file) | decomposed | working tree only; will run again on the new graph | GRAPHVCF-5 |
| `bin/p5_matrix.py:84-104,139-143` (P5) | nolab | repository | GRAPHVCF-6 |
| `bin/p5_sanity.py:160-172` | nolab | repository | GRAPHVCF-6 (report only) |
| `bin/p4_place.py:366-404` (P4 inherited half) | nolab | repository | sound (below) |
| `bin/p4b_place_sv.py:256-330` (P4b inherited half, LV logic) | nolab | repository | sound (below) |
| `bin/sv_intervals.py` via `p5_svgt.sh --catalogue` | decomposed | repository | sound (below) |
| `bin/merge_cohort_vcf.py:456-458` (AA join) | `assets/ancestral.tsv` | repository | sound, but inherits GRAPHVCF-5's missing sites |
| `bin/vcf_to_alignment.py` (chain step 1; also built `cx333.snps.fasta`) | cohort merged VCF / collapsed panel VCF | working tree | sound on its current inputs (below) |

---

## GRAPHVCF-1 (HIGH, already known as Fault A; confirmed): `add_outgroup.py` takes the last duplicate record

- **File:** `/n/netscratch/.../MtbPangenome/assoc/bin/add_outgroup.py:116-121`
  (working tree; there is no repository copy; run by `assoc/bin/cohort_assoc_tail.sh` step 2).
- **What is wrong:**
  ```python
  for rref, ralt, gt in want.get(pos, ()):
      ...
      if rref == ref and "," not in ralt and ralt == alt:
          exact = g          # last duplicate wins, including a '.' (None)
  ```
  The decomposed VCF holds the outgroup's ALT in one of several identical
  (pos,ref,alt) records.
- **Evidence:** I re-implemented the script and reproduced
  `data/trees/<cohort>.og.fasta` exactly at every column (an assert over all
  sites). I then replaced last-wins with "ALT if any exact record is 1". Results:
  - **scale200:** 386 of 46,379 columns go REF->ALT. 1,596 sites have several
    exact records.
  - **gwas1000:** 323 of 102,199 columns go REF->ALT. 1,277 sites have several
    exact records.

  HANDOFF 0j measured 527 and 461 against direct alignment, which also counts
  GRAPHVCF-2's cells and alignment differences.
- **Effect:** `data/trees/{scale200,gwas1000}.og.fasta` -> the rooted cohort trees
  -> `write_event_matrix.py`'s outgroup fallback.
- **Fix:** collect every exact record. Write ALT if any has GT 1, REF if any has
  GT 0 and none has GT 1, and N only if all are missing. Better still, read
  `all_variants.collapsed.vcf.gz` (the nolab file), whose keys are unique. If you
  do, apply the GRAPHVCF-5 trim first, or the 376 padded SNP keys will not match.

## GRAPHVCF-2 (MEDIUM, new): `add_outgroup.py` writes REF where the outgroup has no call, and inside its own deletions

- **File:** `assoc/bin/add_outgroup.py:118, 127-133` (working tree, production).
- **What is wrong, part (a):** missing becomes REF. When every exact record is
  `.`, `exact` stays `None`, so control falls to the `any_rec` branch. That
  branch writes REF (`ref_other_record`). The docstring says "missing -> N".
- **What is wrong, part (b):** "no record means REF". The script also writes REF
  at a site where the outgroup has no record at that position but carries a
  deletion, starting upstream, that spans the site. The outgroup has no base
  there.
- **Evidence:** same reproduction as GRAPHVCF-1.

  | cohort | (a) all exact records missing, written REF | (b) site inside an outgroup deletion, written REF |
  |---|---|---|
  | scale200 | 333 | 108 |
  | gwas1000 | 304 | 233 |

  For (b) I took deletions where `GCF_035581225` has GT >= 1 in the decomposed
  VCF and a span of `pos+len(alt) .. pos+len(ref)-1`.
- **Effect:** same products as GRAPHVCF-1, so about 441 (scale200) and 537
  (gwas1000) further outgroup cells should be N rather than REF.
- **Fix:**
  - write N when all exact records are missing;
  - before "no record -> REF", test whether the position lies in a deletion
    allele the outgroup carries, and write N if so.

## GRAPHVCF-3 (MEDIUM, new): panel polarity table built from the decomposed VCF, loaded last-wins, with per-record allele frequencies

- **Files:**
  - `/n/netscratch/.../MtbPangenome/bin/panel_polarity.py:34-70` (working tree,
    run by hand on 2026-09-25 into `refbias/assets/panel_polarity.tsv`);
  - `/n/netscratch/.../MtbPangenome/assoc/bin/write_event_matrix.py:311-312,
    467-471`, which loads it by default from the chain (step 4 passes no
    `--panel-polarity`, and the default path exists).
- **What is wrong:**
  - The table has one row per **decomposed record**, keyed (pos,ref,alt):
    - `outgroup_gt` is the outgroup's GT in that record only;
    - `panel_af = INFO/AC / INFO/AN` is that record's partial carrier count.
  - The reader does `pol_tab[(pos, ref, alt)] = (r["ancestral"], r["panel_af"])`,
    so the last row wins. This is Fault A again, on every variant AA does not
    resolve (all indels, plus SNPs missing from ancestral.tsv).
  - The `--min-panel-af 0.05` gate then tests a per-record fraction that
    understates the real panel frequency.
- **Evidence:**
  - The table has 156,087 rows over 89,736 keys. 5,807 keys are duplicated, and
    1,423 of them have rows that disagree on `ancestral`. At 903 keys the last
    row's `ancestral` differs from the union over duplicates.
  - I reproduced the production `polarity` column of `assoc/<C>/events/variants.tsv`
    exactly for every `*_outgroup` and `unpolarised` H37Rv row. I then recomputed
    it with the union outgroup GT and the collapsed-key allele frequency
    (script `pol.py`):

    | cohort | REF-ancestral -> ALT-ancestral (polarity inverted) | REF-ancestral -> unpolarised | unpolarised -> ALT-ancestral | rows that use a duplicated key |
    |---|---|---|---|---|
    | gwas1000 | 31 | 8 | 0 | 222 |
    | scale200 | 32 | 9 | 1 | 236 |

    The changed rows are listed in `polarity_changes.<cohort>.tsv`.
- **Effect on current outputs:** 39 (gwas1000) and 42 (scale200) variants have
  inverted or over-confident polarity in `events/variants.tsv`, which affects
  their gains and losses. Of these, 6 and 8 reach `scan.tsv`, and none has
  `q_branch` < 0.05. The number is small now only because most variants are SNPs
  polarised by AA. It applies to every indel.
- **Related, pattern 3/4:**
  - `panel_polarity.tsv` is not a build asset. It is not regenerated by P0 or by
    the chain, and it is not tied to the graph.
  - `write_event_matrix.py` skips it silently when the file is absent
    (`if a.panel_polarity and os.path.exists(...)`). A new graph would then
    either silently use the old graph's table or silently lose all indel
    polarity.
- **Fix:**
  - Build the table from `all_variants.collapsed.vcf.gz` (unique keys, union
    genotypes, correct AC/AN), trimmed as in GRAPHVCF-5. If it stays on
    decomposed, aggregate by key: ancestral = ALT if any outgroup GT is 1, and
    AF = sum(AC)/max(AN).
  - Make it a P0 build asset beside `ancestral.tsv`, and pass it explicitly from
    `cohort_assoc_tail.sh`.
  - Make a missing table FATAL, or at least print a loud warning in the step
    output.

## GRAPHVCF-4 (MEDIUM, new): P1 reference selection keys the panel matrix by position only

- **File:** `bin/t8_select_reference.py:42-67`. The repository copy runs
  (`runroot/bin` -> repository `bin`); the working-tree copy has the same md5.
- **What is wrong:**
  ```python
  if len(alts) != 1 or len(f[3]) != 1 or len(alts[0]) != 1: continue
  pos_idx[int(f[1])] = len(rows)          # a second record at the same POS overwrites
  ...
  called.add(int(f[1]))                    # isolate profile: position only, ALT ignored
  ```
  - `panel_snps.vcf.gz` is split biallelic, so a position with C>A and C>G has
    two rows. `pos_idx` keeps only the last. The earlier row stays in `G` with
    `q = 0` forever, so every panel genome carrying that allele is charged one
    SNP whatever the isolate carries. **473 positions; 481 shadowed rows.**
  - The isolate side is allele-blind: an isolate with C>G at a site where the
    surviving panel row is C>A is scored as matching C>A carriers.
  - The `len(REF)==1` filter drops the **319 SNPs** whose REF was padded by
    GRAPHVCF-5 (for example `5075 CG>TG`), carrying 10,306 carrier cells.
- **Evidence:** I re-implemented the original exactly. Its rank-1 reference and
  distance match every `refbias/<C>/p1/<sample>.candidates.tsv`. I also ran a
  fixed version, keyed on trimmed (pos,ref,alt), with an allele-aware isolate
  profile and 75,587 sites. Both were run over all P1 isolate VCFs
  (`refbias/work/{scale200,gwas1000}_p1/*.h37rv.vcf.gz`; script `t8cmp.py`,
  output `t8cmp.tsv`):

  | cohort | isolates | rank-1 changes | strictly better under the fixed metric (the rest are ties) | production `refmap.tsv` still on the worse reference after the tie-break |
  |---|---|---|---|---|
  | scale200 | 200 | 9 | 3 | 2 (SAMN08795146, SAMN13169763; 1 SNP each) |
  | gwas1000 | 997 | 27 | 16 | 8 (SAMEA112806457, SAMEA1403576, SAMEA1403630, SAMN03647468, SAMN03647614, SAMN06091702, SAMN07660156, SAMN13208027) |

  - The 8 gwas1000 isolates are worse by 1-2 SNPs, except **SAMN07660156**,
    whose fixed choice `GCF_045348265` is better than production's
    `GCF_014899625` by **12 SNPs**.
  - For the other 9 strictly-better changes, the interval tie-break had already
    landed on the fixed choice.
  - The old metric overstates every distance by 21 SNPs on average (scale200)
    and 24 (gwas1000).
- **Effect:**
  - P1's choice for the 10 isolates listed, and everything downstream of P2-P5
    for them: one reference per isolate.
  - The extra distance is a near-constant offset for most isolates. It changes
    the winner only where the top candidates are within a few SNPs, which is
    also where `p1_summary.py`'s tie-break (`--tie-margin 5`) acts.
  - `margin` and `snp_distance` in `refmap.tsv` are inflated for every isolate.
- **Fix:**
  - Key panel rows on the trimmed (pos, ref, alt).
  - Build the isolate profile from (pos, ref, alt) with GT > 0.
  - Accept padded SNPs by trimming, or fix them upstream (GRAPHVCF-5).
  - Also take up CODE_REVIEW 6.11's open point: an uncovered isolate site
    currently counts as REF.

## GRAPHVCF-5 (MEDIUM, new): the collapse round-trip pads REF and loses carriers that share a position

- **Files:**
  - `/n/netscratch/.../MtbPangenome/bin/vcf_decompose.sh:151-153`
  - `bin/vcf_collapse.sh:66-68` (working tree only; this is what the new graph
    build will run)

  ```
  bcftools norm -m +any IN | bcftools norm -m -any | bcftools +fill-tags
  ```
- **What is wrong, part (a): REF padding.**
  - `-m +any` merges every record at a POS into one record whose REF is the
    longest REF there. `-m -any` splits it again but does **not** re-trim. So a
    SNP that shared its POS with a longer REF comes out as, for example,
    `CG>TG` or `GG>CG`.
  - Count: **1,399** records in nolab are padded, **376** of them SNPs (319 in
    `snps.vcf.gz`; the rest went to other classes by size).
  - Every reader that requires `len(REF)==1`, or joins on exact (pos,ref,alt),
    loses them:
    - `vcf_to_alignment.py` dropped them from `cx333.snps.fasta`, so
      **0 of the 376 are in `assets/ancestral.tsv`**;
    - `t8_select_reference.py` drops them (GRAPHVCF-4).
  - In the cohort event matrices, H37Rv-frame SNP variants at these keys have
    no AA:

    | cohort | SNPs without AA | polarised by the outgroup instead | unpolarised |
    |---|---|---|---|
    | gwas1000 | 112 | 92 | 20 |
    | scale200 | 112 | 92 | 20 |

    The outgroup fallback is subject to GRAPHVCF-3.
- **What is wrong, part (b): a haploid merge can hold one allele per sample.**
  - At **64** (pos, sample) cells, over 11 positions, a panel genome carries two
    different ALT alleles at one POS in the decomposed VCF (62 involve a SNP plus
    an indel or insertion anchored on the same base).
  - After the collapse, **61** of those carrier cells read 0. Comparing the union
    over decomposed duplicates with collapsed, per key and sample:
    - 61 cells go 1->0;
    - 7,060 cells go 0->`.`;
    - 1,357 cells go `.`->0.
  - The 0->`.` cells raise missingness in `vcf_to_alignment` (`--max-missing`)
    and in t8 (`G=-1` excluded).
- **Effect on current outputs:**
  - (a) 376 panel SNP sites are missing from the panel alignment, the panel
    tree's input, `ancestral.tsv` and P1.
  - (b) 61 panel genotype cells are wrong in collapsed, nolab and snps; this is
    small.
- **Fix:**
  - Add a normalisation after the split: `bcftools norm -f <H37Rv PanSN fasta>`,
    or a trim pass with `bin/mtb_norm.py`'s trim. Then assert no duplicate keys
    and no padded SNPs.
  - For (b), collapse by grouping identical (pos,ref,alt) records and OR-ing the
    GTs (a 30-line pysam script) instead of `-m +any`, which cannot represent two
    alleles in one haploid sample.
  - Re-split the classes and rebuild `cx333.snps.fasta`/`ancestral.tsv` after the
    new graph's collapse.

## GRAPHVCF-6 (LOW, new): `panel_af` and `h37rv_minor` are looked up by position, not allele

- **Files:**
  - `bin/p5_matrix.py:84-104` (load: `panel_af[pos] = max AF over every record
    at pos`) and `:139-143` (lookup by `h37rv_pos` only). The repository copy
    runs.
  - `bin/p5_sanity.py:160-172`: same pattern, report only.
- **What is wrong:** a key gets the allele frequency of whatever record at that
  POS is most frequent: another SNP allele, or an indel anchored on the same
  base. A key whose allele is not in the panel at all still gets a value.
- **Evidence:** I recomputed with allele-specific AF over trimmed nolab keys
  (`af.py`) against `refbias/<C>/p5/sites.tsv`:

  | cohort | `h37rv_minor` 1 -> 0 | flagged 1, but the allele is not in the panel | `panel_af` filled from a different allele |
  |---|---|---|---|
  | gwas1000 | 64 | 118 | 4,681 |
  | scale200 | 55 | 73 | 2,439 |

- **Effect:** `sites.tsv` columns and the P5 summary count only. `merge_cohort_vcf.py`
  treats both as META and does not carry them into the VCF, and the assoc chain
  does not read them.
- **Fix:** key `panel_af` on trimmed (pos,ref,alt) from collapsed/nolab, and
  compute AF from AC/AN of that record.

## GRAPHVCF-7 (LOW, analysis only): other panel-matrix builders with the same pattern

- `analysis/combined_tree/build_alignment.py:113`: `pan[(rec.pos, ref, alt)] = calls`,
  last wins over the whole 332-genome column. This is already noted in
  HANDOFF 0j (about 200 wrong cells per panel genome in the `<cohort>_cx333`
  trees). Not in the chain.
- `analysis/combined_tree/tree_checks.py:44-48` (`panel_lineages`, also used by
  `itol_annotate.py` and `analysis/pangenome_compare/genome_overlap.py`) reads
  decomposed. It handles duplicates correctly when GT 0 follows a 1, because 0
  keeps the previous value. But a `.` in a later duplicate overwrites an earlier
  ALT with N, and a position with no record is H37Rv's base ("no record means
  REF"). It only affects lineage labels in figures. UNVERIFIED count.
- `analysis/pangenome_compare/compare.py` and `analysis/refeval_cx333/build_backbone.py`
  read decomposed and say they pool records into events. These are
  analysis-only and I did not audit them in depth.

---

## Checked and found sound

- **`bin/sv_intervals.py` (P5svgt catalogue, decomposed):** carriers are unioned
  across duplicate and nested records by clustering (lines 120-141), so
  duplicates are harmless. Its skip of multi-allelic records loses nothing: none
  of the 322 multi-allelic records has a deletion allele of 50 bp or more. A
  reference with `.` is a non-carrier, but `p5_sv_genotype.py:284` then falls
  through to measurement rather than calling REF.
- **`bin/p4_place.py:366-404` (inherited half, nolab):** one row per record where
  the reference's GT > 0. nolab has no duplicate keys. `.` and 0 are skipped,
  which is correct for "R's own differences".
- **`bin/p4b_place_sv.py:256-330` (LV filter):** it skips an LV=0 SV only when an
  LV>0 record overlaps it. In nolab there are 1,595 LV=0 SV records and none has
  an overlapping child, so nothing is dropped (`lv.py`).
- **`bin/vcf_to_alignment.py`:** one column per record, so it is safe only on
  inputs with unique keys. Both its inputs qualify: the cohort `merged.vcf.gz`
  (sites tables have 0 duplicate keys) and the collapsed-lineage panel VCF
  (`cx333.sites.tsv` has 0 duplicate keys). It would be wrong on decomposed;
  consider asserting key uniqueness.
- **`bin/merge_cohort_vcf.py` AA join:** keyed on exact (pos,ref,alt) against
  `ancestral.tsv`, which has unique keys. The PanSN `chrom` in `ancestral.tsv`
  is ignored, which is correct for a single-contig frame. The only loss is
  GRAPHVCF-5(a).
- **Accessory loci (`accessory_loci.tsv` <- `build_accessory_panel.py` <-
  `t2_extract_candidates.py`):** built from a deduplicated VCF (1,643 candidate
  alleles, no duplicate (pos, sequence)), and carriers are taken per allele.
- **PanSN against bare names:**
  - `add_outgroup.py` maps `NC_000962.3` to `GCF_000195955#1#NC_000962.3`
    explicitly, and sends `node_*` sites to N;
  - t8, p4, p4b and p5_matrix ignore CHROM, which is safe because the graph VCFs
    have one contig;
  - `exclude_samples.sh` dropped 0 samples (no lab strains are in CX333), so
    nolab and collapsed share the same 332 columns.
- **`*` alleles:** none occur in any graph VCF. They occur only in cohort merged
  VCFs, which these readers do not treat as panel data.
- **`graph_frame_offsets.tsv`:** derived from FASTAs, not the VCF. Outside this
  area's VCF pattern.

## Open questions

1. GRAPHVCF-1/2: the HANDOFF figures from direct alignment (527 and 461) and my
   code-level figures (386+333+108 and 323+304+233) count different things. I
   did not reconcile them site by site.
2. GRAPHVCF-4: I did not re-run `p1_summary.py`'s interval tie-break on the
   fixed distances, so the final reference for the 10 isolates (and for any
   isolate whose top two move within the 5-SNP margin) may differ from the fixed
   rank-1. UNVERIFIED.
3. Whether the 7,060 0->`.` cells from the collapse (GRAPHVCF-5b) are right or
   wrong depends on what vg's `.` means at each snarl (no traversal against
   conflict). I did not settle this.
4. In the working tree, `bin/annotate_svs.py`, `crosscheck_sv_methods.py`,
   `locus_sv_density.py`, `phage_patterning.py`, `phirv1_figure_data.py`,
   `t16_graph_frame.py`, `t17_clade_specific_sites.py`, `t6_encoding_hazard.py`,
   `tree_coevolution.py`, `stage0_sv_partition.py`, `stage1_selector.py`,
   `is6110/bin/is6110_*`, `is6110/bin/is_junction_*`, `sv2frame/bin/clip_discovery.py`,
   `assoc/bin/caller_graph_overlap.py` and `accessory/bin/{accessory_call_census,merge_catalogues}.py`
   all read graph VCFs. None is called by `refbias_run.sh` (which runs the
   repository's `is6110/bin` and `sv2frame/bin`, which contain no graph-VCF
   readers) or by `cohort_assoc_tail.sh`, so I did not audit them.
