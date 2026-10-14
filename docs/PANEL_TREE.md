# The panel tree, and the other P0 inputs made outside P0

P0 (`bin/p0_prepare.sh`) builds every per-build asset from the graph and its
own files, except for four things that need a cluster job or a choice: the
panel tree, the IS6110 crossmaps, the accessory panel, and the outgroup. This
file says how to make each one for a new graph. Every command is run from the
repository root after `source config/project_env.sh`, with
`B=refbias/build/<build id>` (the build P0 made; `OG=<graph.og> bash
bin/p0_prepare.sh --list` prints it).

No step has a CX333 default any more. A step whose input is missing stops
with an error. It does not fall back to `graphs/CX333.*`, `data/trees/cx333.*`,
`accessory/assets/`, `refbias/panel/`, `refbias/t11/` or
`refbias/build/7713a8d71d8e`. P0 needs `OG` (or `GRAPH_DIR`) and, for step
`assets`, `ACCESSORY_DIR` (section 4).

**The outgroup** is set in one place, `MTB_OUTGROUP` in
`config/project_env.sh` (default `GCF_035581225`, CX333's canettii; empty for
a graph with none; `GCF_022870265`, the L5.1.1 genome, for the P200 panel).
P0 writes it into `$B/build_info.tsv` as `outgroup`, and
`assoc/bin/cohort_assoc_tail.sh` roots each cohort tree and names the event
writer's outgroup from that record. It must be a genome in the graph. The
panel's polarity is no longer that one genome's read: step `ancestral`
reconstructs it by parsimony over the panel tree with all the tree's outgroup
leaves, `ANC_OUTGROUPS` (section 1; decision D44 and the user's decision of
2026-10-13), and refuses a stamped outgroup that is not one of them.

## 1. The panel tree, for step `ancestral`

Step `ancestral` runs `bin/panel_ancestral.py` and writes the build's two
polarity tables, from one Fitch reconstruction over the panel tree (user,
2026-10-13: "yes, use parsimony over all outgroups for every variant class"):

| table | keys | read by |
|---|---|---|
| `$B/assets/ancestral.tsv` | every SNP key of the graph VCF: `chrom pos ref alt AA flag polarity`, flag `TIED:<bases>` or `NODATA` | `bin/p5_finish.sh` → the merged VCF's `AA` tag |
| `$B/assets/panel_polarity.tsv` | every key, every class (SNP, MNP, indel, SV, IS6110, accessory records), with a resolved state: `chrom pos ref alt outgroup_gt ancestral panel_af n_records method` | `assoc/bin/write_event_matrix.py --panel-polarity`, where `AA` did not resolve |

It replaces `bin/ancestral_alleles.py` (SNP alignment only; a tie at the
ancestral node broken by the nearest outgroup) and `bin/panel_polarity.py`
(one outgroup genome's genotype; step `panel_polarity`, which now refuses with
a pointer here). The method, step by step:

1. **Keys.** Every biallelic record of `$B/assets/graph_collapsed.vcf.gz`,
   left-aligned against the build's H37Rv with `bin/mtb_norm.py` exactly as
   the cohort VCF's keys are; copies of one key pooled.
2. **Leaf states per key.** 1 if the genome's GT is 1 in any copy, 0 if 0 in
   some copy and 1 in none, else no call. At a SNP key, GT 0 while carrying a
   sibling allele that changes the same base is no call (TP-7). H37Rv is 0 at
   every key. **Kit genomes** (outside the graph) take SNP states from their
   direct alignment (`analysis/panel_checks/bin/genome_states.py`), and
   indel, MNP and SV states from outgroup kit route A (`bin/outgroup_kit.py`:
   the assembly aligned straight to H37Rv; user, 2026-10-14, option A, until
   OUT-1's pseudo-isolates replace it). Route A's rules, any ambiguity being
   no call:
   - **indel or MNP (< 50 bp):** 1 if the kit has the same normalised allele
     (its `.var` V lines go through `mtb_norm` too; an MNP is never 1, as
     paftools writes single bases); 0 if the key's span, anchor and one
     flanking base included, lies in ONE callable interval (an `.var` R
     line: H37Rv covered by exactly one alignment) and no kit variant
     touches it; else no call.
   - **SV (>= 50 bp):** each kit SV (`.var` indels of >= 50 bp, with their
     sequence; `.sv.tsv` DEL/INS rows between alignments, without) is
     matched to graph SV keys of its type with the start within 50 bp and
     the net length within 20% (the rule of
     `analysis/inherited_matched_check/compare_h37rv.py`). The kit's exact
     allele wins and the other candidates are no call; without one, a
     single candidate is 1 and several are all no call. An unmatched SV key
     is 0 if its span lies in one callable interval with no kit SV, and no
     kit indel of 10 bp or more, within 50 bp; else no call. REPL, INV and
     REARR rows only block a 0.
   Without route A files, a kit genome is no call at every non-SNP key.
3. **Fitch upward pass** over the whole tree.
4. **The ancestral node** is the MRCA of every leaf not in `ANC_OUTGROUPS`:
   for the P200 panel, of the lineage 1-4 and 7 genomes. If an outgroup is
   beneath it the step stops (audit fault B).
5. **Fitch downward pass** from the root to that node: all outgroups at once,
   their order never read; no nearest-first step. One state: that is the
   ancestral allele (`method` `ingroup` if the node's own upward set decided
   it, `outgroups` if only the downward pass did). Both states: `TIED`. Tied
   keys are not written to `panel_polarity.tsv`, so the event writer leaves
   them unpolarised; `ancestral.tsv` keeps them as `AA=.`, `TIED:..`.
6. Unchanged downstream: the event writer's 5% panel-AF guard on an
   ALT-ancestral row (`--min-panel-af`), "presence is derived" for accessory
   sequence absent from H37Rv, and `AA` before the table.
   `outgroup_gt` is now the ancestral node's state (0/1); `panel_af` is as
   before (carriers over genomes called, each genome once).

**a. Alignment and sites, from the build's collapsed graph VCF** (seconds):

```bash
mkdir -p data/trees
"$MTB_PY" bin/vcf_to_alignment.py --vcf "$B/assets/graph_collapsed.vcf.gz" \
    --out data/trees/<panel>.snps.fasta --sites-out data/trees/<panel>.sites.tsv \
    --ref-sample GCF_000195955 --max-missing 0.10
```

- `--ref-sample GCF_000195955` adds H37Rv, the VCF's reference, as an all-REF
  row. Without it H37Rv is missing from its own panel tree, and the step's
  taxa check, which needs exactly the build's genomes, refuses the tree.
- Only biallelic SNP records become columns. A `X,*` record is the SNP REF>X,
  with the `*` samples written N. Use the build's copy of the VCF so that the
  tree and every other asset come from the same file.
- The alignment is used for the tree only. The reconstruction reads the
  graph VCF's genotypes directly, so the `--max-missing` filter does not
  limit which keys get an ancestral allele.

**b. The kit genomes' rows** (seconds). The P200 panel keeps lineage 8
(`GCF_012923765`) and canettii ET1291 (`GCF_035581225`) out of the graph, as
the outgroup kit, but they belong on the tree: canettii roots it, and both
take part in the reconstruction at SNP keys. Their rows come from their
direct-alignment states:

```bash
S=analysis/panel_checks/states
"$MTB_PY" bin/append_kit_rows.py --alignment data/trees/<panel>.snps.fasta \
    --sites data/trees/<panel>.sites.tsv \
    --kit GCF_012923765=$S/GCF_012923765.npz --kit GCF_035581225=$S/GCF_035581225.npz \
    --out data/trees/<panel>.kit.snps.fasta
```

Each kit cell is the genome's base if it is the column's REF or ALT, else N.
The columns are the panel's variable sites, so the kit genomes' private SNPs
are missing and their terminal branches are short; the topology and the root
are what the reconstruction uses. On CX333, the kit rows made this way agree
with the same genomes' graph rows at all but 6 of 72,643 (L8) and 76 of
71,226 (ET1291) cells called by both.

**c. Tree** (a cluster job, about 15 minutes on 24 cores for 333 genomes):

```bash
sbatch bin/build_snp_tree.sh data/trees/<panel>.kit.snps.fasta GCF_035581225 data/trees/<panel>
```

This runs IQ-TREE with GTR+F+ASC+G4. `+ASC` is needed because the alignment
holds only variable sites. The script writes `data/trees/<panel>.rooted.nwk`,
rooted on canettii.

**Choosing the outgroups** (HANDOFF 0j and 0k; decided 2026-10-13):

- **In the graph:** one or two clean lineage 5/6 genomes, the sister clade of
  lineages 1-4,7. For P200: `GCF_022870265` (L5.1.1, carries the tatC-helY
  element) and `GCF_000253355` (L6.1.1). One of them is `MTB_OUTGROUP`.
- **The kit, outside the graph:** lineage 8 `GCF_012923765` and canettii
  `GCF_035581225`.
- **Cross-check:** lineage 1 must be the first split inside the lineage 1-4,7
  clade, and nothing else may lie under that clade's MRCA (the step refuses
  it).

**d. The step:**

```bash
S=analysis/panel_checks/states
MTB_OUTGROUP=GCF_022870265 \
ANC_TREE=data/trees/<panel>.rooted.nwk \
ANC_OUTGROUPS=GCF_022870265,GCF_000253355,GCF_012923765,GCF_035581225 \
ANC_KIT=GCF_012923765=$S/GCF_012923765.npz,GCF_035581225=$S/GCF_035581225.npz \
ANC_KIT_EVENTS=GCF_012923765=analysis/outgroup_kit/out/GCF_012923765,GCF_035581225=analysis/outgroup_kit/out/GCF_035581225 \
    OG=<graph.og> bash bin/p0_prepare.sh --step ancestral
```

- **`ANC_OUTGROUPS`** (required, no default): every outgroup leaf, in-graph
  and kit. Their order does not matter.
- **`ANC_KIT`** (optional): `acc=states.npz` for each kit genome. Each must
  be in `ANC_OUTGROUPS`, and the tree's leaves less these must be exactly the
  build's genomes (`bin/p0_check.py tree-taxa`).
- **`ANC_KIT_EVENTS`** (optional): `acc=prefix` for each kit genome with
  route A files `prefix.var` and `prefix.sv.tsv` (`bin/outgroup_kit.py`,
  `analysis/outgroup_kit/README.md`). Each must be in `ANC_KIT`. Both files
  are copied and checksummed like the rest.
- **The stamped outgroup** (`MTB_OUTGROUP`, in `build_info.tsv`) must be one
  of `ANC_OUTGROUPS` and not a kit genome.
- **Needs** steps `accessions`, `assets` and `refs`.
- **Copies:** the tree and kit files go into `$B/assets/ancestral_inputs/`,
  and the marker records their checksums with the VCF's, H37Rv's and the
  code's. If any changes later, the step refuses to run again.
- `ANC_ALN` and `ANC_SITES` are no longer read.

**Why the kit's event rows matter: the tatC-helY element.** L5.1
(`GCF_022870265`) carries the element and L6.1 (`GCF_000253355`) does not;
they join the lineage 1-4,7 clade at one node, and that clade ties on it.
Without the kit's SV rows the reconstruction calls it TIED, and the event
writer would leave it unpolarised, defaulting to the insertion as derived --
the inverted reading the old single-outgroup table had fixed. L8 and ET1291
both carry it by route A (the exact 68 bp allele), so with `ANC_KIT_EVENTS`
it is ancestral-present under every outgroup order.

**Validation on CX333** (`analysis/parsimony_check/compare.py`, report in
`out/report.txt`):
- At production's ancestral node (the MTBC ancestor, canettii outgroups),
  the new AA equals the build's `ancestral.tsv` at all 73,256 of its sites,
  and adds 3,763 resolved SNP keys the alignment had dropped.
- Route A's event states against the same genomes' graph genotypes (both
  are in CX333's graph), where both are called: L8 indels 99.84% (8,033),
  SVs 99.80% (2,494), MNPs 100% (2,083); ET1291 indels 99.82% (7,158), SVs
  99.70% (2,308), MNPs 100% (1,238). Route A calls about 88% / 79% of indel
  keys and 54% / 50% of SV keys.
- Emulating P200 (lineages 1-4,7 + L5.1 + L6.1 in the graph, L8 and ET1291
  as kit rows with route A events) against the full CX333 panel at the same
  node: SNPs agree at 77,042 keys and differ at 8; indels agree at 9,068 and
  differ at 2, SVs 4,592 and 1; 3 indels, 2 SVs and 1 MNP stay tied that the
  full panel resolves. The events resolve 4 indel, 4 SV and 1 MNP keys that
  were tied without them, the element among them, and change no key that
  was resolved.

## 2. The accessory catalogue: step `catalogue` (in P0, no action needed)

`$B/assets/accessory_catalogue.{tsv,fasta}`, the catalogue that level-1
presence genotyping and the merge read, is now made by P0:

```bash
"$MTB_PY" accessory/bin/merge_catalogues.py --loci "$B/assets/accessory_loci.tsv" \
    --graph-vcf "$B/assets/graph_collapsed.vcf.gz" \
    --clusters "${INSGT_CLUSTERS:-}" --routing "${INSGT_ROUTING:-}" --census "" \
    --out "$B/assets/accessory_catalogue.tsv" --out-fasta "$B/assets/accessory_catalogue.fasta"
```

P0 then indexes the FASTA with `bwa index`.

**CX333's copy.** `accessory/assets/accessory_catalogue.*` was made by this
script with these inputs:

- `--graph-vcf graphs/CX333.s10k.k23.K15/all_variants.decomposed.vcf.gz`;
- `--clusters insgt/assets/insertions.tsv`;
- `--routing insgt/assets/insertion_routing.tsv`;
- `--census accessory/gwas1000_accessory_census.tsv`.

Re-run on 2026-10-06, these inputs reproduced both files byte for byte.
P0 differs from that run in three ways:

- **The collapsed VCF** (audit PGB-9). On CX333, 157 of the 802 rows
  change: 23 of them in the representative allele's length, the rest in
  carrier sets.
- **No census:** its three columns are one cohort's call counts, and nothing
  reads them.
- **insgt's clusters and routing only when given:** they are CX333-specific.
  Without them, the route columns are empty and the merge writes no
  `ACCROUTE`.

**Who reads the build's copy:**

- `accessory/bin/locus_presence_array.sh`;
- `bin/p5_finish.sh`, which passes `--accessory-catalogue` to
  `merge_cohort_vcf.py`.

**Its input,** `$B/assets/accessory_loci.tsv`, is copied by step `assets`
from `ACCESSORY_DIR/panel_manifest.tsv`, made from this graph as in section 4.

## 3. IS6110 crossmaps for every panel genome: step `is6110_intervals`

The P1 tie-break needs every panel genome's IS6110 interval count.
`is6110/bin/p1i_build_matched.sh` now takes the build from `MTB_BUILD_DIR`
and, without `REFMAP`, builds every genome in `$B/assets/accessions.txt`.
Each genome needs its stage-1 GFF first.

```bash
MTB_BUILD_DIR=$B sbatch is6110/bin/p1i_discover_matched.sh   # stage 1, every panel genome
MTB_BUILD_DIR=$B sbatch is6110/bin/p1i_build_matched.sh      # stage 2, after stage 1
OG=<graph.og> bash bin/p0_prepare.sh --step is6110_intervals
```

Both scripts take the build from `MTB_BUILD_DIR` (or the one completed
build), read `$B/refs`, and without `REFMAP` cover every genome in
`$B/assets/accessions.txt`. Both record the sha256 of the reference each
product was made from (`<R>.ref.sha256` in `is6110/assets/matched_gff/` and
`is6110/assets/isclean_matched/`) and redo any reference whose bytes differ,
or that has no record, so a GFF or crossmap made from another build's
sequence is not kept. (The first run after this change re-discovers every
CX333 GFF once, because none has a record; minimap2, seconds each.)

## 4. The accessory panel, for step `assets`

Step `assets` copies `ACCESSORY_DIR/panel_manifest.tsv` into the build as
`assets/accessory_loci.tsv`, the accessory locus table that P3, P4, P5 and
step `catalogue` read, and `ACCESSORY_DIR/accessory_{novel,mosaic}.fasta`,
the panel P3 aligns to. CX333's were `refbias/panel/`, made on 2026-09-14 by
four scripts that lived only in the working tree; they are now in `bin/`:

| step | script | what it does |
|---|---|---|
| a | `bin/t2_extract_candidates.py` | every pure insertion of 500 bp or more in the graph VCF, with its carriers |
| b | `bin/t2_validate_accessory.sh` | blasts every candidate against each panel assembly (one task per genome) |
| c | `bin/t2_summary.py` | carrier vs non-carrier identity and full-length rate per allele |
| d | `bin/build_accessory_panel.py` | one locus per position, representative allele, novelty against H37Rv |

None of them has a default input or output any more (they were the pilot's
`refbias/t2/`, `refbias/T2.accessory_validation.tsv`, `refbias/panel/` and
`loci/hely_tatc/samples.txt`). With `G=<graph dir>` and `P=<panel dir>`, a
new folder:

```bash
mkdir -p $P/t2
"$MTB_BCFTOOLS" query -l $G/all_variants.collapsed.vcf.gz > $P/t2/genomes.txt
"$MTB_PY" bin/t2_extract_candidates.py --vcf $G/all_variants.collapsed.vcf.gz \
    --bcftools "$MTB_BCFTOOLS" --out-fasta $P/t2/candidates.fasta --out-meta $P/t2/candidates.tsv
SAMPLES=$P/t2/genomes.txt QUERY=$P/t2/candidates.fasta OUT=$P/t2/hits \
    sbatch --array=1-$(wc -l < $P/t2/genomes.txt)%40 bin/t2_validate_accessory.sh
"$MTB_PY" bin/t2_summary.py --hits $P/t2/hits --meta $P/t2/candidates.tsv \
    --out $P/t2/accessory_validation.tsv
"$MTB_PY" bin/build_accessory_panel.py --validation $P/t2/accessory_validation.tsv \
    --candidates $P/t2/candidates.tsv --fasta $P/t2/candidates.fasta \
    --genomes $P/t2/genomes.txt --h37rv "$MTB_H37RV" --blastn "$MTB_QC_BIN/blastn" \
    --outdir $P
OG=<graph.og> ACCESSORY_DIR=$P bash bin/p0_prepare.sh --step assets
```

- **The VCF** is the graph's collapsed file, the one step `assets` copies
  into the build. On CX333 the collapsed and the `nolab` file give the same
  1,643 candidates, byte for byte; the decomposed file gives 2,757 (its
  duplicate records).
- **The assemblies** step b blasts are `data/rotated/<genome>.dnaA_rotated.fasta`
  (run from the run root); one task is a few seconds.
- **`--genomes`** is the graph VCF's samples, the genomes step b blasted
  (H37Rv, the VCF's reference, is not one). `carrier_frac` is over them (it
  was over a fixed 332). `build_accessory_panel.py` refuses candidates whose
  carriers are not all in it, and writes it as `$P/panel_manifest.genomes.txt`.
- **P0 checks that record.** Step `assets` refuses an `ACCESSORY_DIR` with no
  `panel_manifest.genomes.txt`, or one whose genomes plus H37Rv are not
  exactly the build's (`bin/p0_check.py accessory-panel`). CX333's
  `refbias/panel/` has no record and is refused; remake it as above.

**CX333 reproduced.** On 2026-10-06 the four steps above, on
`graphs/CX333.s10k.k23.K15/all_variants.collapsed.vcf.gz`, reproduced
`refbias/t2/candidates.{tsv,fasta}`, every `refbias/t2/hits` table,
`refbias/T2.accessory_validation.tsv`, and `refbias/panel/panel_manifest.tsv`,
`accessory_novel.fasta` and `accessory_mosaic.fasta` byte for byte (hits
compared sorted; blastn's row order is not fixed).

**The IS6110 locus list and anchor sets** (`IS6110_LOCI`, `IS6110_ANCHORS`)
are no longer copied by default. They were `refbias/t11/loci_clean.tsv` and
`anchor_sets.tsv`, a pilot test's files; no pass of the chain reads them
(only the working tree's `p1b`/`p1c` benchmarks), and nothing in this
repository makes the locus list (`is6110/bin/is6110_anchor_sets.py` makes the
anchor sets from it). Step `assets` copies them only when both are named.
