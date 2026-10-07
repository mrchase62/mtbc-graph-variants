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
a graph with none). P0 writes it into `$B/build_info.tsv` as `outgroup`, step
`panel_polarity` polarises the panel with it (and refuses another), and
`assoc/bin/cohort_assoc_tail.sh` roots each cohort tree and names the event
writer's outgroup from that record. The panel tree's own outgroup leaves for
step `ancestral` are `ANC_OUTGROUPS` (below; decision D44).

## 1. The panel tree, for step `ancestral`

Step `ancestral` writes `$B/assets/ancestral.tsv`, the ancestral allele at
every panel SNP site, which every cohort VCF's `AA` tag joins. Its inputs are
the panel's SNP alignment, the alignment's sites table and a rooted tree
built from that alignment. CX333's were `data/trees/cx333.{snps.fasta,
sites.tsv,rooted.nwk}`, made on 2026-09-24 the same way as below. Those files
are not inputs for any other graph.

**a. Alignment and sites, from the build's collapsed graph VCF** (seconds):

```bash
mkdir -p data/trees
"$MTB_PY" bin/vcf_to_alignment.py --vcf "$B/assets/graph_collapsed.vcf.gz" \
    --out data/trees/<panel>.snps.fasta --sites-out data/trees/<panel>.sites.tsv \
    --ref-sample GCF_000195955 --max-missing 0.10
```

- `--ref-sample GCF_000195955` adds H37Rv, the VCF's reference, as an all-REF
  row. Without it H37Rv is missing from its own panel tree, and step
  `ancestral`'s taxa check, which needs exactly the build's genomes, refuses
  the alignment.
- Only biallelic SNP records become columns. A `X,*` record is the SNP REF>X,
  with the `*` samples written N. Use the build's copy of the VCF so that the
  tree and every other asset come from the same file.

**b. Tree** (a cluster job, about 15 minutes on 24 cores for 333 genomes):

```bash
sbatch bin/build_snp_tree.sh data/trees/<panel>.snps.fasta <outgroup> data/trees/<panel>
```

This runs IQ-TREE with GTR+F+ASC+G4. `+ASC` is needed because the alignment
holds only variable sites. The script writes `data/trees/<panel>.rooted.nwk`,
rooted on `<outgroup>`.

**Choosing the outgroup** (HANDOFF 0j and 0k):

- **Use a genome that is in the graph,** so that it has a row in the
  alignment.
- **For a lineage 1-4 (+7) graph:** the recommended in-graph outgroup is one
  or two clean lineage 5/6 genomes. They are the sister clade and add about
  3% of SNP sites.
- **Canettii and lineage 8:** keep them out of the graph. They belong to the
  outgroup kit, for questions about the root of the whole complex.
- **Do not reuse CX333's default** `GCF_035581225` (canettii ET1291) unless
  the graph contains it.
- **Cross-check:** lineage 1 must be the first split inside the lineage 1-4,7
  clade.

**c. The step:**

```bash
ANC_TREE=data/trees/<panel>.rooted.nwk ANC_ALN=data/trees/<panel>.snps.fasta \
ANC_SITES=data/trees/<panel>.sites.tsv ANC_OUTGROUPS=<outgroup>[,<second>] \
    bash bin/p0_prepare.sh --step ancestral
```

- **`ANC_OUTGROUPS`:** the tree's outgroup leaves. The ancestral node is the
  MRCA of every other leaf (audit fault B), and the outgroups break ties
  there, nearest first.
- **Unset:** `bin/ancestral_alleles.py` falls back to the two CX333 canettii,
  and refuses to run if they are not leaves of the tree.
- **Copies:** the step copies all three inputs into
  `$B/assets/ancestral_inputs/` and records their checksums. If any input
  changes later, the step refuses to run again.

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
