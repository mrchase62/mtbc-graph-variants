# The panel tree, and the other P0 inputs made outside P0

P0 (`bin/p0_prepare.sh`) builds every per-build asset from the graph and its
own files, except for three things that need a cluster job or a choice. This file says how to make each one for a
new graph. Every command is run from the repository root after
`source config/project_env.sh`, with `B=refbias/build/<build id>` (the build
P0 made; `bash bin/p0_prepare.sh --list` prints it).

No step has a CX333 default any more. A step whose input is missing stops
with an error. It does not fall back to `data/trees/cx333.*`,
`accessory/assets/` or `refbias/build/7713a8d71d8e`.

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

**Not yet produced by P0:** `$B/assets/accessory_loci.tsv` is still copied
from `ACCESSORY_DIR/panel_manifest.tsv` (`refbias/panel/`). Its producer is
outside this repository.

## 3. IS6110 crossmaps for every panel genome: step `is6110_intervals`

The P1 tie-break needs every panel genome's IS6110 interval count.
`is6110/bin/p1i_build_matched.sh` now takes the build from `MTB_BUILD_DIR`
and, without `REFMAP`, builds every genome in `$B/assets/accessions.txt`.
Each genome needs its stage-1 GFF first.

```bash
MTB_BUILD_DIR=$B REFS=$B/refs <stage 1 over every panel genome>   # see below
MTB_BUILD_DIR=$B sbatch is6110/bin/p1i_build_matched.sh
bash bin/p0_prepare.sh --step is6110_intervals
```

**Known gap: stage 1** (`is6110/bin/p1i_discover_matched.sh`):

- It still reads the refmap's references.
- It defaults to the CX333 build's `refs/`.
- It skips any GFF that already exists, even one made from an older build's
  sequence.

Until it is changed to match `p1i_build_matched.sh`, run its loop by hand over
`$B/assets/accessions.txt` with `REFS=$B/refs`, into an empty `OUTDIR`.

`p1i_build_matched.sh` records the sha256 of the reference each clean build
was made from (`<R>.ref.sha256`). It rebuilds any reference whose bytes
differ, so a crossmap made from another build's sequence is not kept.
