# Input specification

What you must supply to run this pipeline from genome assemblies to a merged,
genotyped VCF. Derived from the scripts, not from memory: every requirement
below traces to a `mtb_require_file`, an argument check, or a default in
`config/project_env.sh`.

There are **three stages** and they have different inputs. You do not need the
stage-1 inputs to run stage 3 against an existing build.

---

## Stage 1 — build the pangenome graph

**The stage-1 scripts are not in this repository.** `pggb_build.sh`, the QC
scripts in the table below and `QC_PIPELINE.md` belong to the companion
repository **mtbc-pangenome-graph**, which builds the graph this one consumes.
They are listed here so the whole input contract is in one place.

`bin/pggb_build.sh <input.fasta.gz> <output-name> [pggb args...]`

| input | form | notes |
|---|---|---|
| **complete genome assemblies** | one bgzipped multi-FASTA | **PanSN naming is mandatory**: `SAMPLE#HAP#CONTIG`, e.g. `GCF_000195955#1#NC_000962.3`. Everything downstream parses sample identity out of these names. |
| pggb | container or conda env | `MTB_PGGB_SIF`, or the components `wfmash`, `seqwish`, `smoothxg`, `gfaffix`, `odgi` |

Output: `$MTB_GRAPHS/<output-name>/*.smooth.final.{gfa,og}`.

The 2026 graph used 333 genomes at `-s 10000 -p 95 -k 51 -K 21`. **Those
parameters and the component versions are part of the result** — the graph
filename carries their fingerprint (`f4f5ee2.11fba48.36e68b6`), and rebuilding
with different versions produces a different graph that invalidates every
coordinate projection and node-frame key produced from the old one.

### What the assemblies must be, and the QC that decides it

Complete or near-complete, one contig where possible: fragmented input produces
spurious graph structure that the accessory and structural arms then genotype as
real.

**Correction to an earlier draft of this document, which said nothing checks
this for you. That was wrong.** There is a full per-assembly QC stage, and it is
upstream of `pggb_build.sh` rather than inside it, which is why it did not
appear in the dependency closure computed from the graph builder. It is
documented in `QC_PIPELINE.md` stages 0-2 and it is 17 scripts:

| step | script | what it decides |
|---|---|---|
| 1.1 | `assembly_qc_stats.py` | contig count, length, N content |
| 1.2 | `collinearity_qc.py` | mis-assembly and rearrangement against the reference |
| 1.3 | `assembly_provenance_screen.py` | provenance and insertion content |
| 1.4 | `foreign_insertion_screen.py` | foreign or engineered sequence |
| 1.5 | `tbprofiler_assemblies.sh` -> `tbprofiler_collect.py` | lineage assignment |
| 1.6 | `data/annotation/lab_strains.tsv` | known laboratory and vaccine strains |
| 2.1 | `rotate_to_dnaa.sh` | rotate to a common origin, **and verify** the achieved offset |
| 2.2 | `build_panel.py` | apply exclusions |
| 2.3 | `qc_master_table.py` | join every measurement into one table |
| 2.4 | `make_pansn_fasta.sh` | PanSN naming, without which graph coordinates do not match |
| 2.5 | `snp_nonredundant.py` | optionally collapse clonal redundancy |

`QC_PIPELINE.md` records that steps 1.3 and 1.4 matter most and were the last to
be found, which is worth heeding: a provenance or engineered-sequence problem is
invisible in contig statistics.

**Clonal redundancy is the one to take seriously.** Public complete genomes are
redundant and any analysis treating them as independent is inflated. Lineage 4.9
looked highly conserved for the repeat element until 34 of its 38 members proved
to sit in one cluster within 50 SNPs -- 7.6x redundant. Panel-wide it is a milder
1.6x. Use single-linkage, which never overstates independence, and check
threshold sensitivity: at 200 SNPs the clusters chain and collapse.

---

## Stage 2 — prepare the build assets

`bin/p0_prepare.sh`

| input | default | required? | what it is |
|---|---|---|---|
| graph | none: `OG` (or `GRAPH_DIR`) | **yes** | the `.og` from stage 1 |
| accessory panel | none: `ACCESSORY_DIR` | **yes** | `panel_manifest.tsv`, its genome record and the panel FASTAs, made from this graph (`docs/PANEL_TREE.md` section 4) |
| outgroup | `MTB_OUTGROUP` in `config/project_env.sh` | | stamped into `build_info.tsv`; empty for a graph with none |
| panel SNP VCF | `<graph>/snps.vcf.gz` | **yes** | `vg deconstruct` output over the graph |
| H37Rv reference | `MTB_H37RV`, `MTB_REF_FASTA` | **yes** | plain and PanSN-named copies of the coordinate reference |
| assemblies | `data/assemblies` | yes | the same genomes, unpacked, for panel extraction |
| repeat mask | `data/annotation/H37Rv_repeat_mask.bed` | **yes** | see below |
| element sequence + sites | `data/annotation/IS6110.fasta`, `H37Rv_IS6110.pansn.gff` | for the element arm | |
| gene annotation | `data/annotation/H37Rv_snpeff_dump.txt` | for annotation | |
| known RDs | `data/annotation/known_RDs.bed` | for validation only | `chrom start end NAME|LINEAGE` |
| NCBI summary | `data/ncbi/assembly_summary_refseq.txt` | for naming | |
| bwa, GATK | `MTB_BWA`, `MTB_GATK_SIF` | **yes** | |

Output: `refbias/build/<build-id>/assets/` — the panel genomes with bwa indexes,
the accessory locus catalogue, the repeat mask, the element GFF, and
`graph_frame_offsets.tsv`.

That last one comes from `--step frames`, which runs after `--step refs`. It
measures, for every accession, how the graph's coordinate frame relates to the
deposited sequence in `refs/`. For 110 of 333 accessions they differ by a
rotation, and for 22 also by strand. Every projection reads it.
`bin/refbias_run.sh` passes it to each job as `MTB_GRAPH_FRAMES`.

    export OG=<graph.og> ACCESSORY_DIR=<panel dir>
    bash bin/p0_prepare.sh                      # cheap steps
    bash bin/p0_prepare.sh --step gff
    sbatch --array=1-333 bin/p0_prepare.sh --step refs
    bash bin/p0_prepare.sh --step frames
    bash bin/p0_prepare.sh --step manifest

**The repeat mask is generated, not supplied.** `bin/build_repeat_mask.py`
measures it: 50-mer uniqueness for paralogy and 9-mer recurrence for tandem
repetition, unioned with a curated name-based list. Do not hand-write one — a
name-based mask missed a PE_PGRS-family gene whose repeat tract then produced
the four most widely shared false "novel deletions" in a 200-isolate cohort.

---

## Stage 3 — run a cohort to a VCF

`bin/refbias_run.sh --cohort <name>`

This is the stage you will run most, and it needs **four things**.

### 3.1 A read collection

| variable | required | what |
|---|---|---|
| `MTB_CRAM_ROOT` | **yes** | root of a CRAM or BAM collection of reads aligned to one reference |
| `MTB_CRAM_REF` | **yes** | the reference those files were encoded against, needed to decode CRAM |

**Declared in `config/project_env.sh` with no default, deliberately.** Set them
in `config/site.local.sh`, which is gitignored:

    export MTB_CRAM_ROOT=/path/to/collection
    export MTB_CRAM_REF=/path/to/collection/reference.fasta

The 2026 runs used a collection staged on scratch by a colleague. That was an
opportunistic reuse, it is transient, and it is **one read source among
several** — FASTQ from a new sequencing run will not look like it. Nothing in
the pipeline may assume that particular collection's shape, and any script that
needs the root calls `mtb_require_cram_root` and fails with a clear message
rather than falling back to a path that happened to exist on one cluster.

### 3.2 A cohort table — `sample` and `lineage` are the only columns the pipeline reads

    sample  lineage  deepest  meandepth  coverage  mapping_rate  error_rate  duplicate_rate

| column | used for | required |
|---|---|---|
| `sample` | the key joining every table and output | **yes** |
| `lineage` | the lineage null in the association scan | for the scan |
| the rest | provenance and isolate screening | no |

**Isolate selection is yours, and it matters.** The standing criteria for this
project: paired-end only, at least 60x mean depth, and no mixed samples. Public
sequence archives vary enormously in quality, and a marginal isolate does not
give a marginally worse answer — a low-depth sample produces missing calls that
score as method failures, and a mixed sample carries real minority alleles at
positions a single-strain analysis scores as errors. Both corrupt precision and
recall invisibly. The pipeline does not screen for you.

### 3.3 A read-location table

    sample  cram_relpath

One row per isolate; `cram_relpath` is relative to `MTB_CRAM_ROOT`. Every
sample in the cohort table must appear here or its task fails with
`FATAL: <sample> not in <table>`.

### 3.4 A registry entry

`refbias/cohorts.tsv`, one row:

    cohort  cohort_table  crams_table  outroot  passes  note  workprefix

`passes` is the chain to run; the full graph-to-VCF chain is:

    p1,p2,p3,p4,p4b,p5,p1g,p1i,p1iv,p1is,p5svgt,p5vcf

**`outroot` and `workprefix` are declared rather than derived on purpose.** A
guessed path is how one cohort's calls get scored against another cohort's
reference.

---

## Optional inputs, and what you lose without them

| input | enables | without it |
|---|---|---|
| level-1 presence tables (`accessory/<cohort>/`) | the 802 accessory presence characters, and level-2 conditioning | the accessory arm carries within-insert variation only, with no presence character |
| an outgroup FASTA + sites | resolved ancestral alleles for some variants | more variants fall back to the assumption that ALT is derived |
| a phenotype file (one sample ID per line) | the convergence scan | the chain still produces the VCF and the event matrix |
| `data/annotation/known_RDs.bed` | validation against known regions of difference | no external check on the deletion arm |

---

## Tools

Everything is reached through `config/project_env.sh`; nothing is assumed on
`PATH`. Run `bash bin/show_config.sh` — it prints every resolved path and a
presence check, and tells you what is missing.

| tool | used by |
|---|---|
| bwa, samtools, bcftools, bgzip, tabix, bedtools | throughout |
| GATK | variant calling |
| delly, dysgu | structural variant calling |
| odgi | coordinate projection |
| vg | graph deconstruction to the panel VCF |
| pggb (or wfmash/seqwish/smoothxg/gfaffix) | stage 1 only |
| snpEff | annotation |
| python 3 with numpy, scipy, pysam | throughout |
| IQ-TREE | the tree, for the association stage |

**Containers belong on durable storage, not scratch.** Two were lost to a
scratch purge during development. `bin/stage_in.sh` restores them from an
archive; `bin/show_config.sh` tells you when they are missing.

---

## The shortest real path to a VCF

Assuming a graph and build already exist:

    # 1. site paths, once
    cp config/site.local.sh.example config/site.local.sh   # then edit it

    # 2. check the environment resolves
    bash bin/show_config.sh

    # 3. declare the cohort: two tables and a registry row
    #    refbias/cohort.<name>.tsv   sample + lineage
    #    refbias/<name>.crams.tsv    sample + cram_relpath
    #    refbias/cohorts.tsv         one row

    # 4. run the chain
    bash bin/refbias_run.sh --cohort <name>

    # 5. the deliverable
    #    refbias/<name>/p5/merged.vcf.gz
