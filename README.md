# mtbc-graph-variants

Genotype SNPs, indels, structural variants, repeat-element insertion sites and
accessory-locus presence into **one merged VCF**, using a pangenome graph to
choose a per-isolate reference and to project every call into a common frame.

**Input:** a pangenome graph, plus sequencing reads aligned to a reference.
**Output:** `refbias/<cohort>/p5/merged.vcf.gz`.

See `INPUTS.md` for the full input specification. The graph is built by the
companion repository **mtbc-pangenome-panel**.

## What makes this different from calling against one reference

Every isolate is aligned to the panel genome **nearest to it by SNP distance**,
not to H37Rv. Calls are then projected back into H37Rv coordinates through the
graph. The point is the accessory genome: sequence absent from H37Rv cannot be
genotyped against H37Rv at all, and a matched reference carries it.

Five record classes reach the merged VCF, and their completeness differs a great
deal. Measured on 997 isolates:

| class | records | cells stated |
|---|---:|---:|
| small variants | 173,660 | 98.4% |
| repeat-element insertion sites | 2,250 | 98.4% |
| deletions (catalogued intervals) | 3,658 | 95.1% |
| accessory locus presence | 802 | 92.6% |
| **insertions (caller-derived)** | **7,313** | **0.9%** |

**Read that last row before using the insertion records.** They carry zero
reference calls, because depth in the H37Rv frame cannot genotype an insertion —
there is no interval there to measure coverage over. A non-carrier and an
untested isolate are the same cell. The accessory presence class exists to answer
that question properly for the 802 catalogued loci; the remaining caller
insertion records are presence-only and should be treated as such.

## The chain

    pggb graph  ->  p0_prepare  ->  p1 p2 p3 p4 p4b p5  ->  p5vcf  ->  merged.vcf.gz
                                     p1g p1i p1iv p1is
                                     p5svgt

| pass | script | what it does |
|---|---|---|
| p0 | `bin/p0_prepare.sh` | build assets from the graph: panel genomes, accessory catalogue, repeat mask |
| p1 | `bin/p1_select_reference.sh` | pick each isolate's nearest panel genome by SNP distance |
| p2 | `bin/p2_call.sh` | align to that reference and call, in its coordinates |
| p3 | `bin/p3_accessory.sh` | accessory content from reads the reference could not place |
| p4 | `bin/p4_place.sh` | project small calls into the common frame, routed by region |
| p4b | `bin/p4b_place_sv.sh` | project the structural calls |
| p5 | `bin/p5_merge.sh` | union the key space, then state every isolate at every key |
| p1g/p1i | `is6110/bin/` | element-arm references and junction calls |
| p5svgt | `bin/p5_svgt.sh` | genotype deletion **absence** by depth, so the class is not presence-only |
| acc | `accessory/bin/locus_presence_one.sh` | accessory locus presence, both instruments |
| p5vcf | `bin/p5_finish.sh --merge` | one VCF, last, after every arm that contributes |

Run it with `bash bin/refbias_run.sh --cohort <name>`.

## Setup

    # 1. site paths -- gitignored, required, no defaults
    cat > config/site.local.sh <<'EOF'
    export MTB_CRAM_ROOT=/path/to/read/collection
    export MTB_CRAM_REF=/path/to/reference/the/reads/were/encoded/against.fasta
    EOF

    # 2. check everything resolves
    bash bin/show_config.sh

    # 3. declare the cohort -- see INPUTS.md
    # 4. run
    bash bin/refbias_run.sh --cohort <name>

**The read collection is declared with no default, deliberately.** It is site
configuration, it may be transient, and it is one source among several. A script
that needs it fails with a clear message rather than falling back to a path that
happened to exist on one cluster.

## Record which graph you used

The graph's filename carries a fingerprint of the component versions and
parameters that built it. Because that graph now lives in a separate repository,
**this one must record the fingerprint** — in the build assets and in the VCF
header — or a merged VCF cannot be traced to the graph it came from.
`bin/stamp_build_id.sh` stamps a build id; the graph identity belongs with it.

## Isolate selection is yours, and it is not a formality

The criteria this pipeline was developed against: **paired-end only, at least
60x mean depth, no mixed samples.** Public sequence archives vary enormously.
A marginal isolate does not give a marginally worse answer — a low-depth sample
produces missing calls that score as method failures, and a mixed sample carries
real minority alleles at positions a single-strain analysis scores as errors.
Both corrupt precision and recall invisibly. Nothing here screens for you.

## Provenance

Extracted 2026-09-29 from a larger working repository, as the dependency closure
of `p0_prepare.sh` and `refbias_run.sh`. History before that date is not
included: it carried tracked alignments of 145-162 MB in a 3.2 GiB pack.
Convergence and association analysis is downstream of the VCF and is not here.
