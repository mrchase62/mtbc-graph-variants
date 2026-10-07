# Building a panel: Segment 1, corrected

This replaces RUNBOOK.md Segment 1 and QC_PIPELINE.md Stages 0 to 2 (working
tree) for the next panel. Those documents disagree with each other and with the
code; the 2026-10-05 audit lists the differences (`panel_graph_build.md`
PGB-3, -4, -5, -7, -8, -13). Where this file and those disagree, this file
matches the code in `bin/`.

Every command is run from the repository root after
`source config/project_env.sh`. `Q=data/qc` below. Nothing here chooses the
panel: each screen writes a table, a person records decisions in
`exclusions` / `retained` files, and `build_panel.py` applies them.

## 1a. Candidates from RefSeq

```bash
bash bin/refresh_assemblies.sh            # or --summary <archived summary> --summary-only
```

Selection: `Complete Genome`, `latest`, and species taxid 1773 / 78331 /
1305738 **or** organism name `^Mycobacterium (tuberculosis|canetti|orygis)`.
The old name regex dropped both complete *M. orygis* genomes (GCF_015265495,
GCF_033782915); on the 2026-09-07 summary the new rule selects 516 (the old
514 plus those two). Any complete genome with an MTBC-like name that is not
selected is printed for review. GenBank-only genomes are still not searched.

Rotation (`rotate_to_dnaa.sh`) is unchanged.

```bash
ls data/rotated/*.dnaA_rotated.fasta | xargs -n1 basename \
    | sed 's/\.dnaA_rotated\.fasta$//' > $Q/accessions.txt
```

## 1b. Basic statistics and single-base runs

```bash
$MTB_PY bin/assembly_qc_stats.py --dir data/rotated --suffix .dnaA_rotated.fasta \
    --out $Q/assembly_qc_stats.tsv
```

(RUNBOOK's command pointed `--dir data/rotated` at the default `.fna.gz`
suffix and would read nothing.) New columns `longest_base_run*`,
`base_runs_ge100`, and `flag = single_base_run` for any pure one-base run of
100 bp or more. On CX333's 333 panel sequences it flags exactly
GCF_050259585 (1,015 bp G), GCF_045348265 (484 bp A) and GCF_039770655
(377 bp T), and nothing else.

## 1c. Collinearity

Unchanged: `bin/collinearity_qc.py` as in RUNBOOK 1c.

## 1d. Provenance (the rule that built CX333)

```bash
$MTB_PY bin/assembly_provenance_screen.py --accessions $Q/accessions.txt \
    --assembly-dir data/rotated --ref data/ref/H37Rv.fasta \
    --is6110 data/annotation/IS6110.fasta \
    --assembly-summary data/ncbi/assembly_summary_refseq.txt \
    --minimap2 "$MTB_MINIMAP2" --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" \
    --out $Q/assembly_provenance.tsv
```

`flag`: `REFERENCE_STRUCTURE` when the IS6110 profile is identical to
H37Rv's; else `SHORT_READ_SUSPECT` when the NCBI technology is short-read
only; else `ok`. This is production's rule (the code previously excluded on
fewer than 5 insertions of 50 bp or more, which production never used). The
insertion count is still reported (`insertions_ge50`) and gives
`tech_from_sequence`; `review` lists `unknown_technology`,
`short_read_assembler` (SOAPdenovo and the like on a long-read label) and
metadata/sequence technology conflicts. These notes do not change the flag;
they need a person.

Use `--metadata <previous table>` to rerun on recorded metadata without
querying NCBI. `--no-insertion-count` skips the whole-genome alignment.

Checked on the 484 first-round genomes (recorded metadata,
`--no-insertion-count`): the flag matches production for 483 of 484 (95
REFERENCE_STRUCTURE, 55 SHORT_READ_SUSPECT, 334 ok; production 94 / 56 /
334). The exception, GCF_002886775, moves from SHORT_READ_SUSPECT to
REFERENCE_STRUCTURE; it is excluded either way. `tech_class` matches for
483 (GCF_000195835, "MiSeq/PacBio RSII", is now hybrid; still ok).
`short_read_assembler` fires on exactly GCF_050259585 (SOAPdenovo) and
GCF_000195835 (mummer/soapdenovo).

## 1e. Foreign / engineered DNA

```bash
$MTB_PY bin/foreign_insertion_screen.py --accessions $Q/accessions.txt \
    --assembly-dir data/rotated --ref data/ref/H37Rv.fasta \
    --lineages $Q/lineages.all.tsv --is6110 data/annotation/IS6110.fasta \
    --univec "$MTB_UNIVEC" --blastn "$MTB_BLASTN" --makeblastdb "$MTB_MAKEBLASTDB" \
    --minimap2 "$MTB_MINIMAP2" --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" \
    --threads 8 --workdir $Q/work_foreign --out $Q/foreign_insertions.tsv
```

This is the reconciled command. RUNBOOK.md passed the whole panel FASTA as the
background, which QC_PIPELINE.md 1.4 says must not be done (1,167 false
foreign calls on the external assemblies). `--lineages` builds the background
in the workdir: one genome per sublineage plus every genome with no call,
PanSN-named. A hand-made `--background` must be PanSN-named and at most
`--max-background` (120) genomes. Settings: `-cx asm10 --secondary=yes -N 100
-p 0.05` (QC_PIPELINE.md 1.4; the code had `-N 50 -p 0.1`).

Changes from the CX333 run: inserts between alignment blocks are extracted
(`source = gap`); the score is the best single other genome, not the sum;
an insert carried by only one other genome is `REVIEW_one_homologue`, so two
genomes with the same construct no longer pass each other; every non-native
insert is re-checked with the frequency filter off and against IS6110
(`final_verdict`: `native_is6110`, `native_recheck`). Read `final_verdict`.

Every insert, native or not, is also searched against NCBI's UniVec_Core with
VecScreen's blastn settings (D29; fetch it once with `bash
bin/fetch_univec.sh`, which records the build and checksum). A strong VecScreen
match makes `final_verdict` `VECTOR`; `univec_strong_bp` and `univec_hit` say
how much and what. On CX333's 449 inserts and the external assemblies' 2,254
it flags exactly the two known constructs (GCF_044324775, pJEB; GCF_021535155,
attB vector) and nothing else. Like every screen here it flags for review and
does not exclude by itself (D31).
M. canettii's divergent sequence still reads FOREIGN; review it by hand.

## 1f. Lineage

Unchanged: `tbprofiler_assemblies.sh` then `tbprofiler_collect.py` to
`$Q/lineages.all.tsv`. Mixed lineage within one root is still not detected
(PGB-16).

## 1g. SNP, indel and homopolymer-indel outliers

```bash
$MTB_PY bin/variant_counts.py --accessions $Q/accessions.txt \
    --assembly-dir data/rotated --h37rv data/ref/H37Rv.fasta \
    --minimap2 "$MTB_MINIMAP2" --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" \
    --threads 8 --out $Q/variants
$MTB_PY bin/snp_outlier_screen.py --counts $Q/variants/counts.snp.tsv \
    --accessions $Q/accessions.txt --lineages $Q/lineages.all.tsv \
    --out $Q/snp_outliers.tsv
$MTB_PY bin/snp_outlier_screen.py --counts $Q/variants/counts.indel.tsv \
    --accessions $Q/accessions.txt --lineages $Q/lineages.all.tsv \
    --column indel1 --out $Q/indel_outliers.tsv
$MTB_PY bin/snp_outlier_screen.py --counts $Q/variants/counts.hpindel.tsv \
    --accessions $Q/accessions.txt --lineages $Q/lineages.all.tsv \
    --column indel1_hp_private --out $Q/hpindel_outliers.tsv
bash bin/mash_nearest.sh data/rotated $Q/mash_nearest.tsv
```

`variant_counts.py` is the producer the SNP screen never had (RUNBOOK's
`snp_counts.tsv` did not exist; CX333's counts came from the 2025 graph and
covered 254 of 333). It aligns every candidate directly to H37Rv. Run all
candidates in one call: "private" homopolymer indels are private within the
run. The 1 bp indel and private homopolymer-indel screens are the long-read
error checks of the external-assembly QC. A sublineage with MAD 0 falls back to
its major lineage; if that is also 0, a genome off the median is `NO_SPREAD`.

## 1h. Decide, then write the panel

Record each decision, with its reason, in a TSV (`accession`, `reason`):
`$Q/manual_exclusions.tsv` (SNP excess and the like), further
`--exclusions` files (e.g. `engineered_exclusions.tsv`), and
`$Q/panel_retained_exceptions.tsv` for genomes kept against a flag.

```bash
$MTB_PY bin/build_panel.py --manifest data/ncbi/panel_manifest.tsv \
    --collinearity $Q/collinearity_qc.tsv \
    --summary data/ncbi/assembly_summary_refseq.txt \
    --barcode $Q/lineages.all.tsv \
    --manual-exclusions $Q/manual_exclusions.tsv \
    --exclusions $Q/engineered_exclusions.tsv \
    --provenance $Q/assembly_provenance.tsv \
    --retained $Q/panel_retained_exceptions.tsv \
    --review $Q/assembly_qc_stats.tsv --review $Q/snp_outliers.tsv \
    --review $Q/indel_outliers.tsv --review $Q/hpindel_outliers.tsv \
    --out $Q/panel.<name>.tsv --excluded $Q/panel.<name>.excluded.tsv
bash bin/make_pansn_fasta.sh --panel $Q/panel.<name>.tsv --out mtb.<name>
```

`build_panel.py` writes a list, never the FASTA, and refuses an `--out` that
is FASTA-named or one of its inputs. It stops if a genome reaching the
provenance rule has no provenance row, or if a `--review` flag (a non-empty
`flag` column) has no recorded decision. It writes `<out>.inputs.tsv` with
every input's sha256. `--barcode` is required. `make_pansn_fasta.sh` needs
`--panel` and `--out` and never overwrites an existing `mtb.<name>.fasta.gz`.

FOREIGN or REVIEW inserts in `foreign_insertions.tsv` are reviewed by hand,
and each genome kept or excluded is recorded the same way.

**CX333 from its recorded decisions.** The command above, with the CX333-era
`assembly_provenance.tsv` and without the `--review` tables (they did not
exist), gives 333 genomes identical to `qc_master_table.tsv` `in_panel` and to
the paths of `mtb.complex333.fasta.gz`, and the same exclusion reason for all
181 excluded. `--no-provenance` in place of the provenance, exclusion and
retained inputs gives the first-round 484 of `panel.rebuild.tsv`.

## 1i. Join everything

```bash
$MTB_PY bin/qc_master_table.py --qc-dir $Q --ncbi-dir data/ncbi \
    --summary data/ncbi/assembly_summary_refseq.txt \
    --panel-fasta data/fastas/mtb.<name>.fasta.gz --out $Q/qc_master_table.tsv
```

Run after the panel is built (RUNBOOK ran it before 1i).

## Validation on CX333 (2026-10-05, local, read-only inputs)

| screen | what it reproduces |
|---|---|
| `refresh_assemblies.sh --summary <2026-09-07 summary>` | 516 selected: the 514 of `panel_manifest.tsv` plus GCF_015265495 and GCF_033782915 (M. orygis) |
| `build_panel.py` with the recorded decisions | 333 genomes, identical to `in_panel` and to the panel FASTA's paths; all 181 exclusion reasons identical to `qc_master_table.tsv` |
| `build_panel.py --no-provenance` | the 484 of `panel.rebuild.tsv` |
| `assembly_provenance_screen.py` | flag identical for 483 of 484 (see 1d) |
| `assembly_qc_stats.py` on the 333 panel sequences | `single_base_run` on exactly GCF_039770655, GCF_045348265, GCF_050259585 |
| `foreign_insertion_screen.py`, 7 genomes against a 71-genome PanSN background (the external QC's 69 one-per-sublineage genomes plus both attB-vector genomes) | 121 inserts (37 from alignment gaps). First pass 51 FOREIGN, 2 REVIEW, 68 native; after the re-check 48 of the 51 are IS6110. Left: GCF_044324775 and GCF_021535155 (the attB vectors, REVIEW_one_homologue: each found only in the other; the old code called both native), and the three single-base-run artifacts FOREIGN, including GCF_039770655's 1,144 bp insert at 2,001,528, which only the gap extraction finds |
| `variant_counts.py` on 3 genomes | snps, indel1 and indel1_hp identical to the external QC's counts; private counts differ because privacy depends on the set run |

The long-read checks on CX333's 332 non-reference genomes (plus M. orygis
GCF_015265495, which was in that run), from the external
QC's counts (same algorithm, privacy over 487 genomes) and
`snp_outlier_screen.py` with production lineages: 1 bp indel EXCESS 15,
private homopolymer-indel EXCESS 16 (GCF_002116815 121, GCF_008761675 106,
GCF_039906425 51, GCF_001544705 42, GCF_965124535 32, and 11 more), SNP
EXCESS 6 and DEFICIT 11 (the lineage 4.9 H37Rv-like genomes read as DEFICIT
against major lineage 4). Run on the next candidate set, these are review
flags; none was a CX333 decision.
