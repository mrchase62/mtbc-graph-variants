# Outgroup kit, route A (2026-10-13)

The event states of L8 (GCF_012923765) and M. canettii ET1291
(GCF_035581225), read from each assembly aligned straight to H37Rv. No graph.
Script: `bin/outgroup_kit.py` (method in its docstring). Local, ~6 s per genome.
SNP states already exist (`analysis/panel_checks/states/<g>.npz`).

    O=analysis/outgroup_kit/out; B=runroot/refbias/build/7713a8d71d8e-fix1/assets
    $MTB_PY bin/outgroup_kit.py --genome $MTB_ASSEMBLIES/<g>.fna.gz --name <g> \
      --h37rv $MTB_H37RV --is6110 $MTB_DATA/annotation/IS6110.fasta \
      --accessory $B/accessory_catalogue.fasta $B/accessory_catalogue.tsv --outdir $O

| | ET1291 | L8 |
|---|---:|---:|
| H37Rv callable (one alignment) | 4.18 Mb | 4.38 Mb |
| SNPs / indels <50 / >=50 inside alignments | 10,493 / 389 / 147 | 2,321 / 268 / 119 |
| SV rows (INS/DEL/REPL/REARR) | 114/34/78/1 | 70/52/16/1 |
| IS6110 copies: site H37Rv lacks / shared / unresolved or unplaced | 3 / 0 / 3 | 13 / 1 / 7 |

**Accessory check against the CX333 graph (which still holds both genomes),
fix1 catalogue:**
- novel loci: "present" matches a graph carrier at 21 of 23; "absent" matches
  a non-carrier at 62 of 62.
- mosaic loci (part H37Rv sequence): unreliable; ET1291 "present" at 30 loci
  the graph says it lacks. Treat mosaic as unknown until checked.
- copy_number loci: not callable by alignment (679).
- A first version without the position check called sequences present when
  they sat elsewhere in the genome (e.g. CIPT-assigned loci); fixed.

**Open:** the accessory catalogue is rebuilt with the new graph, so rerun the
accessory step then. No unit tests yet.
