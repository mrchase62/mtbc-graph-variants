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

## PPE38 locus check (2026-10-13, user request)

CDC1551 (GCF_000008585) carries PPE38, two esx genes and PPE71 at
2,628,768-2,632,079. By tblastn the esx genes are esxN-like (2,630,145-429)
and esxJ-like (2,630,444-758), presumably esxX/esxY; H37Rv has only PPE38
there. Each genome's genes were found by blastn of the CDC1551 copies, and
their flanks were placed on H37Rv with minimap2 (MAPQ 60):

- **ET1291:** carries all four (esx 100%, PPE copies 99.7%). The ~3.3 kb block
  sits exactly in place of H37Rv PPE38: flanks end at 2,632,922 and resume at
  2,634,193. Its esxN placed at H37Rv's esxN (2,030,347-2,030,977).
- **L8:** has none of the four. H37Rv 2,627,050-2,637,333 (plcC, plcB, plcA,
  PPE38, PPE39) is replaced by one IS6110 copy; this is the kit's REPL row.
  Its esxN placed at H37Rv's esxN (2,030,347-2,030,977).
- The kit's alignment is right at all of these. A quick nearest-alignment
  projection put L8's esxN 13 kb off, because of L8's 3.1 kb deletion at
  2,025,302. Use flank placement, not `project()`, for exact gene positions.
