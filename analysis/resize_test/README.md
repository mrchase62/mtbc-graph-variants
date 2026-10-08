# Smaller P1/P2 jobs give identical results (2026-10-07)

Analysis only. The user approved cutting P1 and P2 from 8 cores / 16 GB to
4 cores / 8 GB per sample, with three supporting changes:

1. `bwa mem -K 80000000` in `bin/simulate_and_call.sh`. bwa estimates the
   insert-size distribution per batch, and its batch size grows with the
   thread count, so alignments could change with -t. -K fixes the batch at the
   size 8 threads used.
2. P1 passes its core count to bwa and HaplotypeCaller (`MTB_THREADS`), as
   P2 already did. Before, P1 ran them at 8 threads whatever the job had.
3. HaplotypeCaller's Java heap `-Xmx8g` → `-Xmx6g` (`MTB_GATK_XMX`), below the
   job's 8 GB.

**Test:** scale200's largest sample, SAMN12126251 (2,143,912 read pairs),
through P1's H37Rv alignment and calling (`run_arm.sbatch`, `compare.sh`):

| run | job | alignments (md5 of every read record) | calls (md5 of every VCF record) | wall time | peak memory (sacct) |
|---|---|---|---|---|---|
| production P1, 2026-09-23 | 8 cores, 16 GB | a4efe76063e5 | 26b20bc07380 (1,210 records) | | |
| old script, 8 threads | 8 cores, 16 GB | a4efe76063e5 | 26b20bc07380 | 3:24 | 2.2 GB |
| new script, 8 threads | 8 cores, 16 GB | a4efe76063e5 | 26b20bc07380 | 3:13 | 1.4 GB |
| **new script, 4 threads** | **4 cores, 8 GB** | **a4efe76063e5** | **26b20bc07380** | 4:04 | 2.2 GB |

All four are identical, read for read and call for call. The old-script run
stopped at its last step (stamping the build into the VCF header) because the
copied script could not find its helper; that step touches header lines only,
which the comparison excludes.

**Cost:** the 4-core job ran 20% longer, but is billed at 6 units an hour
instead of 12 (one per core, a quarter per GB), so about 40% less per sample.
P1 + P2 for scale200: about 300 → about 180 billing-hours.

**Not tested here:** P2's extra steps (GVCF, delly, dysgu). Its measured
peak on gwas1000 was 6.3 GB, under the new 8 GB. Test cost: under 1
billing-hour.
