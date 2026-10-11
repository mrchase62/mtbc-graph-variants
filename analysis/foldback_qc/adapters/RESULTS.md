# Adapter read-through vs fold-back rate (2026-10-11; report only)

Method: `adapter_count.py`, first 200,000 reads per sample. A read counts if
it contains the first 12 bases of the Nextera/Tn5 adapter (CTGTCTCTTATA) or
the TruSeq adapter (AGATCGGAAGAG), in either orientation.

The call is one adapter at >= 100 reads per million and >= 5 times the
other; "unclear" otherwise. No sample had both at >= 100/M; the two adapter
types separate cleanly.

Comparability:
- Phase B: raw FASTQ.
- Production: CRAMs made from fastp-trimmed reads, so only residual adapter
  is left. Compare rates within a set only.

TruSeq-style adapters are used by both sonication and enzymatic
(nuclease) kits; this does not separate them.
"Unclear" means little read-through: long inserts, or reads trimmed before
upload.

## Phase B (62 raw runs; fold-back v2 on raw reads)
- Tn5: 29. TruSeq: 20, all MtbSetB HiSeq 2500. Unclear: 13, MtbSetA HiSeq 2500.
- All 11 chimeric controls (6.45-64.7 per 1,000) are Tn5 libraries. 10 of
  them have 3,500-142,000 adapter reads per million, i.e. many inserts
  shorter than the read. The exception is RW_TB008 (MiSeq, 230/M).
- Clean groups:
  - TruSeq set B: 0.01-0.07;
  - unclear set A: 0.01-0.39;
  - Tn5 NextSeq set E, about 300/M read-through: 0.33-1.55.
- One Tn5 run with heavy read-through is clean: N1177 (111,000/M, 0.03).
  So short Tn5 inserts do not by themselves make fold-backs.
- In this study, library type and study are the same split.

## Production (500 random samples; old fold-back count per 1,000)

| adapter call | n | <=0.3 | 0.3-1 | 1-2 | 2-5 | >5 | median |
|---|---|---|---|---|---|---|---|
| Tn5 | 276 | 124 | 24 | 18 | 37 | 73 | 0.44 |
| TruSeq | 75 | 36 | 6 | 9 | 7 | 17 | 0.52 |
| unclear | 149 | 88 | 20 | 11 | 10 | 20 | 0.16 |

- High fold-back occurs in both library types (> 5: Tn5 26%, TruSeq 23%), so
  Tn5 alone does not explain it.
- **Within TruSeq, fold-back tracks short inserts.** Spearman with residual
  adapter rate is 0.61; with the fraction of reads trimmed shorter than the
  longest read, 0.75. By adapter-rate tercile, fold-back > 5 is 0/25, 2/25
  and 15/25. This is the pattern expected if fragmentation that makes short
  fragments also makes inverted-repeat chimeras (sonication or nuclease;
  Haile et al. 2019). It is not proven here.
- **Within Tn5**, residual adapter rate does not track fold-back (rho 0.07).
  The trimmed-read-length fraction does weakly (0.29).
- Instruments:
  - Tn5 is the majority on NextSeq 500 (108 of 131), HiSeq 2500 (62) and
    MiSeq (54).
  - TruSeq is spread across HiSeq 2000 (24), NextSeq 2000 (17), MiSeq (15)
    and others.

## Use
- The adapter call and rates are columns of the selfcal calibration table
  (`selfcal/collect.py`), so the false-call curve can be split by library type.
- Short-insert fraction looks like a candidate QC feature next to the
  fold-back rate. It is cheap, but needs the insert size from paired reads to
  be firm. The calibration table records the insert mean and sd from the P2 BAMs.
