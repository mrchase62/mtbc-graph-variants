#!/usr/bin/env python3
"""Library type from adapter read-through (user, 2026-10-11): which samples
were made with Tn5 tagmentation (Nextera adapters) and which with TruSeq-style
adapters (sonication or enzymatic fragmentation; the two are not separated
here). Metadata (ENA library_selection) does not say: nearly all are RANDOM.

A read whose insert is shorter than the read runs into the adapter. The
first 12 bases of each adapter are counted anywhere in the read, both
orientations (reads on the reverse strand of a BAM/CRAM are stored reverse
complemented):
  Nextera / Tn5 mosaic end   CTGTCTCTTATA   (rc TATAAGAGACAG)
  TruSeq                     AGATCGGAAGAG   (rc CTCTTCCGATCT)
Per sample: the first --n reads (FASTQ: read 1 and read 2 alternately; BAM/
CRAM: primary alignments, which in a coordinate-sorted file are the start of
the chromosome; adapter read-through does not depend on position). Rates are
reads per million with each adapter, and the read-length median and the
fraction of reads shorter than the longest read (short inserts after
trimming). Raw FASTQ and trimmed CRAMs differ: production CRAMs were made
from fastp-trimmed reads, so only residual adapter is left there and rates
must be compared within a set, not across sets.

  adapter_count.py --sample S --fastq R1.fq.gz R2.fq.gz --out row.tsv
  adapter_count.py --sample S --cram X.cram --ref H37Rv.fasta --samtools ... --out row.tsv
"""
import argparse
import gzip
import itertools
import statistics
import subprocess

NX, TS = "CTGTCTCTTATA", "AGATCGGAAGAG"
RC = str.maketrans("ACGT", "TGCA")
NXR, TSR = NX.translate(RC)[::-1], TS.translate(RC)[::-1]


def fq(path):
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.rstrip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastq", nargs="+")
    ap.add_argument("--cram")
    ap.add_argument("--ref")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--n", type=int, default=200000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.fastq:
        its = [fq(p) for p in a.fastq]
        seqs = (s for grp in itertools.zip_longest(*its) for s in grp if s)
        src = "fastq"
    else:
        p = subprocess.Popen([a.samtools, "view", "-F", "0x900", "-T", a.ref, a.cram],
                             stdout=subprocess.PIPE, text=True)
        seqs = (line.split("\t", 10)[9] for line in p.stdout)
        src = "cram"
    n = nx = ts = 0
    lens = []
    for s in itertools.islice(seqs, a.n):
        n += 1
        lens.append(len(s))
        nx += (NX in s or NXR in s)
        ts += (TS in s or TSR in s)
    if not a.fastq:
        p.kill()
    mx = max(lens) if lens else 0
    nxr, tsr = (round(1e6 * x / n, 1) if n else 0 for x in (nx, ts))
    # call: one adapter at least 5 times the other and at least 100 per million
    call = ("nextera_tn5" if nxr >= 100 and nxr >= 5 * tsr else
            "truseq" if tsr >= 100 and tsr >= 5 * nxr else "unclear")
    with open(a.out, "w") as fo:
        fo.write("sample\tsource\treads\tnextera_per_M\ttruseq_per_M\tadapter_call\tread_len_median\t"
                 "read_len_max\tfrac_shorter_than_max\n")
        fo.write(f"{a.sample}\t{src}\t{n}\t{nxr}\t{tsr}\t{call}\t{statistics.median(lens) if lens else ''}\t"
                 f"{mx}\t{round(sum(x < mx for x in lens) / n, 4) if n else ''}\n")


if __name__ == "__main__":
    main()
