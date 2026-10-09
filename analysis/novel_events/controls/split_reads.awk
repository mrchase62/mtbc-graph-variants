# Split-read counts from `samtools view -F 0x904` (primary alignments only).
# Inverted local chimera ("fold-back"): the first SA entry is on the same
# contig, the opposite strand, and starts within 1 kb. Other split: same
# contig otherwise. dup_like: same contig, start, CIGAR, strand and mate start
# as an earlier read. Rates per 1,000 primary alignments. Set -v s=<sample>.
{n++; k=$3":"$4":"$6":"and($2,16)":"$8; if(k in seen) d++; seen[k]=1}
/SA:Z:/{split($0,a,"SA:Z:"); split(a[2],b,","); st=and($2,16)?"-":"+"
        if(b[1]==$3 && b[3]!=st && b[2]-$4<1000 && $4-b[2]<1000) inv++; else if(b[1]==$3) oth++}
END{printf "%s\t%d\t%.4f\t%.3f\t%.3f\n", s, n, n?d/n:0, n?1000*inv/n:0, n?1000*oth/n:0}
