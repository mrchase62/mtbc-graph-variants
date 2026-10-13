#!/usr/bin/env python3
"""B vs D (sparse mapping) follow-up: are the calls that differ real genotype
differences, or the same variant written differently, or in masked repeats?
Analysis only. For every genome and every ALT call present in one arm only:
  shifted  - the other arm has, in the same genome, a call of the same class
             and the same length change within 50 bp (representation shift);
  masked   - overlaps the build's repeat mask;
  other    - neither.
    python3 explain_B_D.py B.vcf.gz D.vcf.gz repeat_mask.bed > out
"""
import bisect, collections, gzip, sys
sys.path.insert(0, ".")
from compare_arms import cls, load

def mask_of(path):
    iv = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(path) if l.strip() and not l.startswith(("#", "track")))
    st = [s for s, _ in iv]; reach, m = [], 0
    for _, e in iv:
        m = max(m, e); reach.append(m)
    def hit(p, e):
        i = bisect.bisect_left(st, e) - 1
        while i >= 0 and reach[i] > p:
            if iv[i][1] > p and iv[i][0] < e:
                return True
            i -= 1
        return False
    return hit

b, d, bed = sys.argv[1:4]
hit = mask_of(bed)
sb, altB, missB, _, _ = load(b, "NONE")
sd, altD, missD, _, _ = load(d, "NONE")
tot = collections.Counter(); per = collections.Counter()
for s in sorted(set(sb) & set(sd)):
    for tag, only, other in (("onlyB", altB[s] - altD[s], altD[s]), ("onlyD", altD[s] - altB[s], altB[s])):
        idx = collections.defaultdict(list)
        for k in other:
            idx[(cls(k[1], k[2]), len(k[2]) - len(k[1]))].append(k[0])
        for v in idx.values():
            v.sort()
        for k in only:
            c = cls(k[1], k[2]); dl = len(k[2]) - len(k[1])
            ps = idx.get((c, dl), [])
            i = bisect.bisect_left(ps, k[0] - 50)
            if i < len(ps) and ps[i] <= k[0] + 50:
                why = "shifted"
            elif hit(k[0], k[0] + max(len(k[1]), 1)):
                why = "masked"
            else:
                why = "other"
            tot[tag, c, why] += 1
            if why == "other":
                per[s] += 1
print("class  arm    shifted  masked  other")
for c in ("SNP", "MNP", "indel", "SV"):
    for t in ("onlyB", "onlyD"):
        print(f"{c:<6} {t:<6} {tot[t,c,'shifted']:>7} {tot[t,c,'masked']:>7} {tot[t,c,'other']:>6}")
print("genomes with most 'other' differences: " + ", ".join(f"{s} {n}" for s, n in per.most_common(6)))
