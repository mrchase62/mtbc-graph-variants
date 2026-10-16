#!/usr/bin/env python3
"""Extract SV-like records (CLASS sv / is6110 / accessory_presence, plus
CLASS=small indels >= 50 bp, which OLD kept as small and NEW folds into
CLASS=sv SVSOURCE=small_allele) from a merged cohort VCF into a pickle:
per record a dict of fields, per sample ST and EV codes.
    extract.py <vcf.gz> <out.pkl>
"""
import gzip, pickle, sys
import numpy as np
ST = {"REF": 0, "ALT": 1, "ABSENT": 2, "NOCALL": 3, "DUP": 4}
vcf, out = sys.argv[1:3]
fh = gzip.open(vcf, "rt")
for l in fh:
    if l.startswith("#CHROM"):
        samples = l.rstrip("\n").split("\t")[9:]; break
recs, sts, evs = [], [], []
evcode = {".": 0, "": 0}
for l in fh:
    f = l.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
    cls = info.get("CLASS")
    ref, alt = f[3], f[4]
    if cls == "small":
        if "," in alt or max(len(ref), len(alt)) < 50 or abs(len(ref) - len(alt)) < 50: continue
    fmt = f[8].split(":")
    ist = fmt.index("ST"); iev = fmt.index("EV") if "EV" in fmt else None
    s = np.full(len(samples), 9, np.int8); e = np.zeros(len(samples), np.int16)
    for j, c in enumerate(f[9:]):
        x = c.split(":")
        s[j] = ST.get(x[ist], 9) if len(x) > ist else 9
        if iev is not None and len(x) > iev:
            e[j] = evcode.setdefault(x[iev], len(evcode))
    pos = int(f[1])
    if alt.startswith("<") or ref.startswith("<"):
        seq = None
    else:
        seq = alt[1:] if len(alt) > len(ref) else None
    if cls == "small":
        svtype = "INS" if len(alt) > len(ref) else "DEL"
        svlen = abs(len(alt) - len(ref)); end = pos + len(ref) - 1 if svtype == "DEL" else pos
    else:
        svtype = info.get("SVTYPE", ".")
        svlen = abs(int(info.get("SVLEN", info.get("ACCLEN", "0")) or 0))
        end = int(info["END"]) if "END" in info else (pos + len(ref) - 1 if svtype == "DEL" else pos)
    info.pop("SMALLALLELE", None)
    recs.append(dict(chrom=f[0], pos=pos, end=end, svtype=svtype, svlen=svlen, cls=cls,
                     seq=seq, reflen=len(ref), info=info))
    sts.append(s); evs.append(e)
inv = {v: k for k, v in evcode.items()}; inv[0] = "."
pickle.dump(dict(samples=samples, recs=recs, st=np.array(sts), ev=np.array(evs), evnames=inv),
            open(out, "wb"), protocol=4)
print(len(recs), "records", len(samples), "samples")
