#!/usr/bin/env python3
"""Read check of the sampled OLD/NEW disagreement cells (out/readcheck_sample.tsv)
against each isolate's own P1 alignment to H37Rv
(runroot/refbias/work/new1k_p200_p1/<s>.h37rv.bam). Reads: primary or
supplementary, not duplicate / QC-fail, MAPQ >= 20. Each cell is checked at the
coordinates of the build that called ALT (OLD record for 'old_alt' cells, NEW
record for 'new_alt' cells), so the verdict says whether that ALT call holds.

DEL [b, e] = deleted bases, L = e - b + 1
  depth ratio = mean depth inside / mean depth of the flanks (F = clamp(L, 100, 500))
  alt reads   = reads with a CIGAR deletion >= 0.6 L starting within 50 bp of b,
                or soft-clipped >= 15 bp with the clip point within 15 bp of
                b-1 (right clip) or e+1 (left clip)
  ref reads   = reads crossing the left (b) or right (e) breakpoint by >= 15 bp
                each side with no clip there and no insertion/deletion >= 15 bp
                anywhere in the read (mean of the two breakpoints)
  ALT  if alt >= 3 and alt/(alt+ref) >= 0.5 (breakpoint reads), or flank
       depth >= 5, ratio < 0.2, ref <= 2 and the inside is not a repeat
       (MAPQ>=0 depth > 3x MAPQ>=20 depth) (depth)
  REF  if (alt <= 1 or alt/(alt+ref) < 0.1) and ref >= 8 (and ratio > 0.6 for
       L >= 300), or L >= 300, flank >= 8, ratio > 0.6 and alt <= 2
  unclear otherwise
INS at POS, length L, inserted sequence S where the record carries one
  alt reads   = reads with a CIGAR insertion >= min(0.5 L, 30) within 30 bp of
                POS, or soft-clipped >= 15 bp with the clip point within 30 bp
                (when S is known: only clips whose 15 clipped bases next to
                the clip point occur in S, or in S + 30 bp of reference on
                either side, which covers tandem duplications)
  ref reads   = reads crossing POS-30..POS+30 (and the duplication junction)
                with no clip there and no insertion/deletion >= 15 bp anywhere
  tandem duplications (S >= 90% identical to the L bp beside POS): the
                junction at POS+L / POS-L and insertions of 0.8-1.25 L within
                L+30 bp also count; for L >= 200 the depth of the duplicated
                copy vs flanks (>= 1.6 -> ALT)
  ALT if alt >= 3 and alt/(alt+ref) >= 0.3; REF if alt <= 1 and ref >= 8
IS6110 at POS
  junction reads = soft-clipped >= 15 bp within 20 bp of POS whose 15 clipped
                bases next to the clip point are in IS6110 (either strand), or
                whose SA tag places a part inside an H37Rv IS6110 copy
  ref reads   = reads crossing POS-20..POS+20 unclipped
  ALT if junction >= 3 and junction/(junction+ref) >= 0.2; REF if junction <= 1
  and ref >= 8 (distance to the nearest H37Rv IS6110 copy is reported)
Accessory presence (bonus, ALT<->NOCALL cells)
  reads that are unmapped, soft-clipped >= 20 bp, MAPQ < 20 or NM >= 4 in the
  P1 BAM are re-aligned (bwa mem) to H37Rv + the locus representative; breadth
  (>= 2 reads MAPQ >= 20) over the locus positions that 150-mers tiled from
  the locus place uniquely on it. present if >= 0.7, absent if <= 0.15
    readcheck.py [--acc-threads N]
"""
import argparse, collections, os, pickle, subprocess, sys, tempfile
import numpy as np, pandas as pd, pysam
R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
D = f"{R}/analysis/sv_eval_1k/out"
BAMD = f"{R}/runroot/refbias/work/new1k_p200_p1"
B = f"{R}/runroot/refbias/build/3893a3762ed1"
REF = f"{B}/refs/GCF_000195955.fasta"
CH = "NC_000962.3"
BWA = os.environ.get("MTB_BWA", "bwa")
fa = pysam.FastaFile(REF)
rname = fa.references[0]
G = fa.fetch(rname)
IS = [(int(l.split("\t")[3]), int(l.split("\t")[4])) for l in open(f"{B}/assets/is6110_elements.gff")
      if "\tmobile_genetic_element\t" in l]
def rc(s): return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
E = G[IS[0][0] - 1:IS[0][1]]
K = 15
ISK = {E[i:i + K] for i in range(len(E) - K + 1)} | {rc(E)[i:i + K] for i in range(len(E) - K + 1)}


def ok(r):
    return not (r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail) and r.mapping_quality >= 20


def clips(r):
    """(left_clip_len, right_clip_len) of soft clips."""
    c = r.cigartuples
    return (c[0][1] if c[0][0] == 4 else 0), (c[-1][1] if c[-1][0] == 4 else 0)


def big_ops(r, op, minlen):
    """[(ref_start_1based, len)] of CIGAR ops `op` (1 = I, 2 = D) >= minlen."""
    out, p = [], r.reference_start + 1
    for o, l in r.cigartuples:
        if o == op and l >= minlen: out.append((p, l))
        if o in (0, 2, 3, 7, 8): p += l
    return out


def depth(bam, a, b, q):
    cov = bam.count_coverage(CH, a - 1, b, quality_threshold=0,
                             read_callback=lambda r: not (r.is_unmapped or r.is_secondary or r.is_duplicate
                                                          or r.is_qcfail) and r.mapping_quality >= q)
    return np.sum(cov, axis=0)


def check_del(bam, b, e):
    L = e - b + 1
    F = int(min(500, max(100, L)))
    if L > 30000:
        ins = np.mean([depth(bam, x, x, 20)[0] for x in np.linspace(b, e, 200).astype(int)])
        ins0 = np.mean([depth(bam, x, x, 0)[0] for x in np.linspace(b, e, 200).astype(int)])
    else:
        ins = depth(bam, b, e, 20).mean(); ins0 = depth(bam, b, e, 0).mean()
    fl = np.concatenate([depth(bam, max(1, b - F), b - 1, 20), depth(bam, e + 1, e + F, 20)]).mean()
    alt, refl, refr = set(), 0, 0
    for a_, z_ in ((b - 300, b + 300), (e - 300, e + 300)):
        for r in bam.fetch(CH, max(0, a_), z_):
            if not ok(r): continue
            lc, rcl = clips(r)
            s1, e1 = r.reference_start + 1, r.reference_end
            if any(abs(p - b) <= 50 and l >= 0.6 * L for p, l in big_ops(r, 2, 15)): alt.add(r.query_name)
            if rcl >= 15 and abs(e1 - (b - 1)) <= 15: alt.add(r.query_name)
            if lc >= 15 and abs(s1 - (e + 1)) <= 15: alt.add(r.query_name)
            if r.is_supplementary: continue
            near_d = bool(big_ops(r, 2, 15)) or bool(big_ops(r, 1, 15))
            if not near_d and (a_ == b - 300) and s1 <= b - 15 and e1 >= b + 15 and not (
                    (rcl >= 15 and abs(e1 - b) < 30) or (lc >= 15 and abs(s1 - b) < 30)):
                refl += 1
            if not near_d and (a_ == e - 300) and s1 <= e - 15 and e1 >= e + 15 and not (
                    (rcl >= 15 and abs(e1 - e) < 30) or (lc >= 15 and abs(s1 - e) < 30)):
                refr += 1
    # reads crossing both windows were counted per breakpoint; alt by name
    ref = (refl + refr) / 2
    ratio = ins / fl if fl > 0 else float("nan")
    na = len(alt)
    fr = na / (na + ref) if na + ref else 0
    rep = ins0 > 3 * max(ins, 0.5) and ins0 > 0.2 * fl
    if na >= 3 and fr >= 0.5: v, why = "ALT", "breakpoint_reads"
    elif fl >= 5 and ratio < 0.2 and not rep and ref <= 2: v, why = "ALT", "depth"
    elif (na <= 1 or fr < 0.1) and ref >= 8 and (L < 300 or ratio > 0.6): v, why = "REF", "spanning_reads"
    elif L >= 300 and fl >= 8 and ratio > 0.6 and na <= 2: v, why = "REF", "depth"
    elif fl < 5 and ins0 < 5: v, why = "unclear", "low_depth"
    elif rep: v, why = "unclear", "repeat_mapq0"
    else: v, why = "unclear", "mixed"
    return dict(verdict=v, why=why, L=L, depth_in=round(float(ins), 1), depth_in_mq0=round(float(ins0), 1),
                depth_flank=round(float(fl), 1), ratio=round(float(ratio), 3), alt_reads=na, ref_reads=ref)


def ident(a, b): return sum(x == y for x, y in zip(a, b)) / max(1, len(a))


def check_ins(bam, pos, L, seq):
    """Tandem duplications (S >= 90% identical to the L bases right or left
    of POS) have their junction L bp away from POS and may be written
    anywhere along the repeat: there an insertion of 0.8-1.25 L within L + 30
    bp, clips at POS+L / POS-L, and depth of the duplicated copy (>= 1.6x the
    flanks for L >= 200) also count."""
    w = 30
    ctx = None
    dup = None
    if seq:
        ctx = G[max(0, pos - 31 - L):pos] + seq + G[pos:pos + 30 + L]
        sk = {ctx[i:i + K] for i in range(len(ctx) - K + 1)}
        if ident(seq, G[pos:pos + L]) >= 0.9: dup = "right"
        elif ident(seq, G[pos - L:pos]) >= 0.9: dup = "left"
    W = w + (L if dup else 0)
    minI = min(max(0.5 * L, 15), 30)
    junc = [pos] + ([pos + L] if dup == "right" else [pos - L] if dup == "left" else [])
    alt, ref, clipped_any = set(), 0, set()
    for r in bam.fetch(CH, max(0, pos - 400 - L), pos + 400 + L):
        if not ok(r): continue
        lc, rcl = clips(r)
        s1, e1 = r.reference_start + 1, r.reference_end
        q = r.query_sequence or ""
        bi = big_ops(r, 1, minI)
        if any(abs(p - pos) <= W and (not dup or 0.8 * L <= l <= 1.25 * L) for p, l in bi): alt.add(r.query_name)
        cl = []
        for j in junc:
            if rcl >= 15 and abs(e1 - j) <= w: cl.append(q[len(q) - rcl:len(q) - rcl + K])
            if lc >= 15 and abs(s1 - j - 1) <= w: cl.append(q[lc - K:lc])
        for c in cl:
            clipped_any.add(r.query_name)
            if ctx is None or c in sk: alt.add(r.query_name)
        if r.is_supplementary: continue
        if s1 <= min(junc) - w and e1 >= max(junc) + w and not (lc >= 15 and min(abs(s1 - j) for j in junc) < 2 * w) \
                and not (rcl >= 15 and min(abs(e1 - j) for j in junc) < 2 * w) \
                and not big_ops(r, 1, 15) and not big_ops(r, 2, 15):
            ref += 1
    dr = float("nan")
    if dup and L >= 200:
        a_, b_ = (pos + 1, pos + L) if dup == "right" else (pos - L + 1, pos)
        inside = depth(bam, a_, b_, 20).mean()
        fl = np.concatenate([depth(bam, max(1, a_ - 300), a_ - 1, 20), depth(bam, b_ + 1, b_ + 300, 20)]).mean()
        dr = inside / fl if fl else float("nan")
    na = len(alt)
    fr = na / (na + ref) if na + ref else 0
    if (na >= 3 and fr >= 0.3) or (dr == dr and dr >= 1.6): v = "ALT"
    elif na <= 1 and ref >= 8 and not (dr == dr and dr >= 1.3): v = "REF"
    else: v = "unclear"
    return dict(verdict=v, why="reads" + ("_seq" if seq else "") + (f"_dup{dup}" if dup else ""), L=L, alt_reads=na,
                ref_reads=ref, clipped_any=len(clipped_any), alt_frac=round(fr, 3), dup_depth_ratio=round(dr, 3))


def in_is(p):
    return any(a - 50 <= p <= b + 50 for a, b in IS)


def check_is(bam, pos):
    w = 20
    junc, ref = set(), 0
    for r in bam.fetch(CH, max(0, pos - 400), pos + 400):
        if not ok(r): continue
        lc, rcl = clips(r)
        s1, e1 = r.reference_start + 1, r.reference_end
        q = r.query_sequence or ""
        hit = False
        if rcl >= 15 and abs(e1 - pos) <= w and q[len(q) - rcl:len(q) - rcl + K] in ISK: hit = True
        if lc >= 15 and abs(s1 - pos - 1) <= w and q[lc - K:lc] in ISK: hit = True
        if (lc >= 15 and abs(s1 - pos - 1) <= w) or (rcl >= 15 and abs(e1 - pos) <= w):
            if r.has_tag("SA"):
                for sa in r.get_tag("SA").rstrip(";").split(";"):
                    f = sa.split(",")
                    if f[0] == CH and in_is(int(f[1])): hit = True
        if hit: junc.add(r.query_name)
        if r.is_supplementary: continue
        if s1 <= pos - w and e1 >= pos + w and not (lc >= 15 and abs(s1 - pos) < 2 * w) and not (
                rcl >= 15 and abs(e1 - pos) < 2 * w):
            ref += 1
    nj = len(junc)
    fr = nj / (nj + ref) if nj + ref else 0
    if nj >= 3 and fr >= 0.2: v = "ALT"
    elif nj <= 1 and ref >= 8: v = "REF"
    else: v = "unclear"
    near = min(min(abs(pos - a), abs(pos - b)) for a, b in IS)
    return dict(verdict=v, why="junction", alt_reads=nj, ref_reads=ref, alt_frac=round(fr, 3), dist_h37rv_copy=near)


ACCFA = f"{B}/assets/accessory_catalogue.fasta"


def check_acc(bam_path, locus, threads, tmp):
    lf = pysam.FastaFile(ACCFA)
    ls = lf.fetch(locus)
    d = tempfile.mkdtemp(dir=tmp)
    comb = f"{d}/comb.fa"
    with open(comb, "w") as f:
        f.write(f">{rname}\n{G}\n>{locus}\n{ls}\n")
    subprocess.run([BWA, "index", comb], check=True, capture_output=True)
    # unique positions of the locus: 150-mers tiled every 5 bp
    with open(f"{d}/tiles.fq", "w") as f:
        for i in range(0, max(1, len(ls) - 150 + 1), 5):
            t = ls[i:i + 150]; f.write(f"@t{i}\n{t}\n+\n{'I' * len(t)}\n")
    uniq = np.zeros(len(ls), bool)
    p = subprocess.run([BWA, "mem", "-t", str(threads), comb, f"{d}/tiles.fq"], check=True, capture_output=True, text=True)
    for l in p.stdout.splitlines():
        if l.startswith("@"): continue
        f = l.split("\t")
        if f[2] == locus and int(f[4]) >= 20 and not int(f[1]) & 0x900:
            i = int(f[0][1:]); uniq[i:i + 150] = True
    bam = pysam.AlignmentFile(bam_path)
    n = 0
    with open(f"{d}/reads.fq", "w") as f:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail: continue
            if not r.is_unmapped:
                lc, rcl = clips(r)
                if not (lc >= 20 or rcl >= 20 or r.mapping_quality < 20 or (r.has_tag("NM") and r.get_tag("NM") >= 4)):
                    continue
            s = r.get_forward_sequence(); qq = r.get_forward_qualities()
            if not s: continue
            qs = "".join(chr(x + 33) for x in qq) if qq is not None else "I" * len(s)
            f.write(f"@{r.query_name}_{n}\n{s}\n+\n{qs}\n"); n += 1
    p = subprocess.Popen([BWA, "mem", "-t", str(threads), comb, f"{d}/reads.fq"], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    cov = np.zeros(len(ls), int)
    for l in p.stdout:
        if l.startswith("@"): continue
        f = l.split("\t")
        if f[2] != locus or int(f[4]) < 20 or int(f[1]) & 0x904: continue
        st = int(f[3]) - 1; span = 0
        import re
        for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            if op in "MDN=X": span += int(ln)
        cov[st:st + span] += 1
    p.wait()
    subprocess.run(["rm", "-rf", d])
    nu = int(uniq.sum())
    br = float(((cov >= 2) & uniq).sum() / nu) if nu else float("nan")
    if nu < 100: v = "unclear"; why = "few_unique_positions"
    elif br >= 0.7: v, why = "ALT", "breadth"
    elif br <= 0.15: v, why = "REF", "breadth"
    else: v, why = "unclear", "partial_breadth"
    return dict(verdict=v, why=why, L=len(ls), unique_pos=nu, breadth_unique=round(br, 3),
                reads_realigned=n, depth_unique=round(float(cov[uniq].mean()) if nu else 0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc-threads", type=int, default=4)
    ap.add_argument("--skip-acc", action="store_true")
    ap.add_argument("--sample", default=f"{D}/readcheck_sample.tsv")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    o = pickle.load(open(f"{D}/old.pkl", "rb")); n = pickle.load(open(f"{D}/new.pkl", "rb"))
    s = pd.read_csv(a.sample, sep="\t", keep_default_na=False)
    tmp = os.environ.get("TMPDIR", "/tmp")
    rows = []
    for _, c in s.iterrows():
        oldside = c.stratum in ("old_alt", "old_only")
        side = o if oldside else n
        rec = side["recs"][int(c.old_idx if oldside else c.new_idx)]
        bp = f"{BAMD}/{c['sample']}.h37rv.bam"
        res = dict(c)
        res.pop("seq", None)
        res.update(chk_pos=rec["pos"], chk_end=rec["end"], chk_len=rec["svlen"])
        if c["class"] == "ACC":
            if a.skip_acc: continue
            res.update(check_acc(bp, c.acclocus, a.acc_threads, tmp))
        else:
            bam = pysam.AlignmentFile(bp)
            if c["class"] == "DEL": res.update(check_del(bam, rec["pos"] + 1, rec["end"]))
            elif c["class"] == "INS": res.update(check_ins(bam, rec["pos"], rec["svlen"], rec["seq"]))
            else: res.update(check_is(bam, rec["pos"]))
        # who is right: the ALT side if reads support ALT, the other side if REF
        alt_side = "old" if oldside else "new"
        other = "new" if alt_side == "old" else "old"
        res["right"] = {"ALT": alt_side, "REF": other}.get(res["verdict"], "unclear")
        rows.append(res)
        print(c["class"], c.stratum, c["sample"], rec["pos"], res["verdict"], res.get("why"), file=sys.stderr, flush=True)
    out = a.out or f"{D}/readcheck_results{'_noacc' if a.skip_acc else ''}.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
