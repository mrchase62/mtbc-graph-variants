#!/usr/bin/env python3
"""H37Rv frame on the 25 controls, against assembly truth (user, 2026-10-14:
"Could we explore what causes the 10-66 mixed [SNVs]? What about snps/indels
any loses?").

Truth: each control's own assembly aligned to H37Rv (minimap2 asm5, paftools
call -l 1000 -L 10000), scored in callable positions outside the repeat mask.

1. LOSSES AND FALSE CALLS, per caller and reads (grid_score.score rules):
   HC       production haploid HaplotypeCaller (P1's H37Rv VCF; synthetic:
            refbias/work/fbmix_syn/syn_<s>.vcf.gz)
   default  Mutect2 defaults + min base quality 20 (controls/)
   user     the user's settings (controls_user/h37rv_*)
2. WHAT THE FALSE MIXED CALLS ARE (user settings), each classed by the first
   rule it meets:
   paralog       the 61-mer around the site carrying the ALT occurs in the
                 isolate's own assembly (either strand): a second copy that
                 differs from the first, collapsed onto one H37Rv locus
   near_sv       within 150 bp of a true indel or SV of >= 10 bp
   near_cluster  within 30 bp of >= 2 true SNVs (dense divergence)
   true_nearby   within 5 bp of a true variant of any kind (alignment shift)
   other
   Calls in uncallable truth regions (not covered once by the assembly
   alignment) are counted separately.

    h37rv_explore.py --out out/h37rv_explore
"""
import argparse, bisect, collections, csv, gzip, os, statistics as S, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grid_score as gs  # noqa: E402

RR = gs.RR
FM = gs.FM
H37 = f"{gs.BUILD}/refs/GCF_000195955.fasta"
RCT = str.maketrans("ACGT", "TGCA")


def read_fa(path):
    op = gzip.open if path.endswith(".gz") else open
    return "".join(l.strip().upper() for l in op(path, "rt") if not l.startswith(">"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    work = f"{FM}/controls_grid/truth"
    h37 = read_fa(H37)
    fb = {r[0]: float(r[3]) for r in csv.reader(open(gs.RATES), delimiter="\t") if r[0] != "sample"}
    truth_file = {r["sample"]: r["file"] for r in csv.DictReader(open(gs.TRUTH), delimiter="\t")}
    rows, cause = [], []
    for s in (l.strip() for l in open(gs.CTRL) if l.strip()):
        grp = "chimeric" if fb.get(s, 0) > 4.8 else "clean"
        cal, tv = gs.truth(s, "GCF_000195955", work)
        cs = sorted(cal); starts = [x[0] for x in cs]

        def in_cal(p):
            i = bisect.bisect_right(starts, p - 1) - 1
            return i >= 0 and cs[i][0] <= p - 1 < cs[i][1]
        ok = lambda p: in_cal(p) and not gs.mask[p - 1]
        files = {("HC", "real"): (f"{RR}/work/marinB63_p1/{s}.h37rv.vcf.gz", True),
                 ("HC", "syn"): (f"{RR}/work/fbmix_syn/syn_{s}.vcf.gz", True),
                 ("default", "real"): (f"{FM}/controls/real_{s}.vcf.gz", False),
                 ("default", "syn"): (f"{FM}/controls/syn_{s}.vcf.gz", False),
                 ("user", "real"): (f"{FM}/controls_user/h37rv_real_{s}.vcf.gz", False),
                 ("user", "syn"): (f"{FM}/controls_user/h37rv_syn_{s}.vcf.gz", False)}
        calls = {k: gs.vcf_calls(p, hc) for k, (p, hc) in files.items()}
        for (caller, arm), cl in calls.items():
            n = gs.score(cl, tv, ok)
            rows.append({"sample": s, "group": grp, "caller": caller, "reads": arm, **n})

        # causes of false mixed calls, user settings
        asm = read_fa(f"{gs.ASM}/{truth_file[s]}")
        asm_rc = asm.translate(RCT)[::-1]
        tsnv = sorted(p for p, r, al in tv if len(r) == len(al) == 1)
        tbig = sorted((p, p + max(len(r), 1)) for p, r, al in tv if abs(len(al) - len(r)) >= 10)
        tall = sorted(p for p, _, _ in tv)
        tset = {(p, al) for p, r, al in tv if len(r) == len(al) == 1}
        for arm in ("real", "syn"):
            for c, p, r, al, af in calls[("user", arm)]:
                if not (0.1 <= af < 0.9) or len(r) != 1 or len(al) != 1: continue
                if gs.mask[p - 1]: continue
                if (p, al) in tset: continue          # a true SNV at low AF: a demotion, not false
                if not in_cal(p):
                    k = "uncallable"
                else:
                    w = h37[p - 31:p - 1] + al + h37[p:p + 30]
                    if len(w) == 61 and (w in asm or w in asm_rc):
                        k = "paralog"
                    elif any(a_ - 150 <= p <= b_ + 150 for a_, b_ in tbig):
                        k = "near_sv"
                    elif bisect.bisect_right(tsnv, p + 30) - bisect.bisect_left(tsnv, p - 30) >= 2:
                        k = "near_cluster"
                    elif bisect.bisect_right(tall, p + 5) - bisect.bisect_left(tall, p - 5) >= 1:
                        k = "true_nearby"
                    else:
                        k = "other"
                cause.append({"sample": s, "group": grp, "reads": arm, "pos": p, "ref": r, "alt": al, "af": af, "cause": k})
    keys = ["sample", "group", "caller", "reads"] + sorted({k for r in rows for k in r} - {"sample", "group", "caller", "reads"})
    with open(f"{a.out}/scores.tsv", "w") as fo:
        fo.write("\t".join(keys) + "\n")
        for r in rows: fo.write("\t".join(str(r.get(k, 0)) for k in keys) + "\n")
    with open(f"{a.out}/false_mixed_causes.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(cause[0]), delimiter="\t"); w.writeheader(); w.writerows(cause)

    print("1. H37Rv frame vs assembly truth, median per control: real (synthetic)")
    print("group     caller   type   truth   fixed_tp    fixed_fp  mixed_false  demoted   missed")
    for g in ("clean", "chimeric"):
        for caller in ("HC", "default", "user"):
            for t in ("snv", "indel"):
                def med(arm, k):
                    v = [r.get(f"{t}_{k}", 0) for r in rows if r["group"] == g and r["caller"] == caller and r["reads"] == arm]
                    return f"{S.median(v):g}"
                print(f"{g:9s} {caller:8s} {t:5s} {med('real', 'truth'):>6s}  {med('real', 'fixed_tp'):>5s} ({med('syn', 'fixed_tp')})"
                      f"  {med('real', 'fixed_fp'):>3s} ({med('syn', 'fixed_fp')})  {med('real', 'mixed_false'):>4s} ({med('syn', 'mixed_false')})"
                      f"  {med('real', 'demoted'):>3s} ({med('syn', 'demoted')})  {med('real', 'missed'):>3s} ({med('syn', 'missed')})")
    print("\n2. False mixed SNVs under the user's settings on H37Rv, by cause (all controls pooled)")
    for arm in ("syn", "real"):
        for g in ("clean", "chimeric"):
            x = collections.Counter(c["cause"] for c in cause if c["reads"] == arm and c["group"] == g)
            tot = sum(x.values())
            print(f"  {arm:4s} {g:9s} {tot:5d}: " + ", ".join(f"{k} {v} ({v / tot:.0%})" for k, v in x.most_common()))


if __name__ == "__main__":
    main()
