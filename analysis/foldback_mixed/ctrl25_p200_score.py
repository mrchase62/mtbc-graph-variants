#!/usr/bin/env python3
"""Score the ctrl25_p200 merged VCF (P200 build, all fixes) against each
control's own assembly (user, 2026-10-14). Analysis only.

Truth: the assembly aligned to H37Rv (grid_score.truth: minimap2 asm5,
paftools call -l 1000 -L 10000). Scored only in positions the assembly
alignment calls once and outside the build's repeat mask.

Per sample and class (H37Rv-frame records only):
  small SNP / indel  ALT cell matching a truth variant = TP; ALT cell with no
                     truth variant = FP; truth variant without an ALT cell =
                     FN, split by the cell's state (REF, NOCALL, ABSENT, or
                     no record)
  SV deletion        <DEL> record ALT vs truth deletions >= 50 bp,
                     reciprocal overlap >= 50%
  SV insertion       INS record ALT vs truth insertions >= 50 bp, start within
                     50 bp and length within 20%
Also the EV reasons of ABSENT / NOCALL cells, and inherited cells' states.

    ctrl25_p200_score.py --vcf <merged.vcf.gz> --out out/ctrl25_p200_score.tsv
"""
import argparse, bisect, collections, gzip, os, statistics as S, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grid_score as gs  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    work = f"{gs.FM}/controls_grid/truth"
    fb = {r[0]: float(r[3]) for r in __import__("csv").reader(open(gs.RATES), delimiter="\t") if r[0] != "sample"}
    fh = gzip.open(a.vcf, "rt")
    for l in fh:
        if l.startswith("#CHROM"):
            samples = l.rstrip("\n").split("\t")[9:]; break
    recs = []
    for l in fh:
        f = l.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        if info.get("FRAME", "h37rv") != "h37rv" and info.get("CLASS") != "sv": continue
        fmt = f[8].split(":")
        cells = [dict(zip(fmt, c.split(":"))) for c in f[9:]]
        recs.append((int(f[1]), f[3], f[4], info, cells))
    rows, ev = [], collections.Counter()
    for j, s in enumerate(samples):
        cal, tv = gs.truth(s, "GCF_000195955", work)
        cs = sorted(cal); st = [x[0] for x in cs]
        def ok(p):
            i = bisect.bisect_right(st, p - 1) - 1
            return i >= 0 and cs[i][0] <= p - 1 < cs[i][1] and not gs.mask[p - 1]
        tsnv = {(p, al) for p, r, al in tv if len(r) == len(al) == 1 and ok(p)}
        tind = [(p, len(al) - len(r)) for p, r, al in tv if not (len(r) == len(al) == 1) and abs(len(al) - len(r)) < 50 and ok(p)]
        tdel = [(p, p + len(r) - len(al)) for p, r, al in tv if len(r) - len(al) >= 50 and ok(p)]
        tins = [(p, len(al) - len(r)) for p, r, al in tv if len(al) - len(r) >= 50 and ok(p)]
        n = collections.Counter()
        snv_state, ind_hit = {}, set()
        for pos, ref, alt, info, cells in recs:
            c = cells[j]; stt = c.get("ST", "."); gt = c.get("GT", ".")
            cls = info.get("CLASS")
            if stt in ("ABSENT", "NOCALL"):
                ev[(stt, c.get("EV", ".") or ".")] += 1
            if cls == "small":
                a0 = alt.split(",")[0]
                if not ok(pos): continue
                if len(ref) == 1 and len(a0) == 1:
                    if (pos, a0) in tsnv: snv_state[(pos, a0)] = stt
                    elif gt == "1": n["snv_fp"] += 1
                else:
                    if gt != "1": continue
                    d = len(a0) - len(ref)
                    hit = next((k for k, (p, dl) in enumerate(tind) if dl == d and abs(p - pos) <= 10), None)
                    if hit is None: n["indel_fp"] += 1
                    else: ind_hit.add(hit)
            elif cls == "sv" and info.get("SVTYPE") == "DEL" and "END" in info:
                if gt != "1" or not ok(pos + 1): continue
                e = int(info["END"]); b = pos + 1
                m = any(min(e, te) - max(b, tb) >= 0.5 * max(e - b, te - tb) for tb, te in tdel)
                n["del_tp" if m else "del_fp"] += 1
                if info.get("COMPONENT") == "interval" and c.get("EV", "").startswith("inherited"):
                    n["del_inherited_" + ("tp" if m else "fp")] += 1
            elif cls == "sv" and info.get("SVTYPE") == "INS":
                if gt != "1" or not ok(pos): continue
                L = abs(int(info.get("SVLEN", "0") or 0))
                m = any(abs(p - pos) <= 50 and min(L, tl) >= 0.8 * max(L, tl) for p, tl in tins)
                n["ins_tp" if m else "ins_fp"] += 1
        for k, v in snv_state.items():
            n["snv_tp" if v == "ALT" else f"snv_fn_{v}"] += 1
        n["snv_fn_norecord"] = len(tsnv) - len(snv_state)
        n["indel_tp"] = len(ind_hit); n["indel_fn"] = len(tind) - len(ind_hit)
        n["truth_snv"], n["truth_indel"], n["truth_del"], n["truth_ins"] = len(tsnv), len(tind), len(tdel), len(tins)
        rows.append({"sample": s, "foldback": fb.get(s), **n})
    keys = ["sample", "foldback"] + sorted({k for r in rows for k in r} - {"sample", "foldback"})
    with open(a.out, "w") as fo:
        fo.write("\t".join(keys) + "\n")
        for r in rows: fo.write("\t".join(str(r.get(k, 0)) for k in keys) + "\n")
    tot = collections.Counter()
    for r in rows:
        for k, v in r.items():
            if k not in ("sample", "foldback"): tot[k] += v
    def pr(tp, fp, fn):
        return f"TP {tp:,}  FP {fp:,}  FN {fn:,}  precision {tp / max(1, tp + fp):.4f}  recall {tp / max(1, tp + fn):.4f}"
    print("All 25 controls, outside the repeat mask, callable in the assembly alignment")
    fn_snv = sum(v for k, v in tot.items() if k.startswith("snv_fn_"))
    print("SNP   " + pr(tot["snv_tp"], tot["snv_fp"], fn_snv))
    print("      missed SNPs by cell state: " + ", ".join(f"{k[7:]} {v}" for k, v in tot.items() if k.startswith("snv_fn_") and v))
    print("indel " + pr(tot["indel_tp"], tot["indel_fp"], tot["indel_fn"]))
    print(f"SV DEL ALT cells: {tot['del_tp']:,} match a truth deletion, {tot['del_fp']:,} do not "
          f"(inherited-evidence cells {tot['del_inherited_tp']} / {tot['del_inherited_fp']}); truth deletions {tot['truth_del']}")
    print(f"SV INS ALT cells: {tot['ins_tp']:,} match, {tot['ins_fp']:,} do not; truth insertions {tot['truth_ins']}")
    print("ABSENT/NOCALL evidence (all regions):")
    for (stt, e), v in ev.most_common(15): print(f"   {stt:6s} {e:34s} {v:,}")


if __name__ == "__main__":
    main()
