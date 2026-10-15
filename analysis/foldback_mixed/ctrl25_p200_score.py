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
  insertion          an ALT call of any insertion-like class (SV INS record,
                     IS6110 site, accessory presence; length SVLEN or ACCLEN)
                     vs truth insertions >= 50 bp, start within 50 bp and
                     length within 20%; recall over truth insertions
Also the EV reasons of ABSENT / NOCALL cells, and inherited cells' states.

    ctrl25_p200_score.py --vcf <merged.vcf.gz> --out out/ctrl25_p200_score.tsv
"""
import argparse, bisect, collections, gzip, os, statistics as S, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grid_score as gs  # noqa: E402


CAT = f"{gs.BUILD.replace('7713a8d71d8e-fix1', '3893a3762ed1')}/assets"


def accessory_truth(s):
    """{locus_id: True/False}: the locus representative sequence covered >= 0.8
    (present) or < 0.2 (absent) in the sample's own assembly; partial and
    copy_number loci are left out."""
    import csv, subprocess, gzip as gz
    fa = next(r["file"] for r in csv.DictReader(open(gs.TRUTH), delimiter="\t") if r["sample"] == s)
    lens, name = {}, None
    for l in open(f"{CAT}/accessory_catalogue.fasta"):
        if l.startswith(">"): name = l[1:].split()[0]; lens[name] = 0
        else: lens[name] += len(l.strip())
    nov = {r["locus_id"]: r["novelty"] for r in csv.DictReader(open(f"{CAT}/accessory_catalogue.tsv"), delimiter="\t")}
    r = subprocess.run([gs.MM2, "-c", "-x", "asm10", "--secondary=no", f"{gs.ASM}/{fa}",
                        f"{CAT}/accessory_catalogue.fasta"], capture_output=True, text=True, check=True)
    iv = collections.defaultdict(list)
    for l in r.stdout.splitlines():
        c = l.split("\t"); iv[c[0]].append((int(c[2]), int(c[3])))
    out = {}
    for k, L in lens.items():
        if nov.get(k) == "copy_number" or not L: continue
        n = last = 0
        for a_, b_ in sorted(iv[k]):
            a_ = max(a_, last)
            if b_ > a_: n += b_ - a_
            last = max(last, b_)
        f = n / L
        if f >= 0.8: out[k] = True
        elif f < 0.2: out[k] = False
    return out


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
        # insertion points often sit at the edge of a callable interval (the
        # assembly alignment breaks there), so truth insertions are taken
        # outside the repeat mask whether or not the point itself is callable
        tins = [(p, len(al) - len(r)) for p, r, al in tv if len(al) - len(r) >= 50 and not gs.mask[p - 1]]
        # a call just outside the mask can sit next to a truth insertion whose
        # point is inside it (repeat edges); calls are matched against every
        # truth insertion, the recall denominator stays unmasked ones
        tins_all = [(p, len(al) - len(r)) for p, r, al in tv if len(al) - len(r) >= 50]
        n = collections.Counter()
        snv_state, ind_hit = {}, set()
        calls = []
        acc_truth = accessory_truth(s)
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
            elif cls == "accessory_presence" and info.get("ACCLOCUS") in acc_truth:
                tr = acc_truth[info["ACCLOCUS"]]
                if gt in ("0", "1"):
                    n[f"acc_{'t' if (gt == '1') == tr else 'f'}{'p' if gt == '1' else 'n'}"] += 1
                else:
                    n["acc_nocall"] += 1
            elif cls == "sv" and info.get("SVTYPE") == "DEL" and "END" in info:
                if gt != "1" or not ok(pos + 1): continue
                e = int(info["END"]); b = pos + 1
                m = any(min(e, te) - max(b, tb) >= 0.5 * max(e - b, te - tb) for tb, te in tdel)
                n["del_tp" if m else "del_fp"] += 1
                if info.get("COMPONENT") == "interval" and c.get("EV", "").startswith("inherited"):
                    n["del_inherited_" + ("tp" if m else "fp")] += 1
            elif gt == "1" and ok(pos) and (
                    (cls == "sv" and info.get("SVTYPE") == "INS") or cls == "is6110"):
                L = abs(int(info.get("SVLEN", "0") or 0))
                if L >= 50:
                    calls.append((pos, L, cls))
        # insertions: an ALT call of any class (H37Rv-frame or caller SV
        # insertion, IS6110 site, accessory presence) matches a truth insertion
        # with start within 50 bp and length within 20%; precision per call,
        # recall per truth insertion (matched by a call of any class)
        hit_t = set()
        site_t = set()
        for pos, L, cls in calls:
            m = [k for k, (p, tl) in enumerate(tins_all) if abs(p - pos) <= 50 and min(L, tl) >= 0.8 * max(L, tl)]
            site = [k for k, (p, tl) in enumerate(tins_all) if abs(p - pos) <= 50]
            n[f"ins_{cls}_{'tp' if m else 'site' if site else 'fp'}"] += 1
            idx = {k for k, t in enumerate(tins_all) if t in set(tins)}
            hit_t.update(tins.index(tins_all[k]) for k in m if k in idx)
            site_t.update(tins.index(tins_all[k]) for k in site if k in idx)
        n["ins_truth_site_found"] = len(site_t)
        n["ins_truth_found"] = len(hit_t)
        n["ins_truth_is6110"] = sum(1300 <= tl <= 1420 for _, tl in tins)
        n["ins_truth_is6110_found"] = sum(1300 <= tins[k][1] <= 1420 for k in hit_t)
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
    print(f"INSERTIONS (>= 50 bp): truth {tot['truth_ins']}, found by some ALT call {tot['ins_truth_found']} "
          f"(recall {tot['ins_truth_found'] / max(1, tot['truth_ins']):.3f}); IS6110-sized truth "
          f"{tot['ins_truth_is6110']}, found {tot['ins_truth_is6110_found']}")
    print(f"   site found (any length) {tot['ins_truth_site_found']} ({tot['ins_truth_site_found'] / max(1, tot['truth_ins']):.3f})")
    for c in ("sv", "is6110"):
        tp, si, fp = tot[f"ins_{c}_tp"], tot[f"ins_{c}_site"], tot[f"ins_{c}_fp"]
        print(f"   {c:8s} ALT calls {tp + si + fp:,}: match {tp:,}, right site but length off by > 20% {si:,}, "
              f"no truth insertion within 50 bp {fp:,}")
    tp, fp, tn, fn = tot["acc_tp"], tot["acc_fp"], tot["acc_tn"], tot["acc_fn"]
    print(f"ACCESSORY presence (locus sequence anywhere in the assembly, minimap2 asm10 coverage >= 0.8 present, "
          f"< 0.2 absent; copy-number loci skipped): TP {tp:,} FP {fp:,} TN {tn:,} FN {fn:,} NOCALL {tot['acc_nocall']:,}; "
          f"precision {tp / max(1, tp + fp):.3f} recall {tp / max(1, tp + fn):.3f}")
    print("ABSENT/NOCALL evidence (all regions):")
    for (stt, e), v in ev.most_common(15): print(f"   {stt:6s} {e:34s} {v:,}")


if __name__ == "__main__":
    main()
