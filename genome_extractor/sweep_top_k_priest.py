"""
Top-K sweep for Known Hotspot Case Study — PRIEST scores only (fast, no GPU).

Sweeps K=1..300 across all 7 post-2022 variant contexts and reports
Precision@K, Recall@K, F1@K, and proximity variants per variant.
Use this to find the optimal K before running the heavier all-models sweep.

Run from genome_extractor/:
    /Users/mac/miniconda3/envs/covmutex/bin/python sweep_top_k_priest.py
"""

import os, sys
sys.path.insert(0, "/Users/mac/Documents/ITU_Bioinformatics/CovMutEx-main/genome_extractor")
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "genome_extractor.settings")
import django; django.setup()

from genome.delta_omicron_retrospective import (
    PRECOMPUTED_VARIANT_CONTEXTS, extract_priest_spike_site_score_rows,
    load_precomputed_variant_spike_site_set, rank_spike_site_rows, compute_overlap_metrics,
)

K_VALUES = list(range(1,51)) + list(range(55,101,5)) + list(range(110,201,10)) + [250,300]
TOTAL = 1273
f1 = lambda p,r: 2*p*r/(p+r) if p+r else 0.0

def sweep(node_id, ctx):
    known = list(load_precomputed_variant_spike_site_set(node_id)["positions"])
    site_rows,_ = extract_priest_spike_site_score_rows(node_context=ctx)
    ranked = rank_spike_site_rows(site_rows)
    rows=[]
    for k in K_VALUES:
        if k>len(ranked): break
        top_k=[r["aa_position"] for r in ranked[:k]]
        m=compute_overlap_metrics(top_k, known, proximity_window=3)
        p,r=m["precision_at_k"],m["recall_against_omicron_sites"]
        pp,pr=m["proximity_precision_at_k"],m["proximity_recall"]
        rows.append(dict(k=k,kp=100*k/TOTAL,eo=m["overlap_count"],p=p,r=r,f=f1(p,r),
                         po=m["proximity_overlap_count"],pp=pp,pr=pr,pf=f1(pp,pr)))
    return len(known), len(ranked), rows

variants=sorted(PRECOMPUTED_VARIANT_CONTEXTS.items(),key=lambda x:x[1].get("sort_order",99))
all_be,all_bp=[],[]

print("="*76)
print("  Top-K Sweep — PRIEST scores — all 7 post-2022 variants")
print("="*76)

for nid,ctx in variants:
    nick=ctx.get("variant_nickname","")
    lbl=f"{ctx['variant_label']} ({nick})" if nick else ctx["variant_label"]
    try:
        nk,nr,rows=sweep(nid,ctx)
        be=max(rows,key=lambda r:(r["f"], -r["k"]))
        bp=max(rows,key=lambda r:(r["pf"],-r["k"]))
        all_be.append(be["k"]); all_bp.append(bp["k"])
        print(f"\n{'─'*76}")
        print(f"  {lbl}   known={nk}  ranked={nr}")
        print(f"{'─'*76}")
        print(f"{'K':>5} {'K%':>5} {'EO':>4} {'Prec':>6} {'Rec':>6} {'F1':>6}  {'PO':>4} {'PPrc':>6} {'PRec':>6} {'PF1':>6}")
        for r in rows:
            tag=""
            if r["k"]==be["k"]: tag="  <- best F1"
            elif r["k"]==bp["k"]: tag="  <- best proxF1"
            print(f"{r['k']:>5} {r['kp']:>4.1f}% {r['eo']:>4} {r['p']:>5.1%} {r['r']:>5.1%} {r['f']:>5.1%}  {r['po']:>4} {r['pp']:>5.1%} {r['pr']:>5.1%} {r['pf']:>5.1%}{tag}")
        print(f"\n  Best exact F1={be['f']:.3f} K={be['k']} ({be['kp']:.1f}%) Prec={be['p']:.1%} Rec={be['r']:.1%}")
        print(f"  Best prox  F1={bp['pf']:.3f} K={bp['k']} ({bp['kp']:.1f}%) Prec={bp['pp']:.1%} Rec={bp['pr']:.1%}")
    except Exception as e:
        import traceback; traceback.print_exc()

print("\n"+"="*76)
print("  GLOBAL SUMMARY")
print("="*76)
if all_be:
    se=sorted(all_be); sp=sorted(all_bp)
    print(f"  Best exact-F1 Ks : {all_be}  median={se[len(se)//2]}")
    print(f"  Best prox-F1  Ks : {all_bp}   median={sp[len(sp)//2]}")
print("""
  WHY K=700 IS BAD FOR A PAPER:
    Precision@K falls as K grows -- you dilute the signal.
    Recall trivially rises -- K=700 = calling 55% of all Spike "hotspots".
    Reviewers will flag this as an unfair baseline.

  RECOMMENDED K:
    K ~ # known sites        equal-budget, balanced F1 (see above)
    K <= 10% of Spike (~127) defensible hard cutoff for a methods paper
    Peak F1 K (above)        mathematically optimal
    Report K=25,50,100 side-by-side  standard IR table convention
""")
