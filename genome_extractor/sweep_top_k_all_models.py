"""
Top-K sweep across all CovMutEx models + PRIEST for all 7 post-2022 variants.
Run: /Users/mac/miniconda3/envs/covmutex/bin/python /tmp/sweep_all_models.py
"""
import os, sys, time
sys.path.insert(0, "/Users/mac/Documents/ITU_Bioinformatics/CovMutEx-main/genome_extractor")
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "genome_extractor.settings")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import django; django.setup()

from genome.delta_omicron_retrospective import (
    PRECOMPUTED_VARIANT_CONTEXTS,
    extract_priest_spike_site_score_rows,
    load_precomputed_variant_spike_site_set,
    rank_spike_site_rows, normalize_ranked_site_scores,
    compute_overlap_metrics,
    build_ranked_spike_site_rows,
    build_precomputed_variant_mutation_tuples,
    SPIKE_START,
)
from genome.feature_extractor import construct_variant_genome, predict_mutations, load_legacy_keras_model
from genome.priest_annotations import PROTEIN_REGIONS
from genome.configs import configs

import numpy as np

# ── Config ────────────────────────────────────────────────────────────────────
BASE         = "/Users/mac/Documents/ITU_Bioinformatics/CovMutEx-main/genome_extractor"
MODEL_DIR    = f"{BASE}/covid19_models/models"
GENOME_PATH  = f"{BASE}/genome/genome.txt"
CACHE_PATH   = f"{BASE}/genome/node_features.h5"
CODON_PATH   = f"{BASE}/genome/codon_aa_mapping.json"
SPIKE_REGION = {"S": (PROTEIN_REGIONS["S"][0] - 1, PROTEIN_REGIONS["S"][1] - 1)}
TOTAL        = 1273
PROXIMITY    = 3

MODELS = {
    "PRIEST":                    None,
    "balanced_data_model":       f"{MODEL_DIR}/balanced_data_model_fixed.keras",
    "single_input_ensemble":     f"{MODEL_DIR}/single_input_ensemble_model_best.keras",
    "multi_input_ensemble":      f"{MODEL_DIR}/multi_input_ensemble_model.keras",
}

K_VALUES = (
    list(range(1, 26))
    + list(range(30, 51, 5))
    + list(range(60, 101, 10))
    + [125, 150, 200]
)

f1 = lambda p, r: 2*p*r/(p+r) if p+r else 0.0

def read_genome(path):
    with open(path) as f:
        next(f)
        return "".join(l.strip() for l in f)

def run_model_predictions(model_obj, genome_seq, mutation_tuples, node_id):
    preds = predict_mutations(
        cache_path=CACHE_PATH,
        genome_seq=genome_seq,
        mutations=mutation_tuples,
        codon_mapper=CODON_PATH,
        config_file=configs(),
        node_ids=[node_id],
        elapsed_day=0,
        depth=0,
        protein_regions=SPIKE_REGION,
        model=model_obj,
    )
    return preds

def sweep_k(ranked_rows, known_positions):
    rows = []
    for k in K_VALUES:
        if k > len(ranked_rows): break
        top_k = [r["aa_position"] for r in ranked_rows[:k]]
        m = compute_overlap_metrics(top_k, known_positions, proximity_window=PROXIMITY)
        p,  r  = m["precision_at_k"],        m["recall_against_omicron_sites"]
        pp, pr = m["proximity_precision_at_k"], m["proximity_recall"]
        rows.append(dict(k=k, kp=100*k/TOTAL,
                         eo=m["overlap_count"], p=p, r=r, f=f1(p,r),
                         po=m["proximity_overlap_count"], pp=pp, pr=pr, pf=f1(pp,pr)))
    return rows

def print_sweep(label, n_known, rows, model_name):
    if not rows: return
    be = max(rows, key=lambda r: (r["f"],  -r["k"]))
    bp = max(rows, key=lambda r: (r["pf"], -r["k"]))
    print(f"\n  [{model_name}]  known={n_known}")
    print(f"  {'K':>5} {'K%':>5} {'EO':>4} {'Prec':>6} {'Rec':>6} {'F1':>6}  "
          f"{'PO':>4} {'PPrc':>6} {'PRec':>6} {'PF1':>6}")
    for r in rows:
        tag=""
        if r["k"]==be["k"]: tag="  <- best F1"
        elif r["k"]==bp["k"]: tag="  <- best proxF1"
        print(f"  {r['k']:>5} {r['kp']:>4.1f}% {r['eo']:>4} "
              f"{r['p']:>5.1%} {r['r']:>5.1%} {r['f']:>5.1%}  "
              f"{r['po']:>4} {r['pp']:>5.1%} {r['pr']:>5.1%} {r['pf']:>5.1%}{tag}")
    print(f"  >> Best exact F1={be['f']:.3f} K={be['k']} ({be['kp']:.1f}%) "
          f"Prec={be['p']:.1%} Rec={be['r']:.1%}")
    print(f"  >> Best prox  F1={bp['pf']:.3f} K={bp['k']} ({bp['kp']:.1f}%) "
          f"Prec={bp['pp']:.1%} Rec={bp['pr']:.1%}")
    return be, bp

# ── Main ──────────────────────────────────────────────────────────────────────
print("="*76)
print("  CovMutEx Top-K Sweep -- ALL MODELS x ALL VARIANTS")
print("="*76)

genome_seq = read_genome(GENOME_PATH)
cfg = configs()
variants = sorted(PRECOMPUTED_VARIANT_CONTEXTS.items(), key=lambda x: x[1].get("sort_order",99))

# Summary table accumulation
summary = {}  # model -> list of (variant_label, best_k_exact, best_f1, best_k_prox, best_pf1)

for model_name, model_path in MODELS.items():
    print(f"\n{'#'*76}")
    print(f"  MODEL: {model_name}")
    print(f"{'#'*76}")

    model_obj = None
    if model_path:
        print(f"  Loading model from {model_path} ...")
        t0 = time.time()
        try:
            model_obj = load_legacy_keras_model(model_path)
            print(f"  Loaded in {time.time()-t0:.1f}s")
        except Exception as e:
            print(f"  [ERROR] Could not load model: {e}")
            continue

    summary[model_name] = []

    for node_id, ctx in variants:
        nick  = ctx.get("variant_nickname","")
        label = f"{ctx['variant_label']} ({nick})" if nick else ctx['variant_label']

        print(f"\n{'='*76}")
        print(f"  VARIANT: {label}")
        print(f"{'='*76}")

        try:
            known = list(load_precomputed_variant_spike_site_set(node_id)["positions"])
            n_known = len(known)

            if model_name == "PRIEST":
                site_rows, _ = extract_priest_spike_site_score_rows(node_context=ctx)
                ranked = rank_spike_site_rows(site_rows)
            else:
                mutation_build = build_precomputed_variant_mutation_tuples(
                    reference_genome_sequence=genome_seq,
                    node_context=ctx,
                )
                mut_tuples = mutation_build["mutation_tuples"]
                variant_genome = construct_variant_genome(genome_seq, mut_tuples)
                t0 = time.time()
                preds = run_model_predictions(model_obj, variant_genome, mut_tuples, node_id)
                print(f"  Predictions: shape={preds.shape} min={preds.min():.4f} "
                      f"max={preds.max():.4f} ({time.time()-t0:.1f}s)")

                ranked, _ = build_ranked_spike_site_rows(
                    selected_model=model_name,
                    node_context={"node_id": node_id, **ctx},
                    predictions=preds,
                    reference_genome_sequence=genome_seq,
                    spike_start_genome_position=SPIKE_START,
                    prediction_start_genome_position=SPIKE_START,
                )

            rows = sweep_k(ranked, known)
            result = print_sweep(label, n_known, rows, model_name)
            if result:
                be, bp = result
                summary[model_name].append((label, n_known, be["k"], be["f"], bp["k"], bp["pf"]))

        except Exception as e:
            import traceback; traceback.print_exc()
            print(f"  [ERROR] {label}: {e}")

# ── Global summary table ──────────────────────────────────────────────────────
print("\n\n" + "="*76)
print("  GRAND SUMMARY TABLE — Best-F1 K per model per variant")
print("="*76)
variants_labels = [
    (f"{ctx['variant_label']} ({ctx.get('variant_nickname','')})", ctx.get("sort_order"))
    for _, ctx in sorted(PRECOMPUTED_VARIANT_CONTEXTS.items(), key=lambda x: x[1].get("sort_order",99))
]
header = f"  {'Variant':<30}" + "".join(f"  {m[:14]:>14}" for m in MODELS)
print(header)
print("  " + "─"*74)

for label, _ in variants_labels:
    row = f"  {label:<30}"
    for model_name in MODELS:
        match = next((r for r in summary.get(model_name,[]) if r[0]==label), None)
        if match:
            row += f"  {'K='+str(match[2])+' F1='+f'{match[3]:.2f}':>14}"
        else:
            row += f"  {'—':>14}"
    print(row)

print("\n  PROXIMITY F1 TABLE")
print("  " + "─"*74)
for label, _ in variants_labels:
    row = f"  {label:<30}"
    for model_name in MODELS:
        match = next((r for r in summary.get(model_name,[]) if r[0]==label), None)
        if match:
            row += f"  {'K='+str(match[4])+' pF1='+f'{match[5]:.2f}':>14}"
        else:
            row += f"  {'—':>14}"
    print(row)

print("\n" + "="*76)
