"""
CovMutEx-X Benchmark Engine

Runs multiple models on the same variant/parameters, computes comparison metrics,
and stores results for the frontend dashboard.

Metrics computed:
- Runtime (seconds)
- Per-protein region average mutation probability
- Model agreement (correlation between models)
- Brier score (if ground-truth mutations available)
- AUROC / AUPRC (if ground-truth mutations available)
- Calibration / ECE (Expected Calibration Error)

Ground-truth: Derived from mutations.txt — positions where mutations actually occurred
are labeled as positive (1), all others as negative (0).
"""

import os
import json
import time
import numpy as np
import traceback
from typing import List, Dict, Any, Optional

# Metrics
from sklearn.metrics import (
    roc_auc_score, 
    average_precision_score, 
    brier_score_loss,
    roc_curve,
    precision_recall_curve,
)


# ============================================
# GROUND TRUTH BUILDER
# ============================================

def build_ground_truth(mutations: list, genome_length: int = 29904, region: tuple = None) -> np.ndarray:
    """
    Build binary ground-truth array from mutations list.
    
    A position is labeled 1 if a mutation occurred there, 0 otherwise.
    This represents "did a mutation actually happen at this position?"
    
    Args:
        mutations: List of (nt_position, original, new, aa_position, aa_change)
        genome_length: Total genome length (default: 29904 for SARS-CoV-2)
        region: Optional (start, end) tuple to slice the ground-truth
        
    Returns:
        Binary numpy array of shape (N,) where N = genome_length or region size
    """
    ground_truth = np.zeros(genome_length, dtype=np.float32)
    
    for mutation in mutations:
        pos = mutation[0]  # nt_position (1-based in mutations.txt)
        if 0 < pos <= genome_length:
            ground_truth[pos - 1] = 1.0  # Convert to 0-based
    
    if region:
        start, end = region
        ground_truth = ground_truth[start:end + 1]
    
    return ground_truth


# ============================================
# METRIC COMPUTATION
# ============================================

def compute_total_mutation_probability(predictions: np.ndarray) -> np.ndarray:
    """
    Convert per-nucleotide predictions to total mutation probability per position.
    
    If predictions shape is (N, 1) → single output, use directly.
    If predictions shape is (N, 4) → sum non-reference probabilities.
    If predictions shape is (4, N) → transpose first.
    """
    if predictions.ndim == 1:
        return predictions
    
    if predictions.shape[-1] == 1:
        return predictions.flatten()
    
    # Shape (4, N) → transpose to (N, 4)
    if predictions.shape[0] == 4 and predictions.shape[1] != 4:
        predictions = predictions.T
    
    # For (N, 4): total mutation prob = sum of all nucleotide probs
    # (This is a simplification — ideally we'd exclude the reference nucleotide)
    return np.sum(predictions, axis=1)


def compute_metrics(
    predictions: np.ndarray, 
    ground_truth: np.ndarray,
    model_name: str = "unknown"
) -> Dict[str, Any]:
    """
    Compute all benchmark metrics for a single model.
    
    Args:
        predictions: Model predictions (mutation probabilities per position)
        ground_truth: Binary array (1 = mutation occurred, 0 = no mutation)
        model_name: Name for logging
        
    Returns:
        Dict with all computed metrics
    """
    metrics = {}
    
    # Ensure same length
    min_len = min(len(predictions), len(ground_truth))
    preds = predictions[:min_len]
    gt = ground_truth[:min_len]
    
    # Clip predictions to valid range
    preds = np.clip(preds, 1e-7, 1 - 1e-7)
    
    # --- Brier Score ---
    try:
        metrics['brier_score'] = float(brier_score_loss(gt, preds))
    except Exception as e:
        metrics['brier_score'] = None
        print(f"[BENCH] Brier score failed for {model_name}: {e}")
    
    # --- AUROC ---
    try:
        if len(np.unique(gt)) > 1:  # Need both classes
            metrics['auroc'] = float(roc_auc_score(gt, preds))
            
            # ROC curve data (downsample for frontend)
            fpr, tpr, _ = roc_curve(gt, preds)
            step = max(1, len(fpr) // 200)  # Max 200 points for chart
            metrics['roc_curve'] = {
                'fpr': fpr[::step].tolist(),
                'tpr': tpr[::step].tolist(),
            }
        else:
            metrics['auroc'] = None
            metrics['roc_curve'] = None
            print(f"[BENCH] AUROC not computable for {model_name}: only one class in ground truth")
    except Exception as e:
        metrics['auroc'] = None
        metrics['roc_curve'] = None
        print(f"[BENCH] AUROC failed for {model_name}: {e}")
    
    # --- AUPRC ---
    try:
        if len(np.unique(gt)) > 1:
            metrics['auprc'] = float(average_precision_score(gt, preds))
            
            # PR curve data
            precision, recall, _ = precision_recall_curve(gt, preds)
            step = max(1, len(precision) // 200)
            metrics['pr_curve'] = {
                'precision': precision[::step].tolist(),
                'recall': recall[::step].tolist(),
            }
        else:
            metrics['auprc'] = None
            metrics['pr_curve'] = None
    except Exception as e:
        metrics['auprc'] = None
        metrics['pr_curve'] = None
        print(f"[BENCH] AUPRC failed for {model_name}: {e}")
    
    # --- ECE (Expected Calibration Error) ---
    try:
        metrics['ece'] = float(compute_ece(preds, gt, n_bins=10))
        metrics['calibration_curve'] = compute_calibration_curve(preds, gt, n_bins=10)
    except Exception as e:
        metrics['ece'] = None
        metrics['calibration_curve'] = None
        print(f"[BENCH] ECE failed for {model_name}: {e}")
    
    # --- Summary stats ---
    metrics['mean_prediction'] = float(np.mean(preds))
    metrics['std_prediction'] = float(np.std(preds))
    metrics['max_prediction'] = float(np.max(preds))
    metrics['min_prediction'] = float(np.min(preds))
    metrics['num_positions'] = int(min_len)
    metrics['num_mutations'] = int(np.sum(gt))
    metrics['mutation_rate'] = float(np.mean(gt))
    
    return metrics


def compute_ece(predictions: np.ndarray, ground_truth: np.ndarray, n_bins: int = 10) -> float:
    """
    Expected Calibration Error.
    
    Measures how well the predicted probabilities match actual mutation rates.
    Lower is better (0 = perfectly calibrated).
    """
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    total = len(predictions)
    
    for i in range(n_bins):
        mask = (predictions >= bin_boundaries[i]) & (predictions < bin_boundaries[i + 1])
        if np.sum(mask) == 0:
            continue
        
        bin_conf = np.mean(predictions[mask])
        bin_acc = np.mean(ground_truth[mask])
        bin_size = np.sum(mask)
        
        ece += (bin_size / total) * abs(bin_acc - bin_conf)
    
    return ece


def compute_calibration_curve(predictions: np.ndarray, ground_truth: np.ndarray, n_bins: int = 10) -> Dict:
    """
    Compute calibration curve data for plotting.
    
    Returns bin midpoints (predicted prob) vs actual fraction of positives.
    """
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    bin_midpoints = []
    bin_accuracies = []
    bin_sizes = []
    
    for i in range(n_bins):
        mask = (predictions >= bin_boundaries[i]) & (predictions < bin_boundaries[i + 1])
        count = int(np.sum(mask))
        if count == 0:
            continue
        
        bin_midpoints.append(float(np.mean(predictions[mask])))
        bin_accuracies.append(float(np.mean(ground_truth[mask])))
        bin_sizes.append(count)
    
    return {
        'predicted': bin_midpoints,
        'actual': bin_accuracies,
        'bin_sizes': bin_sizes,
    }


def compute_per_protein_metrics(
    predictions: np.ndarray,
    ground_truth: np.ndarray,
    protein_regions: Dict[str, list],
) -> Dict[str, Dict]:
    """
    Compute metrics per protein region.
    
    Args:
        predictions: Full-genome predictions
        ground_truth: Full-genome ground truth
        protein_regions: Dict mapping protein name to [start, end]
        
    Returns:
        Dict mapping protein name to metrics dict
    """
    per_protein = {}
    
    for protein, (start, end) in protein_regions.items():
        region_preds = predictions[start:end + 1] if end + 1 <= len(predictions) else predictions[start:]
        region_gt = ground_truth[start:end + 1] if end + 1 <= len(ground_truth) else ground_truth[start:]
        
        if len(region_preds) == 0:
            continue
        
        region_metrics = {
            'mean_prediction': float(np.mean(region_preds)),
            'max_prediction': float(np.max(region_preds)),
            'num_positions': int(len(region_preds)),
            'num_mutations': int(np.sum(region_gt)),
            'mutation_rate': float(np.mean(region_gt)),
        }
        
        # Only compute AUROC/AUPRC if both classes present
        if len(np.unique(region_gt)) > 1 and len(region_gt) > 10:
            try:
                region_metrics['auroc'] = float(roc_auc_score(region_gt, region_preds))
                region_metrics['auprc'] = float(average_precision_score(region_gt, region_preds))
                region_metrics['brier_score'] = float(brier_score_loss(region_gt, region_preds))
            except Exception:
                region_metrics['auroc'] = None
                region_metrics['auprc'] = None
                region_metrics['brier_score'] = None
        else:
            region_metrics['auroc'] = None
            region_metrics['auprc'] = None
            region_metrics['brier_score'] = None
        
        per_protein[protein] = region_metrics
    
    return per_protein


# ============================================
# MODEL AGREEMENT / CORRELATION
# ============================================

def compute_model_agreement(all_predictions: Dict[str, np.ndarray]) -> Dict[str, Dict[str, float]]:
    """
    Compute pairwise correlation between models' predictions.
    
    Returns a correlation matrix as nested dict.
    """
    model_names = list(all_predictions.keys())
    agreement = {}
    
    for i, name_a in enumerate(model_names):
        agreement[name_a] = {}
        for j, name_b in enumerate(model_names):
            preds_a = all_predictions[name_a]
            preds_b = all_predictions[name_b]
            
            min_len = min(len(preds_a), len(preds_b))
            if min_len > 0:
                # Check if either prediction is constant (std=0) which causes NaN
                if np.std(preds_a[:min_len]) == 0 or np.std(preds_b[:min_len]) == 0:
                    agreement[name_a][name_b] = None  # Can't compute correlation
                else:
                    corr = float(np.corrcoef(preds_a[:min_len], preds_b[:min_len])[0, 1])
                    agreement[name_a][name_b] = round(corr, 4) if not np.isnan(corr) else None
            else:
                agreement[name_a][name_b] = None
    
    return agreement


# ============================================
# MAIN BENCHMARK RUNNER
# ============================================

def run_benchmark(
    model_configs: List[Dict],
    node_id: str,
    elapsed_day: int,
    mutations: list,
    genome_sequence: str,
    protein_regions: Dict[str, list],
    predict_fn,
    selected_protein_region: str = None,
) -> Dict[str, Any]:
    """
    Run benchmark across multiple models.
    
    Args:
        model_configs: List of dicts, each with:
            - 'name': str (display name)
            - 'model_path': str (path to model file)  
            - 'source': str ('server' or 'uploaded')
            - 'custom_parameters': dict (optional)
        node_id: Variant node ID
        elapsed_day: Days since emergence
        mutations: Parsed mutations list
        genome_sequence: Reference genome string
        protein_regions: Dict of protein name → [start, end]
        predict_fn: Function that takes prediction params and returns predictions array
        selected_protein_region: Optional protein region to focus on
        
    Returns:
        Complete benchmark results dict
    """
    results = {
        'benchmark_id': f"bench_{int(time.time())}",
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'parameters': {
            'node_id': node_id,
            'elapsed_day': elapsed_day,
            'num_mutations': len(mutations),
            'selected_protein_region': selected_protein_region,
        },
        'models': {},
        'model_agreement': {},
        'summary': {},
    }
    
    # Build ground truth
    genome_length = len(genome_sequence) if genome_sequence else 29904
    
    if selected_protein_region and selected_protein_region in protein_regions:
        region = tuple(protein_regions[selected_protein_region])
        ground_truth = build_ground_truth(mutations, genome_length, region=region)
    else:
        ground_truth = build_ground_truth(mutations, genome_length)
    
    results['parameters']['ground_truth_positives'] = int(np.sum(ground_truth))
    results['parameters']['ground_truth_total'] = int(len(ground_truth))
    results['parameters']['mutation_positions'] = [int(i) for i in np.where(ground_truth == 1)[0]]
    
    all_predictions = {}
    
    # Run each model
    for config in model_configs:
        model_name = config['name']
        print(f"\n[BENCHMARK] Running model: {model_name}")
        
        model_result = {
            'name': model_name,
            'source': config.get('source', 'unknown'),
            'status': 'pending',
            'runtime_seconds': None,
            'metrics': {},
            'per_protein': {},
            'error': None,
        }
        
        try:
            start_time = time.time()
            
            # Run prediction using the provided function
            predictions_raw = predict_fn(
                model_path=config['model_path'],
                model_name=model_name,
                source=config.get('source', 'server'),
                node_id=node_id,
                elapsed_day=elapsed_day,
                mutations=mutations,
                genome_sequence=genome_sequence,
                protein_regions=protein_regions,
                selected_protein_region=selected_protein_region,
                custom_parameters=config.get('custom_parameters', {}),
                extractor_path=config.get('extractor_path'),
            )
            
            end_time = time.time()
            runtime = end_time - start_time
            
            # Convert to total mutation probability
            total_probs = compute_total_mutation_probability(predictions_raw)
            
            model_result['runtime_seconds'] = round(runtime, 3)
            model_result['status'] = 'success'
            
            # Compute metrics against ground truth
            model_result['metrics'] = compute_metrics(total_probs, ground_truth, model_name)
            model_result['metrics']['runtime_seconds'] = round(runtime, 3)
            
            # Per-protein breakdown (only for full genome)
            if not selected_protein_region:
                full_gt = build_ground_truth(mutations, genome_length)
                model_result['per_protein'] = compute_per_protein_metrics(
                    total_probs, full_gt, protein_regions
                )
            
            # Store for agreement computation
            all_predictions[model_name] = total_probs
            
            # Store downsampled predictions for overlay chart (every 100th position)
            step = max(1, len(total_probs) // 300)  # ~300 data points
            sampled_indices = list(range(0, len(total_probs), step))
            model_result['prediction_curve'] = {
                'positions': sampled_indices,
                'values': [round(float(total_probs[i]), 6) for i in sampled_indices],
                'step': step,
                'total_positions': len(total_probs),
            }
            
            print(f"[BENCHMARK] {model_name} completed in {runtime:.2f}s")
            
        except Exception as e:
            model_result['status'] = 'error'
            model_result['error'] = str(e)
            print(f"[BENCHMARK] {model_name} FAILED: {e}")
            print(traceback.format_exc())
        
        results['models'][model_name] = model_result
    
    # Compute model agreement
    if len(all_predictions) >= 2:
        results['model_agreement'] = compute_model_agreement(all_predictions)
    
    # Summary: best model per metric
    successful_models = {k: v for k, v in results['models'].items() if v['status'] == 'success'}
    if successful_models:
        summary = {}
        for metric in ['auroc', 'auprc', 'brier_score', 'ece', 'runtime_seconds']:
            values = {}
            for name, result in successful_models.items():
                val = result['metrics'].get(metric)
                if val is not None:
                    values[name] = val
            
            if values:
                if metric in ['brier_score', 'ece', 'runtime_seconds']:
                    # Lower is better
                    best = min(values, key=values.get)
                else:
                    # Higher is better
                    best = max(values, key=values.get)
                summary[metric] = {'best_model': best, 'value': values[best], 'all_values': values}
        
        results['summary'] = summary
    
    return results