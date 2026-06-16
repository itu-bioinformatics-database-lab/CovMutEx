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
import re
import csv
import io
import json
import time
import urllib.parse
import urllib.request
import urllib.error
import numpy as np
import traceback
from typing import List, Dict, Any, Optional, Tuple

# Metrics
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    brier_score_loss,
    roc_curve,
    precision_recall_curve,
    precision_score,
    recall_score,
    f1_score,
    matthews_corrcoef,
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
# CSV-BASED GROUND TRUTH (cov-spectrum.org format)
# ============================================

# Mutation pattern: e.g. "A1C", "T3-", "C23456T"
# Captures: original base (letter), position (digits), new base (letter or '-')
_COV_SPECTRUM_MUTATION_RE = re.compile(r'^([ACGTU])(\d+)([ACGTU-])$', re.IGNORECASE)


def parse_cov_spectrum_mutation(mutation_str: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse a cov-spectrum.org nucleotide mutation string.

    Format: "<orig><pos><new>"
    - "A1C"  -> ('A', 1, 'C')
    - "T3-"  -> ('T', 3, '-')  (deletion, treated as a regular mutation)
    - "C23456T" -> ('C', 23456, 'T')

    Returns (original, position, new) or None if the string is not parseable.
    """
    if not isinstance(mutation_str, str):
        return None
    match = _COV_SPECTRUM_MUTATION_RE.match(mutation_str.strip())
    if not match:
        return None
    orig, pos, new = match.group(1).upper(), int(match.group(2)), match.group(3).upper()
    return orig, pos, new


def build_ground_truth_from_csv(
    csv_content: str,
    genome_length: int = 29904,
    region: tuple = None,
) -> Dict[str, Any]:
    """
    Build probability-based ground-truth array from a cov-spectrum.org CSV.

    Expected CSV columns: mutation, proportion, count, jaccard
    The `proportion` value (in [0, 1]) is used directly as the per-position
    ground-truth probability — enabling probability-vs-probability comparison
    against model predictions rather than binary classification.

    Deletion mutations (e.g., "T3-") are treated as normal mutations.
    Multiple mutations at the same position (different alt bases) are summed
    and clipped to 1.0.

    Args:
        csv_content: Raw CSV text.
        genome_length: Total genome length (default: 29904 for SARS-CoV-2).
        region: Optional (start, end) tuple to slice the ground-truth.

    Returns:
        {
            'ground_truth': np.ndarray of shape (N,) with per-position probabilities,
            'num_mutations_parsed': int,
            'num_mutations_skipped': int,
            'skipped_examples': list of up to 5 mutation strings that failed to parse,
            'mutation_positions': list of 1-based positions with nonzero probability,
        }
    """
    ground_truth = np.zeros(genome_length, dtype=np.float32)
    num_parsed = 0
    num_skipped = 0
    skipped_examples = []

    reader = csv.DictReader(io.StringIO(csv_content))
    # Normalize field names to lowercase for case-insensitive matching
    if reader.fieldnames:
        reader.fieldnames = [f.strip().lower() for f in reader.fieldnames]

    for row in reader:
        mutation_str = (row.get('mutation') or '').strip()
        prop_str = (row.get('proportion') or '').strip()

        if not mutation_str or not prop_str:
            num_skipped += 1
            if len(skipped_examples) < 5:
                skipped_examples.append(mutation_str or '<empty>')
            continue

        parsed = parse_cov_spectrum_mutation(mutation_str)
        if parsed is None:
            num_skipped += 1
            if len(skipped_examples) < 5:
                skipped_examples.append(mutation_str)
            continue

        _, pos, _ = parsed

        try:
            proportion = float(prop_str)
        except ValueError:
            num_skipped += 1
            if len(skipped_examples) < 5:
                skipped_examples.append(mutation_str)
            continue

        if not (0 < pos <= genome_length):
            num_skipped += 1
            continue

        # Sum proportions if multiple alt bases at same position (clip at 1.0)
        idx = pos - 1  # Convert 1-based to 0-based
        ground_truth[idx] = min(1.0, ground_truth[idx] + proportion)
        num_parsed += 1

    mutation_positions_full = [int(i) for i in np.where(ground_truth > 0)[0]]

    if region:
        start, end = region
        ground_truth = ground_truth[start:end + 1]

    return {
        'ground_truth': ground_truth,
        'num_mutations_parsed': num_parsed,
        'num_mutations_skipped': num_skipped,
        'skipped_examples': skipped_examples,
        'mutation_positions': mutation_positions_full,
    }


# ============================================
# COV-SPECTRUM.ORG API (LAPIS)
# ============================================

# LAPIS is the query API powering cov-spectrum.org.
# Public, no auth required. Nucleotide mutations endpoint returns
# { "data": [ { "mutation": "A1C", "count": 86, "coverage": ..., "proportion": 0.0077 }, ... ] }
COV_SPECTRUM_LAPIS_BASE = "https://lapis.cov-spectrum.org/open/v2"
COV_SPECTRUM_NUC_MUTATIONS_ENDPOINT = COV_SPECTRUM_LAPIS_BASE + "/sample/nucleotideMutations"
COV_SPECTRUM_AGGREGATED_ENDPOINT = COV_SPECTRUM_LAPIS_BASE + "/sample/aggregated"

# Module-level cache for the variant list (LAPIS aggregated by pango lineage).
# Key/value so we can invalidate after TTL without a full rebuild.
_VARIANT_CACHE = {'variants': None, 'timestamp': 0.0}
_VARIANT_CACHE_TTL_SECONDS = 3600  # 1 hour


def _fetch_all_lineages(timeout: float = 30.0) -> List[Dict[str, Any]]:
    """
    Fetch the full list of Pango lineages with sample counts from LAPIS.
    Cached in-process for _VARIANT_CACHE_TTL_SECONDS.

    Returns:
        List of {'lineage': str, 'count': int}, sorted by count descending.
    """
    now = time.time()
    cached = _VARIANT_CACHE.get('variants')
    if cached and (now - _VARIANT_CACHE['timestamp']) < _VARIANT_CACHE_TTL_SECONDS:
        return cached

    url = COV_SPECTRUM_AGGREGATED_ENDPOINT + "?" + urllib.parse.urlencode({
        'fields': 'nextcladePangoLineage',
    })
    req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'CovMutEx-X/1.0'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            payload = json.loads(resp.read().decode('utf-8'))
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"cov-spectrum API returned HTTP {e.code}: {e.reason}")
    except urllib.error.URLError as e:
        raise RuntimeError(f"cov-spectrum API request failed: {e.reason}")
    except (ValueError, json.JSONDecodeError) as e:
        raise RuntimeError(f"cov-spectrum API returned invalid JSON: {e}")

    data = payload.get('data', [])
    variants = []
    for row in data:
        lineage = row.get('nextcladePangoLineage') or row.get('pangoLineage')
        count = row.get('count', 0)
        if lineage:  # Skip null lineages
            variants.append({'lineage': str(lineage), 'count': int(count)})

    variants.sort(key=lambda v: -v['count'])
    _VARIANT_CACHE['variants'] = variants
    _VARIANT_CACHE['timestamp'] = now
    return variants


def search_cov_spectrum_variants(query: str = "", limit: int = 50) -> List[Dict[str, Any]]:
    """
    Search Pango lineages by substring match, sorted by sample count.

    Args:
        query: Substring to match against lineage name (case-insensitive). Empty
            string returns the most common variants.
        limit: Maximum number of results.

    Returns:
        List of {'lineage': str, 'count': int}.
    """
    all_variants = _fetch_all_lineages()
    q = (query or "").strip().upper()
    if not q:
        return all_variants[:limit]
    matches = [v for v in all_variants if q in v['lineage'].upper()]
    return matches[:limit]


def fetch_cov_spectrum_mutations(
    lineage: Optional[str] = None,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
    min_proportion: float = 0.0,
    country: Optional[str] = None,
    timeout: float = 30.0,
) -> Dict[str, Any]:
    """
    Fetch nucleotide mutations from cov-spectrum.org LAPIS API.

    Args:
        lineage: Pango lineage (e.g. "BA.1", "XBB.1.5"). Optional.
        date_from: ISO date "YYYY-MM-DD" (inclusive). Optional.
        date_to: ISO date "YYYY-MM-DD" (inclusive). Optional.
        min_proportion: Only return mutations above this proportion (server-side filter).
        country: Optional country filter (e.g. "Turkey", "USA").
        timeout: HTTP timeout in seconds.

    Returns:
        {
            'data': list of {mutation, proportion, count, ...} from the API,
            'query_params': dict of parameters sent,
            'url': full request URL,
        }

    Raises:
        RuntimeError on HTTP/network errors (caller should catch and surface to user).
    """
    params = {}
    if lineage:
        params['nextcladePangoLineage'] = lineage
    if date_from:
        params['dateFrom'] = date_from
    if date_to:
        params['dateTo'] = date_to
    if country:
        params['country'] = country
    if min_proportion and min_proportion > 0:
        params['minProportion'] = str(min_proportion)

    url = COV_SPECTRUM_NUC_MUTATIONS_ENDPOINT
    if params:
        url = url + "?" + urllib.parse.urlencode(params)

    req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'CovMutEx-X/1.0'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read()
            payload = json.loads(raw.decode('utf-8'))
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"cov-spectrum API returned HTTP {e.code}: {e.reason}")
    except urllib.error.URLError as e:
        raise RuntimeError(f"cov-spectrum API request failed: {e.reason}")
    except (ValueError, json.JSONDecodeError) as e:
        raise RuntimeError(f"cov-spectrum API returned invalid JSON: {e}")

    data = payload.get('data', [])
    if not isinstance(data, list):
        raise RuntimeError(f"Unexpected cov-spectrum API response: missing 'data' list")

    return {
        'data': data,
        'query_params': params,
        'url': url,
    }


def build_ground_truth_from_api_data(
    api_data: List[Dict[str, Any]],
    genome_length: int = 29904,
    region: tuple = None,
) -> Dict[str, Any]:
    """
    Build ground-truth array from a cov-spectrum LAPIS API response list.

    Each item is expected to have at least 'mutation' (string) and 'proportion' (float).
    Internally converts to the same CSV format consumed by `build_ground_truth_from_csv`
    so the parsing/aggregation logic stays in one place.
    """
    # Serialize API items to an in-memory CSV string matching cov-spectrum.org's CSV export
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(['mutation', 'proportion', 'count', 'jaccard'])
    for item in api_data:
        mutation = item.get('mutation', '')
        proportion = item.get('proportion', '')
        count = item.get('count', '')
        # cov-spectrum's API doesn't always return jaccard; leave blank
        jaccard = item.get('jaccard', '')
        writer.writerow([mutation, proportion, count, jaccard])
    return build_ground_truth_from_csv(buf.getvalue(), genome_length=genome_length, region=region)


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
    model_name: str = "unknown",
    gt_binarize_threshold: float = 0.5,
    pred_binarize_threshold: float = 0.5,
) -> Dict[str, Any]:
    """
    Compute all benchmark metrics for a single model.

    Args:
        predictions: Model predictions (mutation probabilities per position)
        ground_truth: Array of per-position values.
            - Binary (0/1) when derived from mutations.txt.
            - Continuous [0, 1] probabilities when derived from a cov-spectrum.org
              CSV — in that case Brier/ECE use the probabilities directly, while
              AUROC/AUPRC binarize at `gt_binarize_threshold` (since sklearn
              requires binary y_true).
        model_name: Name for logging
        gt_binarize_threshold: Threshold used to binarize a probability ground
            truth for ROC/PR metrics (positions with value >= threshold are
            treated as positives).

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

    # Detect whether ground truth is already binary (0/1) or continuous probabilities
    unique_vals = np.unique(gt)
    is_binary_gt = len(unique_vals) <= 2 and set(unique_vals.tolist()).issubset({0.0, 1.0})

    # Binarized version used only for metrics that require binary labels
    if is_binary_gt:
        gt_binary = gt.astype(np.float32)
    else:
        gt_binary = (gt >= gt_binarize_threshold).astype(np.float32)

    metrics['ground_truth_type'] = 'binary' if is_binary_gt else 'probability'
    metrics['gt_binarize_threshold'] = None if is_binary_gt else float(gt_binarize_threshold)

    # --- Brier Score ---
    # brier_score_loss requires binary y_true, so for probability GT we compute
    # the Brier score manually as mean((pred - gt)^2) — a natural extension.
    try:
        if is_binary_gt:
            metrics['brier_score'] = float(brier_score_loss(gt, preds))
        else:
            metrics['brier_score'] = float(np.mean((preds - gt) ** 2))
    except Exception as e:
        metrics['brier_score'] = None
        print(f"[BENCH] Brier score failed for {model_name}: {e}")

    # --- AUROC ---
    try:
        if len(np.unique(gt_binary)) > 1:  # Need both classes
            metrics['auroc'] = float(roc_auc_score(gt_binary, preds))

            # ROC curve data (downsample for frontend)
            fpr, tpr, _ = roc_curve(gt_binary, preds)
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
        if len(np.unique(gt_binary)) > 1:
            metrics['auprc'] = float(average_precision_score(gt_binary, preds))

            # PR curve data
            precision, recall, _ = precision_recall_curve(gt_binary, preds)
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

    # --- Precision / Recall / F1 ---
    # All three need binary predictions and binary ground truth. Predictions are
    # binarized at `pred_binarize_threshold`; GT uses gt_binary (binarized at
    # gt_binarize_threshold for probability GT, identity for binary GT).
    # zero_division=0 prevents warnings/NaN when no positives are predicted.
    try:
        preds_binary = (preds >= pred_binarize_threshold).astype(np.int32)
        gt_int = gt_binary.astype(np.int32)
        metrics['precision'] = float(precision_score(gt_int, preds_binary, zero_division=0))
        metrics['recall'] = float(recall_score(gt_int, preds_binary, zero_division=0))
        metrics['f1_score'] = float(f1_score(gt_int, preds_binary, zero_division=0))
        # Matthews Correlation Coefficient: balanced measure (in [-1, 1]) that
        # stays meaningful under heavy class imbalance. matthews_corrcoef returns
        # 0.0 when a class is entirely absent from preds or gt (no warning).
        metrics['mcc'] = float(matthews_corrcoef(gt_int, preds_binary))
        metrics['pred_binarize_threshold'] = float(pred_binarize_threshold)
        # Confusion-matrix-style counts for transparency
        tp = int(np.sum((preds_binary == 1) & (gt_int == 1)))
        fp = int(np.sum((preds_binary == 1) & (gt_int == 0)))
        fn = int(np.sum((preds_binary == 0) & (gt_int == 1)))
        tn = int(np.sum((preds_binary == 0) & (gt_int == 0)))
        metrics['confusion'] = {'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn}
    except Exception as e:
        metrics['precision'] = None
        metrics['recall'] = None
        metrics['f1_score'] = None
        metrics['mcc'] = None
        metrics['confusion'] = None
        print(f"[BENCH] Precision/Recall/F1/MCC failed for {model_name}: {e}")

    # --- ECE (Expected Calibration Error) ---
    # ECE uses ground truth directly (works with both binary and probability GT:
    # bin_acc = mean(gt[mask]) is a rate in both cases).
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
    # For binary GT: count of 1s. For probability GT: count of positions with nonzero proportion.
    metrics['num_mutations'] = int(np.sum(gt > 0))
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
    gt_binarize_threshold: float = 0.5,
    pred_binarize_threshold: float = 0.5,
) -> Dict[str, Dict]:
    """
    Compute metrics per protein region.

    Args:
        predictions: Full-genome predictions
        ground_truth: Full-genome ground truth (binary 0/1 or continuous [0, 1])
        protein_regions: Dict mapping protein name to [start, end]
        gt_binarize_threshold: Threshold for binarizing probability GT for AUROC/AUPRC

    Returns:
        Dict mapping protein name to metrics dict
    """
    per_protein = {}

    for protein, (start, end) in protein_regions.items():
        region_preds = predictions[start:end + 1] if end + 1 <= len(predictions) else predictions[start:]
        region_gt = ground_truth[start:end + 1] if end + 1 <= len(ground_truth) else ground_truth[start:]

        if len(region_preds) == 0:
            continue

        # Detect binary vs probability GT (per region)
        unique_vals = np.unique(region_gt)
        is_binary_gt = len(unique_vals) <= 2 and set(unique_vals.tolist()).issubset({0.0, 1.0})
        region_gt_binary = region_gt.astype(np.float32) if is_binary_gt else (region_gt >= gt_binarize_threshold).astype(np.float32)

        region_metrics = {
            'mean_prediction': float(np.mean(region_preds)),
            'max_prediction': float(np.max(region_preds)),
            'num_positions': int(len(region_preds)),
            'num_mutations': int(np.sum(region_gt > 0)),
            'mutation_rate': float(np.mean(region_gt)),
        }

        # Only compute AUROC/AUPRC if both classes present after binarization
        if len(np.unique(region_gt_binary)) > 1 and len(region_gt_binary) > 10:
            try:
                region_metrics['auroc'] = float(roc_auc_score(region_gt_binary, region_preds))
                region_metrics['auprc'] = float(average_precision_score(region_gt_binary, region_preds))
                if is_binary_gt:
                    region_metrics['brier_score'] = float(brier_score_loss(region_gt, region_preds))
                else:
                    region_metrics['brier_score'] = float(np.mean((region_preds - region_gt) ** 2))
            except Exception:
                region_metrics['auroc'] = None
                region_metrics['auprc'] = None
                region_metrics['brier_score'] = None
        else:
            region_metrics['auroc'] = None
            region_metrics['auprc'] = None
            # Brier is still computable even with one class
            try:
                if is_binary_gt:
                    region_metrics['brier_score'] = float(brier_score_loss(region_gt, region_preds))
                else:
                    region_metrics['brier_score'] = float(np.mean((region_preds - region_gt) ** 2))
            except Exception:
                region_metrics['brier_score'] = None

        # Precision / Recall / F1 / MCC (works with one-class GT; zero_division=0 keeps it safe)
        try:
            region_preds_binary = (region_preds >= pred_binarize_threshold).astype(np.int32)
            region_gt_int = region_gt_binary.astype(np.int32)
            region_metrics['precision'] = float(precision_score(region_gt_int, region_preds_binary, zero_division=0))
            region_metrics['recall'] = float(recall_score(region_gt_int, region_preds_binary, zero_division=0))
            region_metrics['f1_score'] = float(f1_score(region_gt_int, region_preds_binary, zero_division=0))
            region_metrics['mcc'] = float(matthews_corrcoef(region_gt_int, region_preds_binary))
        except Exception:
            region_metrics['precision'] = None
            region_metrics['recall'] = None
            region_metrics['f1_score'] = None
            region_metrics['mcc'] = None

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
    custom_ground_truth: Optional[np.ndarray] = None,
    custom_ground_truth_full: Optional[np.ndarray] = None,
    ground_truth_source: str = "mutations_txt",
    ground_truth_meta: Optional[Dict[str, Any]] = None,
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
        mutations: Parsed mutations list (used for default ground truth if
            `custom_ground_truth` is not provided)
        genome_sequence: Reference genome string
        protein_regions: Dict of protein name → [start, end]
        predict_fn: Function that takes prediction params and returns predictions array
        selected_protein_region: Optional protein region to focus on
        custom_ground_truth: Optional pre-built ground-truth array (already
            sliced to the selected region). When provided, overrides the
            default mutations.txt-based ground truth. May be continuous (e.g.
            from cov-spectrum.org proportions) rather than binary.
        custom_ground_truth_full: Optional full-genome version of the custom
            ground truth (unsliced). Used for per-protein breakdown.
        ground_truth_source: Label for the GT source ("mutations_txt",
            "cov_spectrum_csv", etc.) stored in results metadata.
        ground_truth_meta: Extra metadata about the ground truth (e.g. parsed
            count, skipped examples) stored in results for transparency.

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
            'ground_truth_source': ground_truth_source,
            'ground_truth_meta': ground_truth_meta or {},
        },
        'models': {},
        'model_agreement': {},
        'summary': {},
    }

    # Build ground truth
    genome_length = len(genome_sequence) if genome_sequence else 29904

    if custom_ground_truth is not None:
        ground_truth = np.asarray(custom_ground_truth, dtype=np.float32)
    else:
        if selected_protein_region and selected_protein_region in protein_regions:
            region = tuple(protein_regions[selected_protein_region])
            ground_truth = build_ground_truth(mutations, genome_length, region=region)
        else:
            ground_truth = build_ground_truth(mutations, genome_length)

    results['parameters']['ground_truth_positives'] = int(np.sum(ground_truth > 0))
    results['parameters']['ground_truth_total'] = int(len(ground_truth))
    # Positions with nonzero ground truth (works for both binary and probability GT)
    nonzero_idx = np.where(ground_truth > 0)[0]
    results['parameters']['mutation_positions'] = [int(i) for i in nonzero_idx]
    # Corresponding GT values at those positions (binary -> all 1.0, probability -> proportions)
    results['parameters']['mutation_values'] = [round(float(ground_truth[i]), 6) for i in nonzero_idx]
    # GT type hint for frontend rendering
    unique_gt = np.unique(ground_truth)
    is_binary_gt = len(unique_gt) <= 2 and set(unique_gt.tolist()).issubset({0.0, 1.0})
    results['parameters']['ground_truth_type'] = 'binary' if is_binary_gt else 'probability'
    
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
                if custom_ground_truth_full is not None:
                    full_gt = np.asarray(custom_ground_truth_full, dtype=np.float32)
                else:
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