"""
CovMutEx-X Benchmark Datasets & Reproducibility

FR-3.1: Curated evaluation splits (predefined variant sets for benchmarking)
FR-3.3: Reproducibility (fixed seeds, deterministic preprocessing, environment capture)

Datasets are predefined sets of (nodeId, elapsedDay) pairs grouped by category.
Running a "dataset benchmark" means running ALL models on ALL variants in the set,
then aggregating the metrics.
"""

import os
import json
import time
import random
import platform
import numpy as np
from typing import Dict, List, Any, Optional

# Try importing version info for reproducibility
try:
    import tensorflow as tf
    TF_VERSION = tf.__version__
except ImportError:
    TF_VERSION = "not installed"

try:
    import torch
    TORCH_VERSION = torch.__version__
except ImportError:
    TORCH_VERSION = "not installed"

try:
    import sklearn
    SKLEARN_VERSION = sklearn.__version__
except ImportError:
    SKLEARN_VERSION = "not installed"


# ============================================
# FR-3.1: CURATED BENCHMARK DATASETS
# ============================================

# Each dataset is a list of (nodeId, elapsedDay) pairs
# These represent different evaluation scenarios

BENCHMARK_DATASETS = {
    "alpha_variants": {
        "name": "Alpha Variants (B.1.1.7)",
        "description": "Alpha variant samples from different geographic regions and time points. "
                       "Tests model generalization across Alpha lineage.",
        "category": "lineage",
        "variants": [
            {
                "nodeId": "England/MILK-9E05B3/2020|OV826873.1|2020-12-09",
                "elapsedDay": 30,
                "label": "England Alpha Early",
            },
            {
                "nodeId": "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07",
                "elapsedDay": 60,
                "label": "USA Alpha Mid",
            },
        ],
    },
    "delta_variants": {
        "name": "Delta Variants (B.1.617.2)",
        "description": "Delta variant samples. Tests performance on highly mutated lineage.",
        "category": "lineage",
        "variants": [
            {
                "nodeId": "EGY/CCHE57357_Wave_3_A029/2021|MZ380261.1|2021-05-11",
                "elapsedDay": 110,
                "label": "Egypt Delta",
            },
        ],
    },
    "temporal_progression": {
        "name": "Temporal Progression",
        "description": "Same variant at different elapsed days. Tests how models respond to time progression.",
        "category": "temporal",
        "variants": [
            {
                "nodeId": "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07",
                "elapsedDay": 30,
                "label": "USA 30 days",
            },
            {
                "nodeId": "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07",
                "elapsedDay": 60,
                "label": "USA 60 days",
            },
            {
                "nodeId": "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07",
                "elapsedDay": 120,
                "label": "USA 120 days",
            },
        ],
    },
    "quick_test": {
        "name": "Quick Test (Single Variant)",
        "description": "Single variant for quick smoke testing. Use this for rapid iteration.",
        "category": "test",
        "variants": [
            {
                "nodeId": "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07",
                "elapsedDay": 60,
                "label": "Quick Test",
            },
        ],
    },
}

# Datasets directory for custom user-created datasets
DATASETS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'benchmark_datasets'
)


def get_available_datasets() -> Dict[str, Any]:
    """
    Get all available benchmark datasets (built-in + user-created).
    
    Returns dict of dataset_id -> dataset info (without full variant list for summary)
    """
    datasets = {}
    
    # Built-in datasets
    for key, ds in BENCHMARK_DATASETS.items():
        datasets[key] = {
            "name": ds["name"],
            "description": ds["description"],
            "category": ds["category"],
            "num_variants": len(ds["variants"]),
            "type": "built-in",
        }
    
    # User-created datasets (from filesystem)
    os.makedirs(DATASETS_DIR, exist_ok=True)
    for filename in os.listdir(DATASETS_DIR):
        if filename.endswith('.json'):
            filepath = os.path.join(DATASETS_DIR, filename)
            try:
                with open(filepath, 'r') as f:
                    ds = json.load(f)
                key = filename.replace('.json', '')
                datasets[key] = {
                    "name": ds.get("name", key),
                    "description": ds.get("description", ""),
                    "category": ds.get("category", "custom"),
                    "num_variants": len(ds.get("variants", [])),
                    "type": "custom",
                }
            except Exception:
                continue
    
    return datasets


def get_dataset(dataset_id: str) -> Optional[Dict]:
    """Get a specific dataset by ID (built-in or custom)."""
    # Check built-in first
    if dataset_id in BENCHMARK_DATASETS:
        return BENCHMARK_DATASETS[dataset_id]
    
    # Check custom
    filepath = os.path.join(DATASETS_DIR, f"{dataset_id}.json")
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    
    return None


def save_custom_dataset(dataset_id: str, dataset: Dict) -> str:
    """Save a user-created dataset."""
    os.makedirs(DATASETS_DIR, exist_ok=True)
    filepath = os.path.join(DATASETS_DIR, f"{dataset_id}.json")
    with open(filepath, 'w') as f:
        json.dump(dataset, f, indent=2)
    return filepath


# ============================================
# FR-3.3: REPRODUCIBILITY
# ============================================

DEFAULT_SEED = 42


def set_reproducibility_seed(seed: int = DEFAULT_SEED):
    """
    Set all random seeds for reproducible results.
    
    Sets seeds for: Python random, NumPy, TensorFlow, PyTorch
    """
    random.seed(seed)
    np.random.seed(seed)
    
    try:
        import tensorflow as tf
        tf.random.set_seed(seed)
        # Deterministic ops (may slow down computation)
        os.environ['TF_DETERMINISTIC_OPS'] = '1'
    except ImportError:
        pass
    
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass
    
    return seed


def capture_environment() -> Dict[str, Any]:
    """
    Capture the current execution environment for reproducibility.
    
    Records: Python version, OS, package versions, seeds, hardware info.
    """
    env = {
        "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
        "python_version": platform.python_version(),
        "os": f"{platform.system()} {platform.release()}",
        "architecture": platform.machine(),
        "packages": {
            "tensorflow": TF_VERSION,
            "pytorch": TORCH_VERSION,
            "scikit-learn": SKLEARN_VERSION,
            "numpy": np.__version__,
        },
        "seed": DEFAULT_SEED,
    }
    
    # GPU info
    try:
        import tensorflow as tf
        gpus = tf.config.list_physical_devices('GPU')
        env["gpu_available"] = len(gpus) > 0
        env["gpu_devices"] = [g.name for g in gpus]
    except Exception:
        env["gpu_available"] = False
        env["gpu_devices"] = []
    
    return env


def create_reproducibility_pack(benchmark_results: Dict, seed: int = DEFAULT_SEED) -> Dict:
    """
    Create a complete reproducibility pack that can recreate the benchmark.
    
    Includes: environment, seed, parameters, dataset, model configs.
    """
    pack = {
        "reproducibility_version": "1.0",
        "created_at": time.strftime('%Y-%m-%d %H:%M:%S'),
        "seed": seed,
        "environment": capture_environment(),
        "benchmark_parameters": benchmark_results.get("parameters", {}),
        "models_used": [],
        "instructions": (
            "To reproduce this benchmark:\n"
            "1. Ensure the same package versions listed in 'environment.packages'\n"
            "2. Use seed value from 'seed' field\n"
            "3. Run benchmark with the same parameters and models\n"
            "4. POST to /api/benchmark/run/ with the request body in 'replay_request'"
        ),
    }
    
    # Extract model info
    for model_name, model_data in benchmark_results.get("models", {}).items():
        pack["models_used"].append({
            "name": model_name,
            "source": model_data.get("source", "unknown"),
            "status": model_data.get("status", "unknown"),
        })
    
    # Create replay request
    pack["replay_request"] = {
        "models": [m["name"] if m["source"] == "server" else f"uploaded:{m['name']}" 
                   for m in pack["models_used"]],
        "nodeId": benchmark_results.get("parameters", {}).get("node_id"),
        "elapsedDay": benchmark_results.get("parameters", {}).get("elapsed_day"),
        "selectedProteinRegion": benchmark_results.get("parameters", {}).get("selected_protein_region"),
        "seed": seed,
    }
    
    return pack


# ============================================
# MULTI-VARIANT BENCHMARK AGGREGATION
# ============================================

def aggregate_multi_variant_results(
    per_variant_results: List[Dict],
) -> Dict[str, Any]:
    """
    Aggregate benchmark results across multiple variants.
    
    Takes a list of single-variant benchmark results and produces:
    - Mean/std of each metric across variants
    - Per-model ranking
    - Best model per metric
    """
    if not per_variant_results:
        return {}
    
    # Collect all model names across variants
    all_models = set()
    for vr in per_variant_results:
        all_models.update(vr.get("models", {}).keys())
    
    # Aggregate metrics per model
    aggregated = {}
    metric_keys = ['auroc', 'auprc', 'brier_score', 'ece', 'runtime_seconds', 'mean_prediction']
    
    for model_name in all_models:
        model_metrics = {k: [] for k in metric_keys}
        num_success = 0
        num_fail = 0
        
        for vr in per_variant_results:
            model_data = vr.get("models", {}).get(model_name)
            if not model_data:
                continue
            
            if model_data.get("status") == "success":
                num_success += 1
                for mk in metric_keys:
                    val = model_data.get("metrics", {}).get(mk)
                    if val is not None:
                        model_metrics[mk].append(val)
            else:
                num_fail += 1
        
        # Compute mean and std
        agg_metrics = {}
        for mk in metric_keys:
            values = model_metrics[mk]
            if values:
                agg_metrics[mk] = {
                    "mean": float(np.mean(values)),
                    "std": float(np.std(values)),
                    "min": float(np.min(values)),
                    "max": float(np.max(values)),
                    "n": len(values),
                    "values": [float(v) for v in values],
                }
            else:
                agg_metrics[mk] = None
        
        aggregated[model_name] = {
            "metrics": agg_metrics,
            "num_success": num_success,
            "num_fail": num_fail,
            "success_rate": num_success / max(1, num_success + num_fail),
        }
    
    # Compute rankings per metric
    rankings = {}
    for mk in metric_keys:
        higher_is_better = mk in ['auroc', 'auprc']
        lower_is_better = mk in ['brier_score', 'ece', 'runtime_seconds']
        
        scored = []
        for model_name, data in aggregated.items():
            metric_data = data["metrics"].get(mk)
            if metric_data and metric_data.get("mean") is not None:
                scored.append((model_name, metric_data["mean"]))
        
        if scored:
            if higher_is_better:
                scored.sort(key=lambda x: x[1], reverse=True)
            elif lower_is_better:
                scored.sort(key=lambda x: x[1])
            
            rankings[mk] = [
                {"rank": i + 1, "model": name, "value": round(val, 6)}
                for i, (name, val) in enumerate(scored)
            ]
    
    return {
        "models": aggregated,
        "rankings": rankings,
        "num_variants": len(per_variant_results),
        "total_models": len(all_models),
    }