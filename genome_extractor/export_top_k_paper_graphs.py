"""
Export paper-ready Top-K precision/recall figures for the Known Hotspot Case Study.

Outputs per model:
  - a long-form CSV with exact/proximity sweep metrics across variants
  - a summary CSV with the best exact/proximity F1 operating points
  - a PNG figure per model/variant pair with all Top-K curves

Example:
  /Users/mac/miniconda3/envs/covmutex/bin/python export_top_k_paper_graphs.py \
    --models balanced_data_model PRIEST \
    --max-k 300 \
    --output-dir paper_figures/top_k_sweeps
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, "/Users/mac/Documents/ITU_Bioinformatics/CovMutEx-main/genome_extractor")
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "genome_extractor.settings")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-covmutex")

import django

django.setup()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from genome.configs import configs
from genome.delta_omicron_retrospective import (  # noqa: E402
    PRECOMPUTED_VARIANT_CONTEXTS,
    SPIKE_START,
    build_precomputed_variant_mutation_tuples,
    build_ranked_spike_site_rows,
    build_top_k_metric_sweep,
    extract_priest_spike_site_score_rows,
    load_precomputed_variant_spike_site_set,
    rank_spike_site_rows,
)
from genome.feature_extractor import (  # noqa: E402
    construct_variant_genome,
    load_legacy_keras_model,
    predict_mutations,
)
from genome.priest_annotations import PROTEIN_REGIONS  # noqa: E402


REPO_ROOT = Path("/Users/mac/Documents/ITU_Bioinformatics/CovMutEx-main")
BASE_DIR = REPO_ROOT / "genome_extractor"
MODEL_DIR = BASE_DIR / "covid19_models" / "models"
GENOME_PATH = BASE_DIR / "genome" / "genome.txt"
CACHE_PATH = BASE_DIR / "genome" / "node_features.h5"
CODON_PATH = BASE_DIR / "genome" / "codon_aa_mapping.json"
SPIKE_REGION = {"S": (PROTEIN_REGIONS["S"][0] - 1, PROTEIN_REGIONS["S"][1] - 1)}
DEFAULT_MODELS = ["balanced_data_model", "PRIEST"]
PROXIMITY_WINDOW = 3


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--models",
        nargs="+",
        default=DEFAULT_MODELS,
        help="Model names to export. Supported: balanced_data_model, single_input_ensemble, multi_input_ensemble, PRIEST",
    )
    parser.add_argument(
        "--max-k",
        type=int,
        default=300,
        help="Maximum Top-K value to include in the exported figures and CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        default="paper_figures/top_k_sweeps",
        help="Directory where CSV and PNG outputs will be written.",
    )
    return parser.parse_args()


def read_genome_sequence(path: Path) -> str:
    with path.open() as handle:
        next(handle)
        return "".join(line.strip() for line in handle)


def resolve_model_path(model_name: str) -> Optional[Path]:
    if model_name.upper() == "PRIEST":
        return None
    if model_name == "balanced_data_model":
        fixed_path = MODEL_DIR / "balanced_data_model_fixed.keras"
        if fixed_path.exists():
            return fixed_path
        return MODEL_DIR / "balanced_data_model.keras"
    if model_name == "single_input_ensemble":
        return MODEL_DIR / "single_input_ensemble_model_best.keras"
    if model_name == "multi_input_ensemble":
        return MODEL_DIR / "multi_input_ensemble_model.keras"
    raise ValueError(f"Unsupported model name: {model_name}")


def slugify(value: str) -> str:
    cleaned = [
        character if character.isalnum() else "_"
        for character in (value or "").strip()
    ]
    slug = "".join(cleaned).strip("_")
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug or "model"


def f1_score(precision: float, recall: float) -> float:
    return (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0


def run_model_predictions(model_obj, genome_seq, mutation_tuples, node_id):
    return predict_mutations(
        cache_path=str(CACHE_PATH),
        genome_seq=genome_seq,
        mutations=mutation_tuples,
        codon_mapper=str(CODON_PATH),
        config_file=configs(),
        node_ids=[node_id],
        elapsed_day=0,
        depth=0,
        protein_regions=SPIKE_REGION,
        model=model_obj,
    )


def build_ranked_rows_for_variant(model_name, ctx, genome_seq, model_obj=None):
    node_id = ctx["node_id"]
    if model_name.upper() == "PRIEST":
        site_rows, _ = extract_priest_spike_site_score_rows(node_context=ctx)
        return rank_spike_site_rows(site_rows)

    mutation_build = build_precomputed_variant_mutation_tuples(
        reference_genome_sequence=genome_seq,
        node_context=ctx,
    )
    mutation_tuples = mutation_build["mutation_tuples"]
    variant_genome = construct_variant_genome(genome_seq, mutation_tuples)
    predictions = run_model_predictions(model_obj, variant_genome, mutation_tuples, node_id)
    ranked_rows, _ = build_ranked_spike_site_rows(
        selected_model=model_name,
        node_context=ctx,
        predictions=predictions,
        reference_genome_sequence=genome_seq,
        spike_start_genome_position=SPIKE_START,
        prediction_start_genome_position=SPIKE_START,
    )
    return ranked_rows


def export_csv_rows(path: Path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def pick_tick_values(max_k: int):
    preferred_ticks = [1, 10, 25, 50, 100, 150, 200, 250, 300]
    ticks = [tick for tick in preferred_ticks if tick <= max_k]
    if max_k not in ticks:
        ticks.append(max_k)
    return sorted(set(ticks))


def variant_slug(ctx):
    nickname = ctx.get("variant_nickname")
    lineage = ctx.get("variant_label") or ctx.get("pangolin_lineage") or "variant"
    if nickname:
        return "{}_{}".format(slugify(nickname), slugify(lineage))
    return slugify(lineage)


def render_variant_png(model_name, variant_payload, output_path: Path, max_k: int):
    sweep_rows = [row for row in variant_payload["sweep_rows"] if row["top_k"] <= max_k]
    x_values = [row["top_k"] for row in sweep_rows]
    ticks = pick_tick_values(max(x_values) if x_values else max_k)

    figure, axis = plt.subplots(figsize=(8.6, 5.6), dpi=300, constrained_layout=True)
    figure.patch.set_facecolor("#f8fafc")
    axis.set_facecolor("#ffffff")
    axis.grid(True, which="major", axis="both", color="#e2e8f0", linestyle="--", linewidth=0.7)

    axis.plot(
        x_values,
        [row["precision_at_k"] for row in sweep_rows],
        color="#0284c7",
        linewidth=2.4,
        label="Exact precision",
    )
    axis.plot(
        x_values,
        [row["recall_against_comparison_sites"] for row in sweep_rows],
        color="#e11d48",
        linewidth=2.4,
        label="Exact recall",
    )
    axis.plot(
        x_values,
        [row["proximity_precision_at_k"] for row in sweep_rows],
        color="#0f766e",
        linewidth=2.4,
        label="Proximity precision",
    )
    axis.plot(
        x_values,
        [row["proximity_recall"] for row in sweep_rows],
        color="#d97706",
        linewidth=2.4,
        label="Proximity recall",
    )

    axis.set_xlim(1, max(x_values) if x_values else max_k)
    axis.set_ylim(0, 1)
    axis.set_xticks(ticks)
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))
    axis.tick_params(axis="x", labelsize=10, colors="#475569")
    axis.tick_params(axis="y", labelsize=10, colors="#475569")
    axis.set_xlabel("Top-K", fontsize=11, color="#334155")
    axis.set_ylabel("Metric value", fontsize=11, color="#334155")
    for spine in axis.spines.values():
        spine.set_color("#cbd5e1")

    figure.suptitle(
        "{} — {}".format(model_name, variant_payload["variant_display_name"]),
        fontsize=16,
        fontweight="bold",
        color="#0f172a",
    )
    axis.set_title(
        "{} | known Spike sites={} | max K={}".format(
            variant_payload["emergence_label"],
            variant_payload["known_site_count"],
            max_k,
        ),
        fontsize=10,
        color="#475569",
        pad=10,
    )
    axis.legend(loc="upper right", frameon=False, fontsize=9)
    figure.savefig(output_path, dpi=300, bbox_inches="tight", facecolor=figure.get_facecolor())
    plt.close(figure)


def export_model_graphs(model_name: str, max_k: int, output_dir: Path):
    print(f"[export] model={model_name}")
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_sequence = read_genome_sequence(GENOME_PATH)
    model_path = resolve_model_path(model_name)
    model_obj = None

    if model_path is not None:
        print(f"  loading model from {model_path}")
        model_obj = load_legacy_keras_model(str(model_path))

    long_rows = []
    summary_rows = []
    variant_output_dir = output_dir / "per_variant"
    variant_output_dir.mkdir(parents=True, exist_ok=True)
    variants = sorted(
        PRECOMPUTED_VARIANT_CONTEXTS.items(),
        key=lambda item: item[1].get("sort_order", 999),
    )

    for node_id, ctx in variants:
        variant_display_name = ctx.get("variant_display_name") or ctx["variant_label"]
        print(f"  variant={variant_display_name}")
        known_positions = list(load_precomputed_variant_spike_site_set(node_id)["positions"])
        ranked_rows = build_ranked_rows_for_variant(model_name, ctx, genome_sequence, model_obj=model_obj)
        sweep_rows = build_top_k_metric_sweep(
            ranked_rows,
            known_positions,
            proximity_window=PROXIMITY_WINDOW,
        )
        truncated_rows = [row for row in sweep_rows if row["top_k"] <= max_k]

        for row in truncated_rows:
            exact_f1 = f1_score(
                row["precision_at_k"],
                row["recall_against_comparison_sites"],
            )
            proximity_f1 = f1_score(
                row["proximity_precision_at_k"],
                row["proximity_recall"],
            )
            long_rows.append(
                {
                    "model_name": model_name,
                    "node_id": node_id,
                    "variant_label": ctx["variant_label"],
                    "variant_display_name": variant_display_name,
                    "variant_nickname": ctx.get("variant_nickname") or "",
                    "emergence_label": ctx.get("emergence_label") or "",
                    "top_k": row["top_k"],
                    "exact_overlap_count": row["overlap_count"],
                    "exact_precision_at_k": row["precision_at_k"],
                    "exact_recall_at_k": row["recall_against_comparison_sites"],
                    "exact_f1_at_k": exact_f1,
                    "proximity_overlap_count": row["proximity_overlap_count"],
                    "proximity_precision_at_k": row["proximity_precision_at_k"],
                    "proximity_recall_at_k": row["proximity_recall"],
                    "proximity_f1_at_k": proximity_f1,
                }
            )

        best_exact = max(
            truncated_rows,
            key=lambda row: (
                f1_score(row["precision_at_k"], row["recall_against_comparison_sites"]),
                -row["top_k"],
            ),
        )
        best_proximity = max(
            truncated_rows,
            key=lambda row: (
                f1_score(row["proximity_precision_at_k"], row["proximity_recall"]),
                -row["top_k"],
            ),
        )
        summary_rows.append(
            {
                "model_name": model_name,
                "node_id": node_id,
                "variant_label": ctx["variant_label"],
                "variant_display_name": variant_display_name,
                "emergence_label": ctx.get("emergence_label") or "",
                "known_site_count": len(known_positions),
                "best_exact_k": best_exact["top_k"],
                "best_exact_precision": best_exact["precision_at_k"],
                "best_exact_recall": best_exact["recall_against_comparison_sites"],
                "best_exact_f1": f1_score(
                    best_exact["precision_at_k"],
                    best_exact["recall_against_comparison_sites"],
                ),
                "best_proximity_k": best_proximity["top_k"],
                "best_proximity_precision": best_proximity["proximity_precision_at_k"],
                "best_proximity_recall": best_proximity["proximity_recall"],
                "best_proximity_f1": f1_score(
                    best_proximity["proximity_precision_at_k"],
                    best_proximity["proximity_recall"],
                ),
            }
        )
        variant_payload = {
            "variant_display_name": variant_display_name,
            "variant_slug": variant_slug(ctx),
            "emergence_label": ctx.get("emergence_label") or "",
            "known_site_count": len(known_positions),
            "sweep_rows": truncated_rows,
        }
        render_variant_png(
            model_name=model_name,
            variant_payload=variant_payload,
            output_path=variant_output_dir / "{}_{}.png".format(
                slugify(model_name),
                variant_payload["variant_slug"],
            ),
            max_k=max_k,
        )

    model_slug = slugify(model_name)
    export_csv_rows(
        output_dir / f"{model_slug}_top_k_sweep.csv",
        long_rows,
        [
            "model_name",
            "node_id",
            "variant_label",
            "variant_display_name",
            "variant_nickname",
            "emergence_label",
            "top_k",
            "exact_overlap_count",
            "exact_precision_at_k",
            "exact_recall_at_k",
            "exact_f1_at_k",
            "proximity_overlap_count",
            "proximity_precision_at_k",
            "proximity_recall_at_k",
            "proximity_f1_at_k",
        ],
    )
    export_csv_rows(
        output_dir / f"{model_slug}_summary.csv",
        summary_rows,
        [
            "model_name",
            "node_id",
            "variant_label",
            "variant_display_name",
            "emergence_label",
            "known_site_count",
            "best_exact_k",
            "best_exact_precision",
            "best_exact_recall",
            "best_exact_f1",
            "best_proximity_k",
            "best_proximity_precision",
            "best_proximity_recall",
            "best_proximity_f1",
        ],
    )
    print(f"  wrote outputs to {output_dir}")


def main():
    args = parse_args()
    output_dir = REPO_ROOT / args.output_dir
    for model_name in args.models:
        export_model_graphs(model_name=model_name, max_k=args.max_k, output_dir=output_dir)


if __name__ == "__main__":
    main()
