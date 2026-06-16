import json
import os
import time
import traceback
from dataclasses import dataclass
from typing import Any, Dict, Optional

import logomaker
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from django.http import HttpResponse, JsonResponse
from rest_framework.decorators import api_view

from .configs import configs
from .covmutex_feature_extractors import load_feature_extractor
from .covmutex_models import load_model as load_covmutex_model
from .feature_extractor_updated import (
    add_terminal_padding_nucleotide,
    construct_variant_genome,
    get_sample_depth,
    parse_mutations,
)
from .helpers import (
    calculate_genome_data,
    calculate_protein_region_probabilities,
    measure_time,
    predict_mutations,
    read_genome_sequence,
)
from .cache_paths import CACHE_DIR, NODE_FEATURES_CACHE_PATH
from .plugin_runtime import (
    BundleResolution,
    collect_helper_uploads,
    normalize_elapsed_day,
    parse_custom_parameters,
    resolve_server_model,
    resolve_uploaded_bundle,
    save_uploaded_bundle,
)

matplotlib.use("Agg")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GENOME_PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))
UPLOADED_MODELS_DIR = os.path.join(GENOME_PROJECT_DIR, "uploaded_models")
MODEL_DIRECTORY = os.path.join(GENOME_PROJECT_DIR, "covid19_models", "models")

CODON_MAPPING_PATH = os.path.join(BASE_DIR, "codon_aa_mapping.json")
CACHE_PATH = NODE_FEATURES_CACHE_PATH
DEPTH_FILE = os.path.join(BASE_DIR, "depth_date.json")
GENOME_FILE_PATH = os.path.join(BASE_DIR, "genome.txt")

PROTEIN_REGIONS = {
    "ORF1ab": [266, 21555],
    "S": [21563, 25384],
    "ORF3a": [25393, 26220],
    "E": [26245, 26472],
    "M": [26523, 27191],
    "ORF6": [27202, 27387],
    "ORF7a": [27394, 27759],
    "ORF7b": [27756, 27887],
    "ORF8": [27894, 28259],
    "N": [28274, 29533],
    "ORF10": [29558, 29674],
}


@dataclass
class PredictionContext:
    node_id: Optional[str]
    node_ids: list[str]
    elapsed_day: int
    selected_protein_region: Optional[str]
    selected_protein_regions: Optional[Dict[str, list[int]]]
    depth: float
    bundle: BundleResolution
    custom_parameters: Dict[str, Any]


def _request_data(request):
    return request.data if request.method == "POST" else request.GET


def _request_files(request):
    return getattr(request, "FILES", getattr(request, "_request", request).FILES)


def _selected_region_map(selected_region: Optional[str]) -> Optional[Dict[str, list[int]]]:
    if selected_region and selected_region in PROTEIN_REGIONS:
        return {selected_region: PROTEIN_REGIONS[selected_region]}
    return None


def _list_uploaded_bundle_names() -> list[str]:
    if not os.path.exists(UPLOADED_MODELS_DIR):
        return []
    return [
        entry
        for entry in os.listdir(UPLOADED_MODELS_DIR)
        if os.path.isdir(os.path.join(UPLOADED_MODELS_DIR, entry)) and not entry.startswith("__")
    ]


def prepare_prediction_context(request) -> PredictionContext:
    data = _request_data(request)
    files = _request_files(request)
    os.makedirs(CACHE_DIR, exist_ok=True)

    node_id = data.get("nodeId")
    elapsed_day = normalize_elapsed_day(data.get("elapsedDay")) or 0
    selected_model = data.get("selectedModel")
    selected_protein_region = data.get("selectedProteinRegion")
    custom_parameters = parse_custom_parameters(data.get("customParameters", "{}"))
    helper_uploads = collect_helper_uploads(files, data)

    uploaded_model = files.get("modelFile")
    uploaded_extractor = files.get("extractorFile")

    os.makedirs(UPLOADED_MODELS_DIR, exist_ok=True)

    if uploaded_model:
        bundle = save_uploaded_bundle(
            uploaded_models_dir=UPLOADED_MODELS_DIR,
            uploaded_model=uploaded_model,
            uploaded_extractor=uploaded_extractor,
            helper_uploads=helper_uploads,
            upload_folder_name=data.get("uploadFolderName"),
            custom_parameters=custom_parameters,
        )
    elif selected_model and str(selected_model).startswith("uploaded:"):
        bundle = resolve_uploaded_bundle(
            uploaded_models_dir=UPLOADED_MODELS_DIR,
            selected_model=str(selected_model),
            request_custom_parameters=custom_parameters,
        )
    else:
        bundle = resolve_server_model(
            model_directory=MODEL_DIRECTORY,
            selected_model=selected_model,
        )

    depth = get_sample_depth(DEPTH_FILE, node_id) if node_id else 0

    return PredictionContext(
        node_id=node_id,
        node_ids=[node_id] if node_id else [],
        elapsed_day=elapsed_day,
        selected_protein_region=selected_protein_region,
        selected_protein_regions=_selected_region_map(selected_protein_region),
        depth=depth,
        bundle=bundle,
        custom_parameters=bundle.custom_parameters,
    )


def _load_runtime_components(context: PredictionContext):
    model_wrapper = load_covmutex_model(
        model_path=context.bundle.model_path,
        model_name=context.bundle.model_name,
        description="COVID-19 mutation prediction model",
        source=context.bundle.source,
    )

    if context.bundle.extractor_path and os.path.exists(context.bundle.extractor_path):
        extractor = load_feature_extractor(
            extractor_type="uploaded",
            module_path=context.bundle.extractor_path,
        )
    else:
        extractor = load_feature_extractor(extractor_type="default")

    return model_wrapper, extractor


def _write_prediction_artifact(bundle: BundleResolution, predictions: np.ndarray) -> None:
    if not bundle.bundle_dir:
        return

    output_file = os.path.join(bundle.bundle_dir, "output_predictions.txt")
    with open(output_file, "w", encoding="utf-8") as handle:
        for index, pred in enumerate(predictions):
            if hasattr(pred, "__len__") and len(pred) == 4:
                handle.write(
                    f"Position {index}: "
                    f"A={pred[0]:.6f}, T={pred[1]:.6f}, G={pred[2]:.6f}, C={pred[3]:.6f}\n"
                )
            elif hasattr(pred, "__len__") and len(pred) > 0:
                handle.write(f"{pred[0]:.6f}\n")
            else:
                handle.write(f"{float(pred):.6f}\n")


def run_prediction(context: PredictionContext, request) -> Dict[str, Any]:
    start_time = time.time()
    model_wrapper, extractor = _load_runtime_components(context)
    model_metadata = model_wrapper.metadata()
    extractor_metadata = extractor.get_metadata()

    genome_start = time.time()
    genome_sequence = read_genome_sequence(GENOME_FILE_PATH)
    mutations = parse_mutations(context.node_id) if context.node_id else []
    variant_genome_sequence = construct_variant_genome(genome_sequence, mutations)
    padded_variant_genome_sequence = add_terminal_padding_nucleotide(variant_genome_sequence)
    measure_time("genome_processing", genome_start)

    prediction_params = {
        "cache_path": CACHE_PATH,
        "genome_seq": padded_variant_genome_sequence,
        "mutations": mutations,
        "codon_mapper": CODON_MAPPING_PATH,
        "config_file": configs(),
        "node_ids": context.node_ids,
        "elapsed_day": context.elapsed_day,
        "depth": context.depth,
        "protein_regions": context.selected_protein_regions,
        "selectedModel": context.bundle.selected_model,
        "model_wrapper": model_wrapper,
        "feature_extractor": extractor,
    }
    prediction_params.update(context.custom_parameters)

    predictions = predict_mutations(**prediction_params)
    predictions = np.asarray(predictions)
    if not context.selected_protein_regions and predictions.shape[0] == len(variant_genome_sequence) + 1:
        predictions = predictions[:len(variant_genome_sequence)]

    _write_prediction_artifact(context.bundle, predictions)

    selected_protein_region = (
        tuple(PROTEIN_REGIONS[context.selected_protein_region])
        if context.selected_protein_region in PROTEIN_REGIONS
        else None
    )
    genome_data = calculate_genome_data(
        genome_sequence,
        predictions,
        selected_protein_region=selected_protein_region,
    )

    predictions_offset = (
        PROTEIN_REGIONS[context.selected_protein_region][0]
        if context.selected_protein_region in PROTEIN_REGIONS
        else 0
    )
    protein_mutation_probs = calculate_protein_region_probabilities(
        predictions,
        PROTEIN_REGIONS,
        genome_seq_length=len(genome_sequence),
        predictions_offset=predictions_offset,
    )

    request.session["genome_data"] = genome_data

    response = {
        "nodeId": context.node_id,
        "elapsedDay": context.elapsed_day,
        "selectedModel": context.bundle.selected_model,
        "selectedProteinRegion": context.selected_protein_region,
        "genomeSequence": variant_genome_sequence,
        "genomeData": genome_data,
        "protein_mutation_probs": protein_mutation_probs,
        "proteinRegionPossibilities": PROTEIN_REGIONS,
        "modelType": model_metadata.get("model_type"),
        "model_metadata": model_metadata,
        "extractor_metadata": extractor_metadata,
        "saved_folder": context.bundle.saved_folder,
        "custom_parameters": context.custom_parameters or None,
    }

    measure_time("total_request_handling_internal", start_time)
    return response


@api_view(["GET"])
def get_models(request):
    try:
        return JsonResponse({"available_models": _list_uploaded_bundle_names()}, status=200)
    except Exception as exc:
        return JsonResponse(
            {
                "error": "Failed to retrieve models",
                "details": str(exc),
                "traceback": traceback.format_exc(),
            },
            status=500,
        )


@api_view(["GET"])
def get_model_parameters(request):
    try:
        model_name = request.GET.get("model_name")
        if not model_name:
            return JsonResponse(
                {
                    "error": "model_name parameter is required",
                    "example": "/api/model-parameters/?model_name=my_custom_model",
                },
                status=400,
            )

        if model_name.startswith("uploaded:"):
            model_name = model_name.replace("uploaded:", "")

        model_dir = os.path.join(UPLOADED_MODELS_DIR, model_name)
        params_file = os.path.join(model_dir, "custom_parameters.json")

        if not os.path.exists(model_dir):
            return JsonResponse(
                {
                    "error": f'Model "{model_name}" not found',
                    "available_models": _list_uploaded_bundle_names(),
                },
                status=404,
            )

        if not os.path.exists(params_file):
            return JsonResponse(
                {
                    "model_name": model_name,
                    "hint": "This model was uploaded without custom parameters or the file is missing.",
                    "parameters": {},
                    "status": "success"
                },
                status=404,
            )

        with open(params_file, "r", encoding="utf-8") as handle:
            parameters = json.load(handle)

        return JsonResponse(
            {
                "model_name": model_name,
                "parameters": parameters,
                "status": "success",
            },
            status=200,
        )
    except json.JSONDecodeError as exc:
        return JsonResponse(
            {"error": "Invalid JSON format in custom_parameters.json", "details": str(exc)},
            status=500,
        )
    except Exception as exc:
        return JsonResponse(
            {
                "error": "Failed to retrieve model parameters",
                "details": str(exc),
                "traceback": traceback.format_exc(),
            },
            status=500,
        )


@api_view(["GET", "POST"])
def predict_genome(request):
    if request.method not in ["GET", "POST"]:
        return JsonResponse({"error": "Invalid request method"}, status=400)

    try:
        context = prepare_prediction_context(request)
        return JsonResponse(run_prediction(context, request))
    except ValueError as exc:
        details = str(exc)
        try:
            parsed = json.loads(details)
        except json.JSONDecodeError:
            parsed = None

        if isinstance(parsed, dict) and "scan_results" in parsed:
            return JsonResponse(
                {
                    "error": "Malware detected",
                    "details": parsed.get("message"),
                    "scan_results": parsed.get("scan_results"),
                },
                status=403,
            )
        return JsonResponse({"error": f"An error occurred: {details}"}, status=400)
    except FileNotFoundError as exc:
        return JsonResponse({"error": f"An error occurred: {exc}"}, status=404)
    except Exception as exc:
        traceback.print_exc()
        return JsonResponse({"error": f"An error occurred: {exc}"}, status=500)


def home(request):
    return predict_genome(request)


@api_view(["POST"])
def generate_weblogo(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST method required"}, status=405)

    try:
        data = json.loads(request.body)
        start = int(data.get("start", 1))
        end = int(data.get("end", 25))
        prob_matrix = data.get("probability_matrix", [])
        ref_seq = data.get("reference_sequence", "")
        nuc_order = data.get("nucleotide_order", ["A", "T", "G", "C"])
        confidence_weights = data.get("confidence_weights")

        if not prob_matrix:
            return JsonResponse({"error": "No probability matrix provided"}, status=400)

        prob_array = np.array(prob_matrix)
        if prob_array.ndim != 2 or prob_array.shape[1] != 4:
            return JsonResponse({"error": "Matrix must be Nx4"}, status=400)

        if confidence_weights and len(confidence_weights) == len(prob_matrix):
            weights = np.array(confidence_weights)
            max_weight = weights.max() if weights.max() > 0 else 1.0
            prob_array = prob_array * (weights / max_weight)[:, np.newaxis]

        df_prob = pd.DataFrame(prob_array, columns=nuc_order)

        plt.figure(figsize=(11, 3.5))
        logo = logomaker.Logo(
            df_prob,
            color_scheme="classic",
            font_name="Arial",
            figsize=(11, 3.5),
            stack_order="big_on_top",
        )

        genome_sequence = read_genome_sequence(GENOME_FILE_PATH)
        num_positions = len(prob_matrix)
        positions = range(start, start + num_positions)
        labels = []
        for pos in positions:
            if 0 <= pos < len(genome_sequence):
                labels.append(f"{pos}-{genome_sequence[pos]}")
            else:
                labels.append(str(pos))

        logo.ax.set_xticks(range(len(positions)))
        logo.ax.set_xticklabels(labels)
        logo.ax.xaxis.set_tick_params(rotation=60)
        logo.ax.set_ylabel("Mutation Probability")
        logo.ax.set_xlabel("Position")
        logo.ax.set_title(f"Mutation Profile (Positions {start}-{end})", pad=10)

        if ref_seq and len(ref_seq) == len(prob_matrix):
            for index, nucleotide in enumerate(ref_seq.upper()):
                if nucleotide in nuc_order:
                    logo.highlight_position(p=index, color="lightgray", alpha=0.3)

        response = HttpResponse(content_type="image/png")
        plt.savefig(response, format="png", dpi=150, bbox_inches="tight")
        plt.close()
        return response
    except Exception as exc:
        return JsonResponse({"error": str(exc)}, status=500)
