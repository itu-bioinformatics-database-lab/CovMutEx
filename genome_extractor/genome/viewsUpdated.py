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
from .covmutex_models import validate_prediction_payload
from .covmutex_adapters import load_model_adapter
from .feature_extractor_updated import (
    add_terminal_padding_nucleotide,
    get_sample_depth,
)
from .helpers import (
    calculate_genome_data,
    calculate_protein_region_probabilities,
    measure_time,
    predict_mutations,
    read_genome_sequence,
)
from .cache_paths import CACHE_DIR, NODE_FEATURES_CACHE_PATH
from .organism_registry import (
    apply_variant_sequence,
    get_variant_mutations,
    find_variant_subtype,
    list_all_influenza_variants,
    list_builtin_organisms,
    list_variants,
    read_builtin_organism,
    read_custom_organism_from_bundle,
    read_organism_metadata,
    read_variant,
)
from .plugin_runtime import (
    BundleResolution,
    HelperUpload,
    collect_helper_uploads,
    normalize_elapsed_day,
    parse_custom_parameters,
    read_bundle_metadata,
    resolve_server_model,
    resolve_uploaded_bundle,
    save_uploaded_bundle,
)

matplotlib.use("Agg")


class _BundleSavedOnly(Exception):
    """Raised when a model bundle was saved successfully but no prediction
    should run yet.  Used as a fallback when no smoke-test variant is
    available; the view converts this to a 200 "bundle_saved" response so the
    frontend modal closes and refreshes the model list."""

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GENOME_PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))
UPLOADED_MODELS_DIR = os.path.join(GENOME_PROJECT_DIR, "uploaded_models")
MODEL_DIRECTORY = os.path.join(GENOME_PROJECT_DIR, "covid19_models", "models")

CODON_MAPPING_PATH = os.path.join(BASE_DIR, "codon_aa_mapping.json")
CACHE_PATH = NODE_FEATURES_CACHE_PATH
DEPTH_FILE = os.path.join(BASE_DIR, "depth_date.json")
GENOME_FILE_PATH = os.path.join(BASE_DIR, "genome.txt")

# Canonical names for custom-organism uploads. The genome / protein-region
# files are dropped into dedicated form slots and stored under these fixed
# names inside the bundle, so authors can rely on them (see PLUGIN_CONTRACT.md).
CUSTOM_GENOME_FILENAME = "genome.fasta"
CUSTOM_PROTEIN_REGIONS_FILENAME = "protein_regions.csv"

# The legacy module-level COVID PROTEIN_REGIONS dict is gone. Per-prediction
# protein regions now flow from the organism registry (organisms/<name>/) for
# built-in organisms or from the bundle's protein_regions_file helper for
# custom organisms. The weblogo endpoint still reads the COVID genome from
# GENOME_FILE_PATH directly because it's strictly COVID-coordinate today.


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
    # Organism-aware: resolved once during prepare_prediction_context and
    # reused throughout run_prediction. Avoids any hardcoded COVID-only paths
    # inside the predict pipeline.
    genome_sequence: str = ""
    organism_protein_regions: Dict[str, list[int]] = None  # type: ignore[assignment]
    # Set when this prediction is actually an upload-time smoke test against a
    # default variant; predict_genome rewrites the response to a "bundle saved"
    # marker instead of returning predictions for an arbitrary variant.
    smoke_test_variant: Optional[str] = None


def _request_data(request):
    return request.data if request.method == "POST" else request.GET


def _request_files(request):
    return getattr(request, "FILES", getattr(request, "_request", request).FILES)


def _selected_region_map(
    selected_region: Optional[str],
    protein_regions: Dict[str, list[int]],
) -> Optional[Dict[str, list[int]]]:
    """Pick the requested region out of the organism's regions dict.

    Returns ``None`` when the caller did not select a region or the region
    name isn't in the organism's table (so the prediction falls back to the
    whole genome).
    """
    if selected_region and selected_region in protein_regions:
        return {selected_region: protein_regions[selected_region]}
    return None


def _resolve_organism_data(
    bundle: BundleResolution,
) -> tuple[str, Dict[str, list[int]]]:
    """Resolve genome + protein_regions for a bundle's declared organism.

    Four paths:

    * **custom organism** → bundle's own helper files (genome_file,
      protein_regions_file). No variant catalog applies.
    * **generic influenza (organism="influenza")** → variant is required; look
      up its subtype (H1N1 / H3N2 / H5N1) and pull protein_regions from that
      subtype's organism folder.
    * **subtype-specific influenza + variant** → variant FASTA replaces the
      genome string, protein_regions stay the reference subtype's.
    * **built-in (no variant)** → reference genome + reference protein_regions.
    """
    if bundle.organism == "custom":
        return read_custom_organism_from_bundle(
            bundle.bundle_dir or "",
            bundle.genome_file or "",
            bundle.protein_regions_file,
        )

    if bundle.organism == "influenza":
        # Generic influenza bundles defer organism selection to the variant:
        # the variant's subtype tells us which protein_regions.csv to use.
        if not bundle.variant:
            raise ValueError(
                'organism="influenza" bundles require a variant — pick a strain '
                "at prediction time (one of the cataloged H1N1/H3N2/H5N1 entries)."
            )
        subtype_organism = find_variant_subtype(bundle.variant)
        if not subtype_organism:
            raise ValueError(
                f"Unknown influenza variant: {bundle.variant!r}. "
                f"Catalog: {[v['name'] for v in list_all_influenza_variants()]}"
            )
        _, protein_regions = read_builtin_organism(subtype_organism)
        variant_seq = read_variant(subtype_organism, bundle.variant)
        return variant_seq, protein_regions

    reference_seq, protein_regions = read_builtin_organism(bundle.organism)
    if bundle.variant:
        variant_seq = read_variant(bundle.organism, bundle.variant)
        return variant_seq, protein_regions
    return reference_seq, protein_regions


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
        # Organism dispatch field comes straight from the upload form. Helper
        # files (genome.fasta, protein_regions.csv, model_adapter.py) live on
        # disk under canonical names; their presence is detected at resolve
        # time, so the manifest itself only needs the organism marker.
        form_organism = (data.get("organism") or "").strip().lower() or None
        organism_metadata: Optional[Dict[str, Any]] = None
        if form_organism:
            organism_metadata = {"organism": form_organism}

            uploaded_genome = files.get("genomeFile")
            uploaded_protein_regions = files.get("proteinRegionsFile")
            uploaded_adapter = files.get("adapterFile")
            if uploaded_genome is not None:
                helper_uploads.append(
                    HelperUpload(
                        index="genome",
                        file_obj=uploaded_genome,
                        stored_name=CUSTOM_GENOME_FILENAME,
                    )
                )
            if uploaded_protein_regions is not None:
                helper_uploads.append(
                    HelperUpload(
                        index="protein_regions",
                        file_obj=uploaded_protein_regions,
                        stored_name=CUSTOM_PROTEIN_REGIONS_FILENAME,
                    )
                )
            if uploaded_adapter is not None:
                helper_uploads.append(
                    HelperUpload(
                        index="adapter",
                        file_obj=uploaded_adapter,
                        stored_name="model_adapter.py",
                    )
                )

        bundle = save_uploaded_bundle(
            uploaded_models_dir=UPLOADED_MODELS_DIR,
            uploaded_model=uploaded_model,
            uploaded_extractor=uploaded_extractor,
            helper_uploads=helper_uploads,
            upload_folder_name=data.get("uploadFolderName"),
            custom_parameters=custom_parameters,
            organism_metadata=organism_metadata,
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

    # For influenza bundles uploaded without a variant, run a smoke-test
    # prediction against the first H1N1 strain in the catalog so that model
    # loading errors / shape mismatches surface immediately at upload time
    # (mirrors the COVID and Custom upload flow, which run a trial prediction).
    # The response gets rewritten downstream so the user sees a "bundle saved
    # + smoke test passed" marker instead of a real prediction for a variant
    # they didn't pick.
    smoke_test_variant: Optional[str] = None
    if uploaded_model and bundle.organism == "influenza" and not bundle.variant:
        smoke_test_variant = next(
            (v["name"] for v in list_variants("influenza_h1n1")), None
        )
        if smoke_test_variant:
            bundle.variant = smoke_test_variant
        else:
            # Catalog empty for any reason — fall back to plain bundle_saved.
            raise _BundleSavedOnly()

    depth = get_sample_depth(DEPTH_FILE, node_id) if node_id else 0

    # Prediction-time variant override: a request payload's ``variant`` field
    # supersedes any variant baked into bundle_metadata.json. This is what
    # makes the new "pick model → pick strain" UI work for influenza, where
    # the user chooses the strain at predict time, not upload time. COVID and
    # custom organisms ignore the field (COVID uses nodeId, custom carries
    # its own genome).
    request_variant = (data.get("variant") or "").strip() or None
    if request_variant:
        bundle.variant = request_variant

    # Organism-aware genome + protein_regions. Resolved once here so the rest
    # of the predict pipeline sees uniform values regardless of whether the
    # bundle targets a built-in organism or ships its own genome as a helper.
    genome_sequence, organism_protein_regions = _resolve_organism_data(bundle)

    return PredictionContext(
        node_id=node_id,
        node_ids=[node_id] if node_id else [],
        elapsed_day=elapsed_day,
        selected_protein_region=selected_protein_region,
        selected_protein_regions=_selected_region_map(
            selected_protein_region, organism_protein_regions
        ),
        depth=depth,
        bundle=bundle,
        custom_parameters=bundle.custom_parameters,
        genome_sequence=genome_sequence,
        organism_protein_regions=organism_protein_regions,
        smoke_test_variant=smoke_test_variant,
    )


def _load_runtime_components(context: PredictionContext):
    if context.bundle.adapter_path and os.path.exists(context.bundle.adapter_path):
        model_wrapper = load_model_adapter(
            adapter_path=context.bundle.adapter_path,
            model_path=context.bundle.model_path,
            model_name=context.bundle.model_name,
            description="COVID-19 mutation prediction model",
            source=context.bundle.source,
        )
    else:
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
    # Genome was resolved organism-aware in prepare_prediction_context;
    # no hardcoded COVID file read here anymore.
    genome_sequence = context.genome_sequence
    organism_protein_regions = context.organism_protein_regions or {}
    # Variant overlay: delegate to organisms/<name>/variant_handler.py if present.
    # Each organism owns its own logic; organisms with no handler use raw genome.
    bundle_organism = (context.bundle.organism or "covid").lower()
    mutations = get_variant_mutations(bundle_organism, context.node_id)
    variant_genome_sequence = apply_variant_sequence(bundle_organism, genome_sequence, context.node_id)
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

    payload = predict_mutations(**prediction_params)
    # Inject the organism's protein_regions into the payload so the frontend
    # renders the right protein overlay regardless of which organism this
    # bundle targets. setdefault() preserves any annotation the adapter
    # already chose to declare.
    payload.setdefault("annotations", {}).setdefault(
        "protein_regions", organism_protein_regions
    )

    variant_length = len(variant_genome_sequence)
    domain = payload.setdefault("domain", {})
    region = domain.get("region")

    predictions_array = np.asarray(payload["predictions"]["values"])

    if region is None:
        # Whole-genome prediction: values span the entire genome. Drop the +1
        # terminal-padding artifact the extractor may emit before validating.
        if predictions_array.ndim >= 1 and predictions_array.shape[0] == variant_length + 1:
            predictions_array = predictions_array[:variant_length]
            payload["predictions"]["values"] = predictions_array.tolist()

    # total_length is the full genome the prediction actually ran on. With the
    # organism layer this is authoritative for built-in AND custom organisms
    # (variant_length = length of the resolved organism genome), so we set it for
    # both whole-genome and sub-region payloads. This also corrects the +1 that a
    # padded sequence would otherwise leak into an adapter-declared length.
    domain["total_length"] = variant_length

    validate_prediction_payload(payload)

    _write_prediction_artifact(context.bundle, predictions_array)

    selected_protein_region = (
        tuple(organism_protein_regions[context.selected_protein_region])
        if context.selected_protein_region in organism_protein_regions
        else None
    )

    legacy_genome_data = None
    legacy_protein_probs = None
    task_kind = payload.get("task", {}).get("kind")
    if task_kind == "categorical_per_position" and predictions_array.ndim == 2:
        legacy_genome_data = calculate_genome_data(
            genome_sequence,
            predictions_array,
            selected_protein_region=selected_protein_region,
        )
        predictions_offset = (
            organism_protein_regions[context.selected_protein_region][0]
            if context.selected_protein_region in organism_protein_regions
            else 0
        )
        legacy_protein_probs = calculate_protein_region_probabilities(
            predictions_array,
            organism_protein_regions,
            genome_seq_length=len(genome_sequence),
            predictions_offset=predictions_offset,
        )
        request.session["genome_data"] = legacy_genome_data

    response = {
        "nodeId": context.node_id,
        "elapsedDay": context.elapsed_day,
        "selectedModel": context.bundle.selected_model,
        "selectedProteinRegion": context.selected_protein_region,
        "genomeSequence": variant_genome_sequence,
        "predictionPayload": payload,
        "modelType": model_metadata.get("model_type"),
        "model_metadata": model_metadata,
        "extractor_metadata": extractor_metadata,
        "saved_folder": context.bundle.saved_folder,
        "custom_parameters": context.custom_parameters or None,
    }
    if legacy_genome_data is not None:
        response["genomeData"] = legacy_genome_data
        response["protein_mutation_probs"] = legacy_protein_probs
        response["proteinRegionPossibilities"] = organism_protein_regions

    measure_time("total_request_handling_internal", start_time)
    return response


@api_view(["GET"])
def get_organisms(request):
    """Return the upload-form's organism options + the consolidated influenza
    catalog the prediction-time variant picker draws from.

    Shape:
        {
          "organisms": [
            { "name": "covid",     "display_name": "SARS-CoV-2",     ... },
            { "name": "influenza", "display_name": "Influenza A",    ... },
            { "name": "custom",    "display_name": "Other (custom)", ... }
          ],
          "influenza_variants": [
            { name, display_name, accession, year, subtype, is_reference,
              cds_start, cds_end, organism } * 9
          ],
          "subtype_details": {
            "influenza_h1n1": { display_name, genome_length, reference_accession },
            "influenza_h3n2": { ... },
            "influenza_h5n1": { ... }
          }
        }

    The flat ``influenza_variants`` list is what the prediction screen renders
    as the 9-strain dropdown; the per-subtype details surface display info
    (reference accession, etc.) for the model card.
    """
    try:
        upload_options = [
            {
                "name": "covid",
                "display_name": "SARS-CoV-2",
                "reference_accession": read_organism_metadata("covid").get(
                    "reference_accession"
                ),
                "genome_length": read_organism_metadata("covid").get("genome_length"),
            },
            {
                "name": "influenza",
                "display_name": "Influenza A",
                "reference_accession": None,
                "genome_length": None,
                "note": "Variant (HA strain) is picked at prediction time.",
            },
            {
                "name": "custom",
                "display_name": "Other — upload your own genome",
                "reference_accession": None,
                "genome_length": None,
            },
        ]
        subtype_details = {}
        protein_regions_by_subtype = {}
        for subtype_organism in (
            "influenza_h1n1",
            "influenza_h3n2",
            "influenza_h5n1",
        ):
            meta = read_organism_metadata(subtype_organism)
            subtype_details[subtype_organism] = {
                "display_name": meta.get("display_name", subtype_organism),
                "subtype": meta.get("subtype"),
                "reference_accession": meta.get("reference_accession"),
                "genome_length": meta.get("genome_length"),
            }
            # Per-subtype protein regions so the frontend can swap the picker
            # contents when the user changes influenza variant.
            try:
                _, regions = read_builtin_organism(subtype_organism)
                protein_regions_by_subtype[subtype_organism] = regions
            except FileNotFoundError:
                protein_regions_by_subtype[subtype_organism] = {}

        # COVID regions go alongside, so the frontend doesn't need to import
        # them from a static JS file anymore.
        try:
            _, covid_regions = read_builtin_organism("covid")
        except FileNotFoundError:
            covid_regions = {}

        return JsonResponse(
            {
                "organisms": upload_options,
                "influenza_variants": list_all_influenza_variants(),
                "subtype_details": subtype_details,
                "protein_regions_by_subtype": protein_regions_by_subtype,
                "covid_protein_regions": covid_regions,
            },
            status=200,
        )
    except Exception as exc:
        return JsonResponse(
            {"error": "Failed to retrieve organism registry", "details": str(exc)},
            status=500,
        )


def _server_model_names() -> list[str]:
    """List built-in keras models the platform ships."""
    if not os.path.exists(MODEL_DIRECTORY):
        return []
    return sorted(
        os.path.splitext(name)[0]
        for name in os.listdir(MODEL_DIRECTORY)
        if name.endswith((".keras", ".h5", ".hdf5", ".pb", ".pt", ".pth", ".pkl"))
    )


def _uploaded_bundle_organism(bundle_name: str) -> str:
    """Return the organism declared in a bundle's bundle_metadata.json.

    Defaults to "covid" so legacy bundles uploaded before the organism field
    existed still surface a sensible badge on the prediction-screen dropdown.
    """
    bundle_dir = os.path.join(UPLOADED_MODELS_DIR, bundle_name)
    metadata = read_bundle_metadata(bundle_dir)
    organism = (metadata.get("organism") or "covid").strip().lower()
    return organism if organism else "covid"


def _uploaded_bundle_protein_regions(bundle_name: str) -> Optional[Dict[str, list[int]]]:
    """Return protein regions for a custom-organism bundle, or None.

    Only meaningful for ``organism="custom"`` bundles; covid / influenza
    uploads use the built-in catalog so the frontend already knows the
    regions. Returns ``None`` when no protein_regions.csv was shipped, so
    the frontend can fall through to "whole genome only".
    """
    bundle_dir = os.path.join(UPLOADED_MODELS_DIR, bundle_name)
    regions_path = os.path.join(bundle_dir, CUSTOM_PROTEIN_REGIONS_FILENAME)
    if not os.path.exists(regions_path):
        return None
    try:
        from .organism_registry import _read_protein_regions_csv
        return _read_protein_regions_csv(regions_path)
    except Exception:
        return None


@api_view(["GET"])
def get_models(request):
    """Return the full model catalog with organism tags for badge rendering.

    Response:
        {
          "available_models": [<uploaded names>, ...],   # legacy, kept for back-compat
          "models": [
            {"name": ..., "organism": "covid|influenza|custom", "uploaded": bool},
            ...
          ]
        }

    The frontend uses ``models`` to render organism-colored badges on the
    prediction-screen model dropdown, then morphs the variant picker based on
    the selected model's organism.
    """
    try:
        uploaded = _list_uploaded_bundle_names()
        models: list[dict] = []
        for name in _server_model_names():
            # All built-in server models target the COVID reference genome.
            models.append({"name": name, "organism": "covid", "uploaded": False})
        for name in uploaded:
            organism = _uploaded_bundle_organism(name)
            entry = {
                "name": name,
                "organism": organism,
                "uploaded": True,
            }
            # Only custom bundles ship their own protein regions; built-in
            # organism uploads (covid / influenza) defer to the catalog.
            if organism == "custom":
                entry["protein_regions"] = _uploaded_bundle_protein_regions(name)
            models.append(entry)
        return JsonResponse(
            {"available_models": uploaded, "models": models},
            status=200,
        )
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
            # Bundle exists but ships no custom_parameters.json — that's a
            # totally normal state (most uploads don't declare any). Return
            # 200 with an empty parameters dict so the frontend doesn't have
            # to treat "no params" as a network error.
            return JsonResponse(
                {
                    "model_name": model_name,
                    "hint": "This model was uploaded without custom parameters.",
                    "parameters": {},
                    "status": "success",
                },
                status=200,
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
        payload = run_prediction(context, request)
        if context.smoke_test_variant:
            # Influenza upload smoke test passed — return a success marker
            # instead of the raw prediction (the user didn't ask for this
            # variant, only for the bundle to be saved).
            return JsonResponse(
                {
                    "bundle_saved": True,
                    "smoke_test_passed": True,
                    "tested_variant": context.smoke_test_variant,
                    "message": (
                        f"Model uploaded and verified against {context.smoke_test_variant} "
                        f"(H1N1 reference). Pick a variant when running your first prediction."
                    ),
                },
                status=200,
            )
        return JsonResponse(payload)
    except _BundleSavedOnly:
        # Bundle uploaded successfully; variant will be chosen at predict time.
        return JsonResponse(
            {"bundle_saved": True, "message": "Model uploaded. Select a variant when running your first prediction."},
            status=200,
        )
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
