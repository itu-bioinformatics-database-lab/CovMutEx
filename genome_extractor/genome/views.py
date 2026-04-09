import os
import sys
import json
import time
import traceback
from collections import defaultdict

from weblogo import *
import numpy as np
import base64

import numpy as np
import tensorflow as tf
from rich import print
from django.core.cache import cache
from django.http import JsonResponse
from rest_framework.decorators import api_view
from django.http import HttpResponse

from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .feature_extractor import (
    parse_mutations,
    cache_precomputed_features,
    construct_variant_genome,
    predict_mutations,
    get_sample_depth,
    load_legacy_keras_model,
)
from .priest_annotations import PROTEIN_REGIONS, build_node_priest_annotation
from .configs import configs
from .delta_omicron_retrospective import (
    DEFAULT_ELAPSED_DAY,
    DEFAULT_MODEL_NAME,
    DEFAULT_TOP_K,
    CaseStudyValidationError,
    build_precomputed_variant_mutation_tuples,
    build_known_hotspot_case_study_payload,
    resolve_delta_context,
)

codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')
cache_path = os.path.join(os.path.dirname(__file__), 'node_features.h5')
depth_file = os.path.join(os.path.dirname(__file__), 'depth_date.json')

MODEL_CACHE = {}
SPIKE_ONLY_PREDICTION_REGION = {"S": (PROTEIN_REGIONS["S"][0] - 1, PROTEIN_REGIONS["S"][1] - 1)}

def get_model(path):
    if path not in MODEL_CACHE:
        load_start = time.time()
        print(f"--- Loading model into memory: {path} ---")
        MODEL_CACHE[path] = load_legacy_keras_model(path)
        measure_time("hard_disk_model_load", load_start)
    return MODEL_CACHE[path]

def measure_time(label, start_time):
    elapsed_time = time.time() - start_time
    print(f"[{label}] Elapsed Time: {elapsed_time:.2f} seconds")

def read_genome_sequence(file_path):
    with open(file_path, 'r') as f:
        next(f)
        return ''.join(line.strip() for line in f)


def resolve_model_name_and_path(selected_model):
    model_directory = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), 'covid19_models', 'models'
    )

    if selected_model:
        return selected_model, os.path.join(model_directory, f"{selected_model}.keras")

    fixed_model_path = os.path.join(model_directory, "balanced_data_model_fixed.keras")
    return (
        DEFAULT_MODEL_NAME,
        fixed_model_path
        if os.path.exists(fixed_model_path)
        else os.path.join(model_directory, "balanced_data_model.keras"),
    )


@api_view(["GET", "POST"])
def predict_genome(request):
    if request.method in ['POST', 'GET']:
        return handle_prediction(request)
    return JsonResponse({"error": "Invalid request method"}, status=400)


def home(request):
    return handle_prediction(request)

from django.http import HttpResponse, JsonResponse
from io import BytesIO
import numpy as np
from rest_framework.decorators import api_view

@api_view(['POST'])
def generate_weblogo(request):
    if request.method != 'POST':
        return JsonResponse({"error": "POST method required"}, status=405)
    
    try:
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import logomaker

        data = json.loads(request.body)
        start = int(data.get('start', 1))
        end = int(data.get('end', 25))
        prob_matrix = data.get('probability_matrix', [])
        ref_seq = data.get('reference_sequence', '')
        nuc_order = data.get('nucleotide_order', ['A', 'T', 'G', 'C'])
        confidence_weights = data.get('confidence_weights', None)
        
        if not prob_matrix:
            return JsonResponse({"error": "No probability matrix provided"}, status=400)
            
        prob_array = np.array(prob_matrix)
        if prob_array.ndim != 2 or prob_array.shape[1] != 4:
            return JsonResponse({"error": "Matrix must be Nx4"}, status=400)

        # Scale letter heights by confidence (raw sigmoid sum)
        if confidence_weights and len(confidence_weights) == len(prob_matrix):
            weights = np.array(confidence_weights)
            max_weight = weights.max() if weights.max() > 0 else 1.0
            weights_normalized = weights / max_weight
            prob_array = prob_array * weights_normalized[:, np.newaxis]

        df_prob = pd.DataFrame(prob_array, columns=nuc_order)
        
        plt.figure(figsize=(11, 3.5))
        logo = logomaker.Logo(df_prob,
                            color_scheme='classic',
                            font_name='Arial',
                            figsize=(11, 3.5),
                            stack_order='big_on_top')
        
        base_dir = os.path.dirname(os.path.abspath(__file__))
        genome_file_path = os.path.join(base_dir, "genome.txt")
        genome_sequence = read_genome_sequence(genome_file_path)
        
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
            for i, nuc in enumerate(ref_seq.upper()):
                if nuc in nuc_order:
                    logo.highlight_position(p=i, color='lightgray', alpha=0.3)
        
        response = HttpResponse(content_type='image/png')
        plt.savefig(response, format='png', dpi=150, bbox_inches='tight')
        plt.close()
        return response

    except ImportError as e:
        return JsonResponse(
            {"error": f"Missing optional weblogo dependency: {e}"},
            status=500,
        )
    except Exception as e:
        return JsonResponse({"error": str(e)}, status=500)

def handle_prediction(request):
    try:
        start_time = time.time()
        print("[START] Handling prediction request")

        data = request.data if request.method == 'POST' else request.GET
        nodeId = data.get('nodeId')
        elapsedDay = int(data.get('elapsedDay', 0))
        selectedModel = data.get('selectedModel')
        selectedProteinRegion = data.get('selectedProteinRegion')

        print(f"nodeId={nodeId}, elapsedDay={elapsedDay}, model={selectedModel}, region={selectedProteinRegion}")

        # ----------------------------------------------------------------
        # Get tree depth for this node from depth_date.json
        # ----------------------------------------------------------------
        depth = get_sample_depth(depth_file, nodeId) if nodeId else 0
        print(f"Tree depth for node: {depth}")

        genome_start = time.time()
        base_dir = os.path.dirname(os.path.abspath(__file__))
        genome_file_path = os.path.join(base_dir, "genome.txt")
        genome_sequence = read_genome_sequence(genome_file_path)
        genome_length = len(genome_sequence)
        mutations = parse_mutations(nodeId)
        variant_genome_sequence = construct_variant_genome(genome_sequence, mutations)
        measure_time("genome_processing", genome_start)

        protein_regions = {
            gene: [start, end]
            for gene, (start, end) in PROTEIN_REGIONS.items()
        }

        priest_annotation = build_node_priest_annotation(
            node_id=nodeId,
            mutation_tuples=mutations,
            reference_genome_sequence=genome_sequence,
            variant_genome_sequence=variant_genome_sequence,
        )

        priest_only = (selectedModel or "").upper() == "PRIEST"
        if priest_only:
            response_data = {
                "nodeId": nodeId,
                "elapsedDay": elapsedDay,
                "selectedModel": "PRIEST",
                "selectedProteinRegion": None,
                "genomeSequence": variant_genome_sequence,
                "genomeData": [],
                "protein_mutation_probs": {},
                "proteinRegionPossibilities": protein_regions,
                "modelType": "priest",
                **priest_annotation,
                "model_metadata": None,
            }
            measure_time("total_request_handling_internal", start_time)
            return JsonResponse(response_data)

        features_start = time.time()
        selectedModel, model_path = resolve_model_name_and_path(selectedModel)

        model_load_start = time.time()
        model = get_model(model_path)
        measure_time("model_access_time", model_load_start)
        print("Model output shape:", model.output_shape)

        predictions = predict_mutations(
            cache_path=cache_path,
            genome_seq=genome_sequence,
            mutations=mutations,
            codon_mapper=codon_mapping_path,
            config_file=configs(),
            node_ids=[nodeId],
            elapsed_day=elapsedDay,
            depth=depth,                    # ← FIXED: pass actual tree depth
            protein_regions={selectedProteinRegion: protein_regions[selectedProteinRegion]}
                            if selectedProteinRegion and selectedProteinRegion in protein_regions
                            else None,
            model=model
        )
        measure_time("feature_extraction_and_prediction", features_start)

        def calculate_genome_data(genome_seq, position_predictions, selected_protein_region=None):
            """
            Convert position predictions to genome_data format.
            position_predictions shape: (num_positions, 4) — columns are [A, T, G, C]

            IMPORTANT: These are independent sigmoid scores (not a softmax distribution),
            so we do NOT normalise them. Normalisation destroys the signal and makes
            all positions look uniform.
            """
            print(f"Position predictions shape: {position_predictions.shape}")
            print(f"Raw value range: min={position_predictions.min():.4f}, "
                  f"max={position_predictions.max():.4f}, "
                  f"mean={position_predictions.mean():.4f}")

            # No normalisation — pass raw sigmoid probabilities directly
            genome_data = position_predictions.T.tolist()  # shape: (4, num_positions)

            print(f"Genome data: {len(genome_data)} nucleotides × {len(genome_data[0])} positions")
            for i, nuc in enumerate(['A', 'T', 'G', 'C']):
                print(f"  {nuc} first 3: {genome_data[i][:3]}")

            return genome_data

        selected_protein_region_tuple = (
            tuple(protein_regions[selectedProteinRegion])
            if selectedProteinRegion and selectedProteinRegion in protein_regions
            else None
        )

        genome_data_start = time.time()
        genome_data = calculate_genome_data(
            genome_sequence, predictions,
            selected_protein_region=selected_protein_region_tuple
        )
        measure_time("genome_data_calculation", genome_data_start)


        def calculate_protein_region_probabilities(position_predictions, protein_regions, genome_seq_length):
            protein_mutation_probs = {}
            avg_mutation_per_position = np.mean(position_predictions, axis=1)

            for protein, (protein_start, protein_end) in protein_regions.items():
                protein_start = max(0, protein_start)
                protein_end = min(genome_seq_length, protein_end)

                if protein_start < len(avg_mutation_per_position) and protein_end <= len(avg_mutation_per_position):
                    relevant = avg_mutation_per_position[protein_start:protein_end]
                    avg_prob = float(np.mean(relevant))
                else:
                    avg_prob = 0.0

                protein_mutation_probs[protein] = avg_prob

            return protein_mutation_probs

        protein_probs_start = time.time()
        protein_mutation_probs = calculate_protein_region_probabilities(
            predictions, protein_regions, genome_seq_length=29904
        )
        measure_time("protein_region_probability_calculation", protein_probs_start)

        session_save_start = time.time()
        request.session["genome_data"] = genome_data
        measure_time("django_session_write_time", session_save_start)

        json_construction_start = time.time()
        response_data = {
            "nodeId": nodeId,
            "elapsedDay": elapsedDay,
            "selectedModel": selectedModel,
            "selectedProteinRegion": selectedProteinRegion,
            "genomeSequence": variant_genome_sequence,
            "genomeData": genome_data,
            "protein_mutation_probs": protein_mutation_probs,
            "proteinRegionPossibilities": protein_regions,
            "modelType": "multi-input" if 'multi' in selectedModel.lower() else "single-input",
            **priest_annotation,
            "model_metadata": {
                "output_shape": str(model.output_shape),
                "input_shape": str(model.input_shape),
                "model_type": (
                    "single-output" if model.output_shape[1] == 1 else
                    "dual-output" if model.output_shape[1] == 2 else
                    "multi-class"
                )
            }
        }
        measure_time("json_data_serialization_prep", json_construction_start)
        measure_time("total_request_handling_internal", start_time)
        return JsonResponse(response_data)

    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": f"An error occurred: {str(e)}"}, status=500)


@api_view(["GET", "POST"])
def known_hotspot_case_study(request):
    try:
        data = request.data if request.method == "POST" else request.GET
        selected_model = data.get("selectedModel")
        requested_node_id = data.get("nodeId")
        elapsed_day = int(data.get("elapsedDay", 0))
        top_k = int(data.get("topK", DEFAULT_TOP_K))
        node_context = resolve_delta_context(requested_node_id)
        # If caller did not supply elapsed_day (or sent 0), derive it from the
        # variant's emergence date: elapsed = node_date − 2019-12-01.
        if not elapsed_day:
            elapsed_day = int(node_context.get("elapsed_days", 0))
        node_id = node_context["node_id"]
        selected_model_upper = (selected_model or "").upper()
        spike_predictions = None
        variant_genome_sequence = None
        genome_sequence = None

        if selected_model_upper != "PRIEST":
            selected_model, model_path = resolve_model_name_and_path(selected_model)
            base_dir = os.path.dirname(os.path.abspath(__file__))
            genome_file_path = os.path.join(base_dir, "genome.txt")
            genome_sequence = read_genome_sequence(genome_file_path)
            mutation_build = build_precomputed_variant_mutation_tuples(
                reference_genome_sequence=genome_sequence,
                node_context=node_context,
            )
            mutations = mutation_build["mutation_tuples"]
            variant_genome_sequence = construct_variant_genome(genome_sequence, mutations)
            depth = int(node_context.get("depth") or 0)

            model = get_model(model_path)
            spike_predictions = predict_mutations(
                cache_path=cache_path,
                genome_seq=variant_genome_sequence,
                mutations=mutations,
                codon_mapper=codon_mapping_path,
                config_file=configs(),
                node_ids=[node_id],
                elapsed_day=elapsed_day,
                depth=depth,
                protein_regions=SPIKE_ONLY_PREDICTION_REGION,
                model=model,
            )
            print(
                f"[case-study] model={selected_model} node={node_id} "
                f"predictions shape={spike_predictions.shape} "
                f"min={spike_predictions.min():.6f} max={spike_predictions.max():.6f} "
                f"mean={spike_predictions.mean():.6f} std={spike_predictions.std():.6f} "
                f"mutations={len(mutations)} depth={depth} "
                f"source={mutation_build['source_file']} "
                f"threshold={mutation_build['consensus_threshold']}"
            )
        else:
            selected_model = "PRIEST"
            mutation_build = None

        response_data = build_known_hotspot_case_study_payload(
            predictions=spike_predictions,
            reference_genome_sequence=genome_sequence,
            node_context=node_context,
            selected_model=selected_model,
            elapsed_day=elapsed_day,
            top_k=top_k,
        )
        if variant_genome_sequence is not None:
            response_data["metadata"]["scoring_context"]["variant_genome_length"] = len(
                variant_genome_sequence
            )
        if mutation_build is not None:
            response_data["metadata"]["scoring_context"]["precomputed_variant_summary"] = {
                "mutation_count": mutation_build["mutation_count"],
                "deletion_count": mutation_build["deletion_count"],
                "consensus_threshold": mutation_build["consensus_threshold"],
                "source_file": mutation_build["source_file"],
                "reference_mismatch_count": mutation_build["reference_mismatch_count"],
                "total_rows": mutation_build["total_rows"],
                "invalid_rows": mutation_build["invalid_rows"],
                "below_threshold_count": mutation_build["below_threshold_count"],
                "selected_support_summary": mutation_build["selected_support_summary"],
                "spike_site_count": mutation_build["spike_site_count"],
                "spike_mutation_count": mutation_build["spike_mutation_count"],
                "spike_support_summary": mutation_build["spike_support_summary"],
            }
        return JsonResponse(response_data)
    except CaseStudyValidationError as exc:
        return JsonResponse({"error": str(exc)}, status=400)
    except Exception as exc:
        traceback.print_exc()
        return JsonResponse({"error": f"An error occurred: {str(exc)}"}, status=500)
