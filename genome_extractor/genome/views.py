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

from .feature_extractor import parse_mutations, cache_precomputed_features, construct_variant_genome, predict_mutations
from .configs import configs

# Assuming codon_mapping.json is in the same directory or a subdirectory
codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')

# Specify the cache path relative to the current file's directory
cache_path = os.path.join(os.path.dirname(__file__), 'node_features.h5')


def measure_time(label, start_time):
    """Helper function to log elapsed time."""
    elapsed_time = time.time() - start_time
    print(f"[{label}] Elapsed Time: {elapsed_time:.2f} seconds")  

def load_keras_model(model_path):
    """Load a Keras model from the specified path"""
    try:
        model = load_model(model_path, compile=False)
        return model
    except Exception as e:
        print(f"Error loading model: {e}")
        traceback.print_exc()
        return None



def read_genome_sequence(file_path):
    """Read and process genome sequence from file"""
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        return ''.join(line.strip() for line in f)


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
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Set the backend to Agg
import matplotlib.pyplot as plt
import logomaker
from rest_framework.decorators import api_view

@api_view(['POST'])
def generate_weblogo(request):
    if request.method != 'POST':
        return JsonResponse({"error": "POST method required"}, status=405)
    
    try:
        data = json.loads(request.body)
        start = int(data.get('start', 1))
        end = int(data.get('end', 25))
        prob_matrix = data.get('probability_matrix', [])
        ref_seq = data.get('reference_sequence', '')
        nuc_order = data.get('nucleotide_order', ['A', 'T', 'G', 'C'])
        
        if not prob_matrix:
            return JsonResponse({"error": "No probability matrix provided"}, status=400)
            
        # Validate and convert to numpy array
        prob_array = np.array(prob_matrix)
        if prob_array.ndim != 2 or prob_array.shape[1] != 4:
            return JsonResponse({"error": "Matrix must be Nx4"}, status=400)
        

        
        # Create dataframe with specified nucleotide order
        df_prob = pd.DataFrame(prob_array, columns=nuc_order)
        
        # Generate logo
        plt.figure(figsize=(11, 3.5))  
        logo = logomaker.Logo(df_prob,
                            color_scheme='classic',
                            font_name='Arial',
                            figsize=(11, 3.5),
                            stack_order='big_on_top')
        
        # Style adjustments with real positions
        #positions = range(start, end)
        relative_path = os.path.join("genome.txt")

        # Get the absolute path of the script's directory
        base_dir = os.path.dirname(os.path.abspath(__file__))

        # Combine the base directory and the relative path
        genome_file_path = os.path.join(base_dir, relative_path)
        genome_sequence = read_genome_sequence(genome_file_path)
        
        num_positions = len(prob_matrix) 
        positions = range(start, start+num_positions)  
        print("Positions:", list(positions)) 
        labels = []
        for pos in positions:
            genome_pos = pos
            if 0 <= genome_pos < len(genome_sequence):
                ref_nuc = genome_sequence[genome_pos]
                labels.append(f"{pos}-{ref_nuc}")
            else:
                labels.append(str(pos))

        # Key adjustment: Shift ticks to align with 1-based labels
        logo.ax.set_xticks(range(len(positions)))  
        logo.ax.set_xticklabels(labels)
        #logo.ax.set_xticks(range(len(positions)))
        #logo.ax.set_xticklabels(labels)
        logo.ax.xaxis.set_tick_params(rotation=60)
        logo.ax.set_ylabel("Mutation Probability")
        logo.ax.set_xlabel("Position")
        logo.ax.set_title(f"Mutation Profile (Positions {start}-{end})", pad=10)
        
        # Highlight reference positions if sequence provided
        if ref_seq and len(ref_seq) == len(prob_matrix):
            for i, nuc in enumerate(ref_seq.upper()):
                if nuc in nuc_order:
                    logo.highlight_position(p=i, color='lightgray', alpha=0.3)
        
        # Save to response
        response = HttpResponse(content_type='image/png')
        plt.savefig(response, format='png', dpi=150, bbox_inches='tight')
        plt.close()

        return response  
    except Exception as e:
        return JsonResponse({"error": str(e)}, status=500)


def handle_prediction(request):
    try:
        start_time = time.time()
        print("[START] Handling prediction request")
        
        
        # Load input data
        data = request.data if request.method == 'POST' else request.GET
        nodeId = data.get('nodeId')
        elapsedDay = int(data.get('elapsedDay', 0))
        selectedModel = data.get('selectedModel')
        print("selected model", selectedModel)
        selectedProteinRegion = data.get('selectedProteinRegion')

        print("protein region selected", selectedProteinRegion)
        
        genome_start = time.time()
        # Read genome sequence
        # Define the correct relative path
        relative_path = os.path.join("genome.txt")

        # Get the absolute path of the script's directory
        base_dir = os.path.dirname(os.path.abspath(__file__))

        # Combine the base directory and the relative path
        genome_file_path = os.path.join(base_dir, relative_path)
        genome_sequence = read_genome_sequence(genome_file_path)
        genome_length = len(genome_sequence)
        mutations = parse_mutations(nodeId)
        variant_genome_sequence = construct_variant_genome(genome_sequence, mutations)
        measure_time("genome_processing", genome_start)
        # Extract features
        print("Extracting features...")
        protein_regions = {
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

        features_start = time.time()
        # model_path = os.path.join(r"F:\COVID19 MUTATION\genome_extractor\covid19_models", f"{selectedModel}.h5")
        model_directory = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'covid19_models','models')
  
        
        # Select the model based on your request
        if selectedModel:
            model_path = os.path.join(model_directory, f"{selectedModel}.keras")
        else:
            model_path = os.path.join(model_directory, "balanced_data_model.keras")
            selectedModel = "balanced_data_model"
        print("model", selectedModel)
        
        model = tf.keras.models.load_model(model_path)
        print("Model output shape:", model.output_shape) 
              
        predictions = predict_mutations(
        cache_path=cache_path,
        genome_seq=genome_sequence,
        mutations=mutations,
        codon_mapper=codon_mapping_path,
        config_file=configs(),
        node_ids=[nodeId],
        elapsed_day=elapsedDay,
        protein_regions={selectedProteinRegion: protein_regions[selectedProteinRegion]} 
                        if selectedProteinRegion and selectedProteinRegion in protein_regions else None,
        model=model
    )
        measure_time("feature_extraction_and_prediction", features_start)
        
        
        def calculate_genome_data(genome_seq, position_predictions, selected_protein_region=None):
            """
            Convert position predictions directly to genome_data format with normalization.
            position_predictions shape: (num_positions, 4) where columns are [A, T, G, C]
            Returns: List of 4 arrays, one per nucleotide across all positions
            """
            print(f"Position predictions shape: {position_predictions.shape}")
            
            if selected_protein_region:
                protein_start, protein_end = selected_protein_region
                protein_start = max(0, protein_start)
                protein_end = min(len(genome_seq), protein_end)
                
                # Predictions are already for the protein region only
                relevant_predictions = position_predictions
            else:
                # Use all predictions
                relevant_predictions = position_predictions
            
            # *** NORMALIZATION: Ensure each position sums to 1.0 ***
            row_sums = relevant_predictions.sum(axis=1, keepdims=True)
            # Avoid division by zero - replace zero sums with 1
            row_sums = np.where(row_sums == 0, 1.0, row_sums)
            normalized_predictions = relevant_predictions / row_sums
            
            # Transpose to get (4, num_positions) - one array per nucleotide
            genome_data = normalized_predictions.T.tolist()
            
            print(f"Genome data shape: {len(genome_data)} x {len(genome_data[0]) if genome_data else 0}")
            print(f"Sample data (first 3 positions):")
            for i, nuc in enumerate(['A', 'T', 'G', 'C']):
                print(f"  {nuc}: {genome_data[i][:3]}")
            
            # Verify normalization on first position
            if genome_data and genome_data[0]:
                sample_sum = sum(genome_data[i][0] for i in range(4))
                print(f"✓ Normalization check - Position 0 sum: {sample_sum:.6f} (should be 1.0)")
            
            return genome_data
        # Example usage:
        selected_protein_region = tuple(protein_regions[selectedProteinRegion]) if selectedProteinRegion and selectedProteinRegion in protein_regions else None

        genome_data_start = time.time()
        genome_data = calculate_genome_data(genome_sequence, predictions, selected_protein_region=selected_protein_region)
        measure_time("genome_data_calculation", genome_data_start)

        def calculate_protein_region_probabilities(position_predictions, protein_regions, genome_seq_length):
            """
            Calculate average mutation probability per protein region.
            position_predictions shape: (num_positions, 4) for [A, T, G, C]
            """
            protein_mutation_probs = {}
            
            # Calculate overall mutation probability (average across all 4 nucleotides per position)
            avg_mutation_per_position = np.mean(position_predictions, axis=1)
            
            for protein, (protein_start, protein_end) in protein_regions.items():
                protein_start = max(0, protein_start)
                protein_end = min(genome_seq_length, protein_end)
                
                if protein_start < len(avg_mutation_per_position) and protein_end <= len(avg_mutation_per_position):
                    relevant_predictions = avg_mutation_per_position[protein_start:protein_end]
                    avg_prob = float(np.mean(relevant_predictions))
                else:
                    avg_prob = 0.0
                
                protein_mutation_probs[protein] = avg_prob
            
            return protein_mutation_probs 
        protein_probs_start = time.time()
        protein_mutation_probs = calculate_protein_region_probabilities(predictions, protein_regions, genome_seq_length=29904 )
        measure_time("protein_region_probability_calculation", protein_probs_start)
        
        # Store genome_data in the session
        request.session["genome_data"] = genome_data
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
            "model_metadata": {
                "output_shape": str(model.output_shape),
                "input_shape": str(model.input_shape),
                "model_type": "single-output" if model.output_shape[1] == 1 else 
                             "dual-output" if model.output_shape[1] == 2 else 
                             "multi-class"
            }
        }
        
        measure_time("total_request_handling", start_time)
        end_time = time.time()
        print(f"Total runtime: {end_time - start_time} seconds")

        return JsonResponse(response_data)

    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": f"An error occurred: {str(e)}"}, status=500)
