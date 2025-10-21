import os
import sys
import json
import time
import traceback
import tempfile
import shutil
from collections import defaultdict

from weblogo import *
import numpy as np
import pandas as pd
import base64
import tensorflow as tf
from rich import print
from django.core.cache import cache
from django.http import JsonResponse
from rest_framework.decorators import api_view
from django.http import HttpResponse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logomaker

from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .feature_extractor import parse_mutations, construct_variant_genome
from .configs import configs
from .covmutex_models import load_model as load_covmutex_model, CovMutExKerasModel
from .helpers import (
    measure_time,
    predict_mutations,
    calculate_genome_data,
    calculate_protein_region_probabilities,
    read_genome_sequence
)

# Assuming codon_mapping.json is in the same directory or a subdirectory
codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')

# Specify the cache path relative to the current file's directory
cache_path = os.path.join(os.path.dirname(__file__), 'node_features.h5')

@api_view(["GET", "POST"])
def predict_genome(request):
    if request.method in ['POST', 'GET']:
        return handle_prediction(request)
    return JsonResponse({"error": "Invalid request method"}, status=400)


def home(request):
    return handle_prediction(request)

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

@api_view(['POST'])
def handle_prediction(request):
    try:
        start_time = time.time()
        print("[START] Handling prediction request")
        
        
        # Load input data
        data = request.data if request.method == 'POST' else request.GET
        nodeId = data.get('nodeId')
        elapsedDay = int(data.get('elapsedDay', 0))
        selectedModel = data.get('selectedModel')
        selectedProteinRegion = data.get('selectedProteinRegion')

        print("selected model", selectedModel)
        print("protein region selected", selectedProteinRegion)

        # NEW: Handle uploaded model file
        uploaded_model = request.FILES.get('modelFile')  # Model file from request

        model_directory = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'covid19_models', 'models')

        # NEW: Determine model path
        if uploaded_model:
            # Save uploaded model temporarily
            import tempfile
            temp_dir = tempfile.mkdtemp()
            model_path = os.path.join(temp_dir, uploaded_model.name)
            
            # Write uploaded file to disk
            with open(model_path, 'wb+') as destination:
                for chunk in uploaded_model.chunks():
                    destination.write(chunk)
            
            print(f"Using uploaded model: {uploaded_model.name}")
            
        elif selectedModel:
            # Use existing model from server
            model_path = os.path.join(model_directory, f"{selectedModel}.keras")
            
            # Try alternative extensions if .keras not found
            if not os.path.exists(model_path):
                for ext in ['.h5', '.pb']:
                    alt_path = os.path.join(model_directory, f"{selectedModel}{ext}")
                    if os.path.exists(alt_path):
                        model_path = alt_path
                        break
        else:
            # Default model
            model_path = os.path.join(model_directory, "balanced_data_model.keras")

        # NEW: Validate model file exists
        if not os.path.exists(model_path):
            return JsonResponse({
                "error": f"Model not found: {model_path}"
            }, status=404)
        
        print(f"Loading model from: {model_path}")
        
        # NEW: Load the model using CovMutEx Protocol
        try:
            model_wrapper = load_covmutex_model(
                model_path=model_path,
                model_name=uploaded_model.name if uploaded_model else selectedModel,
                description="COVID-19 mutation prediction model",
                source="uploaded" if uploaded_model else "server"
            )
            
            # Get metadata for logging and response
            model_metadata = model_wrapper.metadata()
            print(f"Model loaded successfully: {model_metadata['name']}")
            print(f"Model type: {model_metadata['model_type']}")
            print(f"Output type: {model_metadata['output_type']}")
            print(f"Source: {model_metadata['source']}")
            
        except Exception as e:
            return JsonResponse({
                "error": f"Failed to load model: {str(e)}"
            }, status=500)
        
        genome_start = time.time()
        # Read genome sequence
        # Define the correct relative path
        relative_path = os.path.join("genome.txt")

        # Get the absolute path of the script's directory
        base_dir = os.path.dirname(os.path.abspath(__file__))

        # Combine the base directory and the relative path
        genome_file_path = os.path.join(base_dir, relative_path)
        genome_sequence = read_genome_sequence(genome_file_path)
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
        
        # Model is already loaded above, just use it
        predictions = predict_mutations(
            cache_path=cache_path, 
            genome_seq=genome_sequence, 
            mutations=mutations, 
            codon_mapper=codon_mapping_path,
            config_file=configs(),
            node_ids=nodeId,
            elapsed_day=elapsedDay,
            protein_regions={selectedProteinRegion: protein_regions[selectedProteinRegion]} if selectedProteinRegion and selectedProteinRegion in protein_regions else {},
            selectedModel=selectedModel,
            model_wrapper=model_wrapper
            )
        
        measure_time("feature_extraction_and_prediction", features_start)
        
        

                

        # Example usage:
        selected_protein_region = tuple(protein_regions[selectedProteinRegion]) if selectedProteinRegion and selectedProteinRegion in protein_regions else None

        genome_data_start = time.time()
        genome_data = calculate_genome_data(genome_sequence, predictions, selected_protein_region=selected_protein_region)
        measure_time("genome_data_calculation", genome_data_start)




            
        protein_probs_start = time.time()
        protein_mutation_probs = calculate_protein_region_probabilities(predictions, protein_regions, genome_seq_length=29904 )
        measure_time("protein_region_probability_calculation", protein_probs_start)
        


        # NEW: Store genome_data in the session
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
            "modelType": model_metadata['model_type'],
            "model_metadata": model_metadata  # Full metadata from Protocol
        }
        
        measure_time("total_request_handling", start_time)
        end_time = time.time()
        print(f"Total runtime: {end_time - start_time} seconds")

        return JsonResponse(response_data)

    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": f"An error occurred: {str(e)}"}, status=500)
    
    # NEW: Clean up temporary files if uploaded model was used    
    finally:
        if 'temp_dir' in dir() and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir, ignore_errors=True)
                print(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as cleanup_error:
                print(f"Warning: Could not clean up temp dir: {cleanup_error}")
