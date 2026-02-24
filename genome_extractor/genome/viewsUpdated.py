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

from .feature_extractor_updated import parse_mutations, construct_variant_genome
from .configs import configs

# NEW: yeni oluşturduğum dosyalardan importlar
from .covmutex_models import load_model as load_covmutex_model
from .covmutex_feature_extractors import load_feature_extractor   
from .helpers import (
    measure_time,
    predict_mutations,
    calculate_genome_data,
    calculate_protein_region_probabilities,
    read_genome_sequence
)
from .security_scanner import scan_folder_for_malware, validate_uploaded_files

# Assuming codon_mapping.json is in the same directory or a subdirectory
codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')

# NEW: buradaki node_features.h5 dosyasının yolu daha farklı
# Specify the cache path relative to the current file's directory
cache_path = os.path.join(os.path.dirname(__file__), 'node_features.h5')

# NEW: Uploaded models base directory constant
UPLOADED_MODELS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), 
    'uploaded_models'
)

@api_view(["GET"])
def get_models():
    """
    GET endpoint to retrieve model names for uploaded models.
    Returns the JSON file containing model names in a list.
    
    Example: GET /api/models/
    """
    try:
        # Check for uploaded models in the UPLOADED_MODELS_DIR
        available_models = []
        if os.path.exists(UPLOADED_MODELS_DIR):
            available_models = [d for d in os.listdir(UPLOADED_MODELS_DIR) 
                                if os.path.isdir(os.path.join(UPLOADED_MODELS_DIR, d))]
            
            return JsonResponse({
                'available_models': available_models
            }, status=200)
        
    except Exception as e:
        return JsonResponse({
            'error': 'Failed to retrieve models',
            'details': str(e),
            'traceback': traceback.format_exc()
        }, status=500)


@api_view(["GET"])
def get_model_parameters(request):
    """
    GET endpoint to retrieve custom_parameters.json for uploaded models.
    Returns the JSON file containing custom parameters with their values.
    
    Query parameters:
    - model_name: Name of the uploaded model folder (required)
    
    Example: GET /api/model-parameters/?model_name=my_custom_model
    """
    try:
        model_name = request.GET.get('model_name')
        
        if not model_name:
            return JsonResponse({
                'error': 'model_name parameter is required',
                'example': '/api/model-parameters/?model_name=my_custom_model'
            }, status=400)
        
        # handle_prediction'da uploaded: prefix'i eklemiştik model_name'e
        if model_name.startswith('uploaded:'):
            model_name = model_name.replace('uploaded:', '')
        
        # Construct path to model's custom_parameters.json file
        model_dir = os.path.join(UPLOADED_MODELS_DIR, model_name)
        params_file = os.path.join(model_dir, 'custom_parameters.json')
        
        # Check if model directory exists
        if not os.path.exists(model_dir):
            available_models = []
            if os.path.exists(UPLOADED_MODELS_DIR):
                available_models = [d for d in os.listdir(UPLOADED_MODELS_DIR) 
                                   if os.path.isdir(os.path.join(UPLOADED_MODELS_DIR, d))]
            
            return JsonResponse({
                'error': f'Model "{model_name}" not found',
                'available_models': available_models
            }, status=404)
        
        # Check if custom_parameters.json exists
        if not os.path.exists(params_file):
            return JsonResponse({
                'error': f'custom_parameters.json not found for model "{model_name}"',
                'model_dir': model_dir,
                'hint': 'This model was uploaded without custom parameters or the file is missing'
            }, status=404)
        
        # Read and return the parameters
        with open(params_file, 'r', encoding='utf-8') as f:
            parameters = json.load(f)
        
        return JsonResponse({
            'model_name': model_name,
            'parameters': parameters,
            'status': 'success'
        }, status=200)
        
    except json.JSONDecodeError as e:
        return JsonResponse({
            'error': 'Invalid JSON format in custom_parameters.json',
            'details': str(e)
        }, status=500)
    except Exception as e:
        return JsonResponse({
            'error': 'Failed to retrieve model parameters',
            'details': str(e),
            'traceback': traceback.format_exc()
        }, status=500)


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

def handle_prediction(request):
    try:
        start_time = time.time()
        print("[START] Handling prediction request")
        
        
        # Load input data
        data = request.data if request.method == 'POST' else request.GET
        nodeId = data.get('nodeId')
        elapsedDay = int(data.get('elapsedDay', 0))
        selectedModel = data.get('selectedModel')
        selectedProteinRegion = data.get('selectedProteinRegion', None)

        print("selected model", selectedModel)
        print("protein region selected", selectedProteinRegion)

        print("[DEBUG] Checking for uploaded files...")
        # NEW: Handle uploaded model file
        # DRF compatibility: access FILES from underlying request if needed
        files = getattr(request, 'FILES', getattr(request, '_request', request).FILES)
        uploaded_model = files.get('modelFile')  # Model file from request
        # NEW: Handle uploaded extractor file
        uploaded_extractor = files.get('extractorFile')
        print(f"[DEBUG] uploaded_model: {uploaded_model}, uploaded_extractor: {uploaded_extractor}")
        # NEW: Handle additional helper files (multiple files with custom names)
        # Format: {'helperFile_0': file, 'helperFile_1': file, ...}
        # and {'helperFileName_0': 'custom_name.ext', 'helperFileName_1': 'another.json', ...}
        # frontend'den helperFile_ ve helperFileName_ prefix'leri ile geliyor
        
        # // FILES (request.FILES):
        # helperFile_0 → <config.json dosyası>
        # helperFile_1 → <weights.h5 dosyası>
        # helperFile_2 → <metadata.yaml dosyası>

        # // DATA (request.data):
        # helperFileName_0 → "my_config.json"
        # helperFileName_1 → "model_weights.h5"
        # helperFileName_2 → "info.yaml
        
        # helper_files = {
        #     '0': <InMemoryUploadedFile: config.json>,
        #     '1': <TemporaryUploadedFile: weights.h5>,
        #     '2': <InMemoryUploadedFile: metadata.yaml>
        # }
        
        # helper_file_names = {
        #     '0': 'my_config.json',
        #     '1': 'model_weights.h5',
        #     '2': 'info.yaml'
        # }
        
        helper_files = {}
        helper_file_names = {}
        print(f"[DEBUG] Checking helper files from FILES keys: {list(files.keys())}")
        for key in files.keys():
            if key.startswith('helperFile_'):
                index = key.replace('helperFile_', '')
                helper_files[index] = request.FILES[key]
                # Get custom name for this file
                name_key = f'helperFileName_{index}'
                # Eğer kullanıcı isim vermemişse otomatik isim üretiliyor, helper_{index}.bin
                helper_file_names[index] = data.get(name_key, f'helper_{index}.bin')
        
        # NEW: Yukarıdaki dosyaları bir klasörde saklamak için isim
        upload_folder_name = data.get('uploadFolderName')  # Kullanıcının verdiği isim
        
        # NEW: Kullanıcının gönderdiği özel parametreler (custom arguments)
        # Frontend'den JSON string olarak geliyor: '{"learning_rate": 0.001, "batch_size": 32}'
        custom_params_str = data.get('customParameters', '{}')
        try:
            custom_parameters = json.loads(custom_params_str) if isinstance(custom_params_str, str) else custom_params_str
        except json.JSONDecodeError:
            custom_parameters = {}
        
        if custom_parameters:
            print(f"[CUSTOM PARAMS] User provided {len(custom_parameters)} custom parameter(s): {list(custom_parameters.keys())}")


        # NEW: Uploaded models base directory (use constant)
        uploaded_models_dir = UPLOADED_MODELS_DIR
        os.makedirs(uploaded_models_dir, exist_ok=True)


        # NEW: Variables for cleanup
        model_path = None
        extractor_path = None
        folder_name = None  # Track which folder was used (for response)
    
        # NEW: bu if bloğu tamamen yeni eklendi
        # SCENARIO 1: New upload of model (and optionally extractor)
        if uploaded_model:
            print(f"[UPLOAD] Model uploaded: {uploaded_model.name}")
            
            # SECURITY: Validate file extensions
            helper_file_names_dict = {helper_file_names.get(idx, helper_files[idx].name) for idx in helper_files.keys()} if helper_files else None
            is_valid, error_msg = validate_uploaded_files(
                model_file=uploaded_model,
                extractor_file=uploaded_extractor,
                helper_files=helper_file_names_dict
            )
            
            if not is_valid:
                return JsonResponse({
                    'error': 'Invalid file extension',
                    'details': error_msg
                }, status=400)
             
            # Create user's folder
            user_folder = os.path.join(uploaded_models_dir, upload_folder_name)
            os.makedirs(user_folder, exist_ok=True)

            print(f"[UPLOAD] Saving to folder: {upload_folder_name}")
                
            # Save model
            model_ext = os.path.splitext(uploaded_model.name)[1]
            model_path = os.path.join(user_folder, f"model{model_ext}")
            
            with open(model_path, 'wb+') as destination:
                for chunk in uploaded_model.chunks():
                    destination.write(chunk)
            
            print(f"[UPLOAD] Model saved: {model_path}")

            # Save extractor if provided
            if uploaded_extractor:
                extractor_path = os.path.join(user_folder, 'feature_extractor.py')
                
                with open(extractor_path, 'wb+') as destination:
                    for chunk in uploaded_extractor.chunks():
                        destination.write(chunk)
                
                print(f"[UPLOAD] Extractor saved: {extractor_path}")
            else:
                extractor_path = None
                print("[UPLOAD] No extractor provided, will use default")
            
            # NEW: Save helper files if provided
            if helper_files:
                print(f"[UPLOAD] Saving {len(helper_files)} helper file(s)")
                for index, helper_file in helper_files.items():
                    custom_name = helper_file_names.get(index, f'helper_{index}.bin')
                    helper_path = os.path.join(user_folder, custom_name)
                    
                    with open(helper_path, 'wb+') as destination:
                        for chunk in helper_file.chunks():
                            destination.write(chunk)
                    
                    print(f"[UPLOAD] Helper file saved: {custom_name} -> {helper_path}")
            else:
                print("[UPLOAD] No helper files provided")
            
            # SECURITY: Scan all uploaded files for malware
            print("[SECURITY] Scanning uploaded files for malware...")
            all_clean, scan_results = scan_folder_for_malware(user_folder)
            
            if not all_clean:
                # Delete entire folder - malware detected
                shutil.rmtree(user_folder, ignore_errors=True)
                
                return JsonResponse({
                    'error': 'Malware detected',
                    'details': 'One or more uploaded files contain malware',
                    'scan_results': scan_results
                }, status=403)
            
            print(f"[SECURITY] ✅ All files clean - no malware detected")
            
            # NEW: Save custom parameters to JSON file
            if custom_parameters:
                params_path = os.path.join(user_folder, 'custom_parameters.json')
                with open(params_path, 'w') as f:
                    json.dump(custom_parameters, f, indent=2)
                print(f"[UPLOAD] Custom parameters saved: {params_path}")
            
            # Set folder_name for response
            folder_name = upload_folder_name

        # SCENARIO 2: Use existing uploaded model by folder name FRONTEND uploaded ismi ile dönmeli!!!!
        elif selectedModel and selectedModel.startswith('uploaded:'):
            folder_name = selectedModel.replace('uploaded:', '')
            user_folder = os.path.join(uploaded_models_dir, folder_name)
            
            print(f"[REUSE] Using uploaded model from: {folder_name}")
            
            # Find model file
            for filename in os.listdir(user_folder):
                if filename.startswith('model.'):
                    model_path = os.path.join(user_folder, filename)
                    break
            
            # Check for extractor
            extractor_path = os.path.join(user_folder, 'feature_extractor.py')
            if not os.path.exists(extractor_path):
                extractor_path = None
            
            # NEW: Load saved custom parameters as defaults
            params_path = os.path.join(user_folder, 'custom_parameters.json')
            saved_parameters = {}
            if os.path.exists(params_path):
                with open(params_path, 'r') as f:
                    saved_parameters = json.load(f)
                print(f"[REUSE] Loaded {len(saved_parameters)} saved parameter(s)")
            
            # NEW: Merge: saved parameters as base, current request can override
            merged_parameters = saved_parameters.copy()
            if custom_parameters:
                merged_parameters.update(custom_parameters)
                print(f"[REUSE] User provided {len(custom_parameters)} new parameter(s), merged with saved")
            
            custom_parameters = merged_parameters
            
            print(f"[REUSE] Model: {model_path}")
            print(f"[REUSE] Extractor: {extractor_path or 'None (will use default)'}")
            print(f"[REUSE] Final parameters: {len(custom_parameters)} parameter(s)")

        # SCENARIO 3: Use server model 
        else:
            model_directory = os.path.join(
                os.path.dirname(os.path.dirname(__file__)), 
                'covid19_models', 
                'models'
            )
            
            if selectedModel:
                model_path = os.path.join(model_directory, f"{selectedModel}.keras")
                
                # Try alternative extensions
                if not os.path.exists(model_path):
                    for ext in ['.h5', '.pb', '.pt', '.pth']:
                        alt_path = os.path.join(model_directory, f"{selectedModel}{ext}")
                        if os.path.exists(alt_path):
                            model_path = alt_path
                            break
            else:
                model_path = os.path.join(model_directory, "balanced_data_model.keras")
            
            if not os.path.exists(model_path):
                return JsonResponse({
                    "error": f"Model not found: {model_path}"
                }, status=404)
            
            extractor_path = None
            print(f"[SERVER] Using server model: {model_path}")


        print(f"Loading model from: {model_path}")
        print(f"[DEBUG] About to load model with load_covmutex_model...")
        
        # NEW: Load the model using CovMutEx Protocol
        try:
            # Determine model name based on source
            if uploaded_model:
                model_name = uploaded_model.name
            elif selectedModel:
                model_name = selectedModel
            else:
                model_name = os.path.basename(model_path)
            
            print(f"[DEBUG] Model name: {model_name}")
            
            model_wrapper = load_covmutex_model(
                model_path=model_path,
                model_name=model_name,
                description="COVID-19 mutation prediction model",
                source="uploaded" if uploaded_model else "server"
            )
            print(f"[DEBUG] Model loaded successfully!")
            
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
        
        # NEW: Load feature extractor
        try:
            if extractor_path and os.path.exists(extractor_path):
                print(f"Loading uploaded extractor: {extractor_path}")
                extractor = load_feature_extractor(
                    extractor_type='uploaded',
                    module_path=extractor_path
                )
            else:
                print("Using default feature extractor")
                extractor = load_feature_extractor(
                    extractor_type='default'
                )
            
            extractor_metadata = extractor.get_metadata()
            print(f"Extractor loaded: {extractor_metadata['name']}")
            
        except Exception as e:
            return JsonResponse({
                "error": f"Failed to load feature extractor: {str(e)}"
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
        
        # NEW: Use Protocol-based model and feature extractor
        # Merge custom parameters with standard parameters
        prediction_params = {
            'cache_path': cache_path,
            'genome_seq': genome_sequence,
            'mutations': mutations,
            'codon_mapper': codon_mapping_path,
            'config_file': configs(),
            'node_ids': nodeId,
            'elapsed_day': elapsedDay,
            'protein_regions': {selectedProteinRegion: protein_regions[selectedProteinRegion]} if selectedProteinRegion and selectedProteinRegion in protein_regions else {},
            'selectedModel': selectedModel,
            'model_wrapper': model_wrapper,
            'feature_extractor': extractor
        }
        
        # NEW: Add custom parameters (user's parameters override if there's a conflict)
        if custom_parameters:
            print(f"[PREDICT] Merging {len(custom_parameters)} custom parameter(s) into prediction")
            prediction_params.update(custom_parameters)
        
        predictions = predict_mutations(**prediction_params)
        
        measure_time("feature_extraction_and_prediction", features_start)
        
        
        print(f"Predictions obtained: {len(predictions)} mutation probabilities")
        print(f"Sample predictions (first 5): {predictions[:5]}")
        
        # NEW: Save predictions to file for verification SILMEYI UNUTMA
        if folder_name:  # Only save for uploaded models
            output_file = os.path.join(uploaded_models_dir, folder_name, 'output_predictions.txt')
            try:
                with open(output_file, 'w') as f:
                    for i, pred in enumerate(predictions):
                        if hasattr(pred, '__len__') and len(pred) > 0:
                            # Check if it's ATGC format (4 values per position)
                            if len(pred) == 4:
                                f.write(f"Position {i}: A={pred[0]:.6f}, T={pred[1]:.6f}, G={pred[2]:.6f}, C={pred[3]:.6f}\n")
                            else:
                                f.write(f"{pred[0]:.6f}\n")
                        else:
                            f.write(f"{pred:.6f}\n")
                print(f"[DEBUG] Predictions saved to: {output_file}")
            except Exception as e:
                print(f"[WARNING] Could not save predictions: {e}")
        # NEW: SILMEYI UNUTMA
                
        # NEW: BURASI POSTPROCESS'E GIRIYOR MU SOR?
        # Example usage:
        selected_protein_region = tuple(protein_regions[selectedProteinRegion]) if selectedProteinRegion and selectedProteinRegion in protein_regions else None

        # NEW: SILMEYI UNUTMA Initialize variables for uploaded models
        genome_data = None
        protein_mutation_probs = None

        # NEW: SILMEYI UNUTMA
        if not uploaded_model:
            genome_data_start = time.time()
            genome_data = calculate_genome_data(genome_sequence, predictions, selected_protein_region=selected_protein_region)
            measure_time("genome_data_calculation", genome_data_start)



        # NEW: SILMEYI UNUTMA
        if not uploaded_model:
            protein_probs_start = time.time()
            protein_mutation_probs = calculate_protein_region_probabilities(predictions, protein_regions, genome_seq_length=29904 )
            measure_time("protein_region_probability_calculation", protein_probs_start)
        


        # NEW: Store genome_data in the session (use _request to access Django's HttpRequest)
        if hasattr(request, '_request'):
            request._request.session["genome_data"] = genome_data if genome_data else [] # NEW: SILMEYI UNUTMA
        response_data = {
            "nodeId": nodeId,
            "elapsedDay": elapsedDay,
            "selectedModel": selectedModel,
            "selectedProteinRegion": selectedProteinRegion,
            "genomeSequence": variant_genome_sequence,
            "genomeData": genome_data if genome_data else [], # NEW: SILMEYI UNUTMA
            "protein_mutation_probs": protein_mutation_probs if protein_mutation_probs else {}, # NEW: SILMEYI UNUTMA
            "proteinRegionPossibilities": protein_regions,
            "modelType": model_metadata['model_type'],
            "model_metadata": model_metadata,  # Full metadata from Protocol
            # NEW: 
            "extractor_metadata": extractor_metadata,  
            "saved_folder": folder_name if uploaded_model else None,
            "custom_parameters": custom_parameters if custom_parameters else None
        }
        
        measure_time("total_request_handling", start_time)
        end_time = time.time()
        print(f"Total runtime: {end_time - start_time} seconds")

        return JsonResponse(response_data)

    except Exception as e:
        import traceback
        print("=" * 80)
        print("ERROR IN handle_prediction:")
        print(traceback.format_exc())
        print("=" * 80)
        return JsonResponse({"error": f"An error occurred: {str(e)}"}, status=500)
    

""" ARTIK TEMİZLEMEK YERİNE KAYDEDİYORUZ
        # NEW: Clean up temporary files if uploaded model was used    
    finally:
        if 'temp_dir' in dir() and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir, ignore_errors=True)
                print(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as cleanup_error:
                print(f"Warning: Could not clean up temp dir: {cleanup_error}")
 """