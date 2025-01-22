import os
import sys
import json
import time
import traceback
from collections import defaultdict

import numpy as np
import tensorflow as tf
from rich import print
from django.core.cache import cache
from django.http import JsonResponse
from rest_framework.decorators import api_view
from django.http import HttpResponse
from keras.models import load_model
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .feature_extractor import parse_mutations, cache_precomputed_features, construct_variant_genome
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
        print("model", selectedModel)
        
        model = tf.keras.models.load_model(model_path)


        def predict_mutations(cache_path, genome_seq, mutations, codon_mapper, config_file, node_ids, elapsed_day=None, protein_regions=None):
            """
            Predicts mutation probabilities for precomputed features.
            Input features shape: (29904, 205) containing all possible transitions
            Returns matrix of shape (29904, 4) where each row contains probabilities for A,T,G,C transitions
            """
            
            # Get precomputed features - shape (29904, 205)
            print("Selected model is ", selectedModel)
            process_time = time.time()
            print("Using model: ",selectedModel )
            features = cache_precomputed_features(
                cache_path=cache_path,
                genome_seq=genome_seq,
                mutations=mutations,
                codon_mapper=codon_mapper,
                config_file=config_file,
                node_ids=node_ids,
                elapsed_day=elapsed_day,
                protein_regions=protein_regions
            )
            process_time_end = time.time()
            print("Shape of feature received: ", features.shape)
            print("Length of Features :", len(features))

            features = features.squeeze()
            # scaler = MinMaxScaler()
            print("Shape of feature before model input: ", features.shape)
            predict_time = time.time()
            if 'multi' in selectedModel.lower():
                predictions = model.predict([features for _ in range(10)], verbose=0)
            else:
                predictions = model.predict(features, verbose=0)

            
            # # Fit and transform the data
            # scaled_data = scaler.fit_transform(features)
            
            # Batch predict all transitions at once
            
            
            predict_time_end = time.time()
            print("Time spent on prediction", predict_time_end - predict_time)
            print("shape of prediction: ", predictions.shape)

            
            print("Time spent processing variant", process_time_end - process_time)
            return predictions

        
        predictions = predict_mutations(cache_path=cache_path, genome_seq=genome_sequence, 
                                mutations= mutations, codon_mapper=codon_mapping_path, 
                                config_file=configs(), 
                                node_ids=nodeId , elapsed_day= elapsedDay,
                                 protein_regions = {selectedProteinRegion: protein_regions[selectedProteinRegion]} if selectedProteinRegion and selectedProteinRegion in protein_regions else None

)
        measure_time("feature_extraction_and_prediction", features_start)
        
        

                
        def calculate_genome_data(genome_seq, position_predictions, selected_protein_region=None):
            """
            Calculate genome data using raw probabilities.
            Returns the probability for each nucleotide (A, T, G, C) at each position.
            Assumes position_predictions is either of shape (length, 1) or (length, 2).
            If a selected_protein_region is provided, only positions within that range are used.
            """
            genome_data = []
            nucleotides = "ATGC"

            # Determine the range of positions to use (either protein region or full genome)
            if selected_protein_region:
                print("position predictions", position_predictions[:2])
                # Extract the start and end positions of the protein region
                protein_start, protein_end = selected_protein_region

                # Ensure the region is within the genome length
                protein_start = max(0, protein_start)  # Clamp start to be at least 0
                protein_end = min(len(genome_seq), protein_end)  # Clamp end to genome length

                # Calculate the corresponding indices in the position_predictions array
                # The selected protein region is in terms of the genome, but position_predictions is relative to that region's size.
                relevant_positions = range(protein_start, protein_end)
                
                # Map the protein region indices to the correct indices for position_predictions
                prediction_start_idx = protein_start - protein_start
                prediction_end_idx = protein_end - protein_start

                # Ensure position_predictions is a numpy array for easy indexing
                relevant_predictions = position_predictions[prediction_start_idx:prediction_end_idx]
                print("Relevant positions", relevant_positions)
                print("Relevant predictions", relevant_predictions[:2])

                # Ensure relevant_predictions is not empty
                if relevant_predictions.size == 0:
                    raise ValueError("Relevant predictions are empty. Cannot calculate genome data.")

                # Iterate through each relevant position in the genome sequence
                for idx, position in enumerate(relevant_positions):
                    # Wrap-around mapping to ensure position stays within bounds
                    valid_pos = position % len(relevant_predictions)

                    if valid_pos >= len(relevant_predictions):
                        raise IndexError(f"Index {valid_pos} out of bounds for relevant_predictions array of size {len(relevant_predictions)}")

                    # Get raw mutation probabilities for this position
                    mutation_probs = relevant_predictions[valid_pos]  # Array with shape (1,) or (2,)

                    current_nucleotide = genome_seq[position]
                    position_probs = {nucleotide: 0.0 for nucleotide in nucleotides}  # Initialize all nucleotides with 0 probability

                    # Handle different lengths of mutation probabilities for each position
                    if mutation_probs.shape[0] == 1:
                        # For shape (length, 1): Assume a single mutation probability for each position
                        mutation_prob = mutation_probs[0]

                        # Create position probabilities for all 4 nucleotides
                        for nucleotide in nucleotides:
                            if nucleotide == current_nucleotide:
                                # Probability of the current nucleotide remaining unchanged
                                position_probs[nucleotide] = 1 - mutation_prob
                            else:
                                # Probability of mutating to this specific nucleotide
                                position_probs[nucleotide] = mutation_prob

                    elif mutation_probs.shape[0] == 2:
                        # For shape (length, 2): Two output probabilities for each position
                        prob_no_mutation = mutation_probs[0]
                        prob_mutation = mutation_probs[1]

                        # Create position probabilities for all 4 nucleotides
                        for nucleotide in nucleotides:
                            if nucleotide == current_nucleotide:
                                # Probability of the current nucleotide remaining unchanged
                                position_probs[nucleotide] = prob_no_mutation
                            else:
                                # Probability of mutating to this specific nucleotide
                                position_probs[nucleotide] = prob_mutation

                    else:
                        # If there are more than 2 probabilities, additional handling is required
                        raise ValueError(f"Unexpected prediction shape for position {position}: {mutation_probs.shape[0]} probabilities")

                    # Append the probabilities for this position
                    genome_data.append([position_probs[nucleotide] for nucleotide in nucleotides])

            else:
                # If no protein region is selected, use the full genome
                relevant_positions = range(len(genome_seq))
                relevant_predictions = position_predictions

                # Ensure relevant_predictions is a numpy array for easy indexing
                relevant_predictions = np.array(relevant_predictions)

                if relevant_predictions.size == 0:
                    raise ValueError("Relevant predictions are empty. Cannot calculate genome data.")

                # Iterate through each relevant position in the genome sequence
                for idx, position in enumerate(relevant_positions):
                    if idx >= len(relevant_predictions):
                        raise IndexError(f"Index {idx} out of bounds for relevant_predictions array of size {len(relevant_predictions)}")

                    # Get raw mutation probabilities for this position
                    mutation_probs = relevant_predictions[idx]  # Array with shape (1,) or (2,)

                    current_nucleotide = genome_seq[position]
                    position_probs = {nucleotide: 0.0 for nucleotide in nucleotides}  # Initialize all nucleotides with 0 probability

                    # Handle different lengths of mutation probabilities for each position
                    if mutation_probs.shape[0] == 1:
                        # For shape (length, 1): Assume a single mutation probability for each position
                        mutation_prob = mutation_probs[0]

                        # Create position probabilities for all 4 nucleotides
                        for nucleotide in nucleotides:
                            if nucleotide == current_nucleotide:
                                # Probability of the current nucleotide remaining unchanged
                                position_probs[nucleotide] = 1 - mutation_prob
                            else:
                                # Probability of mutating to this specific nucleotide
                                position_probs[nucleotide] = mutation_prob

                    elif mutation_probs.shape[0] == 2:
                        # For shape (length, 2): Assuming two output probabilities for each position
                        prob_no_mutation = mutation_probs[0]
                        prob_mutation = mutation_probs[1]

                        # Create position probabilities for all 4 nucleotides
                        for nucleotide in nucleotides:
                            if nucleotide == current_nucleotide:
                                # Probability of the current nucleotide remaining unchanged
                                position_probs[nucleotide] = prob_no_mutation
                            else:
                                # Probability of mutating to this specific nucleotide
                                position_probs[nucleotide] = prob_mutation

                    else:
                        # If there are more than 2 probabilities, additional handling is required
                        raise ValueError(f"Unexpected prediction shape for position {position}: {mutation_probs.shape[0]} probabilities")

                    # Append the probabilities for this position
                    genome_data.append([position_probs[nucleotide] for nucleotide in nucleotides])

            # Convert to numpy array and transpose to get [P(A), P(T), P(G), P(C)] for each position
            return np.array(genome_data).T.tolist()

        # Example usage:
        selected_protein_region = tuple(protein_regions[selectedProteinRegion]) if selectedProteinRegion and selectedProteinRegion in protein_regions else None

        genome_data_start = time.time()
        genome_data = calculate_genome_data(genome_sequence, predictions, selected_protein_region=selected_protein_region)
        measure_time("genome_data_calculation", genome_data_start)






        def calculate_protein_region_probabilities(position_predictions, protein_regions, genome_seq_length):
            """
            Calculate protein region probabilities using raw probabilities without thresholding, 
            with proper mapping of start and end indices.

            Args:
                position_predictions (numpy array or dict): A numpy array or dictionary with mutation probabilities.
                protein_regions (dict): A dictionary with protein names as keys and tuples of (start, end) indices for each region.
                genome_seq_length (int): Length of the genome sequence.

            Returns:
                dict: A dictionary with protein names as keys and the average mutation probability for that region.
            """
            # Flatten the position_predictions array if it is 2D
            if isinstance(position_predictions, np.ndarray) and len(position_predictions.shape) > 1:
                position_predictions = position_predictions.flatten()

            protein_mutation_probs = {}

            for protein, (protein_start, protein_end) in protein_regions.items():
                # Ensure the region is within the genome length
                protein_start = max(0, protein_start)
                protein_end = min(genome_seq_length, protein_end)

                # Map the protein region indices to the correct indices for position_predictions
                prediction_start_idx = protein_start
                prediction_end_idx = protein_end

                # print(f"Processing {protein}: Start Index = {protein_start}, End Index = {protein_end}")
                # print(f"Mapping to Prediction Indices: Start = {prediction_start_idx}, End = {prediction_end_idx}")

                # Extract predictions
                if isinstance(position_predictions, np.ndarray):
                    if prediction_start_idx < len(position_predictions) and prediction_end_idx <= len(position_predictions):
                        relevant_predictions = position_predictions[prediction_start_idx:prediction_end_idx + 1]
                    else:
                        relevant_predictions = []
                    # print(f"Relevant predictions (numpy): {relevant_predictions}")
                elif isinstance(position_predictions, dict):
                    relevant_predictions = [position_predictions.get(pos, 0) for pos in range(prediction_start_idx, prediction_end_idx + 1)]
                    # print(f"Relevant predictions (dict): {relevant_predictions}")
                else:
                    raise ValueError("position_predictions must be a numpy array or dictionary.")

                # Compute average probability for the region
                avg_prob = sum(relevant_predictions) / len(relevant_predictions) if len(relevant_predictions) > 0 else 0
                protein_mutation_probs[protein] = avg_prob

            # print("Protein mutation probs:", protein_mutation_probs)
            return protein_mutation_probs
        protein_probs_start = time.time()
        protein_mutation_probs = calculate_protein_region_probabilities(predictions, protein_regions, genome_seq_length=29904 )
        # print("Protein mutation probs", protein_mutation_probs)
        measure_time("protein_region_probability_calculation", protein_probs_start)
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
        }
        
        measure_time("total_request_handling", start_time)
        end_time = time.time()
        print(f"Total runtime: {end_time - start_time} seconds")

        return JsonResponse(response_data)

    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": f"An error occurred: {str(e)}"}, status=500)
