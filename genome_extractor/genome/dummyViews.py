import numpy as np
import os
import json
import sys
import traceback
from collections import defaultdict
from django.http import JsonResponse
from rest_framework.decorators import api_view
from keras.models import load_model
from sklearn.preprocessing import OneHotEncoder, StandardScaler
import time
from rich import print
from genome_extractor.genome.feature_extractor import predict_mutations, parse_mutations
from configs import configs

sys.path.append(r'C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master')

codon_mapping_path=r'C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\codon_aa_mapping.json'
cache_path = r"F:\COVID19 MUTATION\node_features.h5"


print(f"File exists: {os.path.exists(cache_path)}")
print("got into this file too")
cache_path = "node_features.h5"         
def load_keras_model(model_path):
    """Load a Keras model from the specified path"""
    try:
        model = load_model(model_path, compile=False)
        return model
    except Exception as e:
        print(f"Error loading model: {e}")
        traceback.print_exc()
        return None

def calculate_genome_data(genome_seq, predictions):
    """Generate genome data from predictions."""
    genome_data = []
    nucleotides = "ATGC"
    
    for position in range(len(genome_seq)):
        position_probs = predictions[position]
        genome_data.append({
            nucleotide: position_probs[i] for i, nucleotide in enumerate(nucleotides)
        })
    
    return genome_data

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

def handle_prediction(request):
    try:
        start_time = time.time()
        
        # Load input data
        data = request.data if request.method == 'POST' else request.query_params
        nodeId = data.get('nodeId')
        elapsedDay = int(data.get('elapsedDay', 0))
        selectedModel = data.get('selectedModel')
        selectedProteinRegion = data.get('selectedProteinRegion')
        
        # Read genome sequence
        genome_file_path = r'C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\genome.txt'
        genome_sequence = read_genome_sequence(genome_file_path)
        genome_length = len(genome_sequence)

        # Load model
        model_path = os.path.join(r"C:\Users\DRX\COVID19 MUTATION\covid19_models", f"{selectedModel}.h5")
        model = load_keras_model(model_path)
        if model is None:
            return JsonResponse({"error": "Model could not be loaded"}, status=500)

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
        # features = precompute_feature_vectors(
        #     genome_txt_path=genome_file_path,
        #     phylo_tree_path=r'C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\phylogenetic_tree.nwk',
        #     protein_regions=selectedProteinRegion,
        #     config_file={},  # Placeholder for configs
        #     k=30
        # )

        predictions = predict_mutations(cache_path=cache_path, genome_seq=genome_sequence, 
                                mutations= parse_mutations(nodeId), codon_mapper=codon_mapping_path, 
                                config_file=configs(), 
                                node_ids=nodeId , elapsed_day= elapsedDay, protein_regions={"ORF3a": [25393,26220]})
        
        # features = np.array(features).reshape(genome_length, -1)

        def transform_predictions_to_stackable_format(predictions, genome_length=29904, protein_region=None):
            """
            Transforms the model predictions into the desired shape (4, genome_length), adjusted for protein region if selected.
            
            Arguments:
            - predictions: A list or array of predictions from the model.
                        Each element should be of shape (genome_length, n) where n is the number of outputs (1 or 2).
            - genome_length: Length of the genome sequence (defaults to 29904).
            - protein_region: A tuple (start, end) indicating the positions of the protein region in the genome.
            
            Returns:
            - A 4xgenome_length matrix with nucleotide probabilities, adjusted for the protein region.
            """
            # Initialize a list to hold the probabilities for each nucleotide (A, T, G, C)
            genome_data = []

            # If a protein region is selected, adjust the indices to only consider the protein region
            start, end = protein_region if protein_region else (0, genome_length)

            # Iterate through predictions for each position in the genome (within the protein region if selected)
            for position in range(start, end):
                position_probs = predictions[position]  # Shape could be (1,) or (2,)

                if len(position_probs) == 1:
                    # If we have only one probability (e.g., for one nucleotide type like PA), fill the appropriate row
                    genome_data.append([position_probs[0], 0, 0, 0])  # For 'A', and zeroes for T, G, C
                elif len(position_probs) == 2:
                    # If we have two probabilities (e.g., PA and PT), fill the appropriate rows
                    genome_data.append([position_probs[0], position_probs[1], 0, 0])  # For 'A' and 'T'
                else:
                    raise ValueError("Prediction shape should be (x, 1) or (x, 2)")

            # Convert to numpy array and transpose to get [P(A), P(T), P(G), P(C)] for each position
            return np.array(genome_data).T.tolist()

        def calculate_genome_data(predictions, protein_region=(25393, 26220)):
            """Prepare genome data in a format suitable for front-end visualization."""
            genome_data = transform_predictions_to_stackable_format(predictions, protein_region=protein_region)
            
            # Return the transformed data
            return genome_data

        # Example usage with dummy data and a protein region selected:
        # predictions_example = np.random.rand(29904, 2)  # Example prediction output with (29904, 2)
        # predictions = [predictions_example[i, :] for i in range(29904)]  # Convert to a list of shape (x, 2)

        # Define the protein region, e.g., "ORF3a" from position 25393 to 26220
        protein_region = (25393, 26220)

        genome_data = calculate_genome_data(predictions, protein_region=protein_region)
        print(np.array(genome_data).shape)  # Should print (29904, 4)


        # # Prepare response
        # genome_data = calculate_genome_data(genome_sequence, predictions)
        response_data = {
            "nodeId": nodeId,
            "elapsedDay": elapsedDay,
            "selectedModel": selectedModel,
            "selectedProteinRegion": selectedProteinRegion,
            "genomeSequence": genome_sequence,
            "genomeData": genome_data,
            "modelType": "multi-input" if 'multi' in selectedModel.lower() else "single-input",
        }
        import numpy as np


        end_time = time.time()
        print(f"Total runtime: {end_time - start_time:.2f} seconds")

        return JsonResponse(response_data)

    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": f"An error occurred: {str(e)}"}, status=500)
