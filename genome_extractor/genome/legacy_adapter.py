# /genome_extractor/genome/legacy_adapter.py

import os
import json
import tensorflow as tf
import numpy as np

# We import the new SDK interface we created from the project's root directory.
# For this to work, the project's main directory might need to be in the Python path.
# If you get an import error, we might need to adjust this part.
from covmutex_sdk.model import CovMutExModel 

# We import the functions used in the existing code from their respective files.
from .feature_extractor import parse_mutations, cache_precomputed_features, construct_variant_genome
from .configs import configs

# --- IMPORTANT NOTE ---
# Normally, these helper functions should be moved to a separate 'helpers.py' file.
# For the adapter to be self-sufficient, we are keeping them here for now.

def read_genome_sequence(file_path):
    """Reads the 'genome.txt' file and returns the genome sequence as a single string."""
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        return ''.join(line.strip() for line in f)

def _get_biologically_distributed_probs(p_no_mutation, p_mutation, ref_nuc, ti_tv_ratio=2.0):
    """Takes the raw mutation probability and converts it into a meaningful distribution for 4 nucleotides based on the Ti/Tv ratio."""
    probabilities = {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}
    probabilities[ref_nuc] = p_no_mutation
    
    purines, pyrimidines = {'A', 'G'}, {'T', 'C'}
    transition_target = list(purines - {ref_nuc})[0] if ref_nuc in purines else list(pyrimidines - {ref_nuc})[0]
        
    if p_mutation > 0:
        prob_tv = p_mutation / (ti_tv_ratio + 2)
        prob_ti = prob_tv * ti_tv_ratio
    else:
        prob_tv, prob_ti = 0.0, 0.0
    
    for nuc in probabilities:
        if nuc == ref_nuc: continue
        elif nuc == transition_target: probabilities[nuc] = prob_ti
        else: probabilities[nuc] = prob_tv

    total_sum = sum(probabilities.values())
    if total_sum > 1e-9:
        for nuc in probabilities: probabilities[nuc] /= total_sum
    else:
        for nuc in probabilities: probabilities[nuc] = 1.0 if nuc == ref_nuc else 0.0
            
    return probabilities

# --- ADAPTER CLASS ---

class LegacyModelAdapter(CovMutExModel):
    
    def __init__(self):
        """
        When an object of the class is first created, it loads the model, reference genome,
        and protein region information into memory.
        """
        # We determine the full path to the model file based on the project structure.
        # This file is located in 'covid19_models/models', one level above the 'genome_extractor' directory.
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        model_path = os.path.join(base_dir, 'covid19_models', 'models', 'balanced_data_model.keras')
        
        print(f"LegacyAdapter: Loading model: {model_path}")
        self.model = tf.keras.models.load_model(model_path)
        print("LegacyAdapter: Model loaded successfully.")

        # Let's also store the reference genome within the class
        genome_file_path = os.path.join(os.path.dirname(__file__), "genome.txt")
        self.reference_genome = read_genome_sequence(genome_file_path)

        self.protein_regions = {
            "ORF1ab": [266, 21555], "S": [21563, 25384], "ORF3a": [25393, 26220],
            "E": [26245, 26472], "M": [26523, 27191], "ORF6": [27202, 27387],
            "ORF7a": [27394, 27759], "ORF7b": [27756, 27887], "ORF8": [27894, 28259],
            "N": [28274, 29533], "ORF10": [29558, 29674],
        }

    def metadata(self) -> dict:
        """Returns the model's identification card (metadata)."""
        return {"name": "Balanced Data Model (Legacy)", "version": "1.0.0", "author": "Original Research Team"}

    def input_schema(self) -> dict:
        """Defines what inputs the model expects from the user."""
        return {
            "type": "object",
            "properties": {
                "nodeId": {"type": "string"},
                "elapsedDay": {"type": "integer", "default": 0},
                "selectedModel": {"type": "string"},
                "selectedProteinRegion": {"type": ["string", "null"]},
            },
            "required": ["nodeId", "selectedModel"],
        }

    def preprocess(self, inputs: dict) -> tuple:
        """
        Takes user inputs and creates the 'features' array that the model understands.
        It also returns the 'variant_genome_sequence' to be used in subsequent steps.
        """
        print("LegacyAdapter: Starting preprocess...")
        
        nodeId = inputs.get('nodeId')
        elapsedDay = inputs.get('elapsedDay', 0)
        selectedProteinRegion = inputs.get('selectedProteinRegion')
        
        mutations = parse_mutations(nodeId)
        variant_genome_sequence = construct_variant_genome(self.reference_genome, mutations)

        codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')
        cache_path = os.path.join(os.path.dirname(__file__), 'node_features.h5')
        
        # If a protein region is selected, we only send that region to the feature extractor
        protein_region_to_process = {}
        if selectedProteinRegion and selectedProteinRegion in self.protein_regions:
            protein_region_to_process = {selectedProteinRegion: self.protein_regions[selectedProteinRegion]}

        features = cache_precomputed_features(
            cache_path=cache_path, genome_seq=variant_genome_sequence, mutations=mutations,
            codon_mapper=codon_mapping_path, config_file=configs(), node_ids=[nodeId],
            elapsed_day=elapsedDay, protein_regions=protein_region_to_process
        )
        
        print("LegacyAdapter: Preprocess finished.")
        # We also return additional information needed in subsequent steps as a tuple
        return (features, variant_genome_sequence)

    def predict(self, batch: tuple) -> tuple:
        """Takes the preprocessed data and makes the raw prediction with the model."""
        print("LegacyAdapter: Making prediction...")
        
        features, variant_genome_sequence = batch # We separate the two values coming from preprocess
        
        # Since this model is single-input, we can simplify the 'multi' check
        raw_predictions = self.model.predict(features, verbose=0)
            
        print("LegacyAdapter: Prediction finished.")
        return (raw_predictions, variant_genome_sequence)

    def postprocess(self, prediction_bundle: tuple) -> dict:
        """Takes the raw predictions and calculates both the probability matrix and the protein summaries."""
        print("LegacyAdapter: Starting postprocess...")
        raw_predictions, variant_genome_sequence = prediction_bundle

        # --- 1. Calculate the Main Probability Matrix (genomeData) ---
        genome_data = []
        nucleotides = "ATGC"
        for idx in range(len(raw_predictions)):
            mutation_probs = raw_predictions[idx]
            current_nucleotide = variant_genome_sequence[idx]
            
            p_mutation = mutation_probs[0]
            p_no_mutation = 1.0 - p_mutation
            
            position_probs_dict = _get_biologically_distributed_probs(
                p_no_mutation, p_mutation, current_nucleotide
            )
            genome_data.append([position_probs_dict[nuc] for nuc in nucleotides])
        
        final_matrix = np.array(genome_data).T.tolist()
        final_probabilities_dict = {'A': final_matrix[0], 'T': final_matrix[1], 'G': final_matrix[2], 'C': final_matrix[3]}

        # --- 2. Calculate Protein Region Summaries (protein_mutation_probs) ---
        protein_mutation_probs = {}
        # We take only the mutation probability column from the raw predictions
        mutation_only_predictions = raw_predictions.flatten()
        
        for protein, (start, end) in self.protein_regions.items():
            # Adjust the boundaries according to the genome length
            start_idx = max(0, start)
            end_idx = min(len(mutation_only_predictions), end)
            
            relevant_predictions = mutation_only_predictions[start_idx:end_idx]
            
            avg_prob = float(np.mean(relevant_predictions)) if len(relevant_predictions) > 0 else 0.0
            protein_mutation_probs[protein] = avg_prob
        
        print("LegacyAdapter: Postprocess finished.")
        
        # Return a single result dictionary in the format expected by the SDK
        return {
            "genome_data": final_probabilities_dict,
            "protein_summary": protein_mutation_probs
        }