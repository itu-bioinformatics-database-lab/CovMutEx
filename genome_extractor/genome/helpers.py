"""
Genome Prediction Helper Functions

This module contains business logic for mutation prediction,
genome data calculation, and probability distribution.

These functions are separated from views.py to improve:
- Code organization
- Testability
- Reusability across different parts of the application
"""

import time
import numpy as np
from typing import Optional, Dict, List, Tuple, Any

from .feature_extractor import cache_precomputed_features


def measure_time(label: str, start_time: float) -> None:
    """
    Helper function to log elapsed time.
    
    Args:
        label: Description of the operation being timed
        start_time: Start timestamp from time.time()
    """
    elapsed_time = time.time() - start_time
    print(f"[{label}] Elapsed Time: {elapsed_time:.2f} seconds")


def _get_biologically_distributed_probs(
    p_no_mutation: float, 
    p_mutation: float, 
    ref_nuc: str, 
    ti_tv_ratio: float = 2.0
) -> Dict[str, float]:
    """
    Takes separate no-mutation and mutation probabilities and returns a
    full, NORMALIZED probability distribution dict [A,T,G,C] based on Ti/Tv bias.

    This is the final, publication-standard helper function.
    
    Args:
        p_no_mutation: Probability of no mutation occurring
        p_mutation: Probability of mutation occurring
        ref_nuc: Reference nucleotide ('A', 'T', 'G', or 'C')
        ti_tv_ratio: Transition/Transversion ratio (default: 2.0)
    
    Returns:
        Dictionary with normalized probabilities for each nucleotide {A, T, G, C}
    """
    probabilities = {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}
    
    # Assign the initial (potentially un-normalized) probability
    probabilities[ref_nuc] = p_no_mutation

    # Define the single transition target
    purines = {'A', 'G'}
    pyrimidines = {'T', 'C'}
    if ref_nuc in purines:
        transition_target = list(purines - {ref_nuc})[0]
    else:  # ref_nuc is a pyrimidine
        transition_target = list(pyrimidines - {ref_nuc})[0]
        
    # Distribute the mutation probability based on Ti/Tv ratio
    if p_mutation > 0:
        # prob_tv is the probability for a single transversion event
        prob_tv = p_mutation / (ti_tv_ratio + 2)
        # prob_ti is the probability for the single transition event
        prob_ti = prob_tv * ti_tv_ratio
    else:
        prob_tv = 0.0
        prob_ti = 0.0
    
    for nuc in probabilities:
        if nuc == ref_nuc:
            continue
        elif nuc == transition_target:
            probabilities[nuc] = prob_ti
        else:  # It's a transversion
            probabilities[nuc] = prob_tv

    # --- FINAL NORMALIZATION STEP ---
    # This is crucial because for dual-output models, p_no_mutation + p_mutation
    # may not equal 1.0. This step guarantees the output is a valid distribution.
    total_sum = sum(probabilities.values())
    
    # Using a small tolerance to avoid division-by-zero for floating point numbers
    if total_sum > 1e-9:
        for nuc in probabilities:
            probabilities[nuc] /= total_sum
    else:
        # If the sum is zero, it implies no probability for anything.
        # A safe fallback is to assign 100% probability to the reference nucleotide.
        # This is a non-informative, zero-mutation state.
        for nuc in probabilities:
            probabilities[nuc] = 1.0 if nuc == ref_nuc else 0.0
            
    return probabilities


def predict_mutations(
    cache_path: str,
    genome_seq: str,
    mutations: List,
    codon_mapper: str,
    config_file: Any,
    node_ids: str,
    elapsed_day: Optional[int] = None,
    protein_regions: Optional[Dict] = None,
    selectedModel: str = "balanced_data_model",
    model_wrapper = None
) -> np.ndarray:
    """
    Predicts mutation probabilities for precomputed features using CovMutEx Model Protocol.
    
    Input features shape: (29904, 205) containing all possible transitions
    Returns matrix of shape (29904, N) where N depends on model output type
    
    Args:
        cache_path: Path to feature cache
        genome_seq: Reference genome sequence
        mutations: List of mutations
        codon_mapper: Path to codon mapping file
        config_file: Configuration object
        node_ids: Node identifier
        elapsed_day: Days elapsed (optional)
        protein_regions: Dict of protein regions (optional)
        selectedModel: Model name for logging
        model_wrapper: CovMutExKerasModel instance (Protocol-based wrapper)
    
    Returns:
        numpy array of predictions with shape (29904, N) where N depends on model output type
    """
    # Get model metadata
    metadata = model_wrapper.metadata()
    print(f"Using model: {metadata['name']}")
    print(f"Model type: {metadata['model_type']}")
    print(f"Output type: {metadata['output_type']}")
    print(f"Output shape: {metadata['output_shape']}")

    # Get precomputed features - shape (29904, 205)
    print("Selected model is ", selectedModel)
    
    process_time = time.time()
    print("Using model: ", selectedModel)

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
    print(f"Feature extraction time: {process_time_end - process_time:.2f}s")
    print(f"Feature shape received: {features.shape}")

    # Preprocess using Protocol's preprocess method
    predict_time = time.time()
    preprocessed = model_wrapper.preprocess({
        'features': features,
        'genome_id': node_ids,
        'elapsed_days': elapsed_day,
        'num_inputs': 10  # for multi-input models
    })
    
    # Predict using Protocol's predict method
    raw_predictions = model_wrapper.predict(preprocessed)
    
    # Postprocess using Protocol's postprocess method
    results = model_wrapper.postprocess(raw_predictions)
    
    predict_time_end = time.time()
    print(f"Prediction time: {predict_time_end - predict_time:.2f}s")
    print(f"Prediction shape: {results['shape']}")
    print(f"Interpretation: {results['interpretation']}")
    
    return results['predictions']


def calculate_genome_data(
    genome_seq: str,
    position_predictions: np.ndarray,
    selected_protein_region: Optional[Tuple[int, int]] = None
) -> List[List[float]]:
    """
    Calculate genome data using raw probabilities.
    
    Returns the probability for each nucleotide (A, T, G, C) at each position,
    using a biologically informed Ti/Tv distribution.
    Assumes position_predictions is either of shape (length, 1) or (length, 2).
    If a selected_protein_region is provided, only positions within that range are used.
    
    Args:
        genome_seq: Reference genome sequence string
        position_predictions: numpy array of predictions (N, 1) or (N, 2)
        selected_protein_region: Optional tuple of (start, end) positions
    
    Returns:
        List of 4 lists containing probabilities for A, T, G, C at each position
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
            position_probs_dict = {}

            # Handle different lengths of mutation probabilities for each position
            if mutation_probs.shape[0] == 1:
                # Single output model: prediction is prob_mutation
                p_mutation = mutation_probs[0]
                p_no_mutation = 1.0 - p_mutation
                position_probs_dict = _get_biologically_distributed_probs(p_no_mutation, p_mutation, current_nucleotide)
                
            elif mutation_probs.shape[0] == 2:
                # Dual output model: prediction is [prob_no_mutation, prob_mutation]
                p_no_mutation = mutation_probs[0]
                p_mutation = mutation_probs[1]
                position_probs_dict = _get_biologically_distributed_probs(p_no_mutation, p_mutation, current_nucleotide)

            else:
                raise ValueError(f"Unexpected prediction shape for position {position}: {mutation_probs.shape[0]} probabilities")

            # Append the probabilities for this position in ATGC order
            genome_data.append([position_probs_dict[nuc] for nuc in nucleotides])

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
                print(f"Warning: Predictions array length ({len(relevant_predictions)}) is shorter than genome. Stopping at position {idx-1}.")
                break

            # Get raw mutation probabilities for this position
            mutation_probs = relevant_predictions[idx]  # Array with shape (1,) or (2,)
            current_nucleotide = genome_seq[position]
            position_probs_dict = {}

            # Handle different lengths of mutation probabilities for each position
            if mutation_probs.shape[0] == 1:
                # Single output model: prediction is prob_mutation
                p_mutation = mutation_probs[0]
                p_no_mutation = 1.0 - p_mutation
                position_probs_dict = _get_biologically_distributed_probs(p_no_mutation, p_mutation, current_nucleotide)

            elif mutation_probs.shape[0] == 2:
                # Dual output model: prediction is [prob_no_mutation, prob_mutation]
                p_no_mutation = mutation_probs[0]
                p_mutation = mutation_probs[1]
                position_probs_dict = _get_biologically_distributed_probs(p_no_mutation, p_mutation, current_nucleotide)

            else:
                raise ValueError(f"Unexpected prediction shape for position {position}: {mutation_probs.shape[0]} probabilities")

            # Append the probabilities for this position in ATGC order
            genome_data.append([position_probs_dict[nuc] for nuc in nucleotides])
    
    # Convert to numpy array and transpose to get [P(A), P(T), P(G), P(C)] for each position
    return np.array(genome_data).T.tolist()


def calculate_protein_region_probabilities(
    position_predictions: np.ndarray,
    protein_regions: Dict[str, List[int]],
    genome_seq_length: int
) -> Dict[str, float]:
    """
    Calculate protein region probabilities using raw probabilities without thresholding.
    
    Uses proper mapping of start and end indices.
    Casts numpy floats to standard python floats for JSON serialization.
    
    Args:
        position_predictions: numpy array of predictions
        protein_regions: Dictionary mapping protein names to [start, end] positions
        genome_seq_length: Total length of genome sequence
    
    Returns:
        Dictionary mapping protein names to average mutation probabilities
    """
    # Flatten the position_predictions array if it is 2D
    if isinstance(position_predictions, np.ndarray) and len(position_predictions.shape) > 1:
        # If shape is (n, 2), we take the second column (mutation prob). If (n, 1), we just use it.
        if position_predictions.shape[1] == 2:
            position_predictions = position_predictions[:, 1]
        else:
            position_predictions = position_predictions.flatten()

    protein_mutation_probs = {}

    for protein, (protein_start, protein_end) in protein_regions.items():
        protein_start = max(0, protein_start)
        protein_end = min(genome_seq_length, protein_end)

        prediction_start_idx = protein_start
        prediction_end_idx = protein_end

        if isinstance(position_predictions, np.ndarray):
            if prediction_start_idx < len(position_predictions) and prediction_end_idx <= len(position_predictions):
                relevant_predictions = position_predictions[prediction_start_idx:prediction_end_idx]
            else:
                relevant_predictions = np.array([])
        else:
            # This case for dict is preserved but less likely with current predict_mutations output
            relevant_predictions = [position_predictions.get(pos, 0) for pos in range(prediction_start_idx, prediction_end_idx)]
        
    
        if len(relevant_predictions) > 0:
            avg_prob = float(np.mean(relevant_predictions))
        else:
            avg_prob = 0.0

        protein_mutation_probs[protein] = avg_prob

    return protein_mutation_probs


def read_genome_sequence(file_path: str) -> str:
    """
    Read and process genome sequence from file.
    
    Args:
        file_path: Path to genome sequence file
    
    Returns:
        Concatenated genome sequence string (header line skipped)
    """
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        return ''.join(line.strip() for line in f)
