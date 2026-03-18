import time
import numpy as np
from typing import Optional, Dict, List, Any

from .plugin_runtime import normalize_elapsed_day, normalize_node_ids


def measure_time(label: str, start_time: float) -> None:
    """
    Helper function to log elapsed time.
    
    Args:
        label: Description of the operation being timed
        start_time: Start timestamp from time.time()
    """
    elapsed_time = time.time() - start_time
    print(f"[{label}] Elapsed Time: {elapsed_time:.2f} seconds")


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


def calculate_protein_region_probabilities(position_predictions, protein_regions, genome_seq_length, predictions_offset=0):
    protein_mutation_probs = {}
    avg_mutation_per_position = np.mean(position_predictions, axis=1)
    n = len(avg_mutation_per_position)

    for protein, (protein_start, protein_end) in protein_regions.items():
        protein_start = max(0, protein_start)
        protein_end = min(genome_seq_length, protein_end)

        local_start = protein_start - predictions_offset
        local_end   = protein_end   - predictions_offset

        if 0 <= local_start < n and local_end <= n and local_start < local_end:
            relevant = avg_mutation_per_position[local_start:local_end]
            avg_prob = float(np.mean(relevant))
        else:
            avg_prob = 0.0

        protein_mutation_probs[protein] = avg_prob

    return protein_mutation_probs


def predict_mutations(
    cache_path: str,
    genome_seq: str,
    mutations: List,
    codon_mapper: str,
    config_file: Any,
    node_ids: Any,
    elapsed_day: Optional[int] = None,
    depth: float = 0,
    protein_regions: Optional[Dict] = None,
    selectedModel: str = "balanced_data_model",
    model_wrapper = None,
    feature_extractor = None,
    **custom_params
) -> np.ndarray:
    """
    Predicts mutation probabilities using the Protocol-based flow:
    feature_extractor.extract_features() → model_wrapper.preprocess() → predict() → postprocess()

    DefaultCovMutExFeatureExtractor generates (N*4, 205) features (4 per position),
    and CovMutExKerasModel.postprocess() reshapes output to (N, 4) [A, T, G, C].
    Custom extractors may return any shape; nucleotides_per_position defaults to 1.

    Returns:
        numpy array of predictions
    """
    metadata = model_wrapper.metadata()
    print(f"Using model: {metadata['name']}")
    print(f"Model type: {metadata['model_type']}")
    print(f"Output shape: {metadata.get('output_shape', 'N/A')}")

    extractor_name = feature_extractor.get_metadata()['name'] if feature_extractor else None
    print(f"[predict_mutations] Protocol approach (extractor: {extractor_name})")

    normalized_node_ids = normalize_node_ids(node_ids)
    normalized_elapsed_day = normalize_elapsed_day(elapsed_day)

    process_time = time.time()
    features = feature_extractor.extract_features(
        genome_seq=genome_seq,
        mutations=mutations,
        node_ids=normalized_node_ids,
        elapsed_day=normalized_elapsed_day,
        protein_regions=protein_regions,
        k=30,
        cache_path=cache_path,
        codon_mapper=codon_mapper,
        config_file=config_file,
        depth=depth,
        **custom_params
    )
    print(f"Feature extraction time: {time.time() - process_time:.2f}s")
    features = np.asarray(features)
    print(f"Feature shape received: {features.shape}")

    predict_time = time.time()
    preprocessed = model_wrapper.preprocess({
        'features': features,
        'genome_id': normalized_node_ids[0] if normalized_node_ids else None,
        'elapsed_days': normalized_elapsed_day,
        'num_inputs': metadata.get('num_inputs_expected', 1)
    })

    raw_predictions = model_wrapper.predict(preprocessed)

    # DefaultCovMutExFeatureExtractor ürettiği 4-per-position feature'ları (N*4, 205) → (N, 4) reshape
    npp = getattr(feature_extractor, 'nucleotides_per_position', 1)
    results = model_wrapper.postprocess(raw_predictions, nucleotides_per_position=npp)

    print(f"Prediction time: {time.time() - predict_time:.2f}s")
    print(f"Prediction shape: {results['shape']}")

    return results['predictions']
