import json
import os
import random
from Bio import Phylo
import pandas as pd
from .configs import configs
from sklearn.utils import resample
import time
import h5py
import numpy as np
import csv
import time
import numpy as np
from sklearn.preprocessing import OneHotEncoder, StandardScaler
import cProfile
import tensorflow as tf
from sklearn.preprocessing import MinMaxScaler

# Paths to files
base_dir = os.path.dirname(os.path.abspath(__file__))

# Define relative paths
codon_mapping_path = os.path.join(base_dir, "covid19-genome-feature-extractor-master", "codon_aa_mapping.json")
genome_txt_path = os.path.join(base_dir, "covid19-genome-feature-extractor-master", "genome.txt")
mutations_txt_path = os.path.join(base_dir, "covid19-genome-feature-extractor-master", "mutations.txt")
phylo_tree_path = os.path.join(base_dir, "covid19-genome-feature-extractor-master", "phylogenetic_tree.nwk")
depth_file = os.path.join(base_dir, 'depth_date.json')
genome_extractor_dir = os.path.dirname(base_dir)
ROOT_PATH = os.path.dirname(genome_extractor_dir)
model_path = os.path.join(genome_extractor_dir, 'covid19_models', "model_for_1k.h5")



start_time = time.time()
# profiler = cProfile.Profile()

# profiler.enable()
# Helper functions
def get_sequence():
    with open(genome_txt_path) as file:
        next(file)
        genome = file.read().replace('\n', '')
    genome += 'A'
    return genome

def parse_mutations(nodeId):
    mutations = []
    with open(mutations_txt_path, 'r') as file:
        header = file.readline().strip().split('\t')
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) >= 3 and fields[10] == nodeId:
                nt_position = int(fields[1])
                mutation = fields[2]
                original, new = mutation.split('>')
                aa_position = int(fields[3])
                aa_change = fields[4]
                mutations.append((nt_position, original, new, aa_position, aa_change))
    return mutations

import h5py
import numpy as np
from Bio import Phylo

CACHE_P_FILE = 'phylo_features_cache.h5'

def normalize_features(features):
    min_val = min(features.values())
    max_val = max(features.values())
    if max_val == min_val:
        return {k: 1.0 for k in features}  # Prevent division by zero
    return {k: (v - min_val) / (max_val - min_val) for k, v in features.items()}

def extract_phylogenetic_features(tree_path):
    """
    Extract phylogenetic features from a tree file, with intelligent caching in an HDF5 file.
    
    Args:
        tree_path (str): Path to the Newick format phylogenetic tree file
    
    Returns:
        tuple: A tuple containing normalized distances and diversity metrics
    """
    normalized_distances = None
    diversity_metrics = {}

    try:
        with h5py.File(CACHE_P_FILE, 'a') as f:
            # Create a group for the specific tree if it doesn't exist
            if tree_path not in f:
                print(f"Extracting phylogenetic features for {tree_path}...")
                
                try:
                    # Read the tree
                    tree = Phylo.read(tree_path, 'newick')
                    
                    # Extract clade distances
                    clade_distances = {}
                    for clade in tree.find_clades():
                        if clade.name:
                            clade_distances[clade.name] = clade.branch_length or 0.0
                    
                    # Validate clade distances
                    if not clade_distances:
                        raise ValueError("No valid clade distances found in the tree")
                    
                    # Normalize distances
                    normalized_distances = normalize_features(clade_distances)
                    
                    # Compute diversity metrics
                    diversity_metrics = calculate_phylogenetic_diversity(tree)
                    
                    # Create a group for this specific tree path
                    tree_group = f.create_group(tree_path)
                    
                    # Store clade names, normalized distances, and diversity metrics
                    tree_group.create_dataset('clade_names', data=list(normalized_distances.keys()))
                    tree_group.create_dataset('normalized_distances', data=list(normalized_distances.values()))
                    
                    # Handle potential byte string encoding for diversity metrics
                    diversity_names = list(diversity_metrics.keys())
                    diversity_values = list(diversity_metrics.values())
                    
                    tree_group.create_dataset('diversity_metrics_names', 
                                              data=[name.encode('utf-8') if isinstance(name, str) else name 
                                                    for name in diversity_names])
                    tree_group.create_dataset('diversity_metrics', data=diversity_values)
                
                except Exception as e:
                    print(f"Error extracting features from {tree_path}: {e}")
                    raise
            
            # Retrieve cached data
            tree_group = f[tree_path]
            clade_names = list(tree_group['clade_names'][:])
            normalized_distances = dict(zip(clade_names, tree_group['normalized_distances'][:]))
            
            # Handle potential byte string decoding for diversity metrics names
            diversity_names = [name.decode('utf-8') if isinstance(name, bytes) else name 
                               for name in tree_group['diversity_metrics_names'][:]]
            diversity_values = tree_group['diversity_metrics'][:]
            diversity_metrics = dict(zip(diversity_names, diversity_values))
            
            print(f"Successfully retrieved/cached phylogenetic features for {tree_path}")
    
    except Exception as e:
        print(f"Critical error processing tree {tree_path}: {e}")
        raise
    
    # Final validation
    if normalized_distances is None:
        raise ValueError("Failed to compute or retrieve normalized distances")
    
    return normalized_distances, diversity_metrics
def calculate_phylogenetic_diversity(tree):
    clade_depths = {clade.name: clade.branch_length or 0.0 for clade in tree.find_clades() if clade.branch_length}
    diversity_metrics = {
        'total_branch_length': sum(clade_depths.values()),
        'max_branch_length': max(clade_depths.values(), default=0),
        'min_branch_length': min(clade_depths.values(), default=0),
    }
    return diversity_metrics


def translate_nucleotides_to_amino_acids(nucleotide_sequence, codon_mapper):
    aa_seq = ""
    for i in range(0, len(nucleotide_sequence), 3):
        codon = nucleotide_sequence[i:i + 3]
        if codon in codon_mapper:
            aa_seq += codon_mapper[codon]
    return aa_seq

def find_protein_region(index, protein_regions):
    for protein, (start, end) in protein_regions.items():
        if start <= index <= end:
            return protein
    return "Non-coding"

def compute_nucleotide_frequencies(sequence, window_size=10):
    freq_features = []
    for i in range(len(sequence)):
        start = max(0, i - window_size // 2)
        end = min(len(sequence), i + window_size // 2)
        window = sequence[start:end]
        freqs = {nt: window.count(nt) / len(window) for nt in 'ATGC'}
        freq_features.append(freqs)
    return freq_features

def construct_variant_genome(genome_seq, mutations):
    """Construct variant genome by applying mutations."""
    variant_genome = list(genome_seq)
    applied_mutations = []

    for mutation in mutations:
        position, original, new, _, _  = mutation 
        if variant_genome[position - 1] == original:  # Adjusted position to 0-based index
            variant_genome[position - 1] = new
            applied_mutations.append((position, original, new))

    print(f"Applied {len(applied_mutations)} mutations to the reference genome")
    return "".join(variant_genome)

def is_synonymous(current_aa, new_aa):
    return int(current_aa == new_aa)

def balance_dataset(dataset):
    data_df = pd.DataFrame(dataset)
    majority_class = data_df[data_df.iloc[:, -1] == 0]
    minority_class = data_df[data_df.iloc[:, -1] == 1]
    
    # Oversample minority class
    minority_oversampled = resample(minority_class, 
                                    replace=True, 
                                    n_samples=len(majority_class), 
                                    random_state=42)
    balanced_dataset = pd.concat([majority_class, minority_oversampled])
    return balanced_dataset.values.tolist()



def preprocess_input(features, expected_size=205, is_multi_input=False):
    """
    Preprocess a single feature vector:
    - Dynamically identifies categorical and numerical data
    - One-hot encodes categorical data
    - Standardizes numerical data
    - Ensures consistent feature vector size
    """
    import numpy as np
    from sklearn.preprocessing import OneHotEncoder, StandardScaler

    # Convert features to a numpy array for easier manipulation
    features = np.array(features, dtype=object)

    # Dynamically determine categorical and numerical indices
    categorical_indices = [
        i for i, val in enumerate(features)
        if isinstance(val, str) and val not in {"*", "-"}
    ]
    numerical_indices = [
        i for i, val in enumerate(features)
        if isinstance(val, (int, float)) and val not in {None, "*", "-"}
    ]

    ### Handle categorical features ###
    if categorical_indices:
        # Extract categorical data
        cat_data = features[categorical_indices]
        # Replace invalid placeholders ("-") with "_"
        cat_data = np.array([str(x) if x != "-" else "_" for x in cat_data], dtype=object)

        # Encode categorical features
        encoder = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        one_hot_encoded = encoder.fit_transform(cat_data.reshape(-1, 1)).flatten()
    else:
        one_hot_encoded = np.array([])  # If no categorical features, set to empty array

    ### Handle numerical features ###
    if numerical_indices:
        # Extract numerical data
        num_data = features[numerical_indices]
        # Replace invalid placeholders with NaN
        num_data = np.array([float(x) if x not in {None, "-", "*"} else np.nan for x in num_data])

        # Standardize numerical data
        scaler = StandardScaler()
        standardized_numerical = scaler.fit_transform(num_data.reshape(-1, 1)).flatten()
    else:
        standardized_numerical = np.array([])  # If no numerical features, set to empty array

    ### Combine all features ###
    combined_features = np.concatenate((one_hot_encoded, standardized_numerical))

    ### Ensure feature vector matches expected size ###
    if len(combined_features) < expected_size:
        combined_features = np.pad(combined_features, (0, expected_size - len(combined_features)), mode="constant")
    else:
        combined_features = combined_features[:expected_size]

    # Handle multi-input format
    if is_multi_input:
        return np.array([combined_features for _ in range(10)]).reshape(1, 10, -1)
    else:
        return combined_features.reshape(1, -1)


#FEATURES FOR THE REFERENCE GENOME

# Function to precompute feature vectors
def precompute_feature_vectors(
    genome_txt_path, codon_mapping_path, phylo_tree_path, protein_regions, config_file, k=30
):
    """
    Precompute feature vectors for all positions in the reference genome without any mutations.

    Parameters:
    - genome_txt_path (str): Path to the genome sequence file.
    - codon_mapping_path (str): Path to the codon-to-amino-acid mapping file.
    - phylo_tree_path (str): Path to the phylogenetic tree file in Newick format.
    - protein_regions (dict): Dictionary of protein regions with start and end indices.
    - config_file (dict): Configuration dictionary with nucleotide and amino acid properties.
    - k (int): Size of the sliding window for feature extraction.

    Returns:
    - dataset (list): List of preprocessed feature vectors for all positions in the genome.
    """


    start_time = time.time()
    mid_point = k // 2
    required_vector_length = k + 29  # Adjust based on expected feature count

    # Load configuration
    config_file = configs()
    genome_seq = get_sequence()
    with open(codon_mapping_path) as codon_file:
        codon_mapper = json.load(codon_file)


    phylo_features, phylo_diversity = extract_phylogenetic_features(phylo_tree_path)
    aa_seq = translate_nucleotides_to_amino_acids(genome_seq, codon_mapper)
    nucleotide_frequencies = compute_nucleotide_frequencies(genome_seq)

    # Initialize padded sequence
    padded_seq = "-" * mid_point + genome_seq + "-" * mid_point

    # Extract features for each position
    dataset = []

    for idx in range(len(genome_seq)):
        try:
            window = padded_seq[idx:idx + k]
            center_nucleotide = genome_seq[idx]
            
            aa_idx = idx // 3
            current_aa = aa_seq[aa_idx] if aa_idx < len(aa_seq) else "X"

            # Build the feature vector
            feature_vector = list(window)  # k-mer nucleotides
            feature_vector.extend([
                center_nucleotide,
                idx,
                phylo_features.get(center_nucleotide, 0),  # Normalized
                current_aa
            ])

            # Add nucleotide frequencies
            nucleotide_freq = nucleotide_frequencies[idx]
            feature_vector.extend(list(nucleotide_freq.values()))

            # Add amino acid features
            aa_features = config_file['AA features']
            feature_vector.extend([
                float(aa_features['hydrophobicity'].get(current_aa, 0)),
                float(aa_features['polarity'].get(current_aa, 0)),
                float(aa_features['iso-electric point'].get(current_aa, 0)),
                float(aa_features['volume'].get(current_aa, 0)),
                float(aa_features['molecular weight'].get(current_aa, 0))
            ])

            # Ensure consistent feature vector length
            if len(feature_vector) < required_vector_length:
                feature_vector.extend([0.0] * (required_vector_length - len(feature_vector)))
            elif len(feature_vector) > required_vector_length:
                feature_vector = feature_vector[:required_vector_length]

            # Preprocess the feature vector
            preprocessed_vector = preprocess_input(feature_vector, expected_size=205)

            # Add to dataset
            dataset.append(preprocessed_vector)
        except Exception as e:
            print(f"Error at position {idx}: {e}")
            continue

    end_time = time.time()
    print(f"Feature extraction and preprocessing completed in {end_time - start_time:.2f} seconds")

    return dataset

cache_path = os.path.join(ROOT_PATH, 'node_features.h5')
#Call precompute function and cache result(s)
def precompute_and_cache_ref_features(
    genome_txt_path, codon_mapping_path, phylo_tree_path, protein_regions, config_file, output_h5_path, k=30
):
    if os.path.exists(output_h5_path):
        print(f"Using cached feature vectors from {output_h5_path}")
        with h5py.File(output_h5_path, 'r') as h5f:
            dataset_np = h5f['features'][:]
            print("Loaded cached data from HDF5.")
        return dataset_np
    else:
        print("Cached file not found, computing features...")
    # Call the function to compute features
    dataset = precompute_feature_vectors(
        genome_txt_path, codon_mapping_path, phylo_tree_path, protein_regions, config_file, k
    )
    
    # Convert dataset to a NumPy array for HDF5 compatibility
    dataset_np = np.array(dataset, dtype=np.float32)

    # Save to HDF5 file
    with h5py.File(output_h5_path, 'w') as h5f:
        h5f.create_dataset("features", data=dataset_np, compression="gzip", compression_opts=9)
        h5f.attrs['genome_length'] = len(dataset)
        h5f.attrs['vector_length'] = dataset_np.shape[1]
    
    print(f"Feature vectors cached to {output_h5_path}")

# Example usage:
precompute_and_cache_ref_features(
    genome_txt_path=genome_txt_path,
    codon_mapping_path= codon_mapping_path,
    phylo_tree_path=phylo_tree_path,
    protein_regions={},
    config_file=configs,  # Passing the configs dictionary here
    output_h5_path='features.h5',
    k=30
)

def get_affected_positions(mutations, genome_len, k, protein_regions=None):
    """
    Get positions affected by mutations within the k-mer window and protein region.
    Only returns positions that actually need feature recomputation.
    
    Args:
        mutations (list): List of mutation tuples (position, ref, alt, aa_pos, aa_change)
        genome_len (int): Length of genome sequence
        k (int): Size of k-mer window
        protein_regions (dict): Optional dictionary of protein regions with protein name as key 
                                and (start, end) tuple as value
        
    Returns:
        set: Positions requiring feature recomputation
    """
    affected_positions = set()
    mid_point = k // 2
    
    for position, ref, alt, aa_pos, aa_change in mutations:
        pos = position - 1  # Convert to 0-based indexing
        
        # If protein region specified, skip mutations outside it
        if protein_regions:
            # Iterate over all the protein regions
            for region_start, region_end in protein_regions.values():
                if region_start <= pos <= region_end:
                    break
            else:
                continue  # Skip mutation if it's outside all provided protein regions
        
        # Add positions within k-mer window that need updating
        window_start = max(0, pos - mid_point)
        window_end = min(genome_len, pos + mid_point + 1)
        
        # Only add positions that actually need feature updates
        for idx in range(window_start, window_end):
            # If within protein region constraint
            if not protein_regions or any(region_start <= idx <= region_end for region_start, region_end in protein_regions.values()):
                affected_positions.add(idx)
                
    return sorted(affected_positions)


def get_aa_features(aa, new_aa, config_file):
    """
    Retrieve features for amino acids, only when needed.
    """
    if aa == new_aa:  # No need to fetch features if the amino acid hasn't changed
        return [0] * 16  # Return a placeholder array if no change

    return [
        config_file['AA features']['hydrophobicity'].get(new_aa, 0),
        config_file['AA features']['polarity'].get(aa, 0),
        config_file['AA features']['polarity'].get(new_aa, 0),
        config_file['AA features']['iso-electric point'].get(aa, 0),
        config_file['AA features']['iso-electric point'].get(new_aa, 0),
        config_file['AA features']['volume'].get(aa, 0),
        config_file['AA features']['volume'].get(new_aa, 0),
        config_file['AA features']['molecular weight'].get(aa, 0),
        config_file['AA features']['molecular weight'].get(new_aa, 0),
        config_file['AA features']['pKa'].get(aa, 0),
        config_file['AA features']['pKa'].get(new_aa, 0),
        config_file['AA features']['pKb'].get(aa, 0),
        config_file['AA features']['pKb'].get(new_aa, 0),
        config_file['AA features']['pKx'].get(aa, 0),
        config_file['AA features']['pKx'].get(new_aa, 0),
        config_file['AA features']['pl'].get(aa, 0),
        config_file['AA features']['pl'].get(new_aa, 0),
    ]

def get_sample_depth(depth_file_path, nodeId):
    """
    Extracts the depth value for a given sample from a depth JSON file.

    Parameters:
    - json_file_path (str): Path to the JSON file containing depth information.
    - sample_name (str): The name of the sample to retrieve the depth for.

    Returns:
    - int: Depth value for the sample if found.
    - None: If the sample is not found in the JSON file.
    """
    try:
        # Load the JSON file
        with open(depth_file, 'r') as file:
            depth_data = json.load(file)
        
        # Check if the sample exists in the JSON data
        if nodeId in depth_data:
            return depth_data[nodeId].get("depth", None)
        else:
            print(f"Sample '{nodeId}' not found in the JSON file.")
            return None
    except FileNotFoundError:
        print(f"File '{depth_file}' not found.")
        return None
    except json.JSONDecodeError:
        print(f"Error decoding JSON file '{depth_file}'.")
        return None

# Example usage
# node_id = "Northern_Ireland/NIRE-0048ed/2021|2021-07-03"
# depth_file = 'covid19-genome-feature-extractor-master/depth_date.json'
# depth = get_sample_depth(depth_file, node_id)
# print(f"Depth for {node_id}: {depth}")

def process_variant(
    genome_seq, 
    mutations, 
    precomputed_features, 
    codon_mapper, 
    config_file, 
    protein_regions=None,  # Make protein_regions optional
    k=30,
    affected_positions_set=None,
    elapsed_day=None,  # Make elapsed_day dynamic
    nodeId=None
):
    """
    Process a variant by recalculating features for mutated regions, 
    using precomputed features for non-mutated regions, and narrowing down
    predictions to the selected protein region if provided. 
    If no protein region is given, it only considers the mutation window.
    """
    startPro = time.time()
    mid_point = k // 2
    padded_seq = "-" * mid_point + genome_seq + "-" * mid_point
    updated_features = precomputed_features.copy()
    nucleotides = ['A', 'T', 'C', 'G']
    with open(codon_mapping_path) as codon_file:
        codon_mapper = json.load(codon_file)

    aa_seq = [codon_mapper.get(genome_seq[i:i+3], 'X') for i in range(0, len(genome_seq), 3)]
    print("Started extracting phylo info")
    startPhylo = time.time()
    phylo_features, phylo_diversity = extract_phylogenetic_features(phylo_tree_path)
    endPhylo = time.time()
    endPhylo = endPhylo - startPhylo 
    print("Finished extracting phylo info", endPhylo)
    print("phylo feature", type(phylo_features))
    if 'total_branch_length' in phylo_diversity:
        print(phylo_diversity['total_branch_length'])
    else:
        print("total_branch_length key not found!")


    # Get affected positions based on mutations
    affected_positions_set = get_affected_positions(mutations, len(genome_seq), mid_point, protein_regions=protein_regions)
    
    if protein_regions:
        # Filter affected positions to be within protein region bounds
        for protein, (region_start, region_end) in protein_regions.items():
            affected_positions_set = {idx for idx in affected_positions_set if idx in range(region_start, region_end + 1)}

    phylo_value = phylo_features.get(nodeId, 0) if nodeId else 0

    # List to store mutated positions
    mutated_positions_indices = []

    # Process affected positions
    for idx in affected_positions_set:
        try:
            # If precomputed features exist, and mutation is not at the current position, skip it
            mutated_positions = {m[0] - 1 for m in mutations}
            if idx in precomputed_features and idx not in mutated_positions:
                updated_features[idx] = precomputed_features[idx]
                continue

            # Mark the position as mutated
            if idx in mutated_positions:
                mutated_positions_indices.append(idx)

            # Find the protein region for the current idx if protein_regions is provided
            protein_reg = find_protein_region(idx, protein_regions) if protein_regions else None

            # Sliding window extraction
            window = padded_seq[idx:idx + k]
            aa_idx = idx // 3
            codon_start = aa_idx * 3
            aa_nucleotides = list(genome_seq[codon_start:codon_start + 3])
            original_aa = aa_seq[aa_idx]
            original_aa_features = get_aa_features(original_aa, original_aa, config_file)

            # Ignore windows outside genome boundaries or incomplete codons
            if len(aa_nucleotides) != 3:
                continue
            
            mutation_position = idx % 3  # Nucleotide position within codon
            
            original_codon = list(genome_seq[codon_start:codon_start + 3])
            for nucleotide in nucleotides:
                mutated_codon = original_codon[:]
                mutated_codon[mutation_position] = nucleotide
                new_codon = ''.join(mutated_codon)
                new_aa = codon_mapper.get(new_codon, 'X')
                new_aa_features = (
                    get_aa_features(original_aa, new_aa, config_file)
                    if new_aa != original_aa
                    else original_aa_features
                )

                # Build feature vector for the mutated position
                sample_data = [
                    *list(window),  # 1...k => sliding window
                    window[mid_point],  # k+1 => center of k-mer
                    nucleotide,  # k+2 => nucleotide after mutation
                    idx,  # k+3 => index
                    config_file['nucleotide sub. matrix'][window[mid_point]][nucleotide],  # k+4 => PAM score
                    aa_seq[aa_idx],  # k+5 => original amino acid
                    new_aa,  # k+6 => mutated amino acid
                    config_file['AA PAM matrix'][aa_seq[aa_idx]][new_aa],  # k+7 => AA PAM score
                    int(aa_seq[aa_idx] == new_aa),  # k+8 => synonymous (1) or non-synonymous (0)
                    elapsed_day,  # k+9 => dynamic elapsed day value
                    protein_reg,  # k+10 => protein region (optional)
                    phylo_value,  # Normalized phylogenetic feature
                    phylo_diversity['total_branch_length'],
                    *new_aa_features,
                ]

                # Preprocess and update the feature vector for the mutated position
                updated_features[idx] = preprocess_input(sample_data, expected_size=205)
                
        except Exception as e:
            print(f"Error processing position {idx}: {e}")
            continue
    print("length of mutation positions",  len(mutated_positions_indices),  mutated_positions_indices)
    # Return both the updated features and the list of mutated positions
    endPro = time.time()
    print("Time spent processing Variant: ", endPro-startPro)
    return updated_features

node_ids = ["Northern_Ireland/NIRE-0048ed/2021|2021-07-03"]
def cache_precomputed_features(
    cache_path, genome_seq, mutations, codon_mapper, config_file, node_ids, 
    elapsed_day=None, protein_regions=None
):
    # Load base precomputed features
    with h5py.File('features.h5', 'r') as h5f:
        base_features = h5f['features'][:]
        print(f"Base features shape: {base_features.shape}")
    
    # Function to map limited features back to full space
    def map_features_to_full_space(limited_features, affected_positions, base_features):
        """
        Map limited mutation probabilities back to the full genome feature space
        while maintaining sequential mapping of positions.
        """
        full_features = np.copy(base_features)
        nucleotides = "ATGC"
        limited_features = np.array(limited_features)

        for idx, pos in enumerate(affected_positions):
            if 0 <= pos < len(base_features):
                current_nucleotide = genome_seq[pos]

                # Initialize position probabilities
                position_probs = {nucleotide: 0.0 for nucleotide in nucleotides}

                if limited_features.shape[1] == 2:  # [No mutation, Mutation]
                    prob_no_mutation, prob_mutation = limited_features[idx]
                    position_probs[current_nucleotide] = prob_no_mutation
                    remaining_nucleotides = [n for n in nucleotides if n != current_nucleotide]
                    for nucleotide in remaining_nucleotides:
                        position_probs[nucleotide] = prob_mutation / 3
                else:
                    raise ValueError(f"Unexpected feature shape: {limited_features.shape[1]} probabilities")

                # Update full feature set at the valid position in ATGC order
                full_features[pos] = [position_probs[nucleotide] for nucleotide in nucleotides]
            else:
                raise ValueError(f"Invalid position {pos} outside genome bounds.")

        return full_features




    with h5py.File(cache_path, 'r+') as hdf:

    
        for node_id in node_ids:
        
            node_id_str = str(node_id)
            affected_positions = sorted({m[0] for m in mutations})
            print(f"Affected positions: {affected_positions}")
            
            # Check if we need to process mutations
            if node_id_str not in hdf:
                print(f"Processing new mutations for NodeId {node_id}")
                
                if affected_positions:
                    # Get features for affected positions
                    limited_features = process_variant(
                        genome_seq,
                        mutations,
                        {},  # Empty dict since we're not using precomputed
                        codon_mapper,
                        config_file,
                        nodeId=node_id,
                        affected_positions_set=set(affected_positions)
                    )
                    print(f"Limited features length: {len(limited_features) if not isinstance(limited_features, dict) else len(limited_features.keys())}")
                    
                    # Map limited features to full space
                    combined_features = map_features_to_full_space(limited_features, affected_positions, base_features)
                else:
                    combined_features = np.copy(base_features)
                
                # Cache the new features
                if node_id_str in hdf:
                    del hdf[node_id_str]
                hdf.create_dataset(node_id_str, data=combined_features)
            else:
                print(f"Loading cached features for NodeId {node_id}")
                combined_features = hdf[node_id_str][:]

            # Update elapsed day if provided
            if elapsed_day is not None:
                combined_features[10] = elapsed_day
                print(f"Updated elapsed day")

            # Handle protein regions if provided
            if protein_regions:
                region_features_dict = {}
                for protein, (region_start, region_end) in protein_regions.items():
                    region_positions = set(range(region_start, region_end + 1))
                    affected_region = region_positions.intersection(affected_positions)
                    
                    if not affected_region:
                        print(f"No mutations in {protein} region")
                        # Use base features for unaffected regions
                        region_features_dict[protein] = base_features[region_start:region_end + 1]
                    else:
                        print(f"Processing mutations in {protein} region")
                        region_mutations = [m for m in mutations if m[0] in affected_region]
                        print("Found the regions, starting Variant processing ")
                        limited_region_features = process_variant(
                            genome_seq,
                            region_mutations,
                            {},
                            codon_mapper,
                            config_file,
                            nodeId=node_id,
                            affected_positions_set=affected_region
                        )
                        
                        # Map region features to their positions
                        region_features = map_features_to_full_space(
                            limited_region_features, 
                            affected_region,
                            combined_features[region_start:region_end + 1]
                        )
                        region_features_dict[protein] = region_features

                # Combine all selected regions' features and return
                combined_features = np.concatenate(
                    [region_features_dict[protein] for protein in protein_regions],
                    axis=0
                )
                return combined_features  # Return only the features for the selected protein regions

            # Final processing
            combined_features = np.array(combined_features, dtype=np.float32)
            combined_features = np.nan_to_num(combined_features)
            
            # Update cache
            if node_id_str in hdf:
                del hdf[node_id_str]
            hdf.create_dataset(node_id_str, data=combined_features)
            
            print(f"Final features shape: {combined_features.shape}")
            return combined_features

    return None

def retrieve_features(cache_path, nodeId):
    with h5py.File(cache_path, 'r') as hdf:
        if str(nodeId) in hdf:
            return hdf[str(nodeId)][:]
        else:
            print(f"NodeId {nodeId} not found in cache.")
            return None

def combine_features(precomputed_features, variant_features, affected_positions):
    """
    Combine precomputed reference features with variant-specific features.
    
    Parameters:
    - precomputed_features: numpy array of reference genome features
    - variant_features: numpy array of variant-specific features
    - affected_positions: array of positions affected by mutations
    
    Returns:
    - combined_features: numpy array containing the merged features
    """
    combined_features = precomputed_features.copy()
    
    # Update only the affected positions with variant features
    for idx, pos in enumerate(affected_positions):
        if pos < len(combined_features):
            combined_features[pos] = variant_features[idx]
    
    return combined_features


# features = retrieve_features(cache_path, "USA/OH-CDC-ASC210099476/2021|OK223319.1|2021-09-02")
genome_sequence = get_sequence()
nodeList = [
	"USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07",
]

node_id = [x for x in nodeList]



# Cache precomputed features
# cache_precomputed_features(cache_path, 
#                            genome_seq= genome_sequence,
#                              mutations= mutations, codon_mapper=json.load(open(codon_mapping_path)), 
#                              config_file=configs(), node_ids=node_ids, elapsed_day=67)


with h5py.File('features.h5', 'r') as h5f:
    precomputed_features = h5f['features'][:]

# print(precomputed_features[:3])
# print(precomputed_features[20664])
# print("length of precomputed", precomputed_features.shape)
precomputed_features = precomputed_features.squeeze()

# Test with a cached NodeId
node_ids = ["Northern_Ireland/NIRE-0048ed/2021|2021-07-03"]
mutations = parse_mutations(node_id[0])
depth = get_sample_depth(depth_file, node_id[0])

# features = cache_precomputed_features(
#         cache_path, 
#         genome_seq= construct_variant_genome(genome_sequence, mutations),
#         mutations=parse_mutations(node_id),
#         codon_mapper=json.load(open(codon_mapping_path)),
#         config_file=configs(),
#         node_ids=node_id,
#         elapsed_day=110,
#         protein_regions=None
#     )




# print("length", len(features), type(features), features.shape)

  # Should directly load cached features without calling process_variant


# # Process variant
# updated_features = process_variant(
#     genome_seq=get_sequence(),
#     mutations=mutations,
#     precomputed_features=precomputed_features,
#     codon_mapper=json.load(open(codon_mapping_path)),
#     config_file=configs(),
#     # protein_regions={"S": [345, 2020]},  # Add real protein regions if available
#     k=30,
#     nodeId = "USA/OH-CDC-ASC210099476/2021|OK223319.1|2021-09-02"
# )




# # Retrieve features for a specific nodeId

# print(features)
# # comb_feat = combine_features()
# # # Check the updated feature vector for a specific position
# # print(updated_features)
# print("length of features", len(updated_features))
# # profiler.disable()  # Stop profiling
# # profiler.print_stats(sort="cumtime")  # Print stats sorted by cumulative time

# from django.core.cache import cache

# def get_selected_model():
#     selected_model = cache.get('selected_model')
#     print(f"Cached selected model is: {selected_model}")
#     return selected_model

# selected_model_name = get_selected_model()

# if not selected_model_name:
#     # Option 1: Set a default model name
#     selected_model_name = "model_for_1k"
# else:
#     print("Using model", selected_model_name)
#prev_dir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
#path_to_covid19_models = os.path.join(prev_dir, 'covid19_models'
current_dir = os.getcwd()
genome_extractor_path = os.path.dirname(current_dir)
# model_path = '../covid19_models/model_for_1k.h5'
# print("model", selected_model_name)
model = tf.keras.models.load_model(model_path)


def predict_mutations(cache_path, genome_seq, mutations, codon_mapper, config_file, node_ids, elapsed_day=None, protein_regions=None):
    """
    Predicts mutation probabilities for precomputed features.
    Input features shape: (29904, 205) containing all possible transitions
    Returns matrix of shape (29904, 4) where each row contains probabilities for A,T,G,C transitions
    """
    
    # Get precomputed features - shape (29904, 205)
    # print("Using model: ",selected_model_name )
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

    print("Shape of feature received: ", features.shape)
    print("Length of Features :", len(features))

    features = features.squeeze()
    print("Shape of feature before model input: ", features.shape)
    predictions = model.predict(features, verbose=0)
    
    return predictions




predictions = predict_mutations(cache_path=cache_path, genome_seq=genome_sequence, 
                                mutations=mutations, codon_mapper=codon_mapping_path, 
                                config_file=configs(), 
                                node_ids=node_id[0], elapsed_day=130, protein_regions=None)

end_time2 = time.time()



print("Time spent updating",  end_time2 - start_time )
