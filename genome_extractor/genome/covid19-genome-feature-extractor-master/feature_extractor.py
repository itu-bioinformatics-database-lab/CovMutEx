import json
import os
import random
from Bio import Phylo
import pandas as pd
from configs import configs
from sklearn.utils import resample
import time
import numpy as np
import pickle
import re
from typing import List, Dict, Tuple, Any
import diskcache as dc 


# Paths to files
codon_mapping_path = r'C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\codon_aa_mapping.json'
genome_txt_path = r"C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\genome.txt"
mutations_txt_path = r"C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\mutations.txt"
phylo_tree_path = r"C:\Users\DRX\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\phylogenetic_tree.nwk"

# Initialize cache
cache_dir = 'feature_cache'
if not os.path.exists(cache_dir):
    print("Created cache dir")
    os.makedirs(cache_dir)

# Initialize cache
cache = dc.Cache(cache_dir)

def get_sequence():
    """Read the genome sequence from the text file."""
    with open(genome_txt_path) as file:
        next(file)
        genome = file.read().replace('\n', '')
    genome += 'A'
    return genome

def parse_mutations(mutations_txt_path, nodeId):
    """Parse mutations from mutations.txt file."""
    mutations = []
    with open('mutations.txt', 'r') as file:
        next(file)
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
    
def normalize_features(features):
    """Normalize feature values to a 0-1 scale."""
    min_val = min(features.values())
    max_val = max(features.values())
    if max_val == min_val:
        return {k: 1.0 for k in features}  # Prevent division by zero
    return {k: (v - min_val) / (max_val - min_val) for k, v in features.items()}

def extract_phylogenetic_features(tree_path):
    """Extract phylogenetic features from the tree."""
    tree = Phylo.read(tree_path, 'newick')
    clade_distances = {}
    for clade in tree.find_clades():
        if clade.name:
            clade_distances[clade.name] = clade.branch_length or 0.0
    normalized_distances = normalize_features(clade_distances)
    diversity_metrics = calculate_phylogenetic_diversity(tree)
    return normalized_distances, diversity_metrics

def calculate_phylogenetic_diversity(tree):
    """Calculate phylogenetic diversity metrics."""
    clade_depths = {clade.name: clade.branch_length or 0.0 for clade in tree.find_clades() if clade.branch_length}
    diversity_metrics = {
        'total_branch_length': sum(clade_depths.values()),
        'max_branch_length': max(clade_depths.values(), default=0),
        'min_branch_length': min(clade_depths.values(), default=0),
    }
    return diversity_metrics

def translate_nucleotides_to_amino_acids(nucleotide_sequence, codon_mapper):
    """Translate nucleotide sequence to amino acid sequence."""
    aa_seq = ""
    for i in range(0, len(nucleotide_sequence), 3):
        codon = nucleotide_sequence[i:i + 3]
        if codon in codon_mapper:
            aa_seq += codon_mapper[codon]
    return aa_seq

def find_protein_region(index, protein_regions):
    """Find the protein region for a given index."""
    for protein, (start, end) in protein_regions.items():
        if start <= index <= end:
            return protein
    return "Non-coding"

def compute_nucleotide_frequencies(sequence, window_size=10):
    """Compute nucleotide frequencies in a sliding window."""
    freq_features = []
    for i in range(len(sequence)):
        start = max(0, i - window_size // 2)
        end = min(len(sequence), i + window_size // 2)
        window = sequence[start:end]
        freqs = {nt: window.count(nt) / len(window) for nt in 'ATGC'}
        freq_features.append(freqs)
    return freq_features

def is_synonymous(current_aa, new_aa):
    """Check if mutation is synonymous."""
    return int(current_aa == new_aa)

def build_feature_vector(
    idx: int,
    window: str,
    center_nucleotide: str,
    nucleotide: str,
    aa_idx: int,
    aa_nucleotides: List[str],
    aa_seq: str,
    mutation_dict: Dict[int, List[Any]],
    codon_mapper: Dict[str, str],
    config_file: Dict[str, Any],
    nucleotide_frequencies: List[Dict[str, float]],
    phylo_features: Dict[str, float],
    phylo_diversity: Dict[str, float],
    protein_regions: Dict[str, Tuple[int, int]],
    required_vector_length: int,
    elapsed_day: int,
    node_id: str = None
) -> List[float]:
    """Build feature vector for a specific position and nucleotide."""
    current_aa_nucleotides = aa_nucleotides.copy()
    mutation_position = idx % 3
    current_aa_nucleotides[mutation_position] = nucleotide
    
    try:
        new_aa = codon_mapper.get(''.join(current_aa_nucleotides), "X")
        current_aa = aa_seq[aa_idx] if aa_idx < len(aa_seq) else "X"
        mutation_observed = int(mutation_dict.get(idx, [None, None])[1] == nucleotide)
        
        feature_vector = []
        feature_vector.extend(list(window))  # k-mer nucleotides
        feature_vector.extend([
            center_nucleotide,
            nucleotide,
            idx,
            float(config_file['nucleotide sub. matrix'].get(center_nucleotide, {}).get(nucleotide, 0)),
            current_aa,
            new_aa,
            float(config_file['AA PAM matrix'].get(current_aa, {}).get(new_aa, 0)),
            elapsed_day,
            phylo_features.get(node_id, 0) if node_id else 0,
            phylo_diversity['total_branch_length'],
            is_synonymous(current_aa, new_aa),
            find_protein_region(idx, protein_regions)
        ])
        
        # Add nucleotide frequencies
        nucleotide_freq = nucleotide_frequencies[idx]
        feature_vector.extend(list(nucleotide_freq.values()))
        
        # Add amino acid features
        aa_features = config_file['AA features']
        for feature in ['hydrophobicity', 'polarity', 'iso-electric point', 'volume', 
                       'molecular weight', 'pKa', 'pKb', 'pKx', 'pl']:
            feature_vector.extend([
                float(aa_features[feature].get(current_aa, 0)),
                float(aa_features[feature].get(new_aa, 0))
            ])

        # Ensure correct vector length
        if len(feature_vector) < required_vector_length:
            feature_vector.extend([0.0] * (required_vector_length - len(feature_vector)))
        elif len(feature_vector) > required_vector_length:
            feature_vector = feature_vector[:required_vector_length]

        feature_vector.append(mutation_observed)
        return feature_vector
        
    except Exception as e:
        print(f"Error building feature vector at position {idx}: {e}")
        return None

def construct_variant_genome(genome_seq, mutations):
    """Construct variant genome by applying mutations."""
    variant_genome = list(genome_seq)
    applied_mutations = []

    for mutation in mutations:
        position, original, new = mutation 
        if variant_genome[position - 1] == original:  # Adjusted position to 0-based index
            variant_genome[position - 1] = new
            applied_mutations.append((position, original, new))

    print(f"Applied {len(applied_mutations)} mutations to the reference genome")
    return "".join(variant_genome)

def initialize_feature_vector_cache(
    genome_seq: str, 
    codon_mapping_path: str, 
    phylo_tree_path: str, 
    config_file: Dict[str, Any], 
    k: int = 30
) -> Dict[str, Any]:
    """
    Initialize feature vector caching system for COVID-19 genome variants.
    """
    # Load codon mapping
    with open(codon_mapping_path) as f:
        codon_mapper = json.load(f)
    
    # Precompute static features
    mid_point = k // 2
    required_vector_length = k + 29
    nucleotides = ['A', 'T', 'G', 'C']
    
    reference_aa_seq = translate_nucleotides_to_amino_acids(genome_seq, codon_mapper)
    nucleotide_frequencies = compute_nucleotide_frequencies(genome_seq)
    phylo_features, phylo_diversity = extract_phylogenetic_features(phylo_tree_path)
    
    # Precompute reference feature vectors
    reference_vectors = precompute_reference_vectors(
        genome_seq, 
        codon_mapper, 
        config_file, 
        reference_aa_seq, 
        nucleotide_frequencies, 
        phylo_features, 
        phylo_diversity, 
        k
    )
    
    return {
        'reference_genome': genome_seq,
        'codon_mapper': codon_mapper,
        'reference_aa_seq': reference_aa_seq,
        'nucleotide_frequencies': nucleotide_frequencies,
        'phylo_features': phylo_features,
        'phylo_diversity': phylo_diversity,
        'reference_vectors': reference_vectors,
        'k': k,
        'config_file': config_file
    }

#Cache this as well 

"""
This part is not cached currently so it computes on every call
"""
def precompute_reference_vectors(
    genome_seq: str,
    codon_mapper: Dict[str, str],
    config_file: Dict[str, Any],
    reference_aa_seq: str,
    nucleotide_frequencies: List[Dict[str, float]],
    phylo_features: Dict[str, float],
    phylo_diversity: Dict[str, float],
    k: int
) -> Dict[int, Dict[str, List[float]]]:
    """
    Precompute feature vectors for reference genome.
    """
    mid_point = k // 2
    nucleotides = ['A', 'T', 'G', 'C']
    required_vector_length = k + 29
    
    padded_seq = "-" * mid_point + genome_seq + "-" * mid_point
    feature_vectors = {}
    
    for idx in range(len(genome_seq)):
        window = padded_seq[idx:idx + k]
        center_nucleotide = genome_seq[idx]
        aa_idx = idx // 3
        aa_nucleotides = list(genome_seq[aa_idx * 3: aa_idx * 3 + 3])
        
        position_vectors = {}
        for nucleotide in nucleotides:
            vector = build_feature_vector(
                idx, window, center_nucleotide, nucleotide,
                aa_idx, aa_nucleotides, reference_aa_seq, {},
                codon_mapper, config_file, nucleotide_frequencies,
                phylo_features, phylo_diversity, 
                config_file['protein regions'],
                required_vector_length, 0
            )
            if vector is not None:
                position_vectors[nucleotide] = vector
                
        if position_vectors:
            feature_vectors[idx] = position_vectors
    
    return feature_vectors


def compute_variant_diff_vectors(
    feature_cache: Dict[str, Any],
    mutations: List[Tuple[int, str, str]],
    nodeId: str
) -> Dict[int, Dict[str, List[float]]]:
    """
    Compute and cache variant-specific feature vectors.
    """
    # Normalize nodeId to ensure consistent cache lookup
    nodeId = nodeId.strip().lower()

    # Check if result for this nodeId exists in the cache
    cached_result = cache.get(nodeId)
    if cached_result:
        print(f"Cache hit for nodeId: {nodeId}")
        return cached_result

    print(f"Cache miss for nodeId: {nodeId}. Computing...")

    # Construct the variant genome
    variant_genome = construct_variant_genome(feature_cache['reference_genome'], mutations)

    # Compute diff vectors for each mutation
    diff_vectors = {}
    for mutation in mutations:
        position, original, new = mutation
        if position not in diff_vectors:
            diff_vectors[position] = compute_mutation_impact_vectors(
                variant_genome, position, original, new, feature_cache
            )

    # Save the computed diff vectors to the cache
    cache[nodeId] = diff_vectors
    print(f"Cache updated for nodeId: {nodeId}")
    return diff_vectors


def compute_mutation_impact_vectors(
    variant_genome: str, 
    mutation_pos: int, 
    original: str, 
    new: str, 
    cache: Dict[str, Any]
) -> Dict[int, Dict[str, List[float]]]:
    """
    Compute feature vectors for positions impacted by a specific mutation.
    """
    mid_point = cache['k'] // 2
    padded_seq = "-" * mid_point + variant_genome + "-" * mid_point
    
    # Determine affected region
    start = max(0, mutation_pos - cache['k'])
    end = min(len(variant_genome), mutation_pos + cache['k'] + 1)
    
    mutation_vectors = {}
    for idx in range(start, end):
        window = padded_seq[idx:idx + cache['k']]
        center_nucleotide = variant_genome[idx]
        aa_idx = idx // 3
        aa_nucleotides = list(variant_genome[aa_idx * 3: aa_idx * 3 + 3])
        
        position_vectors = {}
        for nucleotide in cache['reference_vectors'][0].keys():
            ref_vector = cache['reference_vectors'][idx][nucleotide].copy()
            # Only update features for mutated nucleotide position
            if idx == mutation_pos:
                vector = build_feature_vector(
                    idx, window, center_nucleotide, nucleotide,
                    aa_idx, aa_nucleotides, cache['reference_aa_seq'],
                    {mutation_pos: (original, new)},  # Only the mutation for this index
                    cache['codon_mapper'], cache['config_file'], 
                    cache['nucleotide_frequencies'],
                    cache['phylo_features'], cache['phylo_diversity'], 
                    cache['config_file']['protein regions'],
                    cache['k'] + 29, 
                    0  # Elapsed days (placeholder)
                )
                if vector is not None:
                    # Mark as mutation observed
                    vector[-1] = 1
                    position_vectors[nucleotide] = vector
        
        if position_vectors:
            mutation_vectors[idx] = position_vectors
    
    return mutation_vectors

import cProfile
import pstats

def main():
    
    profiler = cProfile.Profile()
    profiler.enable()

    # Example usage
    genome_seq = get_sequence()  # Get genome sequence
    nodeId = "Scotland/QEUH-1F49661/2021|2021-09-20"
    mutations = parse_mutations(mutations_txt_path, nodeId)

    # Initialize feature cache
    feature_cache = initialize_feature_vector_cache(
        genome_seq, 
        codon_mapping_path, 
        phylo_tree_path, 
        configs()
    )

    # Compute or retrieve cached diff vectors
    variant_diff_vectors = compute_variant_diff_vectors(
        feature_cache, 
        mutations, 
        nodeId
    )

    # Optional: Persist diff vectors to a file
    # with open('variant_diff_feature_vectors.pkl', 'wb') as f:
    #     pickle.dump(variant_diff_vectors, f)
    cache = dc.Cache(cache_dir, statistics=True)
    profiler.disable()
    stats = pstats.Stats(profiler)
    stats.sort_stats('cumulative').print_stats(20)
    
if __name__ == "__main__":
    main()