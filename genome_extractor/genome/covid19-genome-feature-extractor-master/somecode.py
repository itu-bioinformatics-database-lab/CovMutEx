
# def extract_features_for_each_position(nodeIds, elapsedDay, selectedProteinRegions):
#     k = 30  # Size of the sliding window
#     mid_point = k // 2
#     nucleotides = ['A', 'T', 'G', 'C']
#     required_vector_length = k + 29

#     config_file = configs()
#     protein_regions = config_file['protein regions']

#     with open(codon_mapping_path) as codon_mapping_file:
#         codon_mapper = json.load(codon_mapping_file)

#     genome_seq = get_sequence()
#     mutations = parse_mutations()  # Use the new mutation parsing function
#     synthetic_mutations = generate_synthetic_mutations(genome_seq)
#     mutations.extend(synthetic_mutations)

#     phylo_features, phylo_diversity = extract_phylogenetic_features(phylo_tree_path)
#     mutation_dict = {mutation[1]: mutation for mutation in mutations}  # Use the nucleotide position as key
    
#     # Precompute feature vectors for reference genome
#     precompute_feature_vectors(genome_seq, 'reference_feature_vectors.pkl')
    
#     # Precompute diff feature vectors for mutations
#     precompute_diff_feature_vectors(genome_seq, mutations, k, 'variant_diff_feature_vectors.pkl')
    
#     # Pad genome sequence with '-' to maintain consistent k-mers
#     padded_seq = "-" * mid_point + genome_seq + "-" * mid_point
#     aa_seq = translate_nucleotides_to_amino_acids(genome_seq, codon_mapper)
#     nucleotide_frequencies = compute_nucleotide_frequencies(genome_seq)

#     dataset = []
    
#     for idx in range(len(genome_seq)):
#         window = padded_seq[idx:idx + k]
#         center_nucleotide = genome_seq[idx]
#         aa_idx = idx // 3
#         aa_nucleotides = list(genome_seq[aa_idx * 3: aa_idx * 3 + 3])

#         for nucleotide in nucleotides:
#             current_aa_nucleotides = aa_nucleotides.copy()
#             mutation_position = idx % 3
#             current_aa_nucleotides[mutation_position] = nucleotide
            
#             try:
#                 new_aa = codon_mapper.get(''.join(current_aa_nucleotides), "X")
#                 current_aa = aa_seq[aa_idx] if aa_idx < len(aa_seq) else "X"
#                 mutation_observed = int(mutation_dict.get(idx, [None, None])[1] == nucleotide)
                
#                 # Build the feature vector
#                 feature_vector = []
#                 feature_vector.extend(list(window))  # 1..k: Nucleotides in the k-mer
#                 feature_vector.extend([
#                     center_nucleotide,
#                     nucleotide,
#                     idx,
#                     float(config_file['nucleotide sub. matrix'].get(center_nucleotide, {}).get(nucleotide, 0)),
#                     current_aa,
#                     new_aa,
#                     float(config_file['AA PAM matrix'].get(current_aa, {}).get(new_aa, 0)),
#                     elapsedDay,
#                     phylo_features.get(nodeIds[min(idx, len(nodeIds) - 1)], 0),  # Normalized
#                     phylo_diversity['total_branch_length'],
#                     is_synonymous(current_aa, new_aa),
#                     find_protein_region(idx, protein_regions)
#                 ])
                
#                 # Add nucleotide frequencies
#                 nucleotide_freq = nucleotide_frequencies[idx]
#                 feature_vector.extend(list(nucleotide_freq.values()))
                
#                 # Add amino acid features
#                 aa_features = config_file['AA features']
#                 feature_vector.extend([
#                     float(aa_features['hydrophobicity'].get(current_aa, 0)),
#                     float(aa_features['hydrophobicity'].get(new_aa, 0)),
#                     float(aa_features['polarity'].get(current_aa, 0)),
#                     float(aa_features['polarity'].get(new_aa, 0)),
#                     float(aa_features['iso-electric point'].get(current_aa, 0)),
#                     float(aa_features['iso-electric point'].get(new_aa, 0)),
#                     float(aa_features['volume'].get(current_aa, 0)),
#                     float(aa_features['volume'].get(new_aa, 0)),
#                     float(aa_features['molecular weight'].get(current_aa, 0)),
#                     float(aa_features['molecular weight'].get(new_aa, 0)),
#                     float(aa_features['pKa'].get(current_aa, 0)),
#                     float(aa_features['pKa'].get(new_aa, 0)),
#                     float(aa_features['pKb'].get(current_aa, 0)),
#                     float(aa_features['pKb'].get(new_aa, 0)),
#                     float(config_file['AA features']['pKx'].get(current_aa)),
#                     float(config_file['AA features']['pKx'].get(new_aa)),
#                     float(config_file['AA features']['pl'].get(current_aa)),
#                     float(config_file['AA features']['pl'].get(new_aa))
#                 ])

#                 # Ensure correct vector length
#                 if len(feature_vector) < required_vector_length:
#                     feature_vector.extend([0.0] * (required_vector_length - len(feature_vector)))
#                 elif len(feature_vector) > required_vector_length:
#                     feature_vector = feature_vector[:required_vector_length]

#                 # Add target label to the feature vector
#                 feature_vector.append(mutation_observed)  # Target (1 for observed, 0 otherwise)

#                 dataset.append(feature_vector)

#             except Exception as e:
#                 print(f"Error at position {idx} with nucleotide {nucleotide}: {e}")
#                 continue

#     print("Extracted dataset size:", len(dataset))
#     pr_poss = selectedProteinRegions if selectedProteinRegions is not None else None
#     end_time = time.time()
#     elapsed_time = end_time - start_time
#     print(f"Runtime for Feature Extraction: {elapsed_time:.2f} seconds")
    
#     return dataset, genome_seq, pr_poss


# def parse_mutations():
#     mutations = []
#     mutation_pattern = re.compile(r'(\S+)\s+(\d+)\s+(\w+)>(\w+)\s+'
#                                    r'(\d+)\s+(\w+)>(\w+)\s+(\S+)\s+'
#                                    r'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)')

#     with open(mutations_txt_path, 'r') as file:
#         next(file)
#         for line in file:
#             match = mutation_pattern.match(line.strip())
#             if match:
#                 # Extract all relevant fields
#                 protein_region = match.group(1)  # Protein or gene name
#                 nt_position = int(match.group(2))  # Nucleotide position
#                 original_nt = match.group(3)  # Original nucleotide
#                 new_nt = match.group(4)  # New nucleotide
#                 aa_position = int(match.group(5))  # Amino acid position
#                 original_aa = match.group(6)  # Original amino acid
#                 new_aa = match.group(7)  # New amino acid
#                 date = match.group(8)  # Date of mutation
#                 country = match.group(9)  # Country
#                 lineage = match.group(10)  # Lineage
#                 variant = match.group(11)  # Variant
#                 accession_info = match.group(12)  # Accession info
#                 species = match.group(13)  # Species
#                 confidence = int(match.group(14))  # Confidence score

#                 mutations.append((protein_region, nt_position, original_nt, new_nt,
#                                   aa_position, original_aa, new_aa, date,
#                                   country, lineage, variant, accession_info,
#                                   species, confidence))

#     return mutations


#LOAD THE PICKLES

# import pickle

# def load_diff_feature_vectors(file_path):
#     """Loads the precomputed diff feature vectors from the given pickle file."""
#     with open(file_path, 'rb') as f:
#         diff_feature_vectors = pickle.load(f)
#     return diff_feature_vectors

# def get_diff_for_variant(mutation, diff_feature_vectors):
#     """
#     Retrieve the differential feature vectors for a given mutation.
#     The mutation is expected to be a tuple (position, original_nucleotide, new_nucleotide).
#     """
#     position, original, new = mutation
    
#     # Check if the position exists in the diff_feature_vectors dictionary
#     if position in diff_feature_vectors:
#         position_vectors = diff_feature_vectors[position]
        
#         # Now, we need to find the nucleotide vector for the given mutation
#         # For example, we can search for the new nucleotide in the vectors.
#         if new in position_vectors:
#             return position_vectors[new]  # Return the feature vector for the mutated nucleotide
#         else:
#             print(f"No feature vectors found for the nucleotide {new} at position {position}")
#             return None
#     else:
#         print(f"No diff vectors found for position {position}")
#         return None

# # Example usage
# mutation = (10, 'A', 'T')  # Example mutation (position 10, original A, new T)
# diff_feature_vectors = load_diff_feature_vectors('variant_diff_feature_vectors.pkl')

# # Retrieve the diff feature vector for the selected mutation
# diff_vector = get_diff_for_variant(mutation, diff_feature_vectors)

# # Print the retrieved diff feature vector
# if diff_vector:
#     print(f"Feature vector for mutation {mutation}: {diff_vector}")

# # Example usage
# node_id = "England/NORW-304C926/2021|2021-09-24"
# depth_file = r"F:\COVID19 MUTATION\genome_extractor\genome\covid19-genome-feature-extractor-master\depth_date.json"
# depth = get_sample_depth(depth_file, node_id)
# # print(f"Depth for {node_id}: {depth}")

# import h5py
# import numpy as np

# CACHE_FILE = "diff_features.h5"
# import h5py

# CACHE_FILE = 'diff_features.h5'

# try:
#     with h5py.File(CACHE_FILE, 'r') as f:
#         print("File opened successfully!")
#         print("Available datasets:", list(f.keys()))
# except OSError as e:
#     print(f"Failed to open file {CACHE_FILE}: {e}")

# def save_to_cache(nodeId, features):
#     """
#     Save features to the cache file for a given node ID.
#     """
#     import h5py

#     # Open or create the cache file
#     with h5py.File(CACHE_FILE, "a") as h5file:
#         if nodeId in h5file:
#             print(f"Overwriting existing cache for nodeId: {nodeId}")
#             del h5file[nodeId]
#         h5file.create_dataset(nodeId, data=features)


# def load_from_cache(node_id):
#     """Load features from the cache."""
#     with h5py.File(CACHE_FILE, "r") as h5file:
#         if node_id in h5file:
#             cached_data = h5file[node_id][()]
#             # Convert back to a dictionary
#             return {int(item[0]): item[1] for item in cached_data}
#         return {}

# def process_variant(
#     genome_seq, 
#     mutations, 
#     precomputed_features, 
#     codon_mapper, 
#     config_file, 
#     protein_regions=None,  
#     k=30,
#     elapsed_day=None,  
#     nodeId="England/NORW-304C926/2021|2021-09-24" 
# ):
#     """
#     Process a variant, dynamically computing features for affected positions,
#     while using cached features for previously computed positions.
#     """
#     mid_point = k // 2
#     padded_seq = "-" * mid_point + genome_seq + "-" * mid_point
#     updated_features = precomputed_features.copy()  # Start with precomputed features
#     nucleotides = ['A', 'T', 'C', 'G']
#     aa_seq = [codon_mapper.get(genome_seq[i:i+3], 'X') for i in range(0, len(genome_seq), 3)]

#     phylo_features, phylo_diversity = extract_phylogenetic_features(phylo_tree_path)

#     # Load cached data if available
#     cached_features = {}
#     if nodeId:
#         cached_features = load_from_cache(nodeId)
#         updated_features.update(cached_features)  # Merge cached features

#     affected_positions_set = set()
#     if protein_regions:
#         for region in protein_regions.values():
#             start, end = region
#             affected_positions_set.update(range(start, end + 1))
#     else:
#         for mutation in mutations:
#             position = mutation[0] - 1
#             codon_start = (position // 3) * 3
#             affected_positions_set.update(range(codon_start, codon_start + 3))
#             affected_positions_set.update(range(max(0, position - mid_point), min(len(genome_seq), position + mid_point + 1)))

#     for idx in sorted(affected_positions_set):
#         if idx in updated_features:
#             continue  # Skip already cached or precomputed features

#         try:
#             # Sliding window extraction
#             window = padded_seq[idx:idx + k]
#             aa_idx = idx // 3
#             codon_start = aa_idx * 3
#             aa_nucleotides = list(genome_seq[codon_start:codon_start + 3])

#             if len(aa_nucleotides) != 3:
#                 continue

#             mutation_position = idx % 3
#             for nucleotide in nucleotides:
#                 aa_nucleotides[mutation_position] = nucleotide
#                 new_codon = ''.join(aa_nucleotides)
#                 new_aa = codon_mapper.get(new_codon, 'X')

#                 # Build feature vector for the mutated position
#                 sample_data = [
#                     *list(window),
#                     window[mid_point],  # Center of k-mer
#                     nucleotide,  # Mutated nucleotide
#                     idx,  # Genome position
#                     config_file['nucleotide sub. matrix'][window[mid_point]][nucleotide],
#                     aa_seq[aa_idx],  # Original amino acid
#                     new_aa,  # Mutated amino acid
#                     config_file['AA PAM matrix'][aa_seq[aa_idx]][new_aa],
#                     int(aa_seq[aa_idx] == new_aa),  # Synonymous mutation
#                     elapsed_day,
#                     find_protein_region(idx, protein_regions) if protein_regions else None,
#                     # Include additional features as needed...
#                     phylo_features,  # Normalized
#                     phylo_diversity['total_branch_length'],
#                     config_file['AA features']['hydrophobicity'].get(new_aa, 0),
#                     config_file['AA features']['polarity'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['polarity'].get(new_aa, 0),
#                     config_file['AA features']['iso-electric point'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['iso-electric point'].get(new_aa, 0),
#                     config_file['AA features']['volume'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['volume'].get(new_aa, 0),
#                     config_file['AA features']['molecular weight'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['molecular weight'].get(new_aa, 0),
#                     config_file['AA features']['pKa'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['pKa'].get(new_aa, 0),
#                     config_file['AA features']['pKb'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['pKb'].get(new_aa, 0),
#                     config_file['AA features']['pKx'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['pKx'].get(new_aa, 0),
#                     config_file['AA features']['pl'].get(aa_seq[aa_idx], 0),
#                     config_file['AA features']['pl'].get(new_aa, 0),
#                 ]
                
#                 updated_features[idx] = preprocess_input(sample_data, expected_size=205)
#         except Exception as e:
#             print(f"Error processing position {idx}: {e}")
#             continue

#     # Save updated features to cache
#     if nodeId:
#         save_to_cache(nodeId, updated_features)

#     return updated_features
