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
from sklearn.preprocessing import OneHotEncoder, StandardScaler
import tensorflow as tf


# Paths to files
base_dir = os.path.dirname(os.path.abspath(__file__))

codon_mapping_path = os.path.join(base_dir, "codon_aa_mapping.json")
genome_txt_path = os.path.join(base_dir, "genome.txt")
mutations_txt_path = os.path.join(base_dir, "mutations.txt")
phylo_tree_path = os.path.join(base_dir, "phylogenetic_tree.nwk")

depth_file = os.path.join(base_dir, 'depth_date.json')
genome_extractor_dir = os.path.dirname(base_dir)
ROOT_PATH = os.path.dirname(genome_extractor_dir)
model_path = os.path.join(genome_extractor_dir, 'covid19_models', "models", "balanced_data_model.keras")

print(f"Base dir: {base_dir}")
print(f"Genome extractor dir: {genome_extractor_dir}")
print(f"ROOT_PATH: {ROOT_PATH}")

start_time = time.time()

#
# k = 30, so total raw features = 30 + 29 = 59 per data point
#
#  Index | Name                          | Type
# -------|-------------------------------|----------
#  0-29  | k-mer window nucleotides      | Cat (5 values each: A,T,G,C,-)
#  30    | center nucleotide (original)  | Cat
#  31    | mutant nucleotide             | Cat
#  32    | position index                | Num
#  33    | nucleotide PAM250 score       | Num
#  34    | original amino acid           | Cat
#  35    | new amino acid                | Cat
#  36    | AA PAM250 score               | Num
#  37    | elapsed days (k+8)            | Num  ← k+8
#  38    | tree depth    (k+9)           | Num  ← k+9
#  39    | synonymous flag (k+10)        | Num  ← k+10  (1=syn, 0=non-syn)
#  40    | protein region / ORF (k+11)   | Cat  ← k+11
#  41    | hydrophobicity (original AA)  | Num  ← k+12
#  42    | hydrophobicity (new AA)       | Num  ← k+13
#  43    | polarity (original AA)        | Num  ← k+14
#  44    | polarity (new AA)             | Num  ← k+15
#  45    | iso-electric pt (original AA) | Num  ← k+16
#  46    | iso-electric pt (new AA)      | Num  ← k+17
#  47    | volume (original AA)          | Num  ← k+18
#  48    | volume (new AA)               | Num  ← k+19
#  49    | weight (original AA)          | Num  ← k+20
#  50    | weight (new AA)               | Num  ← k+21
#  51    | pKa (original AA)             | Num  ← k+22
#  52    | pKa (new AA)                  | Num  ← k+23
#  53    | pKb (original AA)             | Num  ← k+24
#  54    | pKb (new AA)                  | Num  ← k+25
#  55    | pKx (original AA)             | Num  ← k+26
#  56    | pKx (new AA)                  | Num  ← k+27
#  57    | pl (original AA)              | Num  ← k+28
#  58    | pl (new AA)                   | Num  ← k+29
#
# Total raw: 59 features
# After one-hot encoding (Cat → 5 vals each for nucs, ~20 for AAs, ~13 for ORF)
# + standardisation of Num → pads/truncates to 205
#
# CATEGORICAL column indices in raw vector: [0..31, 34, 35, 40]
# NUMERICAL  column indices in raw vector: [32, 33, 36, 37, 38, 39, 41..58]
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_sequence():
    with open(genome_txt_path) as file:
        next(file)
        genome = file.read().replace('\n', '')
    genome += 'A'
    return genome


def parse_mutations(nodeId):
    mutations = []
    with open(mutations_txt_path, 'r') as file:
        file.readline()  # skip header
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


# ---------------------------------------------------------------------------
# Phylogenetic helpers
# ---------------------------------------------------------------------------

CACHE_P_FILE = 'phylo_features_cache.h5'


def normalize_features(features):
    min_val = min(features.values())
    max_val = max(features.values())
    if max_val == min_val:
        return {k: 1.0 for k in features}
    return {k: (v - min_val) / (max_val - min_val) for k, v in features.items()}


def extract_phylogenetic_features(tree_path):
    normalized_distances = None
    diversity_metrics = {}

    try:
        with h5py.File(CACHE_P_FILE, 'a') as f:
            if tree_path not in f:
                print(f"Extracting phylogenetic features for {tree_path}...")
                try:
                    tree = Phylo.read(tree_path, 'newick')
                    clade_distances = {}
                    for clade in tree.find_clades():
                        if clade.name:
                            clade_distances[clade.name] = clade.branch_length or 0.0
                    if not clade_distances:
                        raise ValueError("No valid clade distances found in the tree")
                    normalized_distances = normalize_features(clade_distances)
                    diversity_metrics = calculate_phylogenetic_diversity(tree)
                    tree_group = f.create_group(tree_path)
                    tree_group.create_dataset('clade_names', data=list(normalized_distances.keys()))
                    tree_group.create_dataset('normalized_distances', data=list(normalized_distances.values()))
                    diversity_names = list(diversity_metrics.keys())
                    diversity_values = list(diversity_metrics.values())
                    tree_group.create_dataset(
                        'diversity_metrics_names',
                        data=[name.encode('utf-8') if isinstance(name, str) else name
                              for name in diversity_names]
                    )
                    tree_group.create_dataset('diversity_metrics', data=diversity_values)
                except Exception as e:
                    print(f"Error extracting features from {tree_path}: {e}")
                    raise

            tree_group = f[tree_path]
            clade_names = list(tree_group['clade_names'][:])
            normalized_distances = dict(zip(clade_names, tree_group['normalized_distances'][:]))
            diversity_names = [
                name.decode('utf-8') if isinstance(name, bytes) else name
                for name in tree_group['diversity_metrics_names'][:]
            ]
            diversity_values = tree_group['diversity_metrics'][:]
            diversity_metrics = dict(zip(diversity_names, diversity_values))
            print(f"Successfully retrieved/cached phylogenetic features for {tree_path}")

    except Exception as e:
        print(f"Critical error processing tree {tree_path}: {e}")
        raise

    if normalized_distances is None:
        raise ValueError("Failed to compute or retrieve normalized distances")

    return normalized_distances, diversity_metrics


def calculate_phylogenetic_diversity(tree):
    clade_depths = {
        clade.name: clade.branch_length or 0.0
        for clade in tree.find_clades()
        if clade.branch_length
    }
    return {
        'total_branch_length': sum(clade_depths.values()),
        'max_branch_length': max(clade_depths.values(), default=0),
        'min_branch_length': min(clade_depths.values(), default=0),
    }


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

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


def construct_variant_genome(genome_seq, mutations):
    variant_genome = list(genome_seq)
    applied_mutations = []
    for mutation in mutations:
        position, original, new, _, _ = mutation
        if variant_genome[position - 1] == original:
            variant_genome[position - 1] = new
            applied_mutations.append((position, original, new))
    print(f"Applied {len(applied_mutations)} mutations to the reference genome")
    return "".join(variant_genome)


def balance_dataset(dataset):
    data_df = pd.DataFrame(dataset)
    majority_class = data_df[data_df.iloc[:, -1] == 0]
    minority_class = data_df[data_df.iloc[:, -1] == 1]
    minority_oversampled = resample(
        minority_class, replace=True, n_samples=len(majority_class), random_state=42
    )
    balanced_dataset = pd.concat([majority_class, minority_oversampled])
    return balanced_dataset.values.tolist()


# ---------------------------------------------------------------------------
# AA biochemical features helper  
# ---------------------------------------------------------------------------

def get_aa_features(original_aa, new_aa, config_file):
    """
    Returns 18 numerical AA biochemical features in the order defined by the
    paper (k+12 .. k+29):
        hydrophobicity (orig, new)
        polarity       (orig, new)
        iso-electric   (orig, new)
        volume         (orig, new)
        weight         (orig, new)
        pKa            (orig, new)
        pKb            (orig, new)
        pKx            (orig, new)
        pl             (orig, new)
    """
    props = config_file['AA features']
    return [
        float(props['hydrophobicity'].get(original_aa, 0)),   # k+12
        float(props['hydrophobicity'].get(new_aa, 0)),         # k+13
        float(props['polarity'].get(original_aa, 0)),          # k+14
        float(props['polarity'].get(new_aa, 0)),               # k+15
        float(props['iso-electric point'].get(original_aa, 0)),# k+16
        float(props['iso-electric point'].get(new_aa, 0)),     # k+17
        float(props['volume'].get(original_aa, 0)),            # k+18
        float(props['volume'].get(new_aa, 0)),                 # k+19
        float(props['molecular weight'].get(original_aa, 0)),  # k+20
        float(props['molecular weight'].get(new_aa, 0)),       # k+21
        float(props['pKa'].get(original_aa, 0)),               # k+22
        float(props['pKa'].get(new_aa, 0)),                    # k+23
        float(props['pKb'].get(original_aa, 0)),               # k+24
        float(props['pKb'].get(new_aa, 0)),                    # k+25
        float(props['pKx'].get(original_aa, 0)),               # k+26
        float(props['pKx'].get(new_aa, 0)),                    # k+27
        float(props['pl'].get(original_aa, 0)),                # k+28
        float(props['pl'].get(new_aa, 0)),                     # k+29
    ]


# ---------------------------------------------------------------------------
# Preprocessing 
# ---------------------------------------------------------------------------

def preprocess_matrix(raw_data_list, expected_size=205):
    """
    Bulk preprocessing matching the training pipeline:
      - Categorical columns → One-Hot Encoded
      - Numerical columns   → Z-score Standardised
    All statistics are derived from the ENTIRE batch so every row lives in the
    same feature space.

    Column layout (indices into the 59-element raw vector):
      Categorical: 0-31 (k-mer window + center + mutant), 34, 35, 40 (AAs + ORF)
      Numerical:   32, 33, 36, 37, 38, 39, 41-58
    """
    df = pd.DataFrame(raw_data_list)

    # -----------------------------------------------------------------------
    # Categorical indices:
    #   0-29  : k-mer nucleotides
    #   30    : center nucleotide
    #   31    : mutant nucleotide
    #   34    : original amino acid
    #   35    : new amino acid
    #   40    : protein region / ORF name
    # -----------------------------------------------------------------------
    cat_indices = list(range(32)) + [34, 35, 40]
    known_categories = (
    [['A', 'T', 'G', 'C', '-']] * 32 +
    [['A','C','D','E','F','G','H','I','K','L',
      'M','N','P','Q','R','S','T','V','W','Y','X']] * 2 +
    [['ORF1ab-1','ORF1ab-2','S','ORF3a','E','M',
      'ORF6','ORF7a','ORF7b','ORF8','N','ORF10','Non-coding']]
)

    encoder = OneHotEncoder(
        sparse_output=False,
        handle_unknown='ignore',
        categories=known_categories
    )

    # -----------------------------------------------------------------------
    # Numerical indices:
    #   32    : position index
    #   33    : nucleotide PAM score
    #   36    : AA PAM score
    #   37    : elapsed days      ← k+8
    #   38    : tree depth        ← k+9
    #   39    : synonymous flag   ← k+10
    #   41-58 : AA biochemical features (18 values)
    # -----------------------------------------------------------------------
    num_indices = [32, 33, 36, 37, 38, 39] + list(range(41, 59))

    # encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
    categorical_data = df.iloc[:, cat_indices].astype(str).values
    encoded_cats = encoder.fit_transform(categorical_data)

    scaler = StandardScaler()
    numerical_data = df.iloc[:, num_indices].astype(float).values
    standardized_nums = scaler.fit_transform(numerical_data)

    combined = np.hstack([encoded_cats, standardized_nums])

    if combined.shape[1] < expected_size:
        padding = np.zeros((combined.shape[0], expected_size - combined.shape[1]))
        combined = np.hstack([combined, padding])
    else:
        combined = combined[:, :expected_size]

    return combined.astype(np.float32)


# ---------------------------------------------------------------------------
# Reference genome feature extraction  (stores RAW lists)
# ---------------------------------------------------------------------------

def precompute_feature_vectors(
    genome_txt_path, codon_mapping_path, phylo_tree_path, protein_regions,
    config_file, depth=0, elapsed_day=0, k=30
):
    """
    Build RAW feature vectors for every position in the reference genome.
    Returns a list of raw lists – no encoding / scaling applied here.

    For the reference genome (no variant): mutant = original, elapsed_day = 0,
    depth = 0, synonymous = 1 (no change).
    """
    start = time.time()
    mid_point = k // 2

    config_file = configs()
    genome_seq = get_sequence()
    with open(codon_mapping_path) as codon_file:
        codon_mapper = json.load(codon_file)

    aa_seq = translate_nucleotides_to_amino_acids(genome_seq, codon_mapper)
    padded_seq = "-" * mid_point + genome_seq + "-" * mid_point

    dataset = []

    for idx in range(len(genome_seq)):
        try:
            aa_idx = idx // 3
            if (aa_idx * 3) + 3 > len(genome_seq):
                continue

            window = list(padded_seq[idx:idx + k])
            center_nuc = genome_seq[idx]
            current_aa = aa_seq[aa_idx] if aa_idx < len(aa_seq) else "X"
            protein_reg = find_protein_region(idx, protein_regions) if protein_regions else "Non-coding"

            aa_feats = get_aa_features(current_aa, current_aa, config_file)

            # Nucleotide PAM: self-mutation score
            nuc_pam = config_file['nucleotide sub. matrix'].get(center_nuc, {}).get(center_nuc, 0)
            # AA PAM: self score
            aa_pam = config_file['AA PAM matrix'].get(current_aa, {}).get(current_aa, 0)

            feature_vector = [
                *window,          # 0-29  (Cat)
                center_nuc,       # 30    (Cat)  original
                center_nuc,       # 31    (Cat)  mutant = same (no mutation)
                idx,              # 32    (Num)  position
                nuc_pam,          # 33    (Num)  nucleotide PAM
                current_aa,       # 34    (Cat)  original AA
                current_aa,       # 35    (Cat)  new AA = same
                aa_pam,           # 36    (Num)  AA PAM
                float(elapsed_day),  # 37  (Num)  elapsed days       ← k+8
                float(depth),        # 38  (Num)  tree depth         ← k+9
                1,                   # 39  (Num)  synonymous=1       ← k+10
                protein_reg,         # 40  (Cat)  ORF name           ← k+11
                *aa_feats,           # 41-58 (Num) biochemical (18)
            ]

            dataset.append(feature_vector)

        except Exception as e:
            print(f"Error at position {idx}: {e}")
            continue

    print(f"Raw reference feature extraction: {time.time() - start:.2f}s  "
          f"({len(dataset)} positions)")
    return dataset


def precompute_and_cache_ref_features(
    genome_txt_path, codon_mapping_path, phylo_tree_path,
    protein_regions, config_file, output_h5_path, k=30
):
    """
    Build (or load from cache) the preprocessed reference feature matrix.
    Also saves the RAW feature list so variant processing can merge before
    bulk preprocessing.

    HDF5 datasets:
        'features'     – (N, 205) float32, preprocessed
        'features_raw' – (N,)     object,  JSON-serialised raw lists
    """
    if os.path.exists(output_h5_path):
        print(f"Using cached feature vectors from {output_h5_path}")
        with h5py.File(output_h5_path, 'r') as h5f:
            dataset_np = h5f['features'][:]
        return dataset_np

    print("Cached file not found – computing features...")
    dataset_raw = precompute_feature_vectors(
        genome_txt_path, codon_mapping_path, phylo_tree_path,
        protein_regions, config_file, k=k
    )

    print("Applying bulk preprocessing to reference genome...")
    dataset_np = preprocess_matrix(dataset_raw, expected_size=205)

    print("Saving to HDF5...")
    with h5py.File(output_h5_path, 'w') as h5f:
        h5f.create_dataset("features", data=dataset_np, compression="gzip", compression_opts=9)
        h5f.attrs['genome_length'] = len(dataset_raw)
        h5f.attrs['vector_length'] = dataset_np.shape[1]

        raw_encoded = [json.dumps(row).encode('utf-8') for row in dataset_raw]
        dt = h5py.special_dtype(vlen=bytes)
        raw_ds = h5f.create_dataset("features_raw", (len(raw_encoded),), dtype=dt)
        for i, row in enumerate(raw_encoded):
            raw_ds[i] = row

    print(f"Feature vectors cached to {output_h5_path}")
    return dataset_np


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

cache_path = os.path.join(ROOT_PATH, 'node_features.h5')
os.makedirs(os.path.dirname(cache_path) if os.path.dirname(cache_path) else '.', exist_ok=True)

if not os.path.exists(cache_path):
    print(f"Creating new cache file at {cache_path}")
    with h5py.File(cache_path, 'w') as hdf:
        pass

precompute_and_cache_ref_features(
    genome_txt_path=genome_txt_path,
    codon_mapping_path=codon_mapping_path,
    phylo_tree_path=phylo_tree_path,
    protein_regions={},
    config_file=configs,
    output_h5_path='features.h5',
    k=30
)


# ---------------------------------------------------------------------------
# Variant processing helpers
# ---------------------------------------------------------------------------

def get_affected_positions(mutations, genome_len, k, protein_regions=None):
    affected_positions = set()
    mid_point = k // 2

    for position, ref, alt, aa_pos, aa_change in mutations:
        pos = position - 1
        if protein_regions:
            for region_start, region_end in protein_regions.values():
                if region_start <= pos <= region_end:
                    break
            else:
                continue
        window_start = max(0, pos - mid_point)
        window_end = min(genome_len, pos + mid_point + 1)
        for idx in range(window_start, window_end):
            if not protein_regions or any(
                rs <= idx <= re for rs, re in protein_regions.values()
            ):
                affected_positions.add(idx)

    return sorted(affected_positions)


def get_sample_depth(depth_file_path, nodeId):
    try:
        with open(depth_file, 'r') as file:
            depth_data = json.load(file)
        if nodeId in depth_data:
            return depth_data[nodeId].get("depth", 0)
        else:
            print(f"Sample '{nodeId}' not found in the JSON file.")
            return 0
    except FileNotFoundError:
        print(f"File '{depth_file}' not found.")
        return 0
    except json.JSONDecodeError:
        print(f"Error decoding JSON file '{depth_file}'.")
        return 0


def process_variant_raw(
    genome_seq,
    mutations,
    codon_mapper,
    config_file,
    protein_regions=None,
    k=30,
    affected_positions_set=None,
    elapsed_day=0,
    depth=0,
    nodeId=None,
):
    """
    Build RAW feature vectors for positions affected by mutations.
    Returns a dict: {genome_position (int): raw_feature_list}

    Feature order matches the paper exactly (k+8=elapsed, k+9=depth,
    k+10=synonymous, k+11=ORF).
    """
    start = time.time()
    mid_point = k // 2
    padded_seq = "-" * mid_point + genome_seq + "-" * mid_point

    aa_seq = [codon_mapper.get(genome_seq[i:i + 3], 'X')
              for i in range(0, len(genome_seq), 3)]

    if affected_positions_set is None:
        affected_positions_set = set(get_affected_positions(
            mutations, len(genome_seq), k, protein_regions=protein_regions
        ))

    if protein_regions:
        affected_positions_set = {
            idx for idx in affected_positions_set
            if any(rs <= idx <= re for rs, re in protein_regions.values())
        }

    updated_raw = {}

    for idx in affected_positions_set:
        try:
            protein_reg = find_protein_region(idx, protein_regions) if protein_regions else "Non-coding"
            window = list(padded_seq[idx:idx + k])
            aa_idx = idx // 3
            codon_start = aa_idx * 3
            if codon_start + 3 > len(genome_seq):
                continue

            original_codon = list(genome_seq[codon_start:codon_start + 3])
            original_aa = aa_seq[aa_idx]
            center_nuc = genome_seq[idx]

            mutated_codon = original_codon[:]
            # codon stays as-is (variant already applied in genome_seq)
            new_aa = codon_mapper.get(''.join(mutated_codon), 'X')

            nuc_pam = config_file['nucleotide sub. matrix'].get(center_nuc, {}).get(center_nuc, 0)
            aa_pam  = config_file['AA PAM matrix'].get(original_aa, {}).get(new_aa, 0)
            synonymous = int(original_aa == new_aa)
            aa_feats = get_aa_features(original_aa, new_aa, config_file)

            raw_vector = [
                *window,              # 0-29  (Cat)
                center_nuc,           # 30    (Cat)  original
                center_nuc,           # 31    (Cat)  mutant placeholder
                idx,                  # 32    (Num)
                nuc_pam,              # 33    (Num)
                original_aa,          # 34    (Cat)
                new_aa,               # 35    (Cat)
                aa_pam,               # 36    (Num)
                float(elapsed_day),   # 37    (Num)  ← k+8
                float(depth),         # 38    (Num)  ← k+9
                synonymous,           # 39    (Num)  ← k+10
                protein_reg,          # 40    (Cat)  ← k+11
                *aa_feats,            # 41-58 (Num)  18 values
            ]

            updated_raw[idx] = raw_vector

        except Exception as e:
            print(f"Error processing position {idx}: {e}")
            continue

    print(f"Raw variant feature extraction: {time.time() - start:.2f}s  "
          f"({len(updated_raw)} positions)")
    return updated_raw


# ---------------------------------------------------------------------------
# Caching: merge raw ref + raw variant → preprocess once → save
# ---------------------------------------------------------------------------

def cache_precomputed_features(
    cache_path, genome_seq, mutations, codon_mapper, config_file, node_ids,
    elapsed_day=0, depth=0, protein_regions=None
):
    """
    For each node_id:
      1. Load RAW reference feature list from features.h5
      2. Compute RAW variant features for affected positions only
      3. Merge (overwrite affected positions in the raw ref list)
      4. Preprocess the FULL merged raw matrix ONCE with preprocess_matrix()
      5. Save preprocessed result to node_features.h5
      6. Return the preprocessed combined features
    """
    with h5py.File('features.h5', 'r') as h5f:
        raw_encoded = h5f['features_raw'][:]

    base_raw = [json.loads(row.decode('utf-8') if isinstance(row, bytes) else row)
                for row in raw_encoded]
    print(f"Loaded {len(base_raw)} raw reference feature vectors")

    if isinstance(codon_mapper, str):
        with open(codon_mapper) as f:
            codon_mapper_dict = json.load(f)
    else:
        codon_mapper_dict = codon_mapper

    combined_features = None

    for node_id in node_ids:
        node_id_str = str(node_id)

        with h5py.File(cache_path, 'a') as hdf:
            if node_id_str in hdf:
                print(f"Loading cached features for NodeId {node_id}")
                combined_np = hdf[node_id_str][:]
                if protein_regions:
                    region_slices = []
                    for _, (rs, re) in protein_regions.items():
                        region_slices.append(combined_np[rs:re + 1])
                    combined_np = np.concatenate(region_slices, axis=0)
                combined_features = np.nan_to_num(combined_np.astype(np.float32))
                continue

        print(f"Processing new mutations for NodeId {node_id}")
        affected_positions = sorted({m[0] - 1 for m in mutations})

        if affected_positions:
            variant_raw = process_variant_raw(
                genome_seq=genome_seq,
                mutations=mutations,
                codon_mapper=codon_mapper_dict,
                config_file=config_file,
                protein_regions=protein_regions,
                affected_positions_set=set(affected_positions),
                elapsed_day=elapsed_day,
                depth=depth,
            )

            merged_raw = list(base_raw)
            for pos, raw_vec in variant_raw.items():
                merged_raw[pos] = raw_vec
        else:
            merged_raw = base_raw

        print(f"Bulk preprocessing {len(merged_raw)} merged raw vectors...")
        combined_np = preprocess_matrix(merged_raw, expected_size=205)

        with h5py.File(cache_path, 'a') as hdf:
            if node_id_str in hdf:
                del hdf[node_id_str]
            hdf.create_dataset(node_id_str, data=combined_np,
                               compression="gzip", compression_opts=9)
        print(f"Cached preprocessed features for {node_id_str}")

        if protein_regions:
            region_slices = []
            for _, (rs, re) in protein_regions.items():
                region_slices.append(combined_np[rs:re + 1])
            combined_np = np.concatenate(region_slices, axis=0)

        combined_features = np.nan_to_num(combined_np.astype(np.float32))
        print(f"Final features shape: {combined_features.shape}")

    return combined_features


def retrieve_features(cache_path, nodeId):
    with h5py.File(cache_path, 'r') as hdf:
        if str(nodeId) in hdf:
            return hdf[str(nodeId)][:]
        else:
            print(f"NodeId {nodeId} not found in cache.")
            return None


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------

def load_legacy_keras_model(model_path):
    print(f"Attempting to load model from: {model_path}")

    try:
        import tf_keras
        model = tf_keras.models.load_model(model_path, compile=False)
        print("Model loaded successfully with tf_keras")
        return model
    except ImportError:
        try:
            import subprocess, sys
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tf-keras'])
            import tf_keras
            model = tf_keras.models.load_model(model_path, compile=False)
            print("Model loaded successfully with tf_keras")
            return model
        except Exception as e:
            print(f"Could not install or use tf_keras: {e}")
    except Exception as e:
        print(f"tf_keras loading failed: {e}")

    try:
        model = tf.keras.models.load_model(model_path, compile=False)
        print("Model loaded successfully with standard method")
        return model
    except Exception as e:
        print(f"Standard loading failed: {e}")

    h5_path = model_path.replace('.keras', '.h5')
    if os.path.exists(h5_path):
        try:
            model = tf.keras.models.load_model(h5_path, compile=False)
            print("Model loaded successfully from .h5 format")
            return model
        except Exception as e:
            print(f"Failed to load .h5 format: {e}")

    raise RuntimeError(
        f"Could not load model from {model_path}.\n"
        "Please try: pip install tf-keras\n"
        f"Current TensorFlow version: {tf.__version__}"
    )


model = load_legacy_keras_model(model_path)


# ---------------------------------------------------------------------------
# Prediction 
# ---------------------------------------------------------------------------

def predict_mutations(
    cache_path, genome_seq, mutations, codon_mapper, config_file,
    node_ids, elapsed_day=0, depth=0, protein_regions=None, k=30, model=None
):
    """
    Generates 4 predictions (A/T/G/C) per genome position.

    Feature vector order matches the paper (Ayaz et al. 2025):
        k+8 = elapsed_days,  k+9 = tree_depth,
        k+10 = synonymous,   k+11 = ORF name

    All vectors are bulk-preprocessed ONCE before inference.

    Returns:
        predictions_reshaped: np.ndarray of shape (num_positions, 4)
                              columns = [A, T, G, C] mutation probabilities
    """
    print("=" * 70)
    print("PREDICTING MUTATIONS (correct feature order from Ayaz et al. 2025)")
    print("=" * 70)

    if model is None:
        model = globals().get('model') or load_legacy_keras_model(model_path)

    if isinstance(codon_mapper, str):
        with open(codon_mapper) as f:
            codon_mapper_dict = json.load(f)
    else:
        codon_mapper_dict = codon_mapper

    mid_point = k // 2
    nucleotides = ['A', 'T', 'G', 'C']

    padded_seq = "-" * mid_point + genome_seq + "-" * mid_point
    aa_seq = [codon_mapper_dict.get(genome_seq[i:i + 3], 'X')
              for i in range(0, len(genome_seq), 3)]

    if protein_regions:
        positions = []
        for _, (start, end) in protein_regions.items():
            positions.extend(range(start, min(end + 1, len(genome_seq))))
    else:
        positions = list(range(len(genome_seq)))

    # --- STEP 1: Accumulate ALL raw feature rows (4 per position) ---
    all_raw_data = []
    print(f"Generating raw features: {len(positions)} positions × 4 nucleotides = "
          f"{len(positions) * 4} rows...")

    for idx in positions:
        aa_idx = idx // 3
        codon_start = aa_idx * 3
        if codon_start + 3 > len(genome_seq):
            continue

        original_codon = list(genome_seq[codon_start:codon_start + 3])
        original_aa = aa_seq[aa_idx]
        window = list(padded_seq[idx:idx + k])
        center_nuc = window[mid_point]
        protein_reg = find_protein_region(idx, protein_regions) if protein_regions else "Non-coding"

        for nuc in nucleotides:
            mut_codon = original_codon[:]
            mut_codon[idx % 3] = nuc
            new_aa = codon_mapper_dict.get(''.join(mut_codon), 'X')

            nuc_pam    = config_file['nucleotide sub. matrix'].get(center_nuc, {}).get(nuc, 0)
            aa_pam     = config_file['AA PAM matrix'].get(original_aa, {}).get(new_aa, 0)
            synonymous = int(original_aa == new_aa)
            aa_feats   = get_aa_features(original_aa, new_aa, config_file)

            sample_data = [
                *window,              # 0-29  (Cat)
                center_nuc,           # 30    (Cat)  original nucleotide
                nuc,                  # 31    (Cat)  candidate mutant     ← varies
                idx,                  # 32    (Num)  position
                nuc_pam,              # 33    (Num)  nucleotide PAM score
                original_aa,          # 34    (Cat)  original AA
                new_aa,               # 35    (Cat)  new AA
                aa_pam,               # 36    (Num)  AA PAM score
                float(elapsed_day),   # 37    (Num)  elapsed days  ← k+8 ✓
                float(depth),         # 38    (Num)  tree depth    ← k+9 ✓
                synonymous,           # 39    (Num)  synonymous    ← k+10 ✓
                protein_reg,          # 40    (Cat)  ORF name      ← k+11 ✓
                *aa_feats,            # 41-58 (Num)  18 biochemical features
            ]
            all_raw_data.append(sample_data)

    # --- STEP 2: Bulk preprocess ONCE ---
    print(f"Bulk preprocessing {len(all_raw_data)} raw vectors...")
    X = preprocess_matrix(all_raw_data, expected_size=205)
    print(f"Preprocessed matrix shape: {X.shape}")

    # --- STEP 3: Batch inference ---
    print("Running batch inference...")
    start_pred = time.time()
    preds = model.predict(X, batch_size=4096, verbose=1)
    print(f"Inference completed in {time.time() - start_pred:.2f}s")
    print(f"Raw model output shape: {preds.shape}, sample values: {preds[:5].ravel()}")

    # Model has sigmoid output → shape is (N, 1) with values in [0, 1]
    # Value = probability of mutation at this position for this nucleotide
    if preds.ndim == 2 and preds.shape[1] == 1:
        prob_values = preds[:, 0]          # shape: (N,)
    elif preds.ndim == 2 and preds.shape[1] == 2:
        prob_values = preds[:, 1]          # binary class → take class-1 prob
    else:
        prob_values = preds.ravel()

    # Reshape to (num_positions, 4) — columns = [A, T, G, C]
    # predictions_reshaped = prob_values.reshape(-1, 4)
    # Reshape to (num_positions, 4) — columns = [A, T, G, C]
    predictions_reshaped = prob_values.reshape(-1, 4)
    # row_sums = predictions_reshaped.sum(axis=1, keepdims=True)
    # row_sums = np.where(row_sums == 0, 1.0, row_sums)
    # predictions_reshaped = predictions_reshaped / row_sums

    print(f"Final predictions shape: {predictions_reshaped.shape}  (positions × [A,T,G,C])")
    print(f"Sample predictions (first 5 positions):\n{predictions_reshaped[:5]}")

    np.save('predictions.npy', predictions_reshaped)
    print("Predictions saved to 'predictions.npy'")

    return predictions_reshaped


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

genome_sequence = get_sequence()

node_ids = [
    "EGY/CCHE57357_Wave_3_A029/2021|MZ380261.1|2021-05-11"
]

mutations = parse_mutations(node_ids[0])
depth = get_sample_depth(depth_file, node_ids[0])

features = cache_precomputed_features(
    cache_path=cache_path,
    genome_seq=construct_variant_genome(genome_sequence, mutations),
    mutations=mutations,
    codon_mapper=json.load(open(codon_mapping_path)),
    config_file=configs(),
    node_ids=node_ids,
    elapsed_day=110,
    depth=depth,
    protein_regions=None,
)

predictions = predict_mutations(
    cache_path=cache_path,
    genome_seq=genome_sequence,
    mutations=mutations,
    codon_mapper=codon_mapping_path,
    config_file=configs(),
    node_ids=node_ids[0],
    elapsed_day=130,
    depth=depth,
    protein_regions=None,
)

end_time2 = time.time()
print(f"Total runtime: {end_time2 - start_time:.2f}s")