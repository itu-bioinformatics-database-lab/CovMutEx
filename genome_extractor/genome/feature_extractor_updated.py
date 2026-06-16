import json
import os
import time
import h5py
import hashlib
import numpy as np
import pandas as pd
from .configs import configs
from .cache_paths import (
    CACHE_DIR,
    NODE_FEATURES_CACHE_PATH,
    PHYLO_FEATURES_CACHE_PATH,
    REFERENCE_FEATURES_CACHE_PATH,
)
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
DEFAULT_TEMPLATE_CACHE_FORMAT = "default_atgc_template_v1"
NODE_RAW_CACHE_FORMAT = "default_atgc_node_raw_v1"
NODE_REQUEST_CACHE_FORMAT = "default_atgc_node_request_v1"

os.makedirs(CACHE_DIR, exist_ok=True)

print(f"Base dir: {base_dir}")
print(f"Genome extractor dir: {genome_extractor_dir}")
print(f"ROOT_PATH: {ROOT_PATH}")

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


def add_terminal_padding_nucleotide(genome_seq, nucleotide="A"):
    if not genome_seq:
        return nucleotide
    return genome_seq + nucleotide


def get_reference_protein_regions(config_file=None):
    cfg = config_file() if callable(config_file) else config_file
    if cfg is None:
        cfg = configs()
    return cfg.get("protein regions", {})


def _encode_raw_rows(raw_rows):
    return [json.dumps(row).encode("utf-8") for row in raw_rows]


def _decode_raw_rows(raw_encoded):
    return [json.loads(row.decode("utf-8") if isinstance(row, bytes) else row) for row in raw_encoded]


def _write_raw_rows_dataset(h5_parent, dataset_name, raw_rows):
    raw_encoded = _encode_raw_rows(raw_rows)
    dt = h5py.special_dtype(vlen=bytes)
    raw_ds = h5_parent.create_dataset(dataset_name, (len(raw_encoded),), dtype=dt)
    for i, row in enumerate(raw_encoded):
        raw_ds[i] = row


def _load_reference_template_raw_rows(features_h5_path):
    with h5py.File(features_h5_path, "r") as h5f:
        return _decode_raw_rows(h5f["features_raw"][:])


def _apply_context_to_raw_rows(raw_rows, elapsed_day=0, depth=0):
    for row in raw_rows:
        row[37] = float(elapsed_day or 0)
        row[38] = float(depth or 0)
    return raw_rows


def _slice_preprocessed_features_by_regions(features_np, protein_regions, nucleotides_per_position=4):
    if not protein_regions:
        return np.nan_to_num(features_np.astype(np.float32))

    region_slices = []
    for _, (start, end) in protein_regions.items():
        row_start = start * nucleotides_per_position
        row_end = (end + 1) * nucleotides_per_position
        region_slices.append(features_np[row_start:row_end])

    if not region_slices:
        return np.zeros((0, features_np.shape[1]), dtype=np.float32)
    return np.nan_to_num(np.concatenate(region_slices, axis=0).astype(np.float32))


def _copy_raw_rows(raw_rows):
    return [list(row) for row in raw_rows]


def _node_raw_cache_key(node_id, k=30):
    key_source = f"{node_id or '__reference__'}|{k}"
    digest = hashlib.md5(key_source.encode("utf-8")).hexdigest()[:12]
    safe_prefix = "reference" if not node_id else "node"
    return f"{safe_prefix}_{digest}_k_{k}"


def _node_request_cache_key(elapsed_day=0, depth=0):
    return f"day_{int(elapsed_day or 0)}__depth_{int(depth or 0)}"


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

CACHE_P_FILE = PHYLO_FEATURES_CACHE_PATH

# Hiçbir yerde kullanılmıyor (extract_phylogenetic_features)
def normalize_features(features):
    min_val = min(features.values())
    max_val = max(features.values())
    if max_val == min_val:
        return {k: 1.0 for k in features}
    return {k: (v - min_val) / (max_val - min_val) for k, v in features.items()}

# Hiçbir yerde kullanılmıyor (extract_phylogenetic_features)
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

# Hiçbir yerde kullanılmıyor
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

# Hiçbir yerde kullanılmıyor
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
# NEW Reference genome feature extraction  (stores both RAW and preprocessed lists)
# ---------------------------------------------------------------------------

def precompute_default_feature_template_cache(
    genome_txt_path, codon_mapping_path, config_file, output_h5_path, k=30
):
    """
    Build the canonical padded reference template for the active ATGC extractor.

    Stored artifacts:
        features_raw: raw rows for the padded reference genome
        features:     preprocess_matrix(features_raw)
    """
    genome_seq = get_sequence()
    genome_md5 = hashlib.md5(genome_seq.encode("utf-8")).hexdigest()
    cfg = config_file() if callable(config_file) else config_file

    if os.path.exists(output_h5_path):
        try:
            with h5py.File(output_h5_path, "r") as h5f:
                if (
                    h5f.attrs.get("cache_format") == DEFAULT_TEMPLATE_CACHE_FORMAT
                    and int(h5f.attrs.get("genome_length", -1)) == len(genome_seq)
                    and h5f.attrs.get("genome_md5") == genome_md5
                    and int(h5f.attrs.get("nucleotides_per_position", -1)) == 4
                    and int(h5f.attrs.get("k", -1)) == k
                ):
                    print(f"Using cached reference ATGC template from {output_h5_path}")
                    return h5f["features"][:]
        except Exception:
            pass

    print("Reference ATGC template missing or incompatible - rebuilding...")
    raw_rows = build_all_raw_feature_rows(
        genome_seq=genome_seq,
        codon_mapper=codon_mapping_path,
        config_file=cfg,
        elapsed_day=0,
        depth=0,
        protein_regions=None,
        k=k,
    )
    features_np = preprocess_matrix(raw_rows, expected_size=205)

    with h5py.File(output_h5_path, "w") as h5f:
        h5f.create_dataset("features", data=features_np, compression="gzip", compression_opts=9)
        _write_raw_rows_dataset(h5f, "features_raw", raw_rows)
        h5f.attrs["genome_length"] = len(genome_seq)
        h5f.attrs["vector_length"] = features_np.shape[1]
        h5f.attrs["nucleotides_per_position"] = 4
        h5f.attrs["cache_format"] = DEFAULT_TEMPLATE_CACHE_FORMAT
        h5f.attrs["k"] = k
        h5f.attrs["genome_md5"] = genome_md5
        h5f.attrs["elapsed_day"] = 0
        h5f.attrs["depth"] = 0
        h5f.attrs["is_padded_reference"] = True

    print(f"Reference ATGC template cached to {output_h5_path}")
    return features_np


def _build_variant_rows_for_positions(
    genome_seq,
    positions,
    codon_mapper,
    config_file,
    elapsed_day=0,
    depth=0,
    k=30,
):
    all_rows = build_all_raw_feature_rows(
        genome_seq=genome_seq,
        codon_mapper=codon_mapper,
        config_file=config_file,
        elapsed_day=elapsed_day,
        depth=depth,
        protein_regions=None,
        k=k,
        selected_positions=positions,
    )
    rows_by_position = {}
    for index, pos in enumerate(positions):
        start = index * 4
        rows_by_position[pos] = all_rows[start:start + 4]
    return rows_by_position


def cache_node_atgc_features(
    cache_path,
    node_id,
    genome_seq,
    mutations,
    codon_mapper,
    config_file,
    elapsed_day=0,
    depth=0,
    protein_regions=None,
    k=30,
):
    """
    Cache mutation-aware ATGC features for a node on the full padded genome.

    Disk cache stores a mutation-aware raw base per node with default context
    (elapsed_day=0, depth=0). Request-specific preprocessing is derived from
    that raw base on every request.
    """
    if not os.path.exists(REFERENCE_FEATURES_CACHE_PATH):
        precompute_default_feature_template_cache(
            genome_txt_path=genome_txt_path,
            codon_mapping_path=codon_mapper,
            config_file=config_file,
            output_h5_path=REFERENCE_FEATURES_CACHE_PATH,
            k=k,
        )

    cfg = config_file() if callable(config_file) else config_file
    variant_md5 = hashlib.md5(genome_seq.encode("utf-8")).hexdigest()
    raw_cache_key = _node_raw_cache_key(node_id=node_id, k=k)
    request_cache_key = _node_request_cache_key(elapsed_day=elapsed_day, depth=depth)

    with h5py.File(cache_path, "a") as hdf:
        if raw_cache_key in hdf:
            group = hdf[raw_cache_key]
            if (
                group.attrs.get("cache_format") == NODE_RAW_CACHE_FORMAT
                and group.attrs.get("variant_genome_md5") == variant_md5
                and group.attrs.get("node_id") == (node_id or "__reference__")
            ):
                print(f"Loading cached node raw features: {raw_cache_key}")
                node_raw_base = _decode_raw_rows(group["features_raw"][:])
                requests_group = group.require_group("requests")
                if request_cache_key in requests_group:
                    request_group = requests_group[request_cache_key]
                    if (
                        request_group.attrs.get("cache_format") == NODE_REQUEST_CACHE_FORMAT
                        and request_group.attrs.get("elapsed_day") == float(elapsed_day or 0)
                        and request_group.attrs.get("depth") == float(depth or 0)
                    ):
                        print(f"Loading cached node request features: {raw_cache_key}/{request_cache_key}")
                        features_np = request_group["features"][:]
                        return _slice_preprocessed_features_by_regions(features_np, protein_regions, 4)
            else:
                del hdf[raw_cache_key]
                node_raw_base = None
        else:
            node_raw_base = None

        if node_raw_base is None:
            print(f"Building cached node raw features: {raw_cache_key}")
            base_raw = _load_reference_template_raw_rows(REFERENCE_FEATURES_CACHE_PATH)
            node_raw_base = _apply_context_to_raw_rows(_copy_raw_rows(base_raw), elapsed_day=0, depth=0)

            if mutations:
                affected_positions = get_affected_positions(mutations, len(genome_seq), k, protein_regions=None)
                if affected_positions:
                    variant_rows = _build_variant_rows_for_positions(
                        genome_seq=genome_seq,
                        positions=affected_positions,
                        codon_mapper=codon_mapper,
                        config_file=cfg,
                        elapsed_day=0,
                        depth=0,
                        k=k,
                    )
                    for pos, rows in variant_rows.items():
                        start = pos * 4
                        node_raw_base[start:start + 4] = rows

            group = hdf.create_group(raw_cache_key)
            _write_raw_rows_dataset(group, "features_raw", node_raw_base)
            group.create_group("requests")
            group.attrs["cache_format"] = NODE_RAW_CACHE_FORMAT
            group.attrs["node_id"] = node_id or "__reference__"
            group.attrs["default_elapsed_day"] = 0.0
            group.attrs["default_depth"] = 0.0
            group.attrs["k"] = int(k)
            group.attrs["genome_length"] = len(genome_seq)
            group.attrs["nucleotides_per_position"] = 4
            group.attrs["variant_genome_md5"] = variant_md5

    request_raw = _apply_context_to_raw_rows(_copy_raw_rows(node_raw_base), elapsed_day=elapsed_day, depth=depth)

    print(f"Bulk preprocessing {len(request_raw)} node rows from raw cache key: {raw_cache_key}")
    features_np = np.nan_to_num(preprocess_matrix(request_raw, expected_size=205).astype(np.float32))

    with h5py.File(cache_path, "a") as hdf:
        group = hdf[raw_cache_key]
        requests_group = group.require_group("requests")
        if request_cache_key in requests_group:
            del requests_group[request_cache_key]
        request_group = requests_group.create_group(request_cache_key)
        request_group.create_dataset("features", data=features_np, compression="gzip", compression_opts=9)
        request_group.attrs["cache_format"] = NODE_REQUEST_CACHE_FORMAT
        request_group.attrs["elapsed_day"] = float(elapsed_day or 0)
        request_group.attrs["depth"] = float(depth or 0)
        request_group.attrs["k"] = int(k)
        request_group.attrs["vector_length"] = features_np.shape[1]
        request_group.attrs["genome_length"] = len(genome_seq)
        request_group.attrs["nucleotides_per_position"] = 4

    return _slice_preprocessed_features_by_regions(features_np, protein_regions, 4)

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
    _features_h5 = REFERENCE_FEATURES_CACHE_PATH
    with h5py.File(_features_h5, 'r') as h5f:
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

# Hiçbir yerde kullanılmıyor
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


# ---------------------------------------------------------------------------
# Prediction 
# ---------------------------------------------------------------------------

def build_all_raw_feature_rows(
    genome_seq, codon_mapper, config_file,
    elapsed_day=0, depth=0, protein_regions=None, k=30, selected_positions=None
):
    """
    Generates 4 raw feature rows per genome position (one per A/T/G/C candidate).

    Shared by both predict_mutations() and DefaultCovMutExFeatureExtractor.extract_features()
    so the normal Protocol-based flow uses the exact same features as the legacy path.

    Returns:
        list of raw feature rows, length = num_positions * 4
    """
    cfg = config_file() if callable(config_file) else config_file

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

    protein_region_lookup = get_reference_protein_regions(cfg)

    if selected_positions is not None:
        positions = [idx for idx in selected_positions if 0 <= idx < len(genome_seq)]
    elif protein_regions:
        positions = []
        for _, (start, end) in protein_regions.items():
            positions.extend(range(start, min(end + 1, len(genome_seq))))
    else:
        positions = list(range(len(genome_seq)))

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
        protein_reg = find_protein_region(idx, protein_region_lookup) if protein_region_lookup else "Non-coding"

        for nuc in nucleotides:
            mut_codon = original_codon[:]
            mut_codon[idx % 3] = nuc
            new_aa = codon_mapper_dict.get(''.join(mut_codon), 'X')

            nuc_pam    = cfg['nucleotide sub. matrix'].get(center_nuc, {}).get(nuc, 0)
            aa_pam     = cfg['AA PAM matrix'].get(original_aa, {}).get(new_aa, 0)
            synonymous = int(original_aa == new_aa)
            aa_feats   = get_aa_features(original_aa, new_aa, cfg)

            all_raw_data.append([
                *window,
                center_nuc,
                nuc,
                idx,
                nuc_pam,
                original_aa,
                new_aa,
                aa_pam,
                float(elapsed_day),
                float(depth),
                synonymous,
                protein_reg,
                *aa_feats,
            ])

    return all_raw_data


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

    # --- STEP 1: Accumulate ALL raw feature rows (4 per position) ---
    all_raw_data = build_all_raw_feature_rows(
        genome_seq=genome_seq,
        codon_mapper=codon_mapper,
        config_file=config_file,
        elapsed_day=elapsed_day,
        depth=depth,
        protein_regions=protein_regions,
        k=k,
    )

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

if __name__ == "__main__":
    start_time = time.time()

    cache_path = NODE_FEATURES_CACHE_PATH
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
        output_h5_path=REFERENCE_FEATURES_CACHE_PATH,
        k=30
    )

    model = load_legacy_keras_model(model_path)

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


"""

C:/Users/afanu/Projeler/CovMutEx/                                 <- ROOT_PATH
|
`-- genome_extractor/                                             <- genome_extractor_dir
    |
    |-- covid19_models/
    |   `-- models/
    |       `-- balanced_data_model.keras                         <- model_path
    |
    `-- genome/                                                   <- base_dir
        |
        |-- feature_extractor_updated.py                          <- current file
        |-- covmutex_feature_extractors.py
        |-- build_cache.py
        |-- configs.py
        |-- codon_aa_mapping.json                                 <- codon_mapping_path
        |-- genome.txt                                            <- genome_txt_path
        |-- mutations.txt                                         <- mutations_txt_path
        |-- phylogenetic_tree.nwk                                 <- phylo_tree_path
        |-- depth_date.json                                       <- depth_file
        `-- cache/
            |-- features.h5                                       <- REFERENCE_FEATURES_CACHE_PATH
            |-- node_features.h5                                  <- NODE_FEATURES_CACHE_PATH
            `-- phylo_features_cache.h5                           <- PHYLO_FEATURES_CACHE_PATH


CALL TREE
==============================================================================

External callers
|
|-- viewsUpdated.py
|   |-- parse_mutations()
|   |-- construct_variant_genome()
|   |-- add_terminal_padding_nucleotide()
|   `-- get_sample_depth()
|
|-- helpers.py
|   `-- predict_mutations()                                      <- active runtime path
|
|-- covmutex_feature_extractors.py
|   |-- DefaultCovMutExFeatureExtractor.extract_features()
|   |   `-- cache_node_atgc_features()                          <- active server-model cache path
|   `-- load_feature_extractor()
|
`-- build_cache.py
    `-- precompute_default_feature_template_cache()              <- regenerates canonical features.h5

------------------------------------------------------------------------------

Active internal call tree
|
|-- cache_node_atgc_features()
|   |-- precompute_default_feature_template_cache()              <- only if features.h5 is missing/incompatible
|   |-- _node_raw_cache_key()
|   |-- _node_request_cache_key()
|   |-- _load_reference_template_raw_rows()
|   |-- _copy_raw_rows()
|   |-- _apply_context_to_raw_rows()
|   |-- get_affected_positions()
|   |-- _build_variant_rows_for_positions()
|   |   `-- build_all_raw_feature_rows()
|   `-- _slice_preprocessed_features_by_regions()
|
|-- precompute_default_feature_template_cache()
|   |-- get_sequence()                                           <- padded reference genome (29904)
|   |-- build_all_raw_feature_rows()
|   |   |-- get_reference_protein_regions()
|   |   |-- find_protein_region()
|   |   `-- get_aa_features()
|   `-- preprocess_matrix()
|
|-- build_all_raw_feature_rows()
|   |-- get_reference_protein_regions()
|   |-- find_protein_region()
|   `-- get_aa_features()
|
`-- predict_mutations()                                          <- legacy standalone inference helper
    |-- build_all_raw_feature_rows()
    |-- preprocess_matrix()
    `-- load_legacy_keras_model()                                <- only if model=None

------------------------------------------------------------------------------

Legacy / secondary paths
|
|-- cache_precomputed_features()                                 <- old node-based cache flow, not used by new server runtime
|   |-- process_variant_raw()
|   |   |-- get_affected_positions()
|   |   |-- find_protein_region()
|   |   `-- get_aa_features()
|   `-- preprocess_matrix()
|
`-- precompute_and_cache_ref_features()                          <- old reference-cache builder, kept for compatibility inside this module
    `-- precompute_feature_vectors()
        |-- get_sequence()
        |-- translate_nucleotides_to_amino_acids()
        |-- find_protein_region()
        `-- get_aa_features()

------------------------------------------------------------------------------

Standalone / currently unused
|
|-- extract_phylogenetic_features()
|   |-- normalize_features()
|   `-- calculate_phylogenetic_diversity()
|
|-- balance_dataset()
`-- retrieve_features()


CACHE CONTENTS
==============================================================================

features.h5
|
|-- Purpose
|   `-- Canonical reference template cache for the active server-model flow
|
|-- Scope
|   |-- Full padded reference genome (29904)
|   |-- 4 candidate nucleotides per position (A/T/G/C)
|   `-- Default context only:
|       |-- elapsed_day = 0
|       |-- depth = 0
|       `-- k = 30
|
|-- Datasets
|   |-- features_raw
|   |   `-- Raw biological feature rows (~59 fields each) before encoding
|   `-- features
|       `-- Preprocessed model-ready matrix with shape (29904 * 4, 205)
|
`-- Attrs
    |-- cache_format = default_atgc_template_v1
    |-- genome_length = 29904
    |-- nucleotides_per_position = 4
    |-- vector_length = 205
    |-- elapsed_day = 0
    |-- depth = 0
    |-- k = 30
    `-- is_padded_reference = True

------------------------------------------------------------------------------

node_features.h5
|
|-- Purpose
|   `-- Node-based mutation-aware raw + request-feature cache
|
|-- Scope
|   |-- Full padded variant genome for one node
|   |-- 4 candidate nucleotides per position (A/T/G/C)
|   |-- Default node context only:
|   |   |-- elapsed_day = 0
|   |   |-- depth = 0
|   |   `-- k = 30
|   `-- Canonical full-genome raw cache first; region slicing happens after
|       request-specific preprocessing
|
|-- Group layout
|   `-- <node_raw_cache_key>/
|       |-- features_raw
|       |   `-- Raw mutation-aware rows with default elapsed_day/depth
|       |-- requests/
|       |   `-- <day_depth_key>/
|       |       `-- features
|       |           `-- Preprocessed matrix for one elapsed_day/depth request
|       `-- attrs
|           `-- Metadata for the node-specific raw base
|
`-- Group attrs
    |-- cache_format = default_atgc_node_raw_v1
    |-- node_id
    |-- default_elapsed_day = 0
    |-- default_depth = 0
    |-- k = 30
    |-- genome_length = 29904
    |-- nucleotides_per_position = 4
    `-- variant_genome_md5

Request-group attrs
    |-- cache_format = default_atgc_node_request_v1
    |-- elapsed_day
    |-- depth
    |-- k = 30
    |-- genome_length = 29904
    |-- nucleotides_per_position = 4
    `-- vector_length = 205

------------------------------------------------------------------------------

EXAMPLE SHAPES AND DATA
==============================================================================

1) Reference template example in features.h5

Full padded reference genome length:
    29904

Candidates per position:
    4  (A, T, G, C)

Preprocessed matrix shape:
    (29904 * 4, 205) = (119616, 205)

Raw rows dataset length:
    119616

Conceptual ordering:
    row 0  -> position 0, candidate A
    row 1  -> position 0, candidate T
    row 2  -> position 0, candidate G
    row 3  -> position 0, candidate C
    row 4  -> position 1, candidate A
    row 5  -> position 1, candidate T
    ...

Mini raw-row example (conceptual, shortened):
    [
      ['A','T','G', ... 30-window ..., 'A', 'G', 100, 0.42, 'I', 'V', 0.18,
       0, 0, 0, 'ORF1ab', ... AA numeric features ...],
      ['A','T','G', ... 30-window ..., 'A', 'C', 100, 0.31, 'I', 'L', 0.25,
       0, 0, 0, 'ORF1ab', ... AA numeric features ...]
    ]

Mini preprocessed example (same rows after encoding, shortened):
    [
      [0,1,0,0,0,  1,0,0,0,0,  ...,  0.23,-0.51,1.02, ..., 0.0],
      [0,1,0,0,0,  0,0,1,0,0,  ..., -0.14, 0.77,0.35, ..., 0.0]
    ]

Each final row:
    205 numeric features

------------------------------------------------------------------------------

2) Node-specific example in node_features.h5

Example raw-base cache key:
    node_ff9fbcf71251_k_30

Meaning:
    node_id     = EGY/CCHE57357_Wave_3_A029/2021|MZ380261.1|2021-05-11
    elapsed_day = 0   (stored default raw base)
    depth       = 0   (stored default raw base)
    k           = 30

Stored layout:
    node_features.h5
    `-- node_ff9fbcf71251_k_30/
        |-- features_raw -> length 119616
        `-- requests/
            `-- day_120__depth_19/
                `-- features -> shape (119616, 205)

What changed relative to features.h5:
    - Same overall shape
    - But rows affected by this node's mutations are rebuilt from the padded
      variant genome, then merged back into the full reference template
    - elapsed_day and depth are kept at default 0 in the stored raw base
    - When a request comes in, elapsed_day/depth are applied to a copy of the
      raw base, then the resulting preprocessed matrix is stored under the same
      node group inside requests/

Mini example for one mutated position:
    Reference position 23403 candidate rows:
        position 23403 -> A
        position 23403 -> T
        position 23403 -> G
        position 23403 -> C

    Node-specific cache:
        The same 4 rows exist, but their raw biological context may differ
        because the node's nearby mutations can change:
        - the local nucleotide window
        - the codon
        - the resulting amino acid
        - synonymous / non-synonymous flag
        - substitution scores

------------------------------------------------------------------------------ 

3) Region slicing example after cache load

Example region:
    ORF10 = [29558, 29674]

Region length:
    29674 - 29558 + 1 = 117 positions

Rows sent to model after slicing:
    117 * 4 = 468 rows

Request-specific preprocessed matrix before slice:
    (119616, 205)

Sliced preprocessed matrix shape:
    (468, 205)

Model output for a single-input model:
    (468, 1)

Reshaped position-level predictions:
    (117, 4)

Meaning:
    117 genomic positions
    x
    4 candidate nucleotide scores per position

------------------------------------------------------------------------------

4) Full-genome prediction example

Input feature matrix after request-specific preprocessing:
    (119616, 205)

Single-input model raw output:
    (119616, 1)

Reshaped position-level output:
    (29904, 4)

User-facing output after dropping internal padding base:
    (29903, 4)

Frontend genomeData layout:
    4 x 29903

Meaning:
    genomeData[0] -> A scores across all real genome positions
    genomeData[1] -> T scores
    genomeData[2] -> G scores
    genomeData[3] -> C scores

"""
