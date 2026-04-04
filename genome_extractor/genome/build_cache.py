"""
Build the canonical cache files under genome/cache/.

Usage:
    python -m genome.build_cache
    python genome/build_cache.py
    python -m genome.build_cache --reset
"""

import argparse
import os
import sys

import h5py

if __package__ in (None, ""):
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.dirname(CURRENT_DIR)
    if PROJECT_DIR not in sys.path:
        sys.path.insert(0, PROJECT_DIR)

    from genome.cache_paths import (
        CACHE_DIR,
        NODE_FEATURES_CACHE_PATH,
        PHYLO_FEATURES_CACHE_PATH,
        REFERENCE_FEATURES_CACHE_PATH,
    )
    from genome.configs import configs
    from genome.feature_extractor_updated import (
        codon_mapping_path,
        genome_txt_path,
        precompute_default_feature_template_cache,
    )
else:
    from .cache_paths import (
        CACHE_DIR,
        NODE_FEATURES_CACHE_PATH,
        PHYLO_FEATURES_CACHE_PATH,
        REFERENCE_FEATURES_CACHE_PATH,
    )
    from .configs import configs
    from .feature_extractor_updated import (
        codon_mapping_path,
        genome_txt_path,
        precompute_default_feature_template_cache,
    )


def touch_h5(path: str) -> None:
    with h5py.File(path, "a"):
        pass


def reset_cache_files() -> None:
    for path in (
        REFERENCE_FEATURES_CACHE_PATH,
        NODE_FEATURES_CACHE_PATH,
        PHYLO_FEATURES_CACHE_PATH,
    ):
        if os.path.exists(path):
            os.remove(path)
            print(f"Removed: {path}")


def build_cache(reset: bool = False) -> None:
    os.makedirs(CACHE_DIR, exist_ok=True)

    if reset:
        reset_cache_files()

    print(f"Cache directory: {CACHE_DIR}")
    print("Building reference feature cache...")
    precompute_default_feature_template_cache(
        genome_txt_path=genome_txt_path,
        codon_mapping_path=codon_mapping_path,
        config_file=configs,
        output_h5_path=REFERENCE_FEATURES_CACHE_PATH,
        k=30,
    )

    if not os.path.exists(NODE_FEATURES_CACHE_PATH):
        touch_h5(NODE_FEATURES_CACHE_PATH)
        print(f"Initialized empty node cache: {NODE_FEATURES_CACHE_PATH}")
    else:
        print(f"Node cache already exists: {NODE_FEATURES_CACHE_PATH}")

    if not os.path.exists(PHYLO_FEATURES_CACHE_PATH):
        touch_h5(PHYLO_FEATURES_CACHE_PATH)
        print(f"Initialized empty phylo cache: {PHYLO_FEATURES_CACHE_PATH}")
    else:
        print(f"Phylo cache already exists: {PHYLO_FEATURES_CACHE_PATH}")

    print("Cache build complete.")
    print(f"Reference cache: {REFERENCE_FEATURES_CACHE_PATH}")
    print(f"Node cache: {NODE_FEATURES_CACHE_PATH}")
    print(f"Phylo cache: {PHYLO_FEATURES_CACHE_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build CovMutEx cache files")
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Delete canonical cache files in genome/cache/ before rebuilding.",
    )
    args = parser.parse_args()
    build_cache(reset=args.reset)


if __name__ == "__main__":
    main()
