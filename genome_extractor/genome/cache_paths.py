import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")

NODE_FEATURES_CACHE_PATH = os.path.join(CACHE_DIR, "node_features.h5")
REFERENCE_FEATURES_CACHE_PATH = os.path.join(CACHE_DIR, "features.h5")
PHYLO_FEATURES_CACHE_PATH = os.path.join(CACHE_DIR, "phylo_features_cache.h5")
