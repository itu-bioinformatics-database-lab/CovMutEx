import importlib
import os
import sys
from typing import Any, Dict, List, Optional, Protocol, Tuple

import numpy as np

from .plugin_runtime import PLUGIN_CONTRACT_VERSION, normalize_elapsed_day, normalize_node_ids


class CovMutExFeatureExtractor(Protocol):
    """
    Canonical feature extractor contract used by the plugin runtime.

    Uploaded extractors are expected to implement:
        extract_features(...)
        get_feature_dimension()
        get_feature_description()
        get_metadata()

    Helper files are not part of the contract. They are expected to live in the
    same bundle directory and can be loaded from the provided `base_dir`.
    """

    def extract_features(
        self,
        genome_seq: str,
        mutations: List[Tuple[int, str, str, int, str]],
        node_ids: List[str],
        elapsed_day: Optional[int],
        protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
        k: int = 30,
        **kwargs,
    ) -> np.ndarray:
        ...

    def get_feature_dimension(self) -> int:
        ...

    def get_feature_description(self) -> Dict[str, str]:
        ...

    def get_metadata(self) -> Dict[str, Any]:
        ...


def _normalize_extractor_kwargs(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    normalized = dict(kwargs)
    normalized["node_ids"] = normalize_node_ids(kwargs.get("node_ids"))
    normalized["node_id"] = normalized["node_ids"][0] if normalized["node_ids"] else None
    normalized["elapsed_day"] = normalize_elapsed_day(kwargs.get("elapsed_day"))
    return normalized


class DefaultCovMutExFeatureExtractor:
    def __init__(self) -> None:
        self._feature_dim = 205
        self.nucleotides_per_position = 4

    def extract_features(
        self,
        genome_seq: str,
        mutations: List[Tuple[int, str, str, int, str]],
        node_ids: List[str],
        elapsed_day: Optional[int],
        protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
        k: int = 30,
        **kwargs,
    ) -> np.ndarray:
        del mutations, node_ids

        from . import feature_extractor_updated as feu

        codon_mapper = kwargs.get("codon_mapper")
        config_file = kwargs.get("config_file")

        if codon_mapper is None:
            feu_dir = os.path.dirname(os.path.abspath(feu.__file__))
            codon_mapper = os.path.join(feu_dir, "codon_aa_mapping.json")

        if config_file is None:
            config_file = feu.configs()

        all_raw_data = feu.build_all_raw_feature_rows(
            genome_seq=genome_seq,
            codon_mapper=codon_mapper,
            config_file=config_file,
            elapsed_day=elapsed_day or 0,
            depth=kwargs.get("depth", 0),
            protein_regions=protein_regions,
            k=k,
        )
        return feu.preprocess_matrix(all_raw_data, expected_size=self._feature_dim)

    def get_feature_dimension(self) -> int:
        return self._feature_dim

    def get_feature_description(self) -> Dict[str, str]:
        return {
            "0-29": "K-mer nucleotides (k=30)",
            "30": "Original nucleotide at position",
            "31": "Mutated nucleotide",
            "32": "Position index in genome",
            "33": "PAM250 score (nucleotide level)",
            "34": "Original amino acid",
            "35": "Mutated amino acid",
            "36": "PAM250 score (amino acid level)",
            "37": "Elapsed days",
            "38": "Phylogenetic depth",
            "39": "Synonymous indicator",
            "40": "Protein region/ORF name",
            "41-58": "AA biochemical properties",
            "59-204": "One-hot encoded and standardized features",
        }

    def get_metadata(self) -> Dict[str, Any]:
        return {
            "name": "DefaultCovMutExFeatureExtractor",
            "version": "1.0",
            "contract_version": PLUGIN_CONTRACT_VERSION,
            "feature_dimension": self._feature_dim,
            "nucleotides_per_position": self.nucleotides_per_position,
            "preprocessing": "one-hot encoding + standardization",
        }


class CustomFeatureExtractor:
    """
    Minimal scaffold for future in-repo custom extractors.
    """

    def __init__(self, **config: Any) -> None:
        self.config = config

    def extract_features(
        self,
        genome_seq: str,
        mutations: List[Tuple[int, str, str, int, str]],
        node_ids: List[str],
        elapsed_day: Optional[int],
        protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
        k: int = 30,
        **kwargs,
    ) -> np.ndarray:
        del mutations, node_ids, elapsed_day, protein_regions, k, kwargs
        return np.zeros((len(genome_seq), self.get_feature_dimension()), dtype=np.float32)

    def get_feature_dimension(self) -> int:
        return 205

    def get_feature_description(self) -> Dict[str, str]:
        return {"0-204": "Placeholder custom feature vector"}

    def get_metadata(self) -> Dict[str, Any]:
        return {
            "name": "CustomFeatureExtractor",
            "version": "1.0",
            "contract_version": PLUGIN_CONTRACT_VERSION,
            "feature_dimension": self.get_feature_dimension(),
            "config": self.config,
        }


class UploadedExtractorWrapper:
    def __init__(self, module: Any, module_path: str, base_dir: str) -> None:
        self._module = module
        self._module_path = module_path
        self._base_dir = base_dir
        self.nucleotides_per_position = getattr(module, "nucleotides_per_position", 1)

    def extract_features(
        self,
        genome_seq: str,
        mutations: List[Tuple[int, str, str, int, str]],
        node_ids: List[str],
        elapsed_day: Optional[int],
        protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
        k: int = 30,
        **kwargs,
    ) -> np.ndarray:
        if not hasattr(self._module, "extract_features"):
            raise ValueError("Uploaded extractor must define extract_features(...)")

        normalized_kwargs = _normalize_extractor_kwargs(
            {
                **kwargs,
                "node_ids": node_ids,
                "elapsed_day": elapsed_day,
                "protein_regions": protein_regions,
                "k": k,
                "base_dir": self._base_dir,
            }
        )
        forwarded_kwargs = dict(normalized_kwargs)
        forwarded_kwargs.pop("node_ids", None)
        forwarded_kwargs.pop("elapsed_day", None)
        forwarded_kwargs.pop("protein_regions", None)
        forwarded_kwargs.pop("k", None)
        features = self._module.extract_features(
            genome_seq=genome_seq,
            mutations=mutations,
            node_ids=normalized_kwargs["node_ids"],
            elapsed_day=normalized_kwargs["elapsed_day"],
            protein_regions=protein_regions,
            k=k,
            **forwarded_kwargs,
        )
        return np.asarray(features)

    def get_feature_dimension(self) -> int:
        if hasattr(self._module, "get_feature_dimension"):
            return self._module.get_feature_dimension()
        return 205

    def get_feature_description(self) -> Dict[str, str]:
        if hasattr(self._module, "get_feature_description"):
            return self._module.get_feature_description()
        return {"info": "Uploaded extractor using bundle-local helper files"}

    def get_metadata(self) -> Dict[str, Any]:
        if hasattr(self._module, "get_metadata"):
            metadata = self._module.get_metadata()
        else:
            metadata = {
                "name": "UploadedExtractor",
                "version": "1.0",
                "feature_dimension": self.get_feature_dimension(),
            }
        metadata["module_path"] = self._module_path
        metadata["base_directory"] = self._base_dir
        metadata["contract_version"] = metadata.get("contract_version", PLUGIN_CONTRACT_VERSION)
        return metadata


def _load_uploaded_module(module_path: str) -> Any:
    extractor_dir = os.path.dirname(os.path.abspath(module_path))
    bundle_name = os.path.basename(extractor_dir)
    uploaded_models_dir = os.path.dirname(extractor_dir)
    project_root = os.path.dirname(uploaded_models_dir)

    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    module_name = f"uploaded_models.{bundle_name}.feature_extractor"
    importlib.invalidate_caches()

    if module_name in sys.modules:
        return importlib.reload(sys.modules[module_name])
    return importlib.import_module(module_name)


def load_feature_extractor(extractor_type: str = "default", **kwargs) -> CovMutExFeatureExtractor:
    if extractor_type == "default":
        return DefaultCovMutExFeatureExtractor()

    if extractor_type == "custom":
        return CustomFeatureExtractor(**kwargs)

    if extractor_type == "uploaded":
        module_path = kwargs.get("module_path")
        if not module_path:
            raise ValueError("module_path required for uploaded extractor")
        try:
            module = _load_uploaded_module(module_path)
        except ImportError as exc:
            raise ImportError(
                f"Failed to load feature extractor: {exc}. "
                "Make sure all required helper files are uploaded into the same bundle directory."
            ) from exc

        extractor_dir = os.path.dirname(os.path.abspath(module_path))
        print(f"Loaded uploaded extractor from: {module_path}")
        print(f"Extractor base directory: {extractor_dir}")
        return UploadedExtractorWrapper(module, module_path, extractor_dir)

    raise ValueError(f"Unknown extractor type: {extractor_type}")
