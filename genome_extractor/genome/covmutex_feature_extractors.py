from typing import Protocol, Dict, List, Tuple, Optional
import numpy as np

class CovMutExFeatureExtractor(Protocol):
    """
    Protocol for feature extraction implementations.
    
    Any feature extractor must implement these methods to be compatible
    with CovMutEx prediction pipeline.
    """
    
    def extract_features(
        self,
        genome_seq: str,
        mutations: List[Tuple[int, str, str, int, str]],
        node_id: str,
        elapsed_day: int,
        protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
        k: int = 30,
        **kwargs
    ) -> np.ndarray:
        """
        Extract features for all positions in the genome.
        
        Args:
            genome_seq: Reference or variant genome sequence
            mutations: List of (position, ref, alt, aa_pos, aa_change)
            node_id: Phylogenetic node identifier
            elapsed_day: Days elapsed for temporal prediction
            protein_regions: Optional dict of protein regions to filter
            k: K-mer window size
            **kwargs: Additional extractor-specific parameters
                     (e.g., cache_path, codon_mapping_path, config_file for default extractor)
            
        Returns:
            numpy array of shape (N, feature_dim) where:
            - N = number of positions (29904 or protein region length)
            - feature_dim = number of features per position
        """
        ...
    
    def get_feature_dimension(self) -> int:
        """
        Return the number of features per position.
        
        Returns:
            int: Feature vector dimension (e.g., 205 for default extractor)
        """
        ...
    
    def get_feature_description(self) -> Dict[str, str]:
        """
        Return description of feature composition.
        
        Returns:
            dict: Mapping of feature ranges to their descriptions
            Example:
            {
                "0-29": "K-mer nucleotides (k=30)",
                "30": "Center nucleotide",
                "31": "Mutated nucleotide",
                ...
            }
        """
        ...
    
    def get_metadata(self) -> Dict[str, any]:
        """
        Return metadata about the feature extractor.
        
        Returns:
            dict: Metadata including version, k-mer size, etc.
        """
        ...


class DefaultCovMutExFeatureExtractor:
    """
    Default feature extractor implementation following the paper's methodology.
    
    Feature composition (k=30):
    - Features 0 to k-1 (0-29): K-mer nucleotides
    - Feature k (30): Original nucleotide at position
    - Feature k+1 (31): Mutated nucleotide
    - Feature k+2 (32): Position index
    - Feature k+3 (33): PAM250 score (nucleotide level)
    - Feature k+4 (34): Original amino acid
    - Feature k+5 (35): Mutated amino acid
    - Feature k+6 (36): PAM250 score (amino acid level)
    - Feature k+7 (37): Elapsed days
    - Feature k+8 (38): Phylogenetic depth
    - Feature k+9 (39): Synonymous indicator (1=yes, 0=no)
    - Feature k+10 (40): ORF/Protein region name
    - Features k+11 to k+28 (41-58): AA biochemical properties (18 features)
      - Hydrophobicity (original, mutated)
      - Polarity (original, mutated)
      - Iso-electricity (original, mutated)
      - Volume (original, mutated)
      - Weight (original, mutated)
      - pKa (original, mutated)
      - pKb (original, mutated)
      - pKx (original, mutated)
      - pI (original, mutated)
    
    After preprocessing (one-hot encoding + standardization): 205 dimensions
    """
    
    def __init__(self):
        """
        Initialize the default feature extractor.
        """
        self._feature_dim = 205  # After preprocessing
        
    def extract_features(
        self,
        genome_seq: str,
        mutations: List[Tuple[int, str, str, int, str]],
        node_id: str,
        elapsed_day: int,
        protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
        k: int = 30,
        **kwargs
    ) -> np.ndarray:
        """
        Extract features using the default methodology.
        
        Args:
            genome_seq: Genome sequence
            mutations: List of mutations
            node_id: Node identifier
            elapsed_day: Elapsed days
            protein_regions: Optional protein regions
            k: K-mer size
            **kwargs: Additional parameters:
                cache_path: Optional HDF5 cache path
                codon_mapping_path: Optional path to codon mapping JSON
                config_file: Optional config dict
        """
        import os
        from . import feature_extractor_updated as feu
        
        # Extract kwargs
        cache_path = kwargs.get('cache_path', None)
        codon_mapping_path = kwargs.get('codon_mapping_path', None)
        config_file = kwargs.get('config_file', None)
        
        # Calculate paths dynamically if not provided
        if codon_mapping_path is None:
            feu_dir = os.path.dirname(os.path.abspath(feu.__file__))
            codon_mapping_path = os.path.join(feu_dir, "codon_aa_mapping.json")
        
        if config_file is None:
            config_file = feu.configs()
        
        # Extract features using feature_extractor_updated
        features = feu.cache_precomputed_features(
            cache_path=cache_path,
            genome_seq=genome_seq,
            mutations=mutations,
            codon_mapper=codon_mapping_path,
            config_file=config_file,
            node_ids=[node_id],
            elapsed_day=elapsed_day,
            protein_regions=protein_regions
        )
        
        return features
    
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
            "37": "Elapsed days (temporal feature)",
            "38": "Phylogenetic depth",
            "39": "Synonymous indicator (1=syn, 0=nonsyn)",
            "40": "Protein region/ORF name",
            "41-58": "AA biochemical properties (18 features)",
            "59-204": "Preprocessed features (one-hot + standardized)"
        }
    
    def get_metadata(self) -> Dict[str, any]:
        return {
            "name": "DefaultCovMutExFeatureExtractor",
            "version": "1.0",
            "k_mer_size": 30,
            "feature_dimension": self._feature_dim,
            "preprocessing": "one-hot encoding + standardization",
            "paper_reference": "CovMutEx methodology"
        }
    
# NEW: gerek yok
class CustomFeatureExtractor:
    """
    Custom feature extractor template.
    
    REQUIRED METHODS:
    - extract_features(): Main feature extraction logic
    - get_feature_dimension(): Return feature vector size
    - get_feature_description(): Describe your features
    - get_metadata(): Provide metadata about your extractor
    """
    
    def __init__(self, **config):
        """
        Initialize your feature extractor with custom configuration.
        
        Args:
            **config: Your custom configuration parameters
        """
        self.config = config
        # TODO: Initialize your resources (e.g., load config files, models, etc.)
        
    def extract_features(
        self,
        genome_seq: str,
        mutations: List[Tuple[int, str, str, int, str]],
        node_id: str,
        elapsed_day: int,
        protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
        k: int = 30,
        **kwargs
    ) -> np.ndarray:
        """
        Extract features for genome positions.
        
        Args:
            genome_seq: Genome sequence string (e.g., "ATGCATGC...")
            mutations: List of mutations as (nt_pos, ref, alt, aa_pos, aa_change)
            node_id: Phylogenetic node identifier
            elapsed_day: Days elapsed for temporal prediction
            protein_regions: Optional dict of protein regions {"S": (21563, 25384), ...}
            k: K-mer window size
            **kwargs: Additional custom parameters
            
        Returns:
            numpy array of shape (N, feature_dim) where:
            - N = genome length (29904) or protein region length
            - feature_dim = your feature vector dimension
            
        Example output shape: (29904, 100) for 100-dimensional features
        """
        # TODO: Implement your feature extraction logic
        
        genome_length = len(genome_seq)
        feature_dim = self.get_feature_dimension()
        
        # Example skeleton:
        features = []
        for position in range(genome_length):
            # Extract features for this position
            feature_vector = self._extract_position_features(
                genome_seq, position, mutations, node_id, elapsed_day, k
            )
            features.append(feature_vector)
        
        features = np.array(features, dtype=np.float32)
        
        # Filter by protein region if specified
        if protein_regions:
            # TODO: Filter positions within protein_regions
            pass
        
        return features
    
    def _extract_position_features(
        self, 
        genome_seq: str, 
        position: int,
        mutations: List[Tuple],
        node_id: str,
        elapsed_day: int,
        k: int
    ) -> np.ndarray:
        """
        Extract features for a single position.
        
        TODO: Implement your position-level feature extraction
        
        Returns:
            numpy array of shape (feature_dim,)
        """
        # Example: Extract k-mer
        mid_point = k // 2
        padded_seq = "-" * mid_point + genome_seq + "-" * mid_point
        k_mer = padded_seq[position:position + k]
        
        # TODO: Add more features (amino acids, biochemical properties, etc.)
        
        feature_vector = np.zeros(self.get_feature_dimension())
        # TODO: Populate feature_vector
        
        return feature_vector
    
    def get_feature_dimension(self) -> int:
        """
        Return the number of features per position.
        
        Returns:
            int: Feature vector dimension
            
        Example: return 100 for 100-dimensional features
        """
        # TODO: Return your feature dimension
        return 205  # Replace with your dimension
    
    def get_feature_description(self) -> Dict[str, str]:
        """
        Describe the composition of your feature vector.
        
        Returns:
            dict: Mapping of feature indices/ranges to descriptions
            
        Example:
            {
                "0-29": "K-mer nucleotides",
                "30": "Position-specific feature X",
                "31-50": "Biochemical properties",
                ...
            }
        """
        # TODO: Document your features
        return {
            "0-29": "Your features here",
            # Add more descriptions
        }
    
    def get_metadata(self) -> Dict[str, any]:
        """
        Return metadata about your feature extractor.
        
        Returns:
            dict: Metadata including name, version, parameters, etc.
            
        Example:
            {
                "name": "MyCustomExtractor",
                "version": "1.0",
                "author": "Your Name",
                "feature_dimension": 100,
                "k_mer_size": 30,
                ...
            }
        """
        # TODO: Provide metadata
        return {
            "name": "CustomFeatureExtractor",
            "version": "1.0",
            "feature_dimension": self.get_feature_dimension(),
            "config": self.config
        }


def load_feature_extractor(
    extractor_type: str = "default",
    **kwargs
) -> CovMutExFeatureExtractor:
    """
    Factory function to load a feature extractor.
    
    Args:
        extractor_type: Type of extractor ("default", "custom", "uploaded", etc.)
        **kwargs: Additional arguments for extractor initialization
        
    Returns:
        Feature extractor instance implementing CovMutExFeatureExtractor protocol
        
    Examples:
        >>> # Use default extractor
        >>> extractor = load_feature_extractor("default")
        
        >>> # Use custom extractor with specific config
        >>> extractor = load_feature_extractor("custom", config_path="my_config.json")
        
        >>> # Use uploaded extractor module
        >>> extractor = load_feature_extractor("uploaded", module_path="/path/to/extractor.py")
    """
    if extractor_type == "default":
        return DefaultCovMutExFeatureExtractor(**kwargs)
    
    elif extractor_type == "custom":
        # Kullanıcının kendi extractor'ını yükle
        return CustomFeatureExtractor(**kwargs)
    
    elif extractor_type == "uploaded":
        # Dinamik olarak yüklenmiş modül
        module_path = kwargs.get("module_path")
        if not module_path:
            raise ValueError("module_path required for uploaded extractor")
        
        import importlib.util
        import os
        
        spec = importlib.util.spec_from_file_location("uploaded_extractor", module_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        
        # Extractor'ın bulunduğu klasör (model ile aynı klasör)
        extractor_dir = os.path.dirname(os.path.abspath(module_path))
        
        # Gerekli fonksiyonları kontrol et
        required_functions = ['extract_features', 'get_feature_dimension', 
                            'get_feature_description', 'get_metadata']
        missing = [fn for fn in required_functions if not hasattr(module, fn)]
        
        if missing:
            raise ValueError(
                f"Uploaded extractor missing required functions: {missing}. "
                f"Please provide: {required_functions}"
            )
        
        # Wrapper class oluştur - modülün fonksiyonlarını Protocol'e adapte et
        class UploadedExtractorWrapper:
            def __init__(self, module, base_dir):
                self._module = module
                self._module_path = module_path
                self._base_dir = base_dir  # Model/extractor klasörü
            
            def extract_features(
                self,
                genome_seq: str,
                mutations: List[Tuple[int, str, str, int, str]],
                node_id: str,
                elapsed_day: int,
                protein_regions: Optional[Dict[str, Tuple[int, int]]] = None,
                k: int = 30,
                **kwargs  # Custom parameters to pass to uploaded extractor
            ) -> np.ndarray:
                # Modülün fonksiyonunu çağır - Protocol parametreleri + custom parameters
                # Kullanıcı os.path.dirname(__file__) ile kendi dosyalarına erişebilir
                # kwargs custom parameters olarak iletilir
                return self._module.extract_features(
                    genome_seq=genome_seq,
                    mutations=mutations,
                    node_id=node_id,
                    elapsed_day=elapsed_day,
                    protein_regions=protein_regions,
                    k=k,
                    **kwargs  # Pass custom parameters to user's extractor
                )
            
            def get_feature_dimension(self) -> int:
                return self._module.get_feature_dimension()
            
            def get_feature_description(self) -> Dict[str, str]:
                return self._module.get_feature_description()
            
            def get_metadata(self) -> Dict[str, any]:
                metadata = self._module.get_metadata()
                metadata['module_path'] = self._module_path
                metadata['base_directory'] = self._base_dir
                return metadata
        
        print(f"Loaded uploaded extractor from: {module_path}")
        print(f"Extractor base directory: {extractor_dir}")
        return UploadedExtractorWrapper(module, extractor_dir)
    
    else:
        raise ValueError(f"Unknown extractor type: {extractor_type}")