"""
CovMutEx Model Protocol

This module defines the Protocol for COVID-19 mutation prediction models.
Individual model implementations will be added in the future.
"""

import os
import numpy as np
import tensorflow as tf
from typing import Protocol, Any, runtime_checkable, Optional

# Try to import PyTorch (optional dependency)
try:
    import torch
    PYTORCH_AVAILABLE = True
except ImportError:
    PYTORCH_AVAILABLE = False
    torch = None

# Try to import joblib for sklearn models (optional dependency)
try:
    import joblib
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False
    joblib = None


@runtime_checkable
class CovMutExModel(Protocol):
    """
    Protocol for COVID-19 mutation prediction models.
    
    All future model implementations must implement these 5 methods.
    This ensures a consistent interface across different prediction approaches.
    """
    
    def metadata(self) -> dict:
        """
        Return model metadata.
        
        Returns:
            dict with model information:
                - name: str - Model name
                - description: str - What the model does
                - version: str - Model version
                - model_type: str - Type of model
                - Any other relevant metadata
        """
        ...
    
    def input_schema(self) -> dict:
        """
        Define expected input schema.
        
        Returns:
            dict describing required and optional inputs:
                - genome_id: string (optional) - Variant identifier
                - elapsed_days: int (optional) - Days since reference
                - region: tuple (optional) - Genomic region (start, end)
                - features: np.ndarray - Genomic features
        """
        ...
    
    def preprocess(self, inputs: dict) -> Any:
        """
        Preprocess inputs before prediction.
        
        Args:
            inputs: dict containing:
                - features: np.ndarray - Required genomic features
                - genome_id: str - Optional variant identifier
                - elapsed_days: int - Optional temporal information
                - region: tuple - Optional genomic region
                - Any other model-specific inputs
            
        Returns:
            Preprocessed data ready for model.predict()
            Can be np.ndarray, list of arrays, or any model-specific format
        """
        ...
    
    def predict(self, batch: Any) -> np.ndarray:
        """
        Run model prediction.
        
        Args:
            batch: Preprocessed input from preprocess() method
            
        Returns:
            numpy array of predictions
            Shape depends on model type (e.g., (N, 1), (N, 2), etc.)
        """
        ...
    
    def postprocess(self, raw: Any) -> dict:
        """
        Postprocess raw predictions.
        
        Args:
            raw: Raw predictions from predict() method
            
        Returns:
            dict with processed predictions and metadata:
                - predictions: np.ndarray - Processed predictions
                - shape: tuple - Shape of predictions
                - prediction_type: str - Type of predictions
                - interpretation: str - How to interpret the results
                - num_positions: int - Number of genomic positions
                - Any other relevant information
        """
        ...


class CovMutExKerasModel:
    """
    Base implementation of CovMutExModel Protocol using Keras/TensorFlow.
    
    This class can be used directly or extended for specific model types.
    All future model implementations should either use or extend this class.
    """
    
    def __init__(self, model_path: str, model_name: str = None, description: str = None, source: str = "server"):
        """
        Initialize with a Keras model.
        
        Args:
            model_path: Path to the .keras or .h5 model file
            model_name: Optional custom name for the model
            description: Optional description of the model
            source: Source of the model - "server", "uploaded", or "registry"
        """
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model file not found: {model_path}")
        
        self.model_path = model_path
        self.model_name = model_name or os.path.basename(model_path).replace('.keras', '').replace('.h5', '')
        self.description = description or "COVID-19 mutation prediction model"
        self.source = source
        
        # Keras modeli yüklemek için
        self.keras_model = tf.keras.models.load_model(model_path, compile=False)
        
        # çok girişli mi tek girişli mi olduğunu belirle
        self._is_multi_input = isinstance(self.keras_model.input_shape, list)
        
        # output türünü belirle
        output_shape = self.keras_model.output_shape
        if output_shape[-1] == 1:
            self._output_type = "single-output"
        elif output_shape[-1] == 2:
            self._output_type = "dual-output"
        else:
            self._output_type = "multi-class"
    
    def metadata(self) -> dict:
        """Return model metadata."""
        return {
            "name": self.model_name,
            "description": self.description,
            "version": "1.0",
            "model_type": "multi-input" if self._is_multi_input else "single-input",
            "output_type": self._output_type,
            "input_shape": str(self.keras_model.input_shape),
            "output_shape": str(self.keras_model.output_shape),
            "is_multi_input": self._is_multi_input,
            "framework": "keras/tensorflow",
            "source": self.source
        }
    
    def input_schema(self) -> dict:
        """Define expected input schema."""
        return {
            "required": {
                "features": f"numpy array of shape (N, 205) - genomic features"
            },
            "optional": {
                "genome_id": "string - variant identifier (e.g., 'hCoV-19/USA/CA-123')",
                "elapsed_days": "int - days since reference date",
                "region": "tuple (start, end) - genomic region to analyze",
                "num_inputs": f"int - number of input copies for multi-input models (default: 10)"
            },
            "notes": f"This model is {'multi-input' if self._is_multi_input else 'single-input'}"
        }
    
    def preprocess(self, inputs: dict) -> Any:
        """
        Preprocess inputs for Keras model.
        
        Args:
            inputs: dict with 'features' (required) and optional metadata
            
        Returns:
            Single np.ndarray or list of arrays (for multi-input models)
        """
        features = inputs.get('features')
        if features is None:
            raise ValueError("'features' key is required in inputs dict")
        
        # fazla boyutları temizle
        if isinstance(features, np.ndarray) and features.ndim > 2:
            features = np.squeeze(features)
        
        # çok girişli modeller için girişleri çoğalt
        if self._is_multi_input:
            num_inputs = inputs.get('num_inputs', 10)
            return [features for _ in range(num_inputs)]
        else:
            return features
    
    def predict(self, batch: Any) -> np.ndarray:
        """
        Run Keras model prediction.
        
        Args:
            batch: Preprocessed features (array or list of arrays)
            
        Returns:
            numpy array of predictions
        """
        predictions = self.keras_model.predict(batch, verbose=0)
        
        # Remove batch dimension but keep prediction dimensions
        # E.g., (1, 29904, 1) -> (29904, 1) NOT (29904,)
        if predictions.ndim > 2:
            predictions = np.squeeze(predictions, axis=0)
        
        return predictions
    
    def postprocess(self, raw: Any) -> dict:
        """
        Postprocess Keras model predictions.
        
        Args:
            raw: Raw predictions from predict()
            
        Returns:
            dict with processed predictions and metadata
        """
        if not isinstance(raw, np.ndarray):
            raw = np.array(raw)
        
        # değerlendirme için yorum ekle
        if self._output_type == "single-output":
            interpretation = "Single mutation probability per position"
        elif self._output_type == "dual-output":
            interpretation = "Column 0: P(no mutation), Column 1: P(mutation)"
        else:
            interpretation = "Multi-class probabilities per position"
        
        return {
            "predictions": raw,
            "shape": raw.shape,
            "prediction_type": self._output_type,
            "interpretation": interpretation,
            "num_positions": raw.shape[0] if len(raw.shape) > 0 else 0
        }


class CovMutExPyTorchModel:
    """
    PyTorch implementation of CovMutExModel Protocol.
    
    Supports PyTorch models (.pt, .pth files) for COVID-19 mutation prediction.
    """
    
    def __init__(self, model_path: str, model_name: str = None, description: str = None, source: str = "server"):
        """
        Initialize with a PyTorch model.
        
        Args:
            model_path: Path to the .pt or .pth model file
            model_name: Optional custom name for the model
            description: Optional description of the model
            source: Source of the model - "server", "uploaded", or "registry"
        """
        if not PYTORCH_AVAILABLE:
            raise ImportError("PyTorch is not installed. Install with: pip install torch")
        
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model file not found: {model_path}")
        
        self.model_path = model_path
        self.model_name = model_name or os.path.basename(model_path).replace('.pt', '').replace('.pth', '')
        self.description = description or "COVID-19 mutation prediction model (PyTorch)"
        self.source = source
        
        # Load PyTorch model
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.torch_model = torch.load(model_path, map_location=self.device)
        
        # Set to evaluation mode
        if hasattr(self.torch_model, 'eval'):
            self.torch_model.eval()
        
        # Detect input/output shapes (try to infer from model)
        self._detect_model_properties()
    
    def _detect_model_properties(self):
        """Detect model properties by doing a test inference."""
        try:
            # Create dummy input to detect shapes
            dummy_input = torch.randn(1, 205, device=self.device)
            
            with torch.no_grad():
                dummy_output = self.torch_model(dummy_input)
            
            # Detect if multi-input (this is a simplification)
            self._is_multi_input = False  # PyTorch models typically don't use multi-input like Keras
            
            # Detect output type
            output_dim = dummy_output.shape[-1] if len(dummy_output.shape) > 1 else 1
            if output_dim == 1:
                self._output_type = "single-output"
            elif output_dim == 2:
                self._output_type = "dual-output"
            else:
                self._output_type = "multi-class"
            
            self._input_shape = f"(None, 205)"
            self._output_shape = f"(None, {output_dim})"
            
        except Exception as e:
            print(f"Warning: Could not auto-detect model properties: {e}")
            self._is_multi_input = False
            self._output_type = "unknown"
            self._input_shape = "(None, 205)"
            self._output_shape = "(None, ?)"
    
    def metadata(self) -> dict:
        """Return model metadata."""
        return {
            "name": self.model_name,
            "description": self.description,
            "version": "1.0",
            "model_type": "multi-input" if self._is_multi_input else "single-input",
            "output_type": self._output_type,
            "input_shape": self._input_shape,
            "output_shape": self._output_shape,
            "is_multi_input": self._is_multi_input,
            "framework": "pytorch",
            "source": self.source,
            "device": str(self.device)
        }
    
    def input_schema(self) -> dict:
        """Define expected input schema."""
        return {
            "required": {
                "features": "numpy array of shape (N, 205) - genomic features"
            },
            "optional": {
                "genome_id": "string - variant identifier (e.g., 'hCoV-19/USA/CA-123')",
                "elapsed_days": "int - days since reference date",
                "region": "tuple (start, end) - genomic region to analyze"
            },
            "notes": f"This model is a PyTorch model running on {self.device}"
        }
    
    def preprocess(self, inputs: dict) -> Any:
        """
        Preprocess inputs for PyTorch model.
        
        Args:
            inputs: dict with 'features' (required) and optional metadata
            
        Returns:
            PyTorch tensor ready for inference
        """
        features = inputs.get('features')
        if features is None:
            raise ValueError("'features' key is required in inputs dict")
        
        # Convert to numpy if not already
        if not isinstance(features, np.ndarray):
            features = np.array(features)
        
        # Clean extra dimensions
        if features.ndim > 2:
            features = np.squeeze(features)
        
        # Convert to PyTorch tensor
        tensor = torch.tensor(features, dtype=torch.float32, device=self.device)
        
        return tensor
    
    def predict(self, batch: Any) -> np.ndarray:
        """
        Run PyTorch model prediction.
        
        Args:
            batch: Preprocessed features (PyTorch tensor)
            
        Returns:
            numpy array of predictions
        """
        with torch.no_grad():
            predictions = self.torch_model(batch)
            
            # Convert to numpy
            predictions = predictions.cpu().numpy()
            
            # Clean extra dimensions if needed
            if predictions.ndim > 2:
                predictions = np.squeeze(predictions)
        
        return predictions
    
    def postprocess(self, raw: Any) -> dict:
        """
        Postprocess PyTorch model predictions.
        
        Args:
            raw: Raw predictions from predict()
            
        Returns:
            dict with processed predictions and metadata
        """
        if not isinstance(raw, np.ndarray):
            raw = np.array(raw)
        
        # Determine interpretation based on output type
        if self._output_type == "single-output":
            interpretation = "Single mutation probability per position"
        elif self._output_type == "dual-output":
            interpretation = "Column 0: P(no mutation), Column 1: P(mutation)"
        else:
            interpretation = "Multi-class probabilities per position"
        
        return {
            "predictions": raw,
            "shape": raw.shape,
            "prediction_type": self._output_type,
            "interpretation": interpretation,
            "num_positions": raw.shape[0] if len(raw.shape) > 0 else 0
        }


class CovMutExSklearnModel:
    """
    Scikit-learn implementation of CovMutExModel Protocol.
    
    Supports sklearn models (joblib .pkl files) for sequence mutation prediction.
    Primarily used for Influenza and other non-COVID models.
    """
    
    def __init__(self, model_path: str, model_name: str = None, description: str = None, source: str = "server"):
        """
        Initialize with a Scikit-learn model.
        
        Args:
            model_path: Path to the .pkl model file
            model_name: Optional custom name for the model
            description: Optional description of the model
            source: Source of the model - "server", "uploaded", or "registry"
        """
        if not JOBLIB_AVAILABLE:
            raise ImportError("Joblib is not installed. Install with: pip install joblib")
        
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model file not found: {model_path}")
        
        self.model_path = model_path
        self.model_name = model_name or os.path.basename(model_path).replace('.pkl', '')
        self.description = description or "Sequence mutation prediction model (Sklearn)"
        self.source = source
        
        # Load sklearn model with metadata
        model_data = joblib.load(model_path)
        
        # Handle both formats: direct model or dict with metadata
        if isinstance(model_data, dict):
            self.sklearn_model = model_data.get('model')
            # Legacy fields
            self.max_len = model_data.get('max_len', None)
            self.n_sequences = model_data.get('n_sequences', None)
            # New fields for position-wise mutation models
            self.genome_length = model_data.get('genome_length', None)
            self.encoded_feature_length = model_data.get('encoded_feature_length', None)
            self.mutation_rate = model_data.get('mutation_rate', None)
            self.model_type_detail = model_data.get('model_type', None)
        else:
            self.sklearn_model = model_data
            self.max_len = None
            self.n_sequences = None
            self.genome_length = None
            self.encoded_feature_length = None
            self.mutation_rate = None
            self.model_type_detail = None
        
        # Detect model type
        self._is_multi_input = False  # Sklearn models are typically single-input
        self._output_type = "sequence-prediction"  # Default for sequence models
    
    def metadata(self) -> dict:
        """Return model metadata."""
        metadata = {
            "name": self.model_name,
            "description": self.description,
            "version": "1.0",
            "model_type": "single-input",
            "output_type": self._output_type,
            "is_multi_input": self._is_multi_input,
            "framework": "sklearn",
            "source": self.source,
            "model_class": type(self.sklearn_model).__name__
        }
        
        # Add input/output shapes for compatibility
        if self.encoded_feature_length is not None:
            metadata["input_shape"] = f"(None, {self.encoded_feature_length})"
        else:
            metadata["input_shape"] = "(None, ?)"
            
        if self.genome_length is not None:
            metadata["output_shape"] = f"(None, {self.genome_length}, 1)"
        else:
            metadata["output_shape"] = "(None, ?, 1)"
        
        # Legacy metadata fields
        if self.max_len is not None:
            metadata["max_sequence_length"] = self.max_len
        if self.n_sequences is not None:
            metadata["training_sequences"] = self.n_sequences
        
        # New metadata fields for position-wise mutation models
        if self.genome_length is not None:
            metadata["genome_length"] = self.genome_length
        if self.encoded_feature_length is not None:
            metadata["encoded_feature_length"] = self.encoded_feature_length
        if self.mutation_rate is not None:
            metadata["mutation_rate"] = self.mutation_rate
        if self.model_type_detail is not None:
            metadata["model_type_detail"] = self.model_type_detail
            
        return metadata
    
    def input_schema(self) -> dict:
        """Define expected input schema."""
        return {
            "required": {
                "features": "numpy array - encoded sequence features"
            },
            "optional": {
                "genome_id": "string - sequence identifier",
                "max_len": "int - maximum sequence length for padding"
            },
            "notes": "Sklearn model for sequence-to-sequence prediction"
        }
    
    def preprocess(self, inputs: dict) -> Any:
        """
        Preprocess inputs for Sklearn model.
        
        Args:
            inputs: dict with 'features' (required) and optional metadata
            
        Returns:
            numpy array ready for prediction
        """
        features = inputs.get('features')
        if features is None:
            raise ValueError("'features' key is required in inputs dict")
        
        # Convert to numpy if not already
        if not isinstance(features, np.ndarray):
            features = np.array(features)
        
        # Clean extra dimensions
        if features.ndim > 2:
            features = np.squeeze(features)
        
        # Ensure 2D array for sklearn (add batch dimension if needed)
        if features.ndim == 1:
            features = features.reshape(1, -1)
        
        return features
    
    def predict(self, batch: Any) -> np.ndarray:
        """
        Run Sklearn model prediction.
        
        Args:
            batch: Preprocessed features (numpy array of shape (1, encoded_length))
            
        Returns:
            numpy array of predictions (genome_length, 1) - mutation probabilities per position
        """
        predictions = self.sklearn_model.predict(batch)
        
        # predictions shape: (1, genome_length) for MultiOutputRegressor
        # Reshape to (genome_length, 1) to match COVID model format
        if predictions.ndim == 2 and predictions.shape[0] == 1:
            predictions = predictions.T  # (1, N) → (N, 1)
        elif predictions.ndim == 1:
            predictions = predictions.reshape(-1, 1)
        
        # Clip to valid probability range [0, 1]
        predictions = np.clip(predictions, 0.0, 1.0)
        
        return predictions
    
    def postprocess(self, raw: Any) -> dict:
        """
        Postprocess Sklearn model predictions.
        
        Args:
            raw: Raw predictions from predict()
            
        Returns:
            dict with processed predictions and metadata
        """
        if not isinstance(raw, np.ndarray):
            raw = np.array(raw)
        
        return {
            "predictions": raw,
            "shape": raw.shape,
            "prediction_type": self._output_type,
            "interpretation": "Encoded sequence prediction (requires decoding)",
            "num_positions": raw.shape[1] if len(raw.shape) > 1 else raw.shape[0]
        }


def load_model(model_path: str, model_name: Optional[str] = None, description: Optional[str] = None, source: str = "server"):
    """
    Load a COVID-19 mutation prediction model.
    
    Automatically detects model type based on file extension and loads the appropriate model class.
    
    Args:
        model_path: Path to the model file (.keras, .h5, .pt, .pth, etc.)
        model_name: Optional custom name for the model
        description: Optional description of what the model does
        source: Source of the model - "server" (default), "uploaded", or "registry"
        
    Returns:
        Model wrapper instance implementing CovMutExModel Protocol
        - CovMutExKerasModel for Keras/TensorFlow models
        - CovMutExPyTorchModel for PyTorch models
        
    Raises:
        FileNotFoundError: If model file doesn't exist
        ValueError: If model file format is not supported
        ImportError: If required framework is not installed
        
    Examples:
        # Load Keras model
        model = load_model("models/balanced_data_model.keras")
        
        # Load PyTorch model
        model = load_model("models/pytorch_model.pt", source="server")
        
        # Load uploaded H5 model
        model = load_model("/tmp/uploaded_model.h5", 
                          model_name="User Model",
                          source="uploaded")
    """
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    
    # Get file extension
    file_ext = os.path.splitext(model_path)[1].lower()
    
    # PyTorch models
    if file_ext in ['.pt', '.pth']:
        if not PYTORCH_AVAILABLE:
            raise ImportError(
                f"PyTorch model detected ({file_ext}) but PyTorch is not installed. "
                "Install with: pip install torch"
            )
        return CovMutExPyTorchModel(model_path, model_name, description, source)
    
    # Sklearn models (Joblib pickle)
    elif file_ext in ['.pkl', '.pickle', '.joblib']:
        if not JOBLIB_AVAILABLE:
            raise ImportError(
                f"Sklearn model detected ({file_ext}) but joblib is not installed. "
                "Install with: pip install joblib"
            )
        return CovMutExSklearnModel(model_path, model_name, description, source)
    
    # Keras/TensorFlow models
    elif file_ext in ['.keras', '.h5', '.hdf5', '.pb']:
        return CovMutExKerasModel(model_path, model_name, description, source)
    
    # Unsupported format
    else:
        supported_formats = ['.keras', '.h5', '.hdf5', '.pb']
        if PYTORCH_AVAILABLE:
            supported_formats.extend(['.pt', '.pth'])
        if JOBLIB_AVAILABLE:
            supported_formats.extend(['.pkl', '.pickle', '.joblib'])
        
        raise ValueError(
            f"Unsupported model format: {file_ext}\n"
            f"Supported formats: {supported_formats}"
        )
