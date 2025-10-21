"""
CovMutEx Model Protocol

This module defines the Protocol for COVID-19 mutation prediction models.
Individual model implementations will be added in the future.
"""

import os
import numpy as np
import tensorflow as tf
from typing import Protocol, Any, runtime_checkable, Optional


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
        
        # eğer çıktı fazla boyutluysa sıkıştır
        if predictions.ndim > 2:
            predictions = np.squeeze(predictions)
        
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


def load_model(model_path: str, model_name: Optional[str] = None, description: Optional[str] = None, source: str = "server"):
    """
    Load a COVID-19 mutation prediction model.
    
    Args:
        model_path: Path to the Keras model file (.keras or .h5) or uploaded file path
        model_name: Optional custom name for the model
        description: Optional description of what the model does
        source: Source of the model - "server" (default), "uploaded", or "registry"
        
    Returns:
        CovMutExKerasModel instance implementing CovMutExModel Protocol
        
    Raises:
        FileNotFoundError: If model file doesn't exist
        ValueError: If model file format is not supported
        
    Example:
        # Load server-side model
        model = load_model("models/balanced_data_model.keras")
        
        # Load uploaded model (from temp directory)
        model = load_model("/tmp/uploaded_model.h5", 
                          model_name="User Uploaded Model",
                          source="uploaded")
    """
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    
    # Validate file extension
    valid_extensions = ['.keras', '.h5', '.pb', '.hdf5']
    file_ext = os.path.splitext(model_path)[1].lower()
    
    if file_ext not in valid_extensions:
        raise ValueError(f"Unsupported model format: {file_ext}. Supported formats: {valid_extensions}")

    return CovMutExKerasModel(model_path, model_name, description, source)
