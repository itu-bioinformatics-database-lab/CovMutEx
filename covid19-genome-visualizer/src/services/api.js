/**
 * CovMutEx API Service
 * 
 * Centralized API calls for the CovMutEx frontend.
 * Handles both JSON and FormData requests.
 */

const API_BASE_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

/**
 * API Service object with all endpoint methods
 */
export const api = {
  /**
   * Run mutation prediction
   * 
   * @param {Object} data - Prediction parameters
   * @param {string} data.nodeId - Variant node ID
   * @param {number} data.elapsedDay - Days elapsed
   * @param {string} data.selectedModel - Model name/path
   * @param {string|null} data.selectedProteinRegion - Protein region (optional)
   * @param {boolean} data.isNewUpload - Whether this is a new model upload
   * @param {File} data.modelFile - Model file (for uploads)
   * @param {File} data.extractorFile - Feature extractor file (optional)
   * @param {File[]} data.helperFiles - Additional helper files (optional)
   * @param {Object} data.customParameters - Custom model parameters
   * @param {string} data.uploadName - Name for the uploaded model
   */
  predict: async (data) => {
    if (data.isNewUpload && data.modelFile) {
      // ========== FormData Request (New Upload) ==========
      const formData = new FormData();
      
      // Required fields
      formData.append("nodeId", data.nodeId || "default_node_id");
      formData.append("elapsedDay", String(data.elapsedDay || 60));
      
      // Upload folder name
      const timestamp = new Date().toISOString().slice(0, 10).replace(/-/g, "");
      const folderName = data.uploadName
        ? `${data.uploadName.replace(/\s+/g, "_")}_${timestamp}`
        : `model_${timestamp}`;
      formData.append("uploadFolderName", folderName);
      
      // Model file (required for upload)
      formData.append("modelFile", data.modelFile);
      
      // Feature extractor (optional)
      if (data.extractorFile) {
        formData.append("extractorFile", data.extractorFile);
      }
      
      // Helper files (optional)
      if (data.helperFiles && data.helperFiles.length > 0) {
        data.helperFiles.forEach((file, index) => {
          formData.append(`helperFile_${index}`, file);
          formData.append(`helperFileName_${index}`, file.name);
        });
      }
      
      // Custom parameters
      if (data.customParameters && Object.keys(data.customParameters).length > 0) {
        formData.append("customParameters", JSON.stringify(data.customParameters));
      }
      
      // Protein region (optional)
      if (data.selectedProteinRegion) {
        formData.append("selectedProteinRegion", data.selectedProteinRegion);
      }

      const response = await fetch(`${API_BASE_URL}/api/predict/`, {
        method: "POST",
        body: formData,
        // Note: Don't set Content-Type header - browser will set it with boundary
      });

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.error || `HTTP ${response.status}`);
      }

      return response.json();
      
    } else {
      // ========== JSON Request (Existing Model) ==========
      const payload = {
        nodeId: data.nodeId || "default_node_id",
        elapsedDay: Number(data.elapsedDay) || 60,
        selectedModel: data.selectedModel || "balanced_data_model",
        selectedProteinRegion: data.selectedProteinRegion || null,
      };

      // Add custom parameters if provided
      if (data.customParameters && Object.keys(data.customParameters).length > 0) {
        payload.customParameters = JSON.stringify(data.customParameters);
      }

      const response = await fetch(`${API_BASE_URL}/api/predict/`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(payload),
      });

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.error || `HTTP ${response.status}`);
      }

      return response.json();
    }
  },

  /**
   * Get list of available models
   * GET /api/models/
   */
  getModels: async () => {
    const response = await fetch(`${API_BASE_URL}/api/models/`);
    
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}`);
    }
    
    return response.json();
  },

  /**
   * Get parameters for a specific model
   * GET /api/model-parameters/?model_name=X
   * 
   * @param {string} modelName - Model name (with or without "uploaded:" prefix)
   */
  getModelParams: async (modelName) => {
    // Remove "uploaded:" prefix if present
    const cleanName = modelName.replace("uploaded:", "");
    
    const response = await fetch(
      `${API_BASE_URL}/api/model-parameters/?model_name=${encodeURIComponent(cleanName)}`
    );
    
    if (!response.ok) {
      // Return empty params if not found (not an error)
      if (response.status === 404) {
        return { parameters: {} };
      }
      throw new Error(`HTTP ${response.status}`);
    }
    
    return response.json();
  },

  /**
   * Generate WebLogo image
   * POST /generate-weblogo/
   * 
   * @param {Object} data - WebLogo parameters
   * @param {number} data.start - Start position
   * @param {number} data.end - End position
   * @param {number[][]} data.probability_matrix - Nx4 matrix
   * @param {string} data.reference_sequence - Reference sequence
   * @param {string[]} data.nucleotide_order - Nucleotide order (default: ['A','T','G','C'])
   */
  generateWeblogo: async (data) => {
    const response = await fetch(`${API_BASE_URL}/generate-weblogo/`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(data),
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      throw new Error(errorData.error || `HTTP ${response.status}`);
    }

    // Returns image/png blob
    return response.blob();
  },
};

/**
 * Transform backend genomeData to frontend format
 * 
 * Backend: { genomeData: [[A_probs], [T_probs], [G_probs], [C_probs]], ... }
 * Frontend: [{ pos, nucleotide, mutationPoss: {A, T, G, C} }, ...]
 * 
 * @param {Object} backendData - Raw backend response
 * @returns {Object} Transformed data for frontend
 */
export function transformGenomeData(backendData) {
  const { genomeData, genomeSequence } = backendData;
  
  if (!genomeData || !Array.isArray(genomeData) || genomeData.length !== 4) {
    console.warn("Invalid genomeData format from backend");
    return [];
  }
  
  const numPositions = genomeData[0].length;
  const transformed = [];
  
  for (let i = 0; i < numPositions; i++) {
    transformed.push({
      pos: i,
      nucleotide: genomeSequence ? genomeSequence[i] : "",
      mutationPoss: {
        A: genomeData[0][i] || 0,
        T: genomeData[1][i] || 0,
        G: genomeData[2][i] || 0,
        C: genomeData[3][i] || 0,
      },
    });
  }
  
  return transformed;
}

/**
 * Transform frontend format back to backend format
 * (Useful if you need to send data back)
 * 
 * @param {Array} frontendData - Array of position objects
 * @returns {Array} 4xN array for backend
 */
export function reverseTransformGenomeData(frontendData) {
  if (!frontendData || !Array.isArray(frontendData)) {
    return [[], [], [], []];
  }
  
  const A_probs = [];
  const T_probs = [];
  const G_probs = [];
  const C_probs = [];
  
  frontendData.forEach((item) => {
    A_probs.push(item.mutationPoss?.A || 0);
    T_probs.push(item.mutationPoss?.T || 0);
    G_probs.push(item.mutationPoss?.G || 0);
    C_probs.push(item.mutationPoss?.C || 0);
  });
  
  return [A_probs, T_probs, G_probs, C_probs];
}

export default api;