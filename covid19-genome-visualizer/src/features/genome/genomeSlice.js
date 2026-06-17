import { createSlice, createAsyncThunk } from "@reduxjs/toolkit";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

// ============================================
// ASYNC THUNKS
// ============================================

/**
 * Fetch available models from the server
 * GET /api/models/
 */
export const fetchAvailableModels = createAsyncThunk(
  "genome/fetchAvailableModels",
  async (_, { rejectWithValue }) => {
    try {
      const response = await fetch(`${API_URL}/api/models/`);
      if (!response.ok) {
        throw new Error("Failed to fetch models");
      }
      const data = await response.json();
      return data.available_models || [];
    } catch (error) {
      console.warn("Could not fetch models from API:", error.message);
      return rejectWithValue(error.message);
    }
  }
);

/**
 * Fetch model parameters for a specific model
 * GET /api/model-parameters/?model_name=X
 */
export const fetchModelParameters = createAsyncThunk(
  "genome/fetchModelParameters",
  async (modelName, { rejectWithValue }) => {
    try {
      const cleanName = modelName.replace("uploaded:", "");
      const response = await fetch(
        `${API_URL}/api/model-parameters/?model_name=${encodeURIComponent(cleanName)}`
      );
      if (!response.ok) {
        return { parameters: {} };
      }
      const data = await response.json();
      return data;
    } catch (error) {
      return rejectWithValue(error.message);
    }
  }
);

/**
 * Main prediction thunk
 * POST /api/predict/
 */
export const fetchPrediction = createAsyncThunk(
  "genome/fetchPrediction",
  async (params, { rejectWithValue }) => {
    try {
      const {
        nodeId,
        elapsedDay,
        selectedModel,
        selectedProteinRegion,
        isNewUpload,
        uploadName,
        modelFile,
        extractorFile,
        helperFiles,
        customParameters,
        variant,
      } = params;

      let response;

      if (isNewUpload && modelFile) {
        // ========== FormData Request (New Upload) ==========
        const formData = new FormData();
        
        formData.append("nodeId", nodeId || "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07");
        formData.append("elapsedDay", String(elapsedDay || 60));
        
        // Use upload name as folder name (no timestamp)
        const folderName = uploadName 
          ? uploadName.replace(/\s+/g, "_")
          : `model_${Date.now()}`;
        formData.append("uploadFolderName", folderName);
        
        formData.append("modelFile", modelFile);
        
        if (extractorFile) {
          formData.append("extractorFile", extractorFile);
        }
        
        if (helperFiles && helperFiles.length > 0) {
          helperFiles.forEach((file, index) => {
            formData.append(`helperFile_${index}`, file);
            formData.append(`helperFileName_${index}`, file.name);
          });
        }
        
        if (customParameters && Object.keys(customParameters).length > 0) {
          formData.append("customParameters", JSON.stringify(customParameters));
        }
        
        if (selectedProteinRegion) {
          formData.append("selectedProteinRegion", selectedProteinRegion);
        }

        response = await fetch(`${API_URL}/api/predict/`, {
          method: "POST",
          body: formData,
        });

      } else {
        // ========== JSON Request (Existing Model) ==========
        const payload = {
          nodeId: nodeId || "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07",
          elapsedDay: Number(elapsedDay) || 60,
          selectedModel: selectedModel || "balanced_data_model",
          selectedProteinRegion: selectedProteinRegion || null,
        };

        // Prediction-time variant override for influenza bundles. Backend
        // routes it via _resolve_organism_data → find_variant_subtype, so the
        // payload's HA region annotation matches the chosen strain's subtype.
        if (variant) {
          payload.variant = variant;
        }

        if (customParameters && Object.keys(customParameters).length > 0) {
          payload.customParameters = JSON.stringify(customParameters);
        }

        response = await fetch(`${API_URL}/api/predict/`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(payload),
        });
      }

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.error || `HTTP ${response.status}`);
      }

      const data = await response.json();
      
      return transformBackendResponse(data, selectedProteinRegion);
      
    } catch (error) {
      console.error("Prediction failed:", error);
      return rejectWithValue(error.message);
    }
  }
);

/**
 * Transform backend response to frontend-compatible format.
 * 
 * IMPORTANT: We keep TWO formats:
 * 
 * 1) `genomeDataRaw` (raw [4][N] array format) 
 *    → Used by GenomeChart (Recharts.js) which does genomeData[nucIdx].map(...)
 * 
 * 2) `dataset` (transformed [{pos, nucleotide, mutationPoss: {A,T,G,C}}] format)
 *    → Used by BarChart and BarChart2 (Chart.js) which reads entry.mutationPoss
 */
function transformBackendResponse(data, selectedProteinRegion) {
  const {
    genomeData,
    genomeSequence,
    protein_mutation_probs,
    proteinRegionPossibilities,
    model_metadata,
    extractor_metadata,
    predictionPayload,
    nodeId,
    elapsedDay,
    selectedModel,
  } = data;

  // === FORMAT 1: Raw [4][N] array (for GenomeChart / Recharts) ===
  // Backend sends this only for categorical_per_position payloads.
  const genomeDataRaw = genomeData;

  // === FORMAT 2: Transformed object array (for BarChart / BarChart2) ===
  let dataset = [];

  if (genomeData && Array.isArray(genomeData) && genomeData.length === 4) {
    const numPositions = genomeData[0].length;

    for (let i = 0; i < numPositions; i++) {
      const nucleotide = genomeSequence ? genomeSequence[i] : "";
      dataset.push({
        pos: i,
        nucleotide: nucleotide,
        mutationPoss: {
          A: genomeData[0][i] || 0,
          T: genomeData[1][i] || 0,
          G: genomeData[2][i] || 0,
          C: genomeData[3][i] || 0,
        },
      });
    }
  }

  return {
    genomeDataRaw,          // [4][N] for GenomeChart (categorical only)
    dataset,                // [{pos, nucleotide, mutationPoss}] for BarChart/BarChart2
    genome: genomeSequence || "",
    proteinMutationProbs: protein_mutation_probs || {},
    proteinRegionPossibilities: proteinRegionPossibilities || {},
    modelMetadata: model_metadata || {},
    extractorMetadata: extractor_metadata || {},
    // Self-describing v2.0 PredictionPayload (see backend PLUGIN_CONTRACT.md).
    // Drives the payload-router so binary / scalar / non-COVID models render
    // without the response having to fit the legacy ATGC shape.
    predictionPayload: predictionPayload || null,
    selectedProteinRegion: selectedProteinRegion || null,
    nodeId,
    elapsedDay,
    selectedModel,
  };
}

// ============================================
// INITIAL STATE
// ============================================

const initialState = {
  // Data - TWO legacy formats for different chart components
  genomeDataRaw: [],    // [4][N] format for GenomeChart (Recharts.js) — categorical only
  dataset: [],          // [{pos, nucleotide, mutationPoss}] for BarChart/BarChart2
  genome: "",
  proteinMutationProbs: {},
  proteinRegionPossibilities: {},

  // v2.0 self-describing PredictionPayload (drives the payload-router).
  // null until the first prediction lands or for legacy responses without it.
  predictionPayload: null,
  
  // Model info
  modelMetadata: {},
  extractorMetadata: {},
  availableModels: [],
  modelParameters: {},
  
  // UI State
  chartTitle: "Full Sequence",
  isWholeSequenceSelected: true,
  selectedProteinRegion: null,
  showDoughnut: true,
  isSelected: false,
  
  // Selection state
  nodeId: "",
  elapsedDay: 0,
  model: "",
  
  // Lists (for dropdowns)
  modelList: [],
  nodeList: [],
  
  // Loading states
  loading: false,
  modelsLoading: false,
  parametersLoading: false,
  error: null,
};

// ============================================
// SLICE
// ============================================

export const genomeSlice = createSlice({
  name: "genome",
  initialState,
  reducers: {
    setDataset: (state, action) => {
      const {
        dataset = [],
        genomeDataRaw = [],
        genome = "",
        proteinMutationProbs = {},
        proteinRegionPossibilities = {},
        isSelected = false,
        selectedProteinRegion = null,
      } = action.payload;

      state.dataset = dataset;
      state.genomeDataRaw = genomeDataRaw;
      state.genome = genome;
      state.proteinMutationProbs = proteinMutationProbs;
      state.proteinRegionPossibilities = proteinRegionPossibilities;
      state.isSelected = isSelected;
      state.selectedProteinRegion = selectedProteinRegion;
      state.showDoughnut = selectedProteinRegion === null;
    },

    showProteinRegion: (state, action) => {
      // Organism-aware: look up the region in the current organism's table
      // (set by the most recent prediction). For built-in COVID, this is the
      // same set as the legacy hardcoded import was.
      const region = state.proteinRegionPossibilities?.[action.payload];
      if (!region) {
        console.error(`Protein region "${action.payload}" not found.`);
        return;
      }

      state.chartTitle = action.payload;
      state.isWholeSequenceSelected = false;
      state.selectedProteinRegion = action.payload;
      state.showDoughnut = false;
    },

    resetChart: (state) => {
      state.isWholeSequenceSelected = true;
      state.chartTitle = "Full Sequence";
      state.selectedProteinRegion = null;
      state.showDoughnut = true;
    },

    updateProteinRegion: (state, action) => {
      state.selectedProteinRegion = action.payload || null;
      state.showDoughnut = action.payload === null;
      state.isWholeSequenceSelected = action.payload === null;
      state.chartTitle = action.payload === null ? "Full Sequence" : action.payload;
    },

    resetProteinRegion: (state) => {
      state.selectedProteinRegion = null;
      state.showDoughnut = true;
      state.isWholeSequenceSelected = true;
      state.chartTitle = "Full Sequence";
    },

    selectNode: (state, action) => {
      const [node, elapsedDay, model] = action.payload;
      state.nodeId = node;
      state.elapsedDay = elapsedDay;
      state.model = model;
    },

    loadNodesAndModels: (state, action) => {
      if (action.payload) {
        const [models, nodes] = action.payload;
        state.modelList = models || [];
        state.nodeList = nodes || [];
      }
    },

    setLoading: (state, action) => {
      state.loading = action.payload;
    },

    clearError: (state) => {
      state.error = null;
    },

    resetState: () => initialState,
  },

  extraReducers: (builder) => {
    // ===== fetchPrediction =====
    builder
      .addCase(fetchPrediction.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(fetchPrediction.fulfilled, (state, action) => {
        state.loading = false;
        state.genomeDataRaw = action.payload.genomeDataRaw;  // [4][N] for GenomeChart (categorical only)
        state.dataset = action.payload.dataset;               // [{mutationPoss}] for BarChart
        state.genome = action.payload.genome;
        state.proteinMutationProbs = action.payload.proteinMutationProbs;
        state.proteinRegionPossibilities = action.payload.proteinRegionPossibilities;
        state.modelMetadata = action.payload.modelMetadata;
        state.extractorMetadata = action.payload.extractorMetadata;
        state.predictionPayload = action.payload.predictionPayload;
        state.selectedProteinRegion = action.payload.selectedProteinRegion;
        state.isSelected = true;
        state.showDoughnut = action.payload.selectedProteinRegion === null;
      })
      .addCase(fetchPrediction.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload || "Prediction failed";
      });

    // ===== fetchAvailableModels =====
    builder
      .addCase(fetchAvailableModels.pending, (state) => {
        state.modelsLoading = true;
      })
      .addCase(fetchAvailableModels.fulfilled, (state, action) => {
        state.modelsLoading = false;
        state.availableModels = action.payload;
      })
      .addCase(fetchAvailableModels.rejected, (state, action) => {
        state.modelsLoading = false;
        console.warn("Models API unavailable:", action.payload);
      });

    // ===== fetchModelParameters =====
    builder
      .addCase(fetchModelParameters.pending, (state) => {
        state.parametersLoading = true;
      })
      .addCase(fetchModelParameters.fulfilled, (state, action) => {
        state.parametersLoading = false;
        state.modelParameters = action.payload.parameters || {};
      })
      .addCase(fetchModelParameters.rejected, (state) => {
        state.parametersLoading = false;
        state.modelParameters = {};
      });
  },
});

// Export actions
export const {
  setDataset,
  showProteinRegion,
  resetChart,
  updateProteinRegion,
  resetProteinRegion,
  selectNode,
  loadNodesAndModels,
  setLoading,
  clearError,
  resetState,
} = genomeSlice.actions;

export default genomeSlice.reducer;