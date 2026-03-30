import React, { useState, useEffect, useMemo, useCallback } from "react";
import { nodeIds as nodes } from "../data/nodeIds";
import { modelList as staticModels } from "../data/modelList";
import { useDispatch, useSelector } from "react-redux";
import { Button, Input } from "@material-tailwind/react";
import {
  MdOutlineCreate,
  MdScience,
  MdTimeline,
  MdBiotech,
  MdWarning,
  MdCheckCircle,
  MdCloudUpload,
  MdClose,
  MdAdd,
  MdDelete,
  MdSettings,
  MdInfo,
} from "react-icons/md";
import Select from "react-select";
import LoadingSpinner from "./Spinner";
import { proteinRegions } from "../data/proteinRegions";
import DropDown from "./DropDown";
import {
  updateProteinRegion,
  resetProteinRegion,
  fetchPrediction,
  fetchAvailableModels,
} from "../features/genome/genomeSlice";
import logo from "../CovMutexLogo-removebg-preview.png";

// Charts
import BarChart from "./BarChart";
import BarChart2 from "./BarChart2";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

// ============================================
// STYLES
// ============================================
const customStyles = {
  placeholder: (provided) => ({
    ...provided,
    color: "#9CA3AF",
    fontSize: "14px",
  }),
  control: (provided, state) => ({
    ...provided,
    minHeight: "42px",
    borderRadius: "0.5rem",
    borderColor: state.isFocused ? "#3B82F6" : "#E5E7EB",
    boxShadow: state.isFocused ? "0 0 0 1px #3B82F6" : "none",
    "&:hover": { borderColor: "#3B82F6" },
  }),
  // Uploaded modelleri farklı renkte göster
  option: (provided, state) => ({
    ...provided,
    backgroundColor: state.isSelected
      ? "#3B82F6"
      : state.isFocused
      ? "#EFF6FF"
      : "white",
    color: state.isSelected ? "white" : state.data?.type === "uploaded" ? "#059669" : "#1F2937",
    fontWeight: state.data?.type === "uploaded" ? "600" : "400",
  }),
};

// ============================================
// SUB-COMPONENTS
// ============================================

const SummaryCard = ({ title, value, subtext, icon: Icon, color }) => (
  <div className="bg-white p-4 rounded-xl shadow-sm border border-gray-100 flex items-center gap-4 transition-transform hover:scale-105">
    <div className={`p-3 rounded-full ${color} text-white`}>
      <Icon size={24} />
    </div>
    <div>
      <p className="text-xs font-semibold text-gray-500 uppercase tracking-wide">
        {title}
      </p>
      <h4 className="text-xl font-bold text-gray-800">{value}</h4>
      {subtext && <p className="text-xs text-gray-400">{subtext}</p>}
    </div>
  </div>
);

const EmptyState = () => (
  <div className="h-full flex flex-col items-center justify-center text-center p-10 opacity-60">
    <div className="bg-blue-50 p-6 rounded-full mb-4">
      <MdBiotech size={64} className="text-blue-300" />
    </div>
    <h3 className="text-xl font-bold text-gray-700 mb-2">Ready to Analyze</h3>
    <p className="text-gray-500 max-w-sm">
      Select a prediction model, enter parameters, and run prediction to see results.
    </p>
  </div>
);

// ============================================
// UPLOAD MODAL
// ============================================
const UploadModal = ({ isOpen, onClose, onSuccess }) => {
  const [name, setName] = useState("");
  const [modelFile, setModelFile] = useState(null);
  const [extractorFile, setExtractorFile] = useState(null);
  const [helperFiles, setHelperFiles] = useState([]);
  const [paramsList, setParamsList] = useState([
  ]);
  const [uploading, setUploading] = useState(false);
  const [error, setError] = useState("");
  const [successMessage, setSuccessMessage] = useState("");

  if (!isOpen) return null;

  const addParam = () =>
    setParamsList([...paramsList, { key: "", value: "", required: false }]);

  const removeParam = (index) => {
    const list = [...paramsList];
    list.splice(index, 1);
    setParamsList(list);
  };

  const updateParam = (index, field, val) => {
    const list = [...paramsList];
    list[index][field] = val;
    setParamsList(list);
  };

  const handleHelperFilesChange = (e) => {
    if (e.target.files) {
      setHelperFiles(Array.from(e.target.files));
    }
  };

  const resetForm = () => {
    setName("");
    setModelFile(null);
    setExtractorFile(null);
    setHelperFiles([]);
    setParamsList([{ key: "batch_size", value: "32", required: false }]);
    setError("");
    setSuccessMessage("");
  };

  const handleUpload = async (e) => {
    e.preventDefault();
    setError("");
    setSuccessMessage("");

    if (!modelFile || !name.trim()) {
      setError("Please fill required fields (Name, Model File).");
      return;
    }

    // Validate required parameters have keys
    const invalidParams = paramsList.filter(
      (p) => p.required && !p.key.trim()
    );
    if (invalidParams.length > 0) {
      setError("Required parameters must have a key name.");
      return;
    }

    setUploading(true);

    try {
      const formData = new FormData();

      // Use model name as folder name (no timestamp)
      const folderName = name.trim().replace(/\s+/g, "_");

      // Required fields
      formData.append("uploadFolderName", folderName);
      formData.append("modelFile", modelFile);
      formData.append("nodeId", "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07");
      formData.append("elapsedDay", "60");

      // Optional extractor
      if (extractorFile) {
        formData.append("extractorFile", extractorFile);
      }

      // Helper files
      helperFiles.forEach((file, index) => {
        formData.append(`helperFile_${index}`, file);
        formData.append(`helperFileName_${index}`, file.name);
      });

      // Custom parameters - save with {value, required, default} structure
      const customParamsObj = {};
      paramsList.forEach((item) => {
        if (item.key.trim()) {
          const rawValue = item.value.trim();
          const numericValue = !isNaN(rawValue) && rawValue !== "" ? parseFloat(rawValue) : rawValue;
          customParamsObj[item.key.trim()] = {
            value: numericValue,
            required: item.required || false,
            default: rawValue,
          };
        }
      });
      formData.append("customParameters", JSON.stringify(customParamsObj));

      const response = await fetch(`${API_URL}/api/predict/`, {
        method: "POST",
        body: formData,
      });

      if (response.ok) {
        setSuccessMessage(`Model "${name}" uploaded successfully!`);
        
        // Trigger model list refresh in parent, then close after a short delay
        if (onSuccess) {
          await onSuccess(folderName);  // Pass the folder name so parent can auto-select it
        }
        
        // Close after showing success
        setTimeout(() => {
          resetForm();
          onClose();
        }, 1500);
      } else {
        const err = await response.json().catch(() => ({}));
        setError(err.error || "Upload failed. Please try again.");
      }
    } catch (err) {
      console.error("Upload error:", err);
      setError("Server error during upload.");
    } finally {
      setUploading(false);
    }
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black bg-opacity-50 backdrop-blur-sm p-4 overflow-y-auto">
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-lg p-6 relative my-8 max-h-[90vh] overflow-y-auto">
        {/* Close button */}
        <button
          onClick={onClose}
          className="absolute top-4 right-4 text-gray-400 hover:text-gray-600"
        >
          <MdClose size={24} />
        </button>

        {/* Header */}
        <h2 className="text-xl font-bold text-blue-900 mb-1 flex items-center gap-2">
          <MdCloudUpload /> Upload New Model
        </h2>
        <p className="text-sm text-gray-500 mb-6">
          Add a new prediction model with custom parameters.
        </p>

        {/* Success message */}
        {successMessage && (
          <div className="mb-4 p-3 bg-green-50 border border-green-200 rounded-lg text-green-700 text-sm flex items-center gap-2">
            <MdCheckCircle size={20} />
            {successMessage}
          </div>
        )}

        {/* Error message */}
        {error && (
          <div className="mb-4 p-3 bg-red-50 border border-red-200 rounded-lg text-red-700 text-sm">
            {error}
          </div>
        )}

        <form onSubmit={handleUpload} className="space-y-4">
          {/* Model Name */}
          <div>
            <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
              Model Name *
            </label>
            <Input
              type="text"
              value={name}
              onChange={(e) => setName(e.target.value)}
              placeholder="e.g. LSTM Variant Predictor"
              required
            />
          </div>

          {/* File Uploads */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {/* Model File */}
            <div>
              <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
                Model File (.keras/.h5) *
              </label>
              <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-blue-400 transition-colors">
                <input
                  type="file"
                  accept=".keras,.h5,.pt,.pth"
                  onChange={(e) => setModelFile(e.target.files[0])}
                  className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
                  required
                />
                <div className="text-center">
                  <MdCloudUpload className="mx-auto text-gray-400 text-2xl mb-1" />
                  <p className="text-xs text-gray-500">
                    {modelFile ? (
                      <span className="text-green-600 font-medium">{modelFile.name}</span>
                    ) : (
                      "Click to upload"
                    )}
                  </p>
                </div>
              </div>
            </div>

            {/* Extractor File */}
            <div>
              <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
                Feature Extractor (.py)
              </label>
              <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-green-400 transition-colors">
                <input
                  type="file"
                  accept=".py"
                  onChange={(e) => setExtractorFile(e.target.files[0])}
                  className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
                />
                <div className="text-center">
                  <MdSettings className="mx-auto text-gray-400 text-2xl mb-1" />
                  <p className="text-xs text-gray-500">
                    {extractorFile ? (
                      <span className="text-blue-600 font-medium">{extractorFile.name}</span>
                    ) : (
                      "Optional"
                    )}
                  </p>
                </div>
              </div>
            </div>
          </div>

          {/* Helper Files */}
          <div>
            <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
              Helper Files (Optional)
            </label>
            <input
              type="file"
              multiple
              onChange={handleHelperFilesChange}
              className="block w-full text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:text-sm file:font-semibold file:bg-gray-50 file:text-gray-700 hover:file:bg-gray-100"
            />
            {helperFiles.length > 0 && (
              <p className="text-xs text-gray-500 mt-1">
                {helperFiles.length} file(s) selected
              </p>
            )}
          </div>

          {/* Custom Parameters */}
          <div>
            <div className="flex justify-between items-center mb-2">
              <label className="text-xs font-bold text-gray-700 uppercase">
                Model Parameters
              </label>
              <button
                type="button"
                onClick={addParam}
                className="text-xs text-blue-600 hover:text-blue-800 flex items-center gap-1"
              >
                <MdAdd /> Add Parameter
              </button>
            </div>
            <div className="space-y-2">
              {paramsList.map((item, index) => (
                <div key={index} className="flex gap-2 items-center bg-gray-50 p-2 rounded-lg">
                  {/* Parameter Key */}
                  <input
                    type="text"
                    placeholder="Key"
                    value={item.key}
                    onChange={(e) => updateParam(index, "key", e.target.value)}
                    className="flex-1 border rounded px-2 py-1 text-sm"
                  />
                  {/* Default Value */}
                  <input
                    type="text"
                    placeholder="Default Value"
                    value={item.value}
                    onChange={(e) => updateParam(index, "value", e.target.value)}
                    className="flex-1 border rounded px-2 py-1 text-sm"
                  />
                  {/* Required Checkbox */}
                  <label className="flex items-center gap-1 text-xs text-gray-600 whitespace-nowrap">
                    <input
                      type="checkbox"
                      checked={item.required}
                      onChange={(e) => updateParam(index, "required", e.target.checked)}
                      className="rounded"
                    />
                    Required
                  </label>
                  {paramsList.length > 1 && (
                    <button
                      type="button"
                      onClick={() => removeParam(index)}
                      className="text-red-500 hover:text-red-700 p-1"
                    >
                      <MdDelete />
                    </button>
                  )}
                </div>
              ))}
            </div>
            <p className="text-xs text-gray-400 mt-1">
              Parameters marked as "Required" must be filled before running predictions.
            </p>
          </div>

          {/* Submit Button */}
          <Button
            type="submit"
            color="blue"
            className="w-full"
            disabled={uploading}
          >
            {uploading ? "Uploading & Testing..." : "Upload Model"}
          </Button>
        </form>
      </div>
    </div>
  );
};

// ============================================
// MAIN NAVBAR COMPONENT
// ============================================
function Navbar({ onNodeSelect, onSubmit, isLoading }) {
  const dispatch = useDispatch();

  // Redux state
  const {
    dataset: reduxDataset,
    genome: reduxGenome,
    selectedProteinRegion,
    availableModels,
    loading: reduxLoading,
  } = useSelector((state) => state.genome);

  // Local state
  const [_nodeId, setNodeId] = useState(null);
  const [_elapsedDay, setElapsedDay] = useState("60");
  const [selectedModel, setSelectedModel] = useState(null);
  const [combinedModelList, setCombinedModelList] = useState([]);
  const [isUploadModalOpen, setIsUploadModalOpen] = useState(false);
  const [modelParameters, setModelParameters] = useState([]);
  const [parametersLoading, setParametersLoading] = useState(false);
  const [localGenomeData, setLocalGenomeData] = useState(null);

  const loading = isLoading || reduxLoading;

  // ============================================
  // FETCH MODELS FROM API
  // ============================================
  const fetchModels = useCallback(async (autoSelectFolder = null) => {
    // Format static models
    const formattedStatic = staticModels.map((m) => ({
      label: m.name,
      value: m.path,
      type: "static",
    }));

    try {
      const response = await fetch(`${API_URL}/api/models/`);

      if (response.ok) {
        const data = await response.json();
        
        // Use the new detailed endpoint
        const uploadedModelsList = data.uploaded_models || [];
        
        const formattedAPI = uploadedModelsList.map((m) => ({
          label: `${m.folder_name} (Uploaded)`,
          value: `uploaded:${m.folder_name}`,
          type: "uploaded",
          hasParams: m.has_parameters,
          hasExtractor: m.has_extractor,
          modelFile: m.model_file,
        }));

        const allModels = [...formattedStatic, ...formattedAPI];
        setCombinedModelList(allModels);

        // If a specific folder was just uploaded, auto-select it
        if (autoSelectFolder) {
          const uploadedValue = `uploaded:${autoSelectFolder}`;
          const found = allModels.find((m) => m.value === uploadedValue);
          if (found) {
            setSelectedModel(found.value);
            return;
          }
        }

        // Default selection if nothing selected yet
        if (!selectedModel && allModels.length > 0) {
          setSelectedModel(allModels[0].value);
        }
      } else {
        // Fallback: backward compat with old endpoint
        setCombinedModelList(formattedStatic);
        if (!selectedModel && formattedStatic.length > 0) {
          setSelectedModel(formattedStatic[0].value);
        }
      }
    } catch (error) {
      console.warn("API offline, using static models:", error.message);
      setCombinedModelList(formattedStatic);
      if (!selectedModel && formattedStatic.length > 0) {
        setSelectedModel(formattedStatic[0].value);
      }
    }
  }, [selectedModel]);

  useEffect(() => {
    fetchModels();
    // eslint-disable-next-line
  }, []);

  // ============================================
  // FETCH MODEL PARAMETERS WHEN MODEL CHANGES
  // ============================================
  useEffect(() => {
    const fetchModelParams = async () => {
      if (!selectedModel) {
        setModelParameters([]);
        return;
      }

      // Only fetch parameters for uploaded models
      if (selectedModel.startsWith("uploaded:")) {
        setParametersLoading(true);
        try {
          const modelName = selectedModel.replace("uploaded:", "");
          const response = await fetch(
            `${API_URL}/api/model-parameters/?model_name=${encodeURIComponent(modelName)}`
          );

          if (response.ok) {
            const data = await response.json();
            const params = data.parameters || {};
            
            // Convert to array format for UI rendering
            // Backend now always returns {value, required, default} format
            const paramsArray = Object.entries(params).map(([key, val]) => {
              if (typeof val === "object" && val !== null) {
                return {
                  key,
                  value: String(val.value ?? val.default ?? ""),
                  required: val.required || false,
                  default: val.default !== undefined ? String(val.default) : "",
                };
              } else {
                // Fallback for old format
                return {
                  key,
                  value: String(val),
                  required: false,
                  default: String(val),
                };
              }
            });

            setModelParameters(paramsArray);
          } else {
            setModelParameters([]);
          }
        } catch (error) {
          console.warn("Could not fetch model parameters:", error);
          setModelParameters([]);
        } finally {
          setParametersLoading(false);
        }
      } else {
        // Static/server models don't have custom parameters
        setModelParameters([]);
      }
    };

    fetchModelParams();
  }, [selectedModel]);

  // ============================================
  // HANDLERS
  // ============================================
  const handleProteinRegionChange = (opt) => {
    const val = opt ? opt.value : null;
    if (val) {
      dispatch(updateProteinRegion(val));
    } else {
      dispatch(resetProteinRegion());
    }
  };

  const updateModelParam = (index, newValue) => {
    const updated = [...modelParameters];
    updated[index].value = newValue;
    setModelParameters(updated);
  };

  // Reset a single parameter to its default value
  const resetModelParam = (index) => {
    const updated = [...modelParameters];
    updated[index].value = updated[index].default || "";
    setModelParameters(updated);
  };

  const handleSubmit = async (e) => {
    e.preventDefault();

    if (!_nodeId || !selectedModel) {
      alert("Please select a Variant ID and Model.");
      return;
    }

    if (!_elapsedDay || _elapsedDay === "0") {
      alert("Please enter Elapsed Days (must be greater than 0).");
      return;
    }

    // Validate required parameters
    const missingRequired = modelParameters.filter(
      (p) => p.required && (!p.value || p.value.trim() === "")
    );
    if (missingRequired.length > 0) {
      alert(`Please fill required parameters: ${missingRequired.map((p) => p.key).join(", ")}`);
      return;
    }

    // Build custom parameters object - send actual values
    const customParamsObj = {};
    modelParameters.forEach((item) => {
      if (item.key.trim()) {
        const isNum = !isNaN(item.value) && item.value.trim() !== "";
        customParamsObj[item.key.trim()] = isNum
          ? parseFloat(item.value)
          : item.value;
      }
    });

    const params = {
      nodeId: _nodeId,
      elapsedDay: Number(_elapsedDay) || 60,
      selectedModel: selectedModel,
      selectedProteinRegion: selectedProteinRegion || null,
      isNewUpload: false,
      customParameters: customParamsObj,
    };

    try {
      const result = await dispatch(fetchPrediction(params)).unwrap();
      setLocalGenomeData(result);

      if (onSubmit) {
        await onSubmit(_nodeId, _elapsedDay, selectedModel, selectedProteinRegion);
      }
    } catch (error) {
      console.error("Prediction failed:", error);
      alert("Prediction failed: " + (error.message || error));
    }
  };

  // Upload success handler - refresh list and auto-select the new model
  const handleUploadSuccess = async (folderName) => {
    await fetchModels(folderName);
  };

  // ============================================
  // COMPUTED VALUES
  // ============================================

  const genomeData = reduxDataset?.length > 0 ? reduxDataset : localGenomeData?.dataset;
  const genomeSequence = reduxGenome || localGenomeData?.genome;

  const stats = useMemo(() => {
    if (!genomeData || genomeData.length === 0) return null;

    let totalConfidence = 0;
    let maxRisk = 0;

    genomeData.forEach((item) => {
      if (item.mutationPoss) {
        const values = Object.values(item.mutationPoss);
        const max = Math.max(...values);
        totalConfidence += max;
        const risk = 1 - max;
        if (risk > maxRisk) maxRisk = risk;
      }
    });

    return {
      length: genomeData.length.toLocaleString(),
      confidence: (totalConfidence / genomeData.length).toFixed(3),
      maxRisk: maxRisk.toFixed(3),
    };
  }, [genomeData]);

  // Check if currently selected model is an uploaded model
  const isUploadedModel = selectedModel && selectedModel.startsWith("uploaded:");

  // ============================================
  // RENDER
  // ============================================
  if (loading && !genomeData) {
    return <LoadingSpinner />;
  }

  return (
    <div className="min-h-screen bg-gray-50 flex flex-col p-6 gap-6">
      {/* Header */}
      <header className="flex items-center justify-between bg-white px-8 py-4 rounded-2xl shadow-sm border border-gray-100">
        <div className="flex items-center gap-4">
          <img src={logo} alt="Logo" className="h-12 w-auto object-contain" />
          <div className="hidden md:block w-px h-10 bg-gray-200"></div>
          <h1 className="hidden md:block text-xl font-bold text-blue-900 tracking-tight">
            CovMutEx - Mutation Explorer
          </h1>
        </div>
      </header>

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6 items-start">
        {/* ============================================ */}
        {/* FORM PANEL (Left Side) */}
        {/* ============================================ */}
        <div className="lg:col-span-4 bg-white p-6 rounded-2xl shadow-lg border border-gray-100 sticky top-6 z-10">
          <div className="mb-6 pb-4 border-b border-gray-100">
            <h2 className="text-lg font-bold text-gray-800 flex items-center gap-2">
              <MdScience className="text-blue-600" /> Analysis Parameters
            </h2>
            <p className="text-xs text-gray-400 mt-1">
              Configure settings or upload a new model.
            </p>
          </div>

          <form className="space-y-5" onSubmit={handleSubmit}>
            {/* Model Selection */}
            <div>
              <div className="flex justify-between items-center mb-1">
                <label className="text-xs font-bold text-gray-500 uppercase block ml-1">
                  Prediction Model *
                </label>
                <button
                  type="button"
                  onClick={() => setIsUploadModalOpen(true)}
                  className="text-xs text-blue-600 font-bold hover:text-blue-800 flex items-center gap-1 bg-blue-50 px-2 py-1 rounded transition-colors"
                >
                  <MdCloudUpload /> Upload New
                </button>
              </div>
              <Select
                options={combinedModelList}
                onChange={(opt) => setSelectedModel(opt ? opt.value : null)}
                value={combinedModelList.find((o) => o.value === selectedModel) || null}
                styles={customStyles}
                placeholder="Choose a model..."
                formatOptionLabel={(option) => (
                  <div className="flex items-center gap-2">
                    <span>{option.label}</span>
                    {option.type === "uploaded" && option.hasParams && (
                      <MdSettings className="text-gray-400" size={14} title="Has custom parameters" />
                    )}
                  </div>
                )}
              />
              {isUploadedModel && (
                <p className="text-xs text-green-600 mt-1 ml-1 flex items-center gap-1">
                  <MdCheckCircle size={12} /> Uploaded model selected
                </p>
              )}
            </div>

            {/* Variant ID */}
            <div>
              <label className="text-xs font-bold text-gray-500 uppercase mb-1 block ml-1">
                Variant ID *
              </label>
              <DropDown items={nodes} setNodeId={setNodeId} />
            </div>

            {/* Elapsed Days */}
            <div>
              <label className="text-xs font-bold text-gray-500 uppercase mb-1 block ml-1">
                Elapsed Days *
              </label>
              <Input
                type="number"
                value={_elapsedDay}
                onChange={(e) => setElapsedDay(e.target.value)}
                min={1}
                placeholder="e.g., 60"
                required
                className="!border !border-gray-300 focus:!border-blue-500"
              />
              <p className="text-xs text-gray-400 mt-1 ml-1">
                Days since variant emergence (affects mutation probability)
              </p>
            </div>

            {/* Region */}
            <div>
              <label className="text-xs font-bold text-gray-500 uppercase mb-1 block ml-1">
                Protein Region
              </label>
              <Select
                options={[
                  { label: "Whole Genome", value: "" },
                  ...Object.keys(proteinRegions).map((pr) => ({
                    label: pr,
                    value: pr,
                  })),
                ]}
                onChange={handleProteinRegionChange}
                styles={customStyles}
                placeholder="Optional"
                isClearable
              />
            </div>

            {/* Model Parameters Loading Spinner */}
            {parametersLoading && (
              <div className="text-center py-4">
                <div className="animate-spin rounded-full h-6 w-6 border-2 border-blue-500 border-t-transparent mx-auto"></div>
                <p className="text-xs text-gray-500 mt-2">Loading parameters...</p>
              </div>
            )}

            {/* Model Parameters Section (only for uploaded models with parameters) */}
            {!parametersLoading && modelParameters.length > 0 && (
              <div className="mt-4 p-4 bg-gradient-to-br from-blue-50 to-indigo-50 rounded-xl border border-blue-100">
                <div className="flex items-center gap-2 mb-3">
                  <MdSettings className="text-blue-500" />
                  <label className="text-xs font-bold text-blue-800 uppercase">
                    Model Parameters
                  </label>
                  <MdInfo className="text-blue-400" title="Parameters specific to this uploaded model" />
                </div>
                <div className="space-y-3">
                  {modelParameters.map((item, index) => (
                    <div
                      key={index}
                      className={`bg-white p-3 rounded-lg border ${
                        item.required
                          ? item.value && item.value.trim() !== ""
                            ? "border-green-200"
                            : "border-red-300 bg-red-50/30"
                          : "border-gray-200"
                      } shadow-sm transition-colors`}
                    >
                      <div className="flex items-center justify-between mb-1">
                        <span className="text-xs font-bold text-gray-700 uppercase flex items-center gap-1">
                          {item.key}
                          {item.required && (
                            <span className="text-red-500 text-[10px] font-bold bg-red-50 px-1 rounded">
                              REQUIRED
                            </span>
                          )}
                          {!item.required && (
                            <span className="text-gray-400 text-[10px] font-normal bg-gray-50 px-1 rounded">
                              optional
                            </span>
                          )}
                        </span>
                        <div className="flex items-center gap-2">
                          {item.default !== undefined && item.default !== "" && (
                            <span className="text-xs text-gray-400">
                              Default: {item.default}
                            </span>
                          )}
                          {/* Reset to default button */}
                          {item.default && item.value !== item.default && (
                            <button
                              type="button"
                              onClick={() => resetModelParam(index)}
                              className="text-xs text-blue-500 hover:text-blue-700 underline"
                              title="Reset to default"
                            >
                              Reset
                            </button>
                          )}
                        </div>
                      </div>
                      <input
                        type="text"
                        className={`w-full text-sm font-medium text-gray-800 border rounded-md px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-blue-500 ${
                          item.required && (!item.value || item.value.trim() === "")
                            ? "border-red-300 bg-red-50"
                            : "border-gray-200"
                        }`}
                        value={item.value}
                        onChange={(e) => updateModelParam(index, e.target.value)}
                        placeholder={
                          item.required
                            ? `Required (default: ${item.default || "none"})`
                            : `Optional (default: ${item.default || "none"})`
                        }
                        required={item.required}
                      />
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Uploaded model with no parameters info */}
            {!parametersLoading && isUploadedModel && modelParameters.length === 0 && (
              <div className="p-3 bg-gray-50 rounded-lg border border-gray-200 text-center">
                <p className="text-xs text-gray-500">
                  This uploaded model has no custom parameters.
                </p>
              </div>
            )}

            {/* Submit Button */}
            <Button
              size="lg"
              color="blue"
              type="submit"
              className="w-full flex justify-center items-center gap-2 shadow-blue-500/20 hover:shadow-blue-500/40 mt-4"
              disabled={loading}
            >
              {loading ? "Processing..." : "Run Prediction"}
              <MdOutlineCreate className="text-lg" />
            </Button>
          </form>
        </div>

        {/* ============================================ */}
        {/* RESULTS PANEL (Right Side) */}
        {/* ============================================ */}
        <div className="lg:col-span-8 space-y-6">
          {/* Empty State */}
          {!genomeData && (
            <div className="bg-white rounded-2xl shadow-sm border border-gray-100 h-[500px] flex items-center justify-center">
              <EmptyState />
            </div>
          )}

          {/* Results */}
          {genomeData && genomeData.length > 0 && stats && (
            <div className="space-y-6">
              {/* Summary Cards */}
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                <SummaryCard
                  title="Genome Length"
                  value={stats.length}
                  subtext="Base pairs"
                  icon={MdTimeline}
                  color="bg-blue-500"
                />
                <SummaryCard
                  title="Avg Confidence"
                  value={stats.confidence}
                  subtext="Certainty"
                  icon={MdCheckCircle}
                  color="bg-green-500"
                />
                <SummaryCard
                  title="Max Risk"
                  value={stats.maxRisk}
                  subtext="Mutation prob"
                  icon={MdWarning}
                  color="bg-orange-500"
                />
              </div>

              {/* Charts */}
              <div className="bg-white p-4 rounded-2xl shadow-sm border border-gray-100">
                <h3 className="text-sm font-bold text-gray-700 mb-2">
                  Mutation Probability Distribution
                </h3>
                <BarChart data={genomeData} />
              </div>

              <div className="bg-white p-4 rounded-2xl shadow-sm border border-gray-100">
                <h3 className="text-sm font-bold text-gray-700 mb-2">
                  Detailed Sequence View
                </h3>
                <BarChart2 data={genomeData} seq={genomeSequence || ""} />
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Upload Modal */}
      <UploadModal
        isOpen={isUploadModalOpen}
        onClose={() => setIsUploadModalOpen(false)}
        onSuccess={handleUploadSuccess}
      />
    </div>
  );
}

export default Navbar;