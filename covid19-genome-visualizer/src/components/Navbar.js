import React, { useState, useEffect } from "react";
import { nodeIds as nodes } from "../data/nodeIds";
import { modelList as staticModels } from "../data/modelList";
import { useDispatch, useSelector } from "react-redux";
import { Button, Input } from "@material-tailwind/react";
import {
  MdOutlineCreate,
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
};

// ============================================
// PARAMETER PRESETS (shared between UploadModal and main form)
// ============================================
const UPLOAD_PARAM_PRESETS = [
  { key: "batch_size", label: "Batch Size", defaultValue: "32", options: ["8", "16", "32", "64", "128", "256"], hint: "Samples per batch" },
  { key: "learning_rate", label: "Learning Rate", defaultValue: "0.001", placeholder: "e.g., 0.001", hint: "Optimizer step size" },
  { key: "epochs", label: "Epochs", defaultValue: "100", options: ["10", "25", "50", "100", "200", "500"], hint: "Training iterations" },
  { key: "optimizer", label: "Optimizer", defaultValue: "adam", options: ["adam", "sgd", "rmsprop", "adamw", "adagrad"], hint: "Optimization algorithm" },
  { key: "dropout_rate", label: "Dropout Rate", defaultValue: "0.2", placeholder: "0.0 - 1.0", hint: "Neuron drop fraction" },
  { key: "sequence_length", label: "Sequence Length", defaultValue: "29903", placeholder: "e.g., 29903", hint: "Genome input length" },
  { key: "hidden_units", label: "Hidden Units", defaultValue: "128", options: ["32", "64", "128", "256", "512"], hint: "Neurons in hidden layers" },
  { key: "activation", label: "Activation", defaultValue: "relu", options: ["relu", "sigmoid", "tanh", "softmax", "leaky_relu", "gelu"], hint: "Activation function" },
];

// ============================================
// UPLOAD MODAL
// ============================================
const UploadModal = ({ isOpen, onClose, onSuccess }) => {
  const [name, setName] = useState("");
  const [modelFile, setModelFile] = useState(null);
  const [extractorFile, setExtractorFile] = useState(null);
  const [helperFiles, setHelperFiles] = useState([]);
  const [paramsList, setParamsList] = useState([
    { key: "batch_size", value: "32", required: false },
  ]);
  const [uploading, setUploading] = useState(false);
  const [error, setError] = useState("");

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
    if (field === "required") {
      list[index][field] = val;
    } else {
      list[index][field] = val;
    }
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
  };

  const handleUpload = async (e) => {
    e.preventDefault();
    setError("");

    if (!modelFile || !name.trim()) {
      setError("Please fill required fields (Name, Model File).");
      return;
    }

    setUploading(true);

    try {
      const formData = new FormData();

      // Use model name as folder name (clean, no timestamp)
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

      // Custom parameters - include required flag
      const customParamsObj = {};
      paramsList.forEach((item) => {
        if (item.key.trim()) {
          customParamsObj[item.key.trim()] = {
            value: isNaN(item.value) ? item.value : parseFloat(item.value),
            required: item.required || false,
            default: item.value,
          };
        }
      });
      formData.append("customParameters", JSON.stringify(customParamsObj));

      const response = await fetch(`${API_URL}/api/predict/`, {
        method: "POST",
        body: formData,
      });

      if (response.ok) {
        resetForm();
        onSuccess();  // refresh model list
        onClose();
      } else {
        const err = await response.json().catch(() => ({}));
        // Model might have been saved even if prediction failed
        // Refresh the model list anyway
        onSuccess();
        setError(err.error || "Upload completed but prediction may have failed. Your model has been saved.");
      }
    } catch (err) {
      console.error("Upload error:", err);
      // Even on network error, try refreshing the list
      onSuccess();
      setError("Server error during upload. If the model was saved, it will appear in the list.");
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

          {/* Custom Parameters - Dropdown Style */}
          <div className="bg-gradient-to-br from-blue-50 to-indigo-50 rounded-xl border border-blue-100 overflow-hidden">
            <div className="bg-gradient-to-r from-blue-600 to-indigo-600 px-4 py-3 flex justify-between items-center">
              <div className="flex items-center gap-2 text-white">
                <MdSettings className="text-lg" />
                <span className="text-sm font-bold">Model Parameters</span>
              </div>
              <button
                type="button"
                onClick={addParam}
                className="flex items-center gap-1 text-xs text-white bg-white/20 hover:bg-white/30 px-3 py-1.5 rounded-lg transition-colors border border-white/30"
              >
                <MdAdd /> Add
              </button>
            </div>
            <div className="p-3 space-y-2">
              {paramsList.map((item, index) => {
                const isKnownParam = UPLOAD_PARAM_PRESETS.some(p => p.key === item.key);
                const presetInfo = UPLOAD_PARAM_PRESETS.find(p => p.key === item.key);
                return (
                  <div key={index} className="bg-white rounded-xl border border-gray-200 overflow-hidden shadow-sm">
                    <div className="p-3 space-y-2">
                      {/* Parameter Name */}
                      <div>
                        <label className="text-[10px] font-bold text-gray-400 uppercase tracking-wider mb-1 block">Parameter</label>
                        <select
                          value={isKnownParam ? item.key : "__custom__"}
                          onChange={(e) => {
                            const val = e.target.value;
                            if (val === "__custom__") {
                              updateParam(index, "key", "");
                              updateParam(index, "value", "");
                            } else {
                              const p = UPLOAD_PARAM_PRESETS.find(pr => pr.key === val);
                              updateParam(index, "key", val);
                              if (p?.defaultValue) updateParam(index, "value", p.defaultValue);
                            }
                          }}
                          className="w-full border-2 border-gray-200 rounded-lg px-3 py-2 text-sm bg-white hover:border-blue-300 focus:border-blue-500 focus:ring-2 focus:ring-blue-200 transition-all cursor-pointer font-medium"
                        >
                          <option value="" disabled>-- Select parameter --</option>
                          {UPLOAD_PARAM_PRESETS.map(p => (
                            <option key={p.key} value={p.key}>{p.label}</option>
                          ))}
                          <option value="__custom__">Custom...</option>
                        </select>
                        {!isKnownParam && (
                          <input
                            type="text"
                            placeholder="Custom parameter name"
                            value={item.key}
                            onChange={(e) => updateParam(index, "key", e.target.value)}
                            className="w-full border-2 border-dashed border-gray-300 rounded-lg px-3 py-2 text-sm mt-1.5 focus:border-blue-400"
                          />
                        )}
                      </div>
                      {/* Value */}
                      <div>
                        <label className="text-[10px] font-bold text-gray-400 uppercase tracking-wider mb-1 block">
                          Value {presetInfo?.hint && <span className="normal-case font-normal text-gray-300">({presetInfo.hint})</span>}
                        </label>
                        {presetInfo?.options ? (
                          <div className="flex flex-wrap gap-1.5">
                            {presetInfo.options.map(opt => (
                              <button
                                key={opt}
                                type="button"
                                onClick={() => updateParam(index, "value", opt)}
                                className={`px-3 py-1.5 rounded-lg text-xs font-medium border transition-all ${
                                  item.value === opt
                                    ? "border-blue-500 bg-blue-50 text-blue-700"
                                    : "border-gray-200 bg-gray-50 text-gray-600 hover:border-blue-300"
                                }`}
                              >
                                {opt}
                              </button>
                            ))}
                          </div>
                        ) : (
                          <input
                            type="text"
                            placeholder={presetInfo?.placeholder || "Value"}
                            value={item.value}
                            onChange={(e) => updateParam(index, "value", e.target.value)}
                            className="w-full border-2 border-gray-200 rounded-lg px-3 py-2 text-sm focus:border-blue-500 focus:ring-2 focus:ring-blue-200"
                          />
                        )}
                      </div>
                      {/* Required + Delete row */}
                      <div className="flex items-center justify-between pt-1">
                        <label className="flex items-center gap-1.5 text-xs text-gray-500 cursor-pointer">
                          <input
                            type="checkbox"
                            checked={item.required}
                            onChange={(e) => updateParam(index, "required", e.target.checked)}
                            className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                          />
                          Required
                        </label>
                        {paramsList.length > 1 && (
                          <button
                            type="button"
                            onClick={() => removeParam(index)}
                            className="text-xs text-red-400 hover:text-red-600 flex items-center gap-1"
                          >
                            <MdDelete /> Remove
                          </button>
                        )}
                      </div>
                    </div>
                  </div>
                );
              })}
            </div>
          </div>

          {/* Submit Button */}
          <Button
            type="submit"
            color="blue"
            className="w-full"
            disabled={uploading}
          >
            {uploading ? "Uploading..." : "Upload Model"}
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
    selectedProteinRegion,
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
  const loading = isLoading || reduxLoading;

  // ============================================
  // FETCH MODELS ON MOUNT
  // ============================================
  const fetchModels = async () => {
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
        const uploadedModels = data.available_models || [];

        const formattedAPI = uploadedModels.map((m) => ({
          label: `${m} (Uploaded)`,
          value: `uploaded:${m}`,
          type: "uploaded",
        }));

        const allModels = [...formattedStatic, ...formattedAPI];
        setCombinedModelList(allModels);

        if (!selectedModel && allModels.length > 0) {
          setSelectedModel(allModels[0].value);
        }
      } else {
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
  };

  useEffect(() => {
    fetchModels();
    // eslint-disable-next-line
  }, []);

  // Also refresh models when component becomes visible (e.g., after navigating back from upload page)
  useEffect(() => {
    const handleFocus = () => fetchModels();
    window.addEventListener("focus", handleFocus);
    return () => window.removeEventListener("focus", handleFocus);
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

      // Only fetch for uploaded models
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
            
            // Convert to array format for display
            const paramsArray = Object.entries(params).map(([key, val]) => {
              // Handle both old format (just value) and new format (object with value, required, default)
              if (typeof val === 'object' && val !== null) {
                return {
                  key,
                  value: String(val.value ?? val.default ?? ""),
                  required: val.required || false,
                  default: val.default,
                };
              } else {
                return {
                  key,
                  value: String(val),
                  required: false,
                  default: val,
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
        // Static models don't have custom parameters
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

    // Build custom parameters object
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
      await dispatch(fetchPrediction(params)).unwrap();

      if (onSubmit) {
        await onSubmit(_nodeId, _elapsedDay, selectedModel, selectedProteinRegion);
      }
    } catch (error) {
      console.error("Prediction failed:", error);
      alert("Prediction failed: " + (error.message || error));
    }
  };

  // ============================================
  // RENDER
  // ============================================
  if (loading) {
    return <LoadingSpinner />;
  }

  return (
    <div className="min-h-[calc(100vh-3.5rem)] p-4 flex flex-col justify-center items-center bg-gradient-to-br from-gray-50 via-blue-50/30 to-gray-50 dark:from-gray-950 dark:via-gray-900 dark:to-gray-950 transition-colors">
      {/* Centered Logo */}
      <div className="flex justify-center items-center w-full">
        <img
          src={logo}
          alt="CovMutEx Logo"
          className="w-[18rem] h-[9rem] sm:w-64 sm:h-[9rem] md:w-80 md:h-[10rem] lg:w-[22rem] lg:h-[14rem] xl:w-[32rem] xl:h-[20rem] object-contain drop-shadow-sm"
        />
      </div>

      {/* Form Card */}
      <div className="w-full max-w-xl">
        <div className="bg-white dark:bg-gray-900 rounded-2xl shadow-xl dark:shadow-gray-900/50 border border-gray-200 dark:border-gray-800 p-6 animate-fade-in">
          {/* Step Indicator */}
          <div className="flex items-center justify-center gap-0 mb-6">
            {[
              { num: 1, label: "Select Model" },
              { num: 2, label: "Configure" },
              { num: 3, label: "Run" },
            ].map((step, idx) => (
              <div key={step.num} className="flex items-center">
                <div className="flex flex-col items-center">
                  <div className={`w-8 h-8 rounded-full flex items-center justify-center text-xs font-bold transition-colors ${
                    (step.num === 1 && selectedModel) || (step.num === 2 && _nodeId) || step.num === 3
                      ? "bg-blue-600 text-white dark:bg-blue-500"
                      : "bg-gray-200 dark:bg-gray-700 text-gray-500 dark:text-gray-400"
                  }`}>
                    {step.num}
                  </div>
                  <span className="text-[10px] text-gray-400 dark:text-gray-500 mt-1 font-medium">{step.label}</span>
                </div>
                {idx < 2 && (
                  <div className="w-16 h-px bg-gray-200 dark:bg-gray-700 mx-2 mb-4" />
                )}
              </div>
            ))}
          </div>

          <form className="space-y-4" onSubmit={handleSubmit}>
            {/* Model Selection */}
            <div>
              <div className="flex justify-between items-center mb-1">
                <label className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase block ml-1">
                  Prediction Model *
                </label>
                <button
                  type="button"
                  onClick={() => setIsUploadModalOpen(true)}
                  className="text-xs text-blue-600 dark:text-blue-400 font-bold hover:text-blue-800 dark:hover:text-blue-300 flex items-center gap-1 bg-blue-50 dark:bg-blue-900/30 px-2 py-1 rounded-lg transition-colors"
                >
                  <MdCloudUpload /> Upload New
                </button>
              </div>
              <Select
                options={combinedModelList}
                onChange={(opt) => setSelectedModel(opt.value)}
                value={combinedModelList.find((o) => o.value === selectedModel)}
                styles={customStyles}
                placeholder="Choose a model..."
                formatOptionLabel={(option) => (
                  <div className="flex items-center justify-between">
                    <span className="text-sm">{option.label}</span>
                    {option.type === "uploaded" && (
                      <span className="text-[10px] bg-green-100 text-green-700 font-bold px-1.5 py-0.5 rounded-full ml-2">
                        Uploaded
                      </span>
                    )}
                  </div>
                )}
              />
            </div>

            {/* Variant ID */}
            <div>
              <label className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block ml-1">
                Variant ID *
              </label>
              <DropDown items={nodes} setNodeId={setNodeId} />
            </div>

            {/* Elapsed Days */}
            <div>
              <label className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block ml-1">
                Elapsed Days *
              </label>
              <Input
                type="number"
                value={_elapsedDay}
                onChange={(e) => setElapsedDay(e.target.value)}
                min={1}
                placeholder="e.g., 60"
                required
                className="!border !border-gray-300 dark:!border-gray-600 focus:!border-blue-500 dark:!bg-gray-800 dark:!text-gray-200"
              />
              <p className="text-xs text-gray-400 dark:text-gray-500 mt-1 ml-1">
                Days since variant emergence (affects mutation probability)
              </p>
            </div>

            {/* Region */}
            <div>
              <label className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block ml-1">
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
                placeholder="Optional — defaults to whole genome"
                isClearable
              />
            </div>

            {/* Model Parameters (for uploaded models) */}
            {parametersLoading && (
              <div className="text-center py-4">
                <div className="animate-spin rounded-full h-6 w-6 border-2 border-blue-500 border-t-transparent mx-auto"></div>
                <p className="text-xs text-gray-500 mt-2">Loading parameters...</p>
              </div>
            )}

            {!parametersLoading && modelParameters.length > 0 && (
              <div className="mt-4 bg-gradient-to-br from-blue-50 to-indigo-50 rounded-xl border border-blue-100 overflow-hidden">
                <div className="bg-gradient-to-r from-blue-600 to-indigo-600 px-4 py-3 flex items-center gap-2">
                  <MdSettings className="text-white text-lg" />
                  <span className="text-sm font-bold text-white">Model Parameters</span>
                  <MdInfo className="text-blue-200" title="Parameters specific to this model" />
                </div>
                <div className="p-3 space-y-2">
                  {modelParameters.map((item, index) => {
                    const presetInfo = UPLOAD_PARAM_PRESETS.find(p => p.key === item.key);
                    return (
                      <div
                        key={index}
                        className={`bg-white rounded-xl border shadow-sm overflow-hidden ${
                          item.required ? "border-orange-200" : "border-gray-200"
                        }`}
                      >
                        <div className="px-3 py-2 bg-gray-50 border-b border-gray-100 flex items-center justify-between">
                          <div className="flex items-center gap-2">
                            <span className="text-xs font-bold text-gray-700 uppercase">
                              {presetInfo?.label || item.key}
                            </span>
                            {item.required && (
                              <span className="text-[10px] bg-red-100 text-red-600 px-1.5 py-0.5 rounded-full font-bold">Required</span>
                            )}
                          </div>
                          {item.default !== undefined && (
                            <span className="text-xs text-gray-400">
                              Default: {item.default}
                            </span>
                          )}
                        </div>
                        <div className="p-3">
                          {presetInfo?.options ? (
                            <div className="flex flex-wrap gap-1.5">
                              {presetInfo.options.map(opt => (
                                <button
                                  key={opt}
                                  type="button"
                                  onClick={() => updateModelParam(index, opt)}
                                  className={`px-3 py-1.5 rounded-lg text-xs font-medium border transition-all ${
                                    item.value === opt
                                      ? "border-blue-500 bg-blue-50 text-blue-700"
                                      : "border-gray-200 bg-gray-50 text-gray-600 hover:border-blue-300"
                                  }`}
                                >
                                  {opt}
                                </button>
                              ))}
                            </div>
                          ) : (
                            <input
                              type="text"
                              className={`w-full text-sm font-medium text-gray-800 border-2 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-blue-500 transition-all ${
                                item.required && !item.value ? "border-red-300 bg-red-50" : "border-gray-200"
                              }`}
                              value={item.value}
                              onChange={(e) => updateModelParam(index, e.target.value)}
                              placeholder={item.required ? "Required" : "Optional"}
                              required={item.required}
                            />
                          )}
                          {presetInfo?.hint && (
                            <p className="text-xs text-blue-500 mt-1.5">{presetInfo.hint}</p>
                          )}
                        </div>
                      </div>
                    );
                  })}
                </div>
              </div>
            )}

            {/* Submit Button */}
            <Button
              size="lg"
              color="blue"
              type="submit"
              className="w-full flex justify-center items-center gap-2 shadow-lg shadow-blue-500/20 hover:shadow-blue-500/40 mt-6 rounded-xl"
              disabled={loading}
            >
              {loading ? "Processing..." : "Run Prediction"}
              <MdOutlineCreate className="text-lg" />
            </Button>
          </form>
        </div>
      </div>

      {/* Upload Modal */}
      <UploadModal
        isOpen={isUploadModalOpen}
        onClose={() => setIsUploadModalOpen(false)}
        onSuccess={fetchModels}
      />
    </div>
  );
}

export default Navbar;