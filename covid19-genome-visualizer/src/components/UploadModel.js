import React, { useState, useEffect, useCallback } from "react";
import { useDispatch, useSelector } from "react-redux";
import { useNavigate } from "react-router-dom";
import { nodeIds as staticNodes } from "../data/nodeIds";
import { modelList as staticModels } from "../data/modelList";
import {
  fetchPrediction,
  fetchAvailableModels,
} from "../features/genome/genomeSlice";
import {
  MdCloudUpload,
  MdAdd,
  MdDelete,
  MdArrowBack,
  MdPlayArrow,
  MdSettings,
  MdCheckCircle,
  MdInfo,
} from "react-icons/md";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";
const DEFAULT_NODE_ID = "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07";

/**
 * UploadModel Component
 * 
 * Standalone page for uploading new models and running predictions.
 * Supports both new uploads and existing model selection.
 * Parameters support: required flag, default value, and editable value.
 */
const UploadModel = () => {
  const dispatch = useDispatch();
  const navigate = useNavigate();

  // Redux state
  const { loading, error } = useSelector((state) => state.genome);

  // Tab state
  const [activeTab, setActiveTab] = useState("upload");

  // Common fields
  const [nodeId, setNodeId] = useState(DEFAULT_NODE_ID);
  const [elapsedDay, setElapsedDay] = useState("60");
  const [selectedProteinRegion, setSelectedProteinRegion] = useState("");

  // Existing model fields
  const [selectedModel, setSelectedModel] = useState("balanced_data_model");
  const [existingModelList, setExistingModelList] = useState([]);
  const [existingModelParams, setExistingModelParams] = useState([]);
  const [paramsLoading, setParamsLoading] = useState(false);

  // Upload fields
  const [uploadName, setUploadName] = useState("");
  const [modelFile, setModelFile] = useState(null);
  const [extractorFile, setExtractorFile] = useState(null);
  const [helperFiles, setHelperFiles] = useState([]);
  const [customParams, setCustomParams] = useState([
    { key: "batch_size", value: "32", required: false },
  ]);

  // UI state
  const [uploadError, setUploadError] = useState("");
  const [uploadSuccess, setUploadSuccess] = useState("");

  // ============================================
  // FETCH MODELS LIST
  // ============================================
  const fetchModelsList = useCallback(async () => {
    try {
      const response = await fetch(`${API_URL}/api/models/`);
      if (response.ok) {
        const data = await response.json();
        
        // Server models
        const serverModels = (data.server_models || []).map((name) => ({
          name,
          value: name,
          type: "server",
        }));

        // Uploaded models
        const uploadedModels = (data.uploaded_models || []).map((m) => ({
          name: `${m.folder_name} (Uploaded)`,
          value: `uploaded:${m.folder_name}`,
          type: "uploaded",
          hasParams: m.has_parameters,
          hasExtractor: m.has_extractor,
        }));

        setExistingModelList([...serverModels, ...uploadedModels]);
      }
    } catch (err) {
      console.warn("Could not fetch models:", err.message);
      // Fallback to static models
      setExistingModelList(
        staticModels.map((m) => ({
          name: m.name,
          value: m.path,
          type: "server",
        }))
      );
    }
  }, []);

  useEffect(() => {
    fetchModelsList();
    dispatch(fetchAvailableModels());
  }, [dispatch, fetchModelsList]);

  // ============================================
  // FETCH PARAMETERS WHEN EXISTING MODEL CHANGES
  // ============================================
  useEffect(() => {
    const fetchParams = async () => {
      if (!selectedModel || !selectedModel.startsWith("uploaded:")) {
        setExistingModelParams([]);
        return;
      }

      setParamsLoading(true);
      try {
        const modelName = selectedModel.replace("uploaded:", "");
        const response = await fetch(
          `${API_URL}/api/model-parameters/?model_name=${encodeURIComponent(modelName)}`
        );

        if (response.ok) {
          const data = await response.json();
          const params = data.parameters || {};

          const paramsArray = Object.entries(params).map(([key, val]) => {
            if (typeof val === "object" && val !== null) {
              return {
                key,
                value: String(val.value ?? val.default ?? ""),
                required: val.required || false,
                default: val.default !== undefined ? String(val.default) : "",
              };
            } else {
              return {
                key,
                value: String(val),
                required: false,
                default: String(val),
              };
            }
          });

          setExistingModelParams(paramsArray);
        } else {
          setExistingModelParams([]);
        }
      } catch (err) {
        console.warn("Could not fetch model parameters:", err);
        setExistingModelParams([]);
      } finally {
        setParamsLoading(false);
      }
    };

    if (activeTab === "existing") {
      fetchParams();
    }
  }, [selectedModel, activeTab]);

  // ============================================
  // UPLOAD PARAMETER HANDLERS
  // ============================================
  const addParam = () => {
    setCustomParams([...customParams, { key: "", value: "", required: false }]);
  };

  const removeParam = (index) => {
    const list = [...customParams];
    list.splice(index, 1);
    setCustomParams(list);
  };

  const updateParam = (index, field, value) => {
    const list = [...customParams];
    list[index][field] = value;
    setCustomParams(list);
  };

  // ============================================
  // EXISTING MODEL PARAMETER HANDLERS
  // ============================================
  const updateExistingParam = (index, newValue) => {
    const updated = [...existingModelParams];
    updated[index].value = newValue;
    setExistingModelParams(updated);
  };

  const resetExistingParam = (index) => {
    const updated = [...existingModelParams];
    updated[index].value = updated[index].default || "";
    setExistingModelParams(updated);
  };

  const handleHelperFilesChange = (e) => {
    if (e.target.files) {
      setHelperFiles(Array.from(e.target.files));
    }
  };

  // ============================================
  // SUBMIT HANDLER
  // ============================================
  const handleSubmit = async (e) => {
    e.preventDefault();
    setUploadError("");
    setUploadSuccess("");

    let params = {
      nodeId: nodeId || DEFAULT_NODE_ID,
      elapsedDay: Number(elapsedDay) || 60,
      selectedProteinRegion: selectedProteinRegion || null,
    };

    if (activeTab === "existing") {
      // Validate required parameters for existing model
      const missingRequired = existingModelParams.filter(
        (p) => p.required && (!p.value || p.value.trim() === "")
      );
      if (missingRequired.length > 0) {
        setUploadError(
          `Please fill required parameters: ${missingRequired.map((p) => p.key).join(", ")}`
        );
        return;
      }

      // Build custom parameters object
      const paramsObj = {};
      existingModelParams.forEach((item) => {
        if (item.key.trim()) {
          const isNum = !isNaN(item.value) && item.value.trim() !== "";
          paramsObj[item.key.trim()] = isNum ? parseFloat(item.value) : item.value;
        }
      });

      params.selectedModel = selectedModel;
      params.isNewUpload = false;
      params.customParameters = paramsObj;
    } else {
      // Upload tab
      if (!uploadName.trim() || !modelFile) {
        setUploadError("Please provide a Model Name and Model File.");
        return;
      }

      // Build custom parameters with {value, required, default} structure
      const paramsObj = {};
      customParams.forEach((item) => {
        if (item.key.trim()) {
          const rawValue = item.value.trim();
          const numericValue = !isNaN(rawValue) && rawValue !== "" ? parseFloat(rawValue) : rawValue;
          paramsObj[item.key.trim()] = {
            value: numericValue,
            required: item.required || false,
            default: rawValue,
          };
        }
      });

      params.isNewUpload = true;
      params.uploadName = uploadName.trim();
      params.modelFile = modelFile;
      params.extractorFile = extractorFile;
      params.helperFiles = helperFiles;
      params.customParameters = paramsObj;
    }

    try {
      await dispatch(fetchPrediction(params)).unwrap();
      
      if (activeTab === "upload") {
        setUploadSuccess(`Model "${uploadName}" uploaded successfully!`);
        // Refresh model list
        await fetchModelsList();
      }
      
      navigate("/genome-mutation-visualization");
    } catch (err) {
      console.error("Prediction/Upload failed:", err);
      setUploadError(err.message || "Operation failed. Please try again.");
    }
  };

  // ============================================
  // RENDER
  // ============================================

  return (
    <div className="min-h-screen bg-gray-50 py-8 px-4">
      <div className="max-w-4xl mx-auto">
        {/* Header */}
        <div className="mb-8">
          <button
            onClick={() => navigate("/")}
            className="flex items-center gap-2 text-gray-600 hover:text-blue-600 mb-4 transition-colors"
          >
            <MdArrowBack /> Back to Home
          </button>
          <h1 className="text-3xl font-bold text-gray-800">
            Model Management & Prediction
          </h1>
          <p className="text-gray-500 mt-2">
            Upload a new model or run predictions with existing models.
          </p>
        </div>

        {/* Main Card */}
        <div className="bg-white rounded-2xl shadow-lg border border-gray-100 overflow-hidden">
          {/* Tabs */}
          <div className="flex border-b border-gray-200">
            <button
              className={`flex-1 py-4 px-6 text-center font-semibold transition-colors ${
                activeTab === "upload"
                  ? "border-b-2 border-green-500 text-green-600 bg-green-50/50"
                  : "text-gray-500 hover:text-green-600 hover:bg-gray-50"
              }`}
              onClick={() => setActiveTab("upload")}
            >
              <MdCloudUpload className="inline mr-2" />
              Upload New Model
            </button>
            <button
              className={`flex-1 py-4 px-6 text-center font-semibold transition-colors ${
                activeTab === "existing"
                  ? "border-b-2 border-blue-500 text-blue-600 bg-blue-50/50"
                  : "text-gray-500 hover:text-blue-600 hover:bg-gray-50"
              }`}
              onClick={() => setActiveTab("existing")}
            >
              <MdSettings className="inline mr-2" />
              Use Existing Model
            </button>
          </div>

          {/* Form */}
          <form onSubmit={handleSubmit} className="p-6 space-y-6">
            {/* Messages */}
            {uploadSuccess && (
              <div className="p-4 bg-green-50 border border-green-200 rounded-lg text-green-700 flex items-center gap-2">
                <MdCheckCircle size={20} />
                {uploadSuccess}
              </div>
            )}
            {uploadError && (
              <div className="p-4 bg-red-50 border border-red-200 rounded-lg text-red-700">
                {uploadError}
              </div>
            )}

            {/* ============================================ */}
            {/* COMMON FIELDS */}
            {/* ============================================ */}
            <div className="bg-gray-50 p-4 rounded-xl space-y-4">
              <h3 className="font-semibold text-gray-700">Prediction Parameters</h3>
              
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {/* Node ID */}
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Variant Node ID
                  </label>
                  <select
                    value={nodeId}
                    onChange={(e) => setNodeId(e.target.value)}
                    className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
                  >
                    {staticNodes && staticNodes.map((node, idx) => (
                      <option key={idx} value={node}>
                        {node.substring(0, 50)}...
                      </option>
                    ))}
                  </select>
                </div>

                {/* Elapsed Day */}
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Elapsed Days
                  </label>
                  <input
                    type="number"
                    value={elapsedDay}
                    onChange={(e) => setElapsedDay(e.target.value)}
                    min={0}
                    className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
                    placeholder="e.g., 60"
                  />
                </div>
              </div>
            </div>

            {/* ============================================ */}
            {/* UPLOAD TAB CONTENT */}
            {/* ============================================ */}
            {activeTab === "upload" && (
              <div className="space-y-4">
                {/* Model Name */}
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Model Name (Folder Name) *
                  </label>
                  <input
                    type="text"
                    value={uploadName}
                    onChange={(e) => setUploadName(e.target.value)}
                    className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-green-500 focus:border-green-500"
                    placeholder="e.g., MyCustomModel_v1"
                    required={activeTab === "upload"}
                  />
                </div>

                {/* File Uploads */}
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  {/* Model File */}
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Model File (.keras, .h5, .pt) *
                    </label>
                    <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-6 hover:border-green-400 transition-colors cursor-pointer">
                      <input
                        type="file"
                        accept=".keras,.h5,.pt,.pth"
                        onChange={(e) => setModelFile(e.target.files[0])}
                        className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
                        required={activeTab === "upload"}
                      />
                      <div className="text-center">
                        <MdCloudUpload className="mx-auto text-4xl text-gray-400 mb-2" />
                        <p className="text-sm text-gray-500">
                          {modelFile ? (
                            <span className="text-green-600 font-medium">{modelFile.name}</span>
                          ) : (
                            "Click or drag to upload"
                          )}
                        </p>
                      </div>
                    </div>
                  </div>

                  {/* Extractor File */}
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Feature Extractor (.py) - Optional
                    </label>
                    <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-6 hover:border-blue-400 transition-colors cursor-pointer">
                      <input
                        type="file"
                        accept=".py"
                        onChange={(e) => setExtractorFile(e.target.files[0])}
                        className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
                      />
                      <div className="text-center">
                        <MdSettings className="mx-auto text-4xl text-gray-400 mb-2" />
                        <p className="text-sm text-gray-500">
                          {extractorFile ? (
                            <span className="text-blue-600 font-medium">{extractorFile.name}</span>
                          ) : (
                            "Custom feature extractor"
                          )}
                        </p>
                      </div>
                    </div>
                  </div>
                </div>

                {/* Helper Files */}
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Additional Helper Files - Optional
                  </label>
                  <input
                    type="file"
                    multiple
                    onChange={handleHelperFilesChange}
                    className="block w-full text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:text-sm file:font-semibold file:bg-gray-100 file:text-gray-700 hover:file:bg-gray-200"
                  />
                  {helperFiles.length > 0 && (
                    <p className="text-xs text-gray-500 mt-1">
                      {helperFiles.length} file(s) selected: {helperFiles.map(f => f.name).join(", ")}
                    </p>
                  )}
                </div>

                {/* Custom Parameters for Upload */}
                <div>
                  <div className="flex justify-between items-center mb-2">
                    <label className="text-sm font-medium text-gray-700 flex items-center gap-1">
                      Custom Parameters
                      <MdInfo className="text-gray-400" title="Define parameters users must fill when using this model" />
                    </label>
                    <button
                      type="button"
                      onClick={addParam}
                      className="text-sm text-green-600 hover:text-green-800 flex items-center gap-1"
                    >
                      <MdAdd /> Add
                    </button>
                  </div>
                  <div className="space-y-2">
                    {customParams.map((item, index) => (
                      <div key={index} className="flex gap-2 items-center bg-gray-50 p-3 rounded-lg">
                        {/* Key */}
                        <input
                          type="text"
                          placeholder="Key (e.g., batch_size)"
                          value={item.key}
                          onChange={(e) => updateParam(index, "key", e.target.value)}
                          className="flex-1 border border-gray-300 rounded-lg px-3 py-2 text-sm"
                        />
                        {/* Default Value */}
                        <input
                          type="text"
                          placeholder="Default Value"
                          value={item.value}
                          onChange={(e) => updateParam(index, "value", e.target.value)}
                          className="flex-1 border border-gray-300 rounded-lg px-3 py-2 text-sm"
                        />
                        {/* Required Toggle */}
                        <label className="flex items-center gap-1.5 text-xs whitespace-nowrap cursor-pointer select-none">
                          <input
                            type="checkbox"
                            checked={item.required}
                            onChange={(e) => updateParam(index, "required", e.target.checked)}
                            className="rounded border-gray-300 text-green-600 focus:ring-green-500"
                          />
                          <span className={item.required ? "text-red-600 font-bold" : "text-gray-500"}>
                            Required
                          </span>
                        </label>
                        {customParams.length > 1 && (
                          <button
                            type="button"
                            onClick={() => removeParam(index)}
                            className="p-2 text-red-500 hover:text-red-700 hover:bg-red-50 rounded-lg transition-colors"
                          >
                            <MdDelete />
                          </button>
                        )}
                      </div>
                    ))}
                  </div>
                  <p className="text-xs text-gray-400 mt-1">
                    Parameters marked "Required" must be filled before running predictions with this model.
                  </p>
                </div>
              </div>
            )}

            {/* ============================================ */}
            {/* EXISTING MODEL TAB CONTENT */}
            {/* ============================================ */}
            {activeTab === "existing" && (
              <div className="space-y-4">
                {/* Model Selection */}
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Select Model
                  </label>
                  <select
                    value={selectedModel}
                    onChange={(e) => setSelectedModel(e.target.value)}
                    className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
                  >
                    {existingModelList.map((model, idx) => (
                      <option key={idx} value={model.value}>
                        {model.name}
                        {model.type === "uploaded" && model.hasParams ? " ⚙️" : ""}
                      </option>
                    ))}
                  </select>
                  {selectedModel && selectedModel.startsWith("uploaded:") && (
                    <p className="text-xs text-green-600 mt-1 flex items-center gap-1">
                      <MdCheckCircle size={12} /> Uploaded model selected
                    </p>
                  )}
                </div>

                {/* Parameters Loading */}
                {paramsLoading && (
                  <div className="text-center py-4">
                    <div className="animate-spin rounded-full h-6 w-6 border-2 border-blue-500 border-t-transparent mx-auto"></div>
                    <p className="text-xs text-gray-500 mt-2">Loading model parameters...</p>
                  </div>
                )}

                {/* Model Parameters for Existing Model */}
                {!paramsLoading && existingModelParams.length > 0 && (
                  <div className="p-4 bg-gradient-to-br from-blue-50 to-indigo-50 rounded-xl border border-blue-100">
                    <div className="flex items-center gap-2 mb-3">
                      <MdSettings className="text-blue-500" />
                      <label className="text-sm font-bold text-blue-800 uppercase">
                        Model Parameters
                      </label>
                    </div>
                    <div className="space-y-3">
                      {existingModelParams.map((item, index) => (
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
                              {item.required ? (
                                <span className="text-red-500 text-[10px] font-bold bg-red-50 px-1 rounded">
                                  REQUIRED
                                </span>
                              ) : (
                                <span className="text-gray-400 text-[10px] font-normal bg-gray-50 px-1 rounded">
                                  optional
                                </span>
                              )}
                            </span>
                            <div className="flex items-center gap-2">
                              {item.default && (
                                <span className="text-xs text-gray-400">
                                  Default: {item.default}
                                </span>
                              )}
                              {item.default && item.value !== item.default && (
                                <button
                                  type="button"
                                  onClick={() => resetExistingParam(index)}
                                  className="text-xs text-blue-500 hover:text-blue-700 underline"
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
                            onChange={(e) => updateExistingParam(index, e.target.value)}
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

                {/* No parameters info */}
                {!paramsLoading &&
                  selectedModel &&
                  selectedModel.startsWith("uploaded:") &&
                  existingModelParams.length === 0 && (
                    <div className="p-3 bg-gray-50 rounded-lg border border-gray-200 text-center">
                      <p className="text-xs text-gray-500">
                        This model has no custom parameters.
                      </p>
                    </div>
                  )}
              </div>
            )}

            {/* ============================================ */}
            {/* SUBMIT BUTTON */}
            {/* ============================================ */}
            <button
              type="submit"
              disabled={loading}
              className={`w-full py-3 px-6 rounded-xl font-semibold text-white flex items-center justify-center gap-2 transition-all ${
                activeTab === "upload"
                  ? "bg-green-600 hover:bg-green-700 disabled:bg-green-300"
                  : "bg-blue-600 hover:bg-blue-700 disabled:bg-blue-300"
              } disabled:cursor-not-allowed`}
            >
              {loading ? (
                <>
                  <div className="animate-spin rounded-full h-5 w-5 border-2 border-white border-t-transparent" />
                  Processing...
                </>
              ) : (
                <>
                  <MdPlayArrow className="text-xl" />
                  {activeTab === "upload" ? "Upload & Predict" : "Run Prediction"}
                </>
              )}
            </button>
          </form>
        </div>
      </div>
    </div>
  );
};

export default UploadModel;