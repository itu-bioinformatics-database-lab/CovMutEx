import React, { useState, useEffect } from "react";
import { useDispatch, useSelector } from "react-redux";
import { useNavigate } from "react-router-dom";
import { nodeIds as staticNodes } from "../data/nodeIds";
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
} from "react-icons/md";

const DEFAULT_NODE_ID = "USA/UT-UPHL-210820924226/2021|OK040008.1|2021-08-07";
const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

/**
 * UploadModel Component
 * 
 * Standalone page for uploading new models and running predictions.
 * Supports both new uploads and existing model selection.
 */
const UploadModel = () => {
  const dispatch = useDispatch();
  const navigate = useNavigate();

  // Redux state
  const { loading, availableModels, error } = useSelector((state) => state.genome);

  // Tab state
  const [activeTab, setActiveTab] = useState("upload");

  // Common fields
  const [nodeId, setNodeId] = useState(DEFAULT_NODE_ID);
  const [elapsedDay, setElapsedDay] = useState("60");
  const [selectedProteinRegion, setSelectedProteinRegion] = useState("");

  // Existing model fields
  const [selectedModel, setSelectedModel] = useState("balanced_data_model");

  // Upload fields
  const [uploadName, setUploadName] = useState("");
  const [modelFile, setModelFile] = useState(null);
  const [extractorFile, setExtractorFile] = useState(null);
  const [helperFiles, setHelperFiles] = useState([]);
  const [customParams, setCustomParams] = useState([
    { key: "batch_size", value: "32" },
  ]);

  // UI state
  const [uploadError, setUploadError] = useState("");

  // Fetch available models on mount
  useEffect(() => {
    dispatch(fetchAvailableModels());
  }, [dispatch]);

  // ============================================
  // HANDLERS
  // ============================================

  const addParam = () => {
    setCustomParams([...customParams, { key: "", value: "" }]);
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

  const handleHelperFilesChange = (e) => {
    if (e.target.files) {
      setHelperFiles(Array.from(e.target.files));
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    setUploadError("");

    // Build parameters object
    const paramsObj = {};
    customParams.forEach((item) => {
      if (item.key.trim()) {
        const isNum = !isNaN(item.value) && item.value.trim() !== "";
        paramsObj[item.key.trim()] = isNum ? parseFloat(item.value) : item.value;
      }
    });

    let params = {
      nodeId: nodeId || DEFAULT_NODE_ID,
      elapsedDay: Number(elapsedDay) || 60,
      selectedProteinRegion: selectedProteinRegion || null,
      customParameters: paramsObj,
    };

    if (activeTab === "existing") {
      // Using existing model
      params.selectedModel = selectedModel;
      params.isNewUpload = false;
    } else {
      // Uploading new model
      if (!uploadName.trim() || !modelFile) {
        setUploadError("Please provide a Model Name and Model File.");
        return;
      }

      params.isNewUpload = true;
      params.uploadName = uploadName.trim();
      params.modelFile = modelFile;
      params.extractorFile = extractorFile;
      params.helperFiles = helperFiles;
    }

    try {
      if (activeTab === "existing") {
        // Existing model → run prediction and show chart
        await dispatch(fetchPrediction(params)).unwrap();
        navigate("/genome-mutation-visualization");
      } else {
        // New upload → save model via predict endpoint, ignore prediction result, go home
        const formData = new FormData();
        formData.append("uploadFolderName", params.uploadName.replace(/\s+/g, "_"));
        formData.append("modelFile", params.modelFile);
        formData.append("nodeId", params.nodeId);
        formData.append("elapsedDay", String(params.elapsedDay));
        if (params.extractorFile) formData.append("extractorFile", params.extractorFile);
        if (params.helperFiles) {
          params.helperFiles.forEach((file, i) => {
            formData.append(`helperFile_${i}`, file);
            formData.append(`helperFileName_${i}`, file.name);
          });
        }
        if (params.customParameters && Object.keys(params.customParameters).length > 0) {
          formData.append("customParameters", JSON.stringify(params.customParameters));
        }

        const response = await fetch(`${API_URL}/api/predict/`, {
          method: "POST",
          body: formData,
        });

        // Whether prediction succeeded or failed, model is saved
        // Navigate home without showing prediction chart
        navigate("/");
      }
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
            {/* Error Message */}
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

                {/* Custom Parameters */}
                <div>
                  <div className="flex justify-between items-center mb-2">
                    <label className="text-sm font-medium text-gray-700">
                      Custom Parameters
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
                      <div key={index} className="flex gap-2">
                        <input
                          type="text"
                          placeholder="Key"
                          value={item.key}
                          onChange={(e) => updateParam(index, "key", e.target.value)}
                          className="flex-1 border border-gray-300 rounded-lg px-3 py-2 text-sm"
                        />
                        <input
                          type="text"
                          placeholder="Value"
                          value={item.value}
                          onChange={(e) => updateParam(index, "value", e.target.value)}
                          className="flex-1 border border-gray-300 rounded-lg px-3 py-2 text-sm"
                        />
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
                </div>
              </div>
            )}

            {/* ============================================ */}
            {/* EXISTING MODEL TAB CONTENT */}
            {/* ============================================ */}
            {activeTab === "existing" && (
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Select Model
                  </label>
                  <select
                    value={selectedModel}
                    onChange={(e) => setSelectedModel(e.target.value)}
                    className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
                  >
                    <option value="balanced_data_model">balanced_data_model (Default)</option>
                    <option value="model_for_1k">model_for_1k</option>
                    <option value="model_for_5k">model_for_5k</option>
                    {/* Uploaded models from API */}
                    {availableModels && availableModels.map((model, idx) => (
                      <option key={idx} value={`uploaded:${model}`}>
                        {model} (Uploaded)
                      </option>
                    ))}
                  </select>
                  <p className="text-xs text-gray-500 mt-1">
                    For uploaded models, use format: <code className="bg-gray-100 px-1 rounded">uploaded:FolderName</code>
                  </p>
                </div>
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