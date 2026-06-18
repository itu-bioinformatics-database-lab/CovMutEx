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

/** Predefined parameter presets for dropdown selection */
const PARAM_PRESETS = [
  {
    key: "batch_size",
    label: "Batch Size",
    type: "number",
    defaultValue: "32",
    placeholder: "e.g., 16, 32, 64",
    options: ["8", "16", "32", "64", "128", "256"],
    hint: "Number of samples per training batch",
  },
  {
    key: "learning_rate",
    label: "Learning Rate",
    type: "number",
    defaultValue: "0.001",
    placeholder: "e.g., 0.001",
    step: "0.0001",
    min: "0",
    hint: "Step size for optimizer",
  },
  {
    key: "epochs",
    label: "Epochs",
    type: "number",
    defaultValue: "100",
    placeholder: "e.g., 50, 100, 200",
    options: ["10", "25", "50", "100", "200", "500"],
    hint: "Number of training iterations",
  },
  {
    key: "optimizer",
    label: "Optimizer",
    type: "select",
    defaultValue: "adam",
    options: ["adam", "sgd", "rmsprop", "adamw", "adagrad"],
    hint: "Optimization algorithm",
  },
  {
    key: "dropout_rate",
    label: "Dropout Rate",
    type: "number",
    defaultValue: "0.2",
    placeholder: "0.0 - 1.0",
    step: "0.05",
    min: "0",
    max: "1",
    hint: "Fraction of neurons to drop during training",
  },
  {
    key: "sequence_length",
    label: "Sequence Length",
    type: "number",
    defaultValue: "29903",
    placeholder: "e.g., 29903",
    hint: "Input genome sequence length",
  },
  {
    key: "hidden_units",
    label: "Hidden Units",
    type: "number",
    defaultValue: "128",
    placeholder: "e.g., 64, 128, 256",
    options: ["32", "64", "128", "256", "512"],
    hint: "Number of neurons in hidden layers",
  },
  {
    key: "activation",
    label: "Activation Function",
    type: "select",
    defaultValue: "relu",
    options: ["relu", "sigmoid", "tanh", "softmax", "leaky_relu", "gelu"],
    hint: "Non-linear activation function",
  },
];

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

  // Organism targeting — drives how the backend resolves the reference genome
  // and protein_regions for this bundle. "custom" means the bundle ships its
  // own genome + (optionally) protein-region CSV as helper files.
  const [organism, setOrganism] = useState("covid");
  const [genomeFile, setGenomeFile] = useState("");
  const [proteinRegionsFile, setProteinRegionsFile] = useState("");

  // UI state
  const [uploadError, setUploadError] = useState("");
  // Track which param card has its dropdown open
  const [openDropdownIndex, setOpenDropdownIndex] = useState(null);

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
    setCustomParams(prev => {
      const list = prev.map((item, i) => i === index ? { ...item, [field]: value } : item);
      return list;
    });
  };

  const updateParamBoth = (index, key, value) => {
    setCustomParams(prev => {
      const list = prev.map((item, i) => i === index ? { ...item, key, value } : item);
      return list;
    });
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

      // Custom organism requires a genome helper file reference.
      if (organism === "custom" && !genomeFile.trim()) {
        setUploadError(
          'For a custom organism, enter the Genome File name (a helper file you upload below, e.g. "my_genome.fasta").'
        );
        return;
      }

      params.isNewUpload = true;
      params.uploadName = uploadName.trim();
      params.modelFile = modelFile;
      params.extractorFile = extractorFile;
      params.helperFiles = helperFiles;
      params.organism = organism;
      params.genomeFile = genomeFile.trim();
      params.proteinRegionsFile = proteinRegionsFile.trim();
      // No variant at upload time anymore — the user picks the strain on
      // the prediction screen, so a single influenza bundle can be re-run
      // against any of the 9 cataloged HA strains.
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
        // Organism dispatch — written into bundle_metadata.json server-side.
        formData.append("organism", params.organism);
        if (params.genomeFile) formData.append("genome_file", params.genomeFile);
        if (params.proteinRegionsFile)
          formData.append("protein_regions_file", params.proteinRegionsFile);
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

                {/* Organism selector */}
                <div className="bg-indigo-50/50 border border-indigo-100 p-4 rounded-xl space-y-3">
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Target Organism
                    </label>
                    <select
                      value={organism}
                      onChange={(e) => setOrganism(e.target.value)}
                      className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
                    >
                      <option value="covid">SARS-CoV-2 (built-in)</option>
                      <option value="influenza">Influenza A — HA Segment (built-in)</option>
                      <option value="custom">Other — I'll upload my own genome</option>
                    </select>
                    <p className="text-xs text-gray-500 mt-1">
                      Built-in organisms use our reference genome &amp; protein regions.
                      Choose <span className="font-medium">Other</span> to predict on a virus we don't ship.
                      For Influenza A you'll pick the specific HA strain (PR/8/34, Cal/07,
                      cattle/Texas, etc.) on the prediction screen — not here.
                    </p>
                  </div>

                  {/* Custom-organism helper-file references */}
                  {organism === "custom" && (
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-3 pt-1">
                      <div>
                        <label className="block text-sm font-medium text-gray-700 mb-1">
                          Genome File Name *
                        </label>
                        <input
                          type="text"
                          value={genomeFile}
                          onChange={(e) => setGenomeFile(e.target.value)}
                          className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
                          placeholder="my_genome.fasta"
                        />
                      </div>
                      <div>
                        <label className="block text-sm font-medium text-gray-700 mb-1">
                          Protein Regions File — Optional
                        </label>
                        <input
                          type="text"
                          value={proteinRegionsFile}
                          onChange={(e) => setProteinRegionsFile(e.target.value)}
                          className="w-full border border-gray-300 rounded-lg px-3 py-2 focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
                          placeholder="my_protein_regions.csv"
                        />
                      </div>
                      <p className="text-xs text-gray-500 md:col-span-2">
                        Upload these as <span className="font-medium">Helper Files</span> below. The
                        names here must match the uploaded file names exactly. Protein-regions CSV
                        format: <span className="font-mono">name,start,end</span> (1-based inclusive),
                        one protein per line.
                      </p>
                    </div>
                  )}
                </div>

                {/* File Uploads */}
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  {/* Model File */}
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      Model File *
                      <span className="ml-1 text-xs text-gray-400 font-normal">.pt .pth .keras .h5 .onnx</span>
                    </label>
                    <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-6 hover:border-green-400 transition-colors cursor-pointer">
                      <input
                        type="file"
                        accept=".keras,.h5,.pt,.pth,.onnx,.bin,.safetensors"
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

                {/* ============================================ */}
                {/* CUSTOM PARAMETERS - Custom Dropdown UI */}
                {/* ============================================ */}
                <div className="bg-white border-2 border-green-200 rounded-2xl overflow-hidden">
                  {/* Section Header */}
                  <div className="bg-gradient-to-r from-green-600 to-emerald-600 px-5 py-4 flex justify-between items-center">
                    <div className="text-white">
                      <h3 className="font-bold text-base flex items-center gap-2">
                        <MdSettings className="text-xl" />
                        Model Hyperparameters
                      </h3>
                      <p className="text-green-100 text-xs mt-0.5">Configure training parameters for your model</p>
                    </div>
                    <button
                      type="button"
                      onClick={addParam}
                      className="flex items-center gap-1.5 bg-white/20 hover:bg-white/30 backdrop-blur text-white text-sm font-medium px-4 py-2 rounded-lg transition-colors border border-white/30"
                    >
                      <MdAdd className="text-lg" /> Add
                    </button>
                  </div>

                  <div className="p-5 space-y-4">
                    {customParams.map((item, index) => {
                      const isPreset = PARAM_PRESETS.some(p => p.key === item.key);
                      const preset = PARAM_PRESETS.find(p => p.key === item.key);
                      const isDropdownOpen = openDropdownIndex === index;

                      return (
                        <div key={index} className="border border-gray-200 rounded-xl bg-gray-50 overflow-hidden">
                          {/* Card Header */}
                          <div className="px-4 py-2 bg-white border-b border-gray-100 flex items-center justify-between">
                            <div className="flex items-center gap-2">
                              <span className="w-6 h-6 rounded-full bg-green-100 text-green-700 text-xs font-bold flex items-center justify-center">
                                {index + 1}
                              </span>
                              <span className="text-sm font-medium text-gray-600">
                                {isPreset ? preset.label : (item.key || "New Parameter")}
                              </span>
                              {isPreset && (
                                <span className="text-xs bg-green-100 text-green-700 px-2 py-0.5 rounded-full font-medium">
                                  Preset
                                </span>
                              )}
                            </div>
                            {customParams.length > 1 && (
                              <button
                                type="button"
                                onClick={() => removeParam(index)}
                                className="flex items-center gap-1 text-xs text-red-400 hover:text-red-600 hover:bg-red-50 px-2 py-1 rounded transition-colors"
                              >
                                <MdDelete /> Remove
                              </button>
                            )}
                          </div>

                          <div className="p-4 space-y-3">
                            {/* Parameter Name - Custom Dropdown */}
                            <div>
                              <label className="block text-xs font-bold text-gray-500 uppercase tracking-wider mb-1.5">
                                Parameter Name
                              </label>
                              <div className="relative">
                                {/* Dropdown Trigger Button */}
                                <button
                                  type="button"
                                  onClick={() => setOpenDropdownIndex(isDropdownOpen ? null : index)}
                                  className={`w-full text-left border-2 rounded-xl px-4 py-3 text-sm font-medium transition-all flex items-center justify-between ${
                                    isDropdownOpen
                                      ? "border-green-500 ring-2 ring-green-200 bg-white"
                                      : "border-gray-200 bg-white hover:border-green-300"
                                  }`}
                                >
                                  <span className={item.key ? "text-gray-800" : "text-gray-400"}>
                                    {isPreset ? preset.label : (item.key || "Click to select a parameter...")}
                                  </span>
                                  <svg className={`w-5 h-5 text-gray-400 transition-transform ${isDropdownOpen ? "rotate-180" : ""}`} fill="none" viewBox="0 0 24 24" stroke="currentColor">
                                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                                  </svg>
                                </button>

                                {/* Dropdown Panel */}
                                {isDropdownOpen && (
                                  <div className="absolute z-50 mt-1 w-full bg-white border-2 border-green-200 rounded-xl shadow-xl max-h-72 overflow-y-auto">
                                    {/* Preset Options */}
                                    <div className="px-3 py-2 bg-gray-50 border-b border-gray-100">
                                      <span className="text-xs font-bold text-gray-500 uppercase tracking-wider">Common Presets</span>
                                    </div>
                                    {PARAM_PRESETS.map(p => (
                                      <button
                                        key={p.key}
                                        type="button"
                                        onClick={() => {
                                          updateParamBoth(index, p.key, p.defaultValue || "");
                                          setOpenDropdownIndex(null);
                                        }}
                                        className={`w-full text-left px-4 py-3 hover:bg-green-50 transition-colors border-b border-gray-50 flex items-center justify-between group ${
                                          item.key === p.key ? "bg-green-50" : ""
                                        }`}
                                      >
                                        <div>
                                          <div className="text-sm font-semibold text-gray-800 group-hover:text-green-700">
                                            {p.label}
                                          </div>
                                          <div className="text-xs text-gray-400 mt-0.5">{p.hint}</div>
                                        </div>
                                        {item.key === p.key && (
                                          <span className="text-green-600 text-lg">✓</span>
                                        )}
                                        {p.defaultValue && item.key !== p.key && (
                                          <span className="text-xs text-gray-300 bg-gray-100 px-2 py-0.5 rounded">
                                            default: {p.defaultValue}
                                          </span>
                                        )}
                                      </button>
                                    ))}

                                    {/* Custom Option */}
                                    <div className="px-3 py-2 bg-gray-50 border-t border-gray-100">
                                      <span className="text-xs font-bold text-gray-500 uppercase tracking-wider">Custom</span>
                                    </div>
                                    <button
                                      type="button"
                                      onClick={() => {
                                        updateParamBoth(index, "", "");
                                        setOpenDropdownIndex(null);
                                      }}
                                      className="w-full text-left px-4 py-3 hover:bg-blue-50 transition-colors flex items-center gap-2"
                                    >
                                      <MdAdd className="text-blue-500" />
                                      <span className="text-sm font-medium text-blue-600">Write custom parameter name...</span>
                                    </button>
                                  </div>
                                )}
                              </div>

                              {/* Custom parameter name input - only if not a preset */}
                              {!isPreset && !isDropdownOpen && item.key !== undefined && (
                                <input
                                  type="text"
                                  placeholder="Type your custom parameter name..."
                                  value={item.key}
                                  onChange={(e) => updateParam(index, "key", e.target.value)}
                                  className="w-full border-2 border-dashed border-gray-300 rounded-xl px-4 py-2.5 text-sm mt-2 focus:border-green-400 focus:ring-2 focus:ring-green-200 transition-all bg-white"
                                />
                              )}
                            </div>

                            {/* Value Input */}
                            <div>
                              <label className="block text-xs font-bold text-gray-500 uppercase tracking-wider mb-1.5">
                                Value
                              </label>
                              {preset?.options ? (
                                <div className="grid grid-cols-3 sm:grid-cols-4 md:grid-cols-6 gap-2">
                                  {preset.options.map(opt => (
                                    <button
                                      key={opt}
                                      type="button"
                                      onClick={() => updateParam(index, "value", opt)}
                                      className={`px-3 py-2 rounded-lg text-sm font-medium border-2 transition-all ${
                                        item.value === opt
                                          ? "border-green-500 bg-green-50 text-green-700 shadow-sm"
                                          : "border-gray-200 bg-white text-gray-600 hover:border-green-300 hover:bg-green-50/50"
                                      }`}
                                    >
                                      {opt}
                                    </button>
                                  ))}
                                </div>
                              ) : (
                                <input
                                  type={preset?.type === "number" ? "number" : "text"}
                                  placeholder={preset?.placeholder || "Enter value..."}
                                  value={item.value}
                                  onChange={(e) => updateParam(index, "value", e.target.value)}
                                  step={preset?.step}
                                  min={preset?.min}
                                  max={preset?.max}
                                  className="w-full border-2 border-gray-200 rounded-xl px-4 py-2.5 text-sm focus:border-green-500 focus:ring-2 focus:ring-green-200 transition-all bg-white"
                                />
                              )}
                            </div>
                          </div>
                        </div>
                      );
                    })}
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