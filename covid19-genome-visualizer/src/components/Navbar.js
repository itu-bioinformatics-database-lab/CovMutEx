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
  const [adapterFile, setAdapterFile] = useState(null);
  const [helperFiles, setHelperFiles] = useState([]);
  const [organism, setOrganism] = useState("covid");
  // Custom-organism uploads. The user drops the actual files; the backend
  // stores them under canonical names (genome.fasta / protein_regions.csv), so
  // there's no fragile "typed name must match a helper file" step.
  const [genomeUpload, setGenomeUpload] = useState(null);
  const [proteinRegionsUpload, setProteinRegionsUpload] = useState(null);
  const [paramsList, setParamsList] = useState([]);
  const [uploading, setUploading] = useState(false);
  const [error, setError] = useState("");

  if (!isOpen) return null;

  const addParam = () =>
    setParamsList([...paramsList, { key: "", value: "", required: false, type: "text", options: "" }]);

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
      const incoming = Array.from(e.target.files);
      setHelperFiles((prev) => {
        const existingNames = new Set(prev.map((f) => f.name));
        const deduped = incoming.filter((f) => !existingNames.has(f.name));
        return [...prev, ...deduped];
      });
      e.target.value = "";
    }
  };

  const removeHelperFile = (name) => {
    setHelperFiles((prev) => prev.filter((f) => f.name !== name));
  };

  const resetForm = () => {
    setName("");
    setModelFile(null);
    setExtractorFile(null);
    setAdapterFile(null);
    setHelperFiles([]);
    setOrganism("covid");
    setGenomeUpload(null);
    setProteinRegionsUpload(null);
    setParamsList([]);
    setError("");
  };

  const handleUpload = async (e) => {
    e.preventDefault();
    setError("");

    if (!modelFile || !name.trim()) {
      setError("Please fill required fields (Name, Model File).");
      return;
    }

    if (organism === "custom" && !genomeUpload) {
      setError("For a custom organism, drop your genome file (FASTA) into the Genome slot.");
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

      // Optional extractor + adapter
      if (extractorFile) {
        formData.append("extractorFile", extractorFile);
      }
      if (adapterFile) {
        formData.append("adapterFile", adapterFile);
      }

      // Organism dispatch — written into bundle_metadata.json server-side.
      // For custom organisms we send the actual files; the backend stores them
      // under canonical names and records those names in the metadata.
      formData.append("organism", organism);
      if (organism === "custom") {
        if (genomeUpload) formData.append("genomeFile", genomeUpload);
        if (proteinRegionsUpload)
          formData.append("proteinRegionsFile", proteinRegionsUpload);
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
          const paramEntry = {
            value: isNaN(item.value) ? item.value : parseFloat(item.value),
            required: item.required || false,
            default: item.value,
            type: item.type || "text",
          };
          if ((item.type === "radio" || item.type === "dropdown") && item.options) {
            paramEntry.options = item.options.split(",").map(o => o.trim()).filter(Boolean);
          }
          customParamsObj[item.key.trim()] = paramEntry;
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

          {/* Organism selector */}
          <div className="bg-indigo-50/60 border border-indigo-100 rounded-lg p-3 space-y-3">
            <div>
              <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
                Target Organism
              </label>
              <select
                value={organism}
                onChange={(e) => setOrganism(e.target.value)}
                className="w-full border border-gray-300 rounded-lg px-3 py-2 text-sm focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500"
              >
                <option value="covid">SARS-CoV-2 (built-in)</option>
                <option value="influenza">Influenza A — HA Segment (built-in)</option>
                <option value="custom">Other — I'll upload my own genome</option>
              </select>
              <p className="text-[11px] text-gray-500 mt-1 leading-snug">
                Built-in organisms use our reference genome &amp; protein regions — no genome
                upload needed. Choose <span className="font-medium">Other</span> to predict on a
                virus we don't ship.
              </p>
            </div>

            {organism === "custom" && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                {/* Genome file drop slot */}
                <div>
                  <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
                    Genome File (FASTA) *
                  </label>
                  <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-3 hover:border-indigo-400 transition-colors">
                    <input
                      type="file"
                      accept=".fasta,.fa,.txt"
                      onChange={(e) => setGenomeUpload(e.target.files[0] || null)}
                      className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
                    />
                    <div className="text-center">
                      <MdCloudUpload className="mx-auto text-gray-400 text-xl mb-1" />
                      <p className="text-[11px] text-gray-500 break-all">
                        {genomeUpload ? (
                          <span className="text-indigo-600 font-medium">{genomeUpload.name}</span>
                        ) : (
                          "Drop genome"
                        )}
                      </p>
                    </div>
                  </div>
                </div>
                {/* Protein-regions file drop slot */}
                <div>
                  <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
                    Protein Regions (CSV)
                  </label>
                  <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-3 hover:border-indigo-400 transition-colors">
                    <input
                      type="file"
                      accept=".csv"
                      onChange={(e) => setProteinRegionsUpload(e.target.files[0] || null)}
                      className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
                    />
                    <div className="text-center">
                      <MdCloudUpload className="mx-auto text-gray-400 text-xl mb-1" />
                      <p className="text-[11px] text-gray-500 break-all">
                        {proteinRegionsUpload ? (
                          <span className="text-indigo-600 font-medium">{proteinRegionsUpload.name}</span>
                        ) : (
                          "Optional"
                        )}
                      </p>
                    </div>
                  </div>
                </div>
                <p className="text-[11px] text-gray-500 md:col-span-2 leading-snug">
                  We store these under fixed names (<span className="font-mono">genome.fasta</span>,{" "}
                  <span className="font-mono">protein_regions.csv</span>) in your bundle, next to the
                  model &amp; feature extractor. Protein-regions CSV format:{" "}
                  <span className="font-mono">name,start,end</span> (1-based inclusive), one per line.
                </p>
              </div>
            )}
          </div>

          {/* File Uploads */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {/* Model File */}
            <div>
              <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
                Model File (.keras/.h5/.pt/.pth/.bin) *
              </label>
              <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-blue-400 transition-colors">
                <input
                  type="file"
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

          {/* Model Adapter File — postprocess output to a different payload shape */}
          <div>
            <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
              Model Adapter (.py)
            </label>
            <div className="relative border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-purple-400 transition-colors">
              <input
                type="file"
                accept=".py"
                onChange={(e) => setAdapterFile(e.target.files[0])}
                className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
              />
              <div className="text-center">
                <MdSettings className="mx-auto text-gray-400 text-2xl mb-1" />
                <p className="text-xs text-gray-500">
                  {adapterFile ? (
                    <span className="text-purple-600 font-medium">{adapterFile.name}</span>
                  ) : (
                    "Optional — exposes preprocess / predict / postprocess to override the default Keras model wrapper"
                  )}
                </p>
              </div>
            </div>
          </div>

          {/* Helper Files */}
          <div>
            <label className="text-xs font-bold text-gray-700 uppercase mb-1 block">
              Helper Files (Optional)
            </label>
            <div
              className="relative border-2 border-dashed border-gray-300 rounded-lg p-3 hover:border-blue-400 transition-colors text-center"
              onDragOver={(e) => { e.preventDefault(); e.currentTarget.classList.add("border-blue-400"); }}
              onDragLeave={(e) => e.currentTarget.classList.remove("border-blue-400")}
              onDrop={(e) => {
                e.preventDefault();
                e.currentTarget.classList.remove("border-blue-400");
                const dropped = Array.from(e.dataTransfer.files);
                if (dropped.length) {
                  setHelperFiles((prev) => {
                    const existingNames = new Set(prev.map((f) => f.name));
                    return [...prev, ...dropped.filter((f) => !existingNames.has(f.name))];
                  });
                }
              }}
            >
              <input
                type="file"
                multiple
                onChange={handleHelperFilesChange}
                className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
              />
              <p className="text-xs text-gray-500 pointer-events-none">
                Click or drag & drop files here
              </p>
            </div>
            {helperFiles.length > 0 && (
              <ul className="mt-2 space-y-1">
                {helperFiles.map((f) => (
                  <li key={f.name} className="flex items-center justify-between text-xs bg-gray-50 rounded px-2 py-1">
                    <span className="text-gray-700 truncate mr-2">{f.name}</span>
                    <button
                      type="button"
                      onClick={() => removeHelperFile(f.name)}
                      className="text-gray-400 hover:text-red-500 shrink-0"
                    >
                      ✕
                    </button>
                  </li>
                ))}
              </ul>
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
              {paramsList.length === 0 && (
                <p className="text-xs text-gray-400 text-center py-4">
                  No parameters added. Click "Add" to define model parameters.
                </p>
              )}
              {paramsList.map((item, index) => (
                <div key={index} className="bg-white rounded-xl border border-gray-200 overflow-hidden shadow-sm">
                  <div className="p-3 space-y-2">
                    {/* Parameter Name */}
                    <div>
                      <label className="text-[10px] font-bold text-gray-400 uppercase tracking-wider mb-1 block">Parameter Name</label>
                      <input
                        type="text"
                        placeholder="e.g. batch_size, learning_rate..."
                        value={item.key}
                        onChange={(e) => updateParam(index, "key", e.target.value)}
                        className="w-full border-2 border-gray-200 rounded-lg px-3 py-2 text-sm focus:border-blue-500 focus:ring-2 focus:ring-blue-200"
                      />
                    </div>
                    {/* Default Value */}
                    <div>
                      <label className="text-[10px] font-bold text-gray-400 uppercase tracking-wider mb-1 block">Default Value</label>
                      <input
                        type="text"
                        placeholder="e.g. 32, 0.001, adam..."
                        value={item.value}
                        onChange={(e) => updateParam(index, "value", e.target.value)}
                        className="w-full border-2 border-gray-200 rounded-lg px-3 py-2 text-sm focus:border-blue-500 focus:ring-2 focus:ring-blue-200"
                      />
                    </div>
                    {/* Input Type */}
                    <div>
                      <label className="text-[10px] font-bold text-gray-400 uppercase tracking-wider mb-1 block">Input Type</label>
                      <div className="flex gap-1.5">
                        {[
                          { value: "text", label: "Text" },
                          { value: "radio", label: "Radio Button" },
                          { value: "dropdown", label: "Dropdown" },
                        ].map((opt) => (
                          <button
                            key={opt.value}
                            type="button"
                            onClick={() => updateParam(index, "type", opt.value)}
                            className={`px-3 py-1.5 rounded-lg text-xs font-medium border transition-all ${
                              (item.type || "text") === opt.value
                                ? "border-blue-500 bg-blue-50 text-blue-700"
                                : "border-gray-200 bg-gray-50 text-gray-600 hover:border-blue-300"
                            }`}
                          >
                            {opt.label}
                          </button>
                        ))}
                      </div>
                    </div>
                    {/* Options (only for radio/dropdown) */}
                    {(item.type === "radio" || item.type === "dropdown") && (
                      <div>
                        <label className="text-[10px] font-bold text-gray-400 uppercase tracking-wider mb-1 block">
                          Options <span className="normal-case font-normal text-gray-300">(comma separated)</span>
                        </label>
                        <input
                          type="text"
                          placeholder="e.g. adam, sgd, rmsprop"
                          value={item.options || ""}
                          onChange={(e) => updateParam(index, "options", e.target.value)}
                          className="w-full border-2 border-gray-200 rounded-lg px-3 py-2 text-sm focus:border-blue-500 focus:ring-2 focus:ring-blue-200"
                        />
                      </div>
                    )}
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
                      <button
                        type="button"
                        onClick={() => removeParam(index)}
                        className="text-xs text-red-400 hover:text-red-600 flex items-center gap-1"
                      >
                        <MdDelete /> Remove
                      </button>
                    </div>
                  </div>
                </div>
              ))}
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
  // Faz 1.6: organism-aware variant picker. The selected model's organism
  // (covid / influenza / custom) decides which picker we render below the
  // model dropdown. _variant is the influenza HA strain choice; nodeId
  // (state above) stays the COVID isolate choice. Catalog comes from
  // /api/organisms/.
  const [_variant, setVariant] = useState("");
  const [influenzaVariants, setInfluenzaVariants] = useState([]);
  // modelName → organism, populated alongside combinedModelList so we can
  // tag dropdown options and pick the right variant UI without an extra
  // round-trip.
  const [modelOrganismMap, setModelOrganismMap] = useState({});
  // bundleName → {region: [start, end], ...}, for custom-organism uploads
  // whose protein_regions.csv was shipped in the bundle.
  const [modelProteinRegionsMap, setModelProteinRegionsMap] = useState({});
  // Per-subtype influenza protein regions ({influenza_h1n1: {...}, ...}).
  // Filled from GET /api/organisms/; used to swap the protein-region dropdown
  // contents when the user picks a different HA variant.
  const [influenzaProteinRegions, setInfluenzaProteinRegions] = useState({});
  // COVID protein regions sourced from the backend organism registry,
  // overriding the static src/data/proteinRegions.js when available.
  const [covidProteinRegions, setCovidProteinRegions] = useState(null);
  const loading = isLoading || reduxLoading;

  // Derive the selected model's organism for picker dispatch.
  const selectedOrganism = (() => {
    if (!selectedModel) return null;
    const key = selectedModel.startsWith("uploaded:")
      ? selectedModel.replace("uploaded:", "")
      : selectedModel;
    // Static (data/modelList.js) models all target COVID by convention.
    return modelOrganismMap[key] || "covid";
  })();

  // Resolve the protein-region catalog that should populate the dropdown for
  // the current selection. Three sources:
  //   - covid: backend's organisms/covid/protein_regions.csv (fallback: static
  //     src/data/proteinRegions.js so the picker still renders if /api/
  //     organisms/ hasn't loaded yet)
  //   - influenza: depends on the picked variant's subtype — swaps when the
  //     user changes variant
  //   - custom: ships with the bundle's protein_regions.csv (may be empty)
  const activeProteinRegions = (() => {
    if (selectedOrganism === "covid") {
      return covidProteinRegions || proteinRegions;
    }
    if (selectedOrganism === "influenza") {
      const variantEntry = influenzaVariants.find((v) => v.name === _variant);
      const subtypeOrg = variantEntry?.organism;
      return subtypeOrg ? (influenzaProteinRegions[subtypeOrg] || {}) : {};
    }
    if (selectedOrganism === "custom") {
      const key = selectedModel.startsWith("uploaded:")
        ? selectedModel.replace("uploaded:", "")
        : selectedModel;
      return modelProteinRegionsMap[key] || {};
    }
    return {};
  })();

  // Some organisms / models don't take an elapsed-day input. Only COVID's
  // built-in extractor consumes it — for influenza and custom uploads we
  // hide the field and send 0 server-side.
  const showElapsedDay = selectedOrganism === "covid";

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
        const richModels = data.models || [];
        // Build modelName → organism map so the variant picker can dispatch
        // the moment the user changes their selection.
        const orgMap = {};
        const regionsMap = {};
        richModels.forEach((m) => {
          orgMap[m.name] = m.organism;
          if (m.protein_regions) {
            regionsMap[m.name] = m.protein_regions;
          }
        });
        setModelOrganismMap(orgMap);
        setModelProteinRegionsMap(regionsMap);

        const formattedAPI = richModels
          .filter((m) => m.uploaded)
          .map((m) => ({
            label: `${m.name} (Uploaded)`,
            value: `uploaded:${m.name}`,
            type: "uploaded",
            organism: m.organism,
          }));

        // Tag static models too — they're all COVID by convention.
        const formattedStaticTagged = formattedStatic.map((m) => ({
          ...m,
          organism: "covid",
        }));

        const allModels = [...formattedStaticTagged, ...formattedAPI];
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

  // Fetch the consolidated influenza variant catalog (9 strains across H1N1
  // / H3N2 / H5N1). Only the prediction-screen variant picker uses it.
  useEffect(() => {
    fetch(`${API_URL}/api/organisms/`)
      .then((response) => (response.ok ? response.json() : { influenza_variants: [] }))
      .then((data) => {
        setInfluenzaVariants(data.influenza_variants || []);
        setInfluenzaProteinRegions(data.protein_regions_by_subtype || {});
        if (data.covid_protein_regions) {
          setCovidProteinRegions(data.covid_protein_regions);
        }
      })
      .catch(() => {
        setInfluenzaVariants([]);
        setInfluenzaProteinRegions({});
      });
  }, []);

  // Reset the variant when the selected model changes — variant catalog and
  // nodeId are organism-scoped, so a stale selection from a previous model
  // would silently corrupt the next prediction.
  useEffect(() => {
    setVariant("");
  }, [selectedModel]);

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

    if (!selectedModel) {
      alert("Please select a Model.");
      return;
    }

    // Organism-aware variant validation: each organism has its own picker
    // (COVID = nodeId from mutations.txt, Influenza = HA strain from catalog,
    // Custom = none, bundle ships its own genome).
    if (selectedOrganism === "covid" && !_nodeId) {
      alert("Please pick a COVID Variant ID.");
      return;
    }
    if (selectedOrganism === "influenza" && !_variant) {
      alert("Please pick an Influenza HA Strain (one of the 9 cataloged).");
      return;
    }

    // Elapsed days is only a COVID-extractor feature. Influenza and custom
    // organisms skip the field entirely (server treats it as 0).
    if (showElapsedDay && (!_elapsedDay || _elapsedDay === "0")) {
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
      elapsedDay: showElapsedDay ? (Number(_elapsedDay) || 60) : 0,
      selectedModel: selectedModel,
      selectedProteinRegion: selectedProteinRegion || null,
      isNewUpload: false,
      customParameters: customParamsObj,
      // Influenza HA strain override — only sent for influenza models. The
      // backend ignores it for COVID (uses nodeId) and custom (uses the
      // bundle's own genome).
      variant: selectedOrganism === "influenza" ? _variant : undefined,
    };

    try {
      if (onSubmit) {
        await onSubmit(params);
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
                formatOptionLabel={(option) => {
                  // Pick a badge color per organism so the dropdown signals
                  // at a glance what kind of variant picker comes next.
                  const orgBadgeStyle = {
                    covid: "bg-blue-100 text-blue-700",
                    influenza: "bg-purple-100 text-purple-700",
                    custom: "bg-orange-100 text-orange-700",
                  }[option.organism] || "bg-gray-100 text-gray-700";
                  const orgLabel = {
                    covid: "COVID",
                    influenza: "Influenza",
                    custom: "Custom",
                  }[option.organism] || (option.organism || "?");
                  return (
                    <div className="flex items-center justify-between">
                      <span className="text-sm">{option.label}</span>
                      <span className="flex items-center gap-1 ml-2">
                        {option.type === "uploaded" && (
                          <span className="text-[10px] bg-green-100 text-green-700 font-bold px-1.5 py-0.5 rounded-full">
                            Uploaded
                          </span>
                        )}
                        <span className={`text-[10px] font-bold px-1.5 py-0.5 rounded-full ${orgBadgeStyle}`}>
                          {orgLabel}
                        </span>
                      </span>
                    </div>
                  );
                }}
              />
            </div>

            {/* Variant picker — morphs by the selected model's organism. */}
            {selectedOrganism === "covid" && (
              <div>
                <label className="text-xs font-bold text-gray-500 uppercase mb-1 block ml-1">
                  COVID Variant ID *
                </label>
                <DropDown items={nodes} setNodeId={setNodeId} />
              </div>
            )}

            {selectedOrganism === "influenza" && (
              <div>
                <label className="text-xs font-bold text-gray-500 uppercase mb-1 block ml-1">
                  Influenza HA Strain *
                </label>
                <Select
                  options={influenzaVariants.map((v) => ({
                    label: `${v.display_name}${v.year ? ` — ${v.year}` : ""}${
                      v.accession ? ` (${v.accession})` : ""
                    }`,
                    value: v.name,
                    subtype: v.subtype,
                    isReference: String(v.is_reference).toLowerCase() === "true",
                  }))}
                  onChange={(opt) => setVariant(opt ? opt.value : "")}
                  value={
                    influenzaVariants
                      .map((v) => ({
                        label: `${v.display_name}${v.year ? ` — ${v.year}` : ""}${
                          v.accession ? ` (${v.accession})` : ""
                        }`,
                        value: v.name,
                      }))
                      .find((o) => o.value === _variant) || null
                  }
                  styles={customStyles}
                  placeholder="Pick an HA strain (PR/8/34, Cal/07, cattle/Texas, ...)"
                  formatOptionLabel={(option) => (
                    <div className="flex items-center justify-between">
                      <span className="text-sm">{option.label}</span>
                      <span className="flex items-center gap-1 ml-2">
                        {option.subtype && (
                          <span className="text-[10px] bg-indigo-100 text-indigo-700 font-bold px-1.5 py-0.5 rounded-full">
                            {option.subtype}
                          </span>
                        )}
                        {option.isReference && (
                          <span className="text-[10px] bg-amber-100 text-amber-700 font-bold px-1.5 py-0.5 rounded-full">
                            Reference
                          </span>
                        )}
                      </span>
                    </div>
                  )}
                />
                <p className="text-xs text-gray-400 mt-1 ml-1">
                  9 cataloged HA strains across H1N1 / H3N2 / H5N1, including the 2024 US dairy
                  cattle isolate (clade 2.3.4.4b).
                </p>
              </div>
            )}

            {selectedOrganism === "custom" && (
              <div className="rounded-xl bg-orange-50 border border-orange-100 p-3 text-xs text-orange-800">
                This bundle ships its own genome and protein regions. No variant
                selection is needed — predictions run on the genome the bundle uploaded.
              </div>
            )}

            {/* Elapsed Days — COVID-only feature */}
            {showElapsedDay && (
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
            )}

            {/* Region — populated from the selected model/variant's catalog */}
            <div>
              <label className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block ml-1">
                Protein Region
              </label>
              <Select
                options={[
                  { label: "Whole Genome", value: "" },
                  ...Object.keys(activeProteinRegions).map((pr) => ({
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
