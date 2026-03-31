import React, { useState, useEffect, useCallback, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import Chart from "chart.js/auto";
import { useRef } from "react";
import { modelList as staticModels } from "../data/modelList";
import {
  MdCompareArrows,
  MdArrowBack,
  MdPlayArrow,
  MdCheckCircle,
  MdError,
  MdTimer,
  MdBarChart,
  MdScience,
  MdExpandMore,
  MdExpandLess,
  MdInfo,
} from "react-icons/md";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

// Color palette for models
const MODEL_COLORS = [
  { bg: "rgba(59, 130, 246, 0.7)", border: "#3B82F6", light: "#EFF6FF" },
  { bg: "rgba(239, 68, 68, 0.7)", border: "#EF4444", light: "#FEF2F2" },
  { bg: "rgba(16, 185, 129, 0.7)", border: "#10B981", light: "#ECFDF5" },
  { bg: "rgba(245, 158, 11, 0.7)", border: "#F59E0B", light: "#FFFBEB" },
  { bg: "rgba(139, 92, 246, 0.7)", border: "#8B5CF6", light: "#F5F3FF" },
];

// ============================================
// COLLAPSIBLE SECTION
// ============================================
const Section = ({ title, icon, children, defaultOpen = true, badge }) => {
  const [open, setOpen] = useState(defaultOpen);
  return (
    <div className="bg-white rounded-xl shadow-sm border border-gray-100 overflow-hidden">
      <button
        onClick={() => setOpen(!open)}
        className="w-full flex items-center justify-between px-5 py-3 hover:bg-gray-50 transition-colors"
      >
        <div className="flex items-center gap-2">
          {icon}
          <h3 className="font-bold text-gray-800 text-sm">{title}</h3>
          {badge && (
            <span className="text-xs bg-blue-100 text-blue-700 px-2 py-0.5 rounded-full font-semibold">
              {badge}
            </span>
          )}
        </div>
        {open ? <MdExpandLess className="text-gray-400" /> : <MdExpandMore className="text-gray-400" />}
      </button>
      {open && <div className="px-5 pb-4">{children}</div>}
    </div>
  );
};

// ============================================
// CALIBRATION CHART
// ============================================
const CalibrationChart = ({ models, results }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current || !results) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const datasets = [];
    models.forEach((modelName, idx) => {
      const modelData = results.models?.[modelName];
      const cal = modelData?.metrics?.calibration_curve;
      if (!cal) return;
      const color = MODEL_COLORS[idx % MODEL_COLORS.length];
      datasets.push({
        label: modelName,
        data: cal.predicted.map((p, i) => ({ x: p, y: cal.actual[i] })),
        borderColor: color.border,
        backgroundColor: color.bg,
        pointRadius: 5,
        showLine: true,
        tension: 0.1,
      });
    });
    datasets.push({
      label: "Perfect",
      data: [{ x: 0, y: 0 }, { x: 1, y: 1 }],
      borderColor: "#9CA3AF",
      borderDash: [5, 5],
      pointRadius: 0,
      showLine: true,
    });

    chartInstance.current = new Chart(chartRef.current, {
      type: "scatter",
      data: { datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          x: { title: { display: true, text: "Predicted Probability" }, min: 0, max: 1 },
          y: { title: { display: true, text: "Actual Fraction" }, min: 0, max: 1 },
        },
        plugins: {
          legend: { position: "bottom", labels: { font: { size: 11 } } },
        },
      },
    });
    return () => { if (chartInstance.current) chartInstance.current.destroy(); };
  }, [results, models]);

  return <canvas ref={chartRef} />;
};

// ============================================
// PER-PROTEIN BAR CHART
// ============================================
const PerProteinChart = ({ models, results, metric = "mean_prediction" }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current || !results) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const firstModel = models.find((m) => results.models?.[m]?.per_protein);
    if (!firstModel) return;
    const proteinNames = Object.keys(results.models[firstModel].per_protein);

    const datasets = models.map((modelName, idx) => {
      const perProtein = results.models?.[modelName]?.per_protein || {};
      const color = MODEL_COLORS[idx % MODEL_COLORS.length];
      return {
        label: modelName,
        data: proteinNames.map((p) => perProtein[p]?.[metric] ?? 0),
        backgroundColor: color.bg,
        borderColor: color.border,
        borderWidth: 1,
      };
    });

    chartInstance.current = new Chart(chartRef.current, {
      type: "bar",
      data: { labels: proteinNames, datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          x: { ticks: { font: { size: 10 }, maxRotation: 45 } },
          y: {
            beginAtZero: true,
            title: {
              display: true,
              text: metric === "mean_prediction" ? "Mean Probability" :
                    metric === "auroc" ? "AUROC" : metric,
              font: { size: 11 },
            },
          },
        },
        plugins: {
          legend: { position: "bottom", labels: { font: { size: 11 } } },
        },
      },
    });
    return () => { if (chartInstance.current) chartInstance.current.destroy(); };
  }, [results, models, metric]);

  return <canvas ref={chartRef} />;
};

// ============================================
// PREDICTION OVERLAY CHART
// ============================================
const PredictionOverlayChart = ({ models, results, useLogScale }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current || !results?.models) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const datasets = [];
    models.forEach((modelName, idx) => {
      const curve = results.models?.[modelName]?.prediction_curve;
      if (!curve) return;
      const color = MODEL_COLORS[idx % MODEL_COLORS.length];
      datasets.push({
        label: modelName,
        data: curve.positions.map((pos, i) => ({ x: pos, y: curve.values[i] })),
        borderColor: color.border,
        backgroundColor: "transparent",
        borderWidth: 1.5,
        pointRadius: 0,
        showLine: true,
        tension: 0.2,
        order: 2,
      });
    });

    const mutPositions = results.parameters?.mutation_positions || [];
    if (mutPositions.length > 0) {
      datasets.push({
        label: "Actual Mutations",
        data: mutPositions.map((pos) => ({ x: pos, y: 1.0 })),
        borderColor: "#DC2626",
        backgroundColor: "#DC262680",
        pointRadius: 5,
        pointStyle: "triangle",
        showLine: false,
        order: 1,
      });
    }

    chartInstance.current = new Chart(chartRef.current, {
      type: "scatter",
      data: { datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false,
        interaction: { mode: "index", intersect: false },
        scales: {
          x: {
            title: { display: true, text: "Genome Position", font: { size: 11 } },
            min: 0,
            max: results.parameters?.ground_truth_total || 29903,
          },
          y: {
            type: useLogScale ? "logarithmic" : "linear",
            title: {
              display: true,
              text: useLogScale ? "Mutation Probability (log)" : "Mutation Probability",
              font: { size: 11 },
            },
            ...(useLogScale ? { min: 0.0001 } : { min: 0, max: 1 }),
            ticks: {
              font: { size: 10 },
              ...(useLogScale ? {
                callback: (val) => {
                  if ([0.0001, 0.001, 0.01, 0.1, 1].includes(val)) return val;
                  return null;
                }
              } : {}),
            },
          },
        },
        plugins: {
          legend: { position: "bottom", labels: { font: { size: 11 } } },
          tooltip: {
            backgroundColor: "rgba(0,0,0,0.85)",
            callbacks: {
              label: (ctx) => `${ctx.dataset.label}: pos ${Math.round(ctx.parsed.x)}, prob ${ctx.parsed.y.toFixed(4)}`,
            },
          },
        },
      },
    });

    return () => { if (chartInstance.current) chartInstance.current.destroy(); };
  }, [results, models, useLogScale]);

  return <canvas ref={chartRef} />;
};

// ============================================
// MAIN BENCHMARK DASHBOARD
// ============================================
const BenchmarkDashboard = () => {
  const navigate = useNavigate();

  // State
  const [availableModels, setAvailableModels] = useState([]);
  const [selectedModels, setSelectedModels] = useState([]);
  const [nodeId, setNodeId] = useState("EGY/CCHE57357_Wave_3_A029/2021|MZ380261.1|2021-05-11");
  const [elapsedDay, setElapsedDay] = useState("110");
  const [selectedRegion, setSelectedRegion] = useState("");
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState("");
  const [pastBenchmarks, setPastBenchmarks] = useState([]);
  const [showPast, setShowPast] = useState(false);
  const [perProteinMetric, setPerProteinMetric] = useState("mean_prediction");
  const [benchmarkMode, setBenchmarkMode] = useState("single");
  const [availableDatasets, setAvailableDatasets] = useState({});
  const [selectedDataset, setSelectedDataset] = useState("");
  const [seed, setSeed] = useState("42");
  const [useLogScale, setUseLogScale] = useState(false);

  // Delete uploaded model
  const handleDeleteModel = async (modelValue, modelName) => {
    const name = modelValue.replace('uploaded:', '');
    if (!window.confirm(`Are you sure you want to delete "${name}"? This cannot be undone.`)) return;
    try {
      const res = await fetch(`${API_URL}/api/models/delete/?model_name=${encodeURIComponent(name)}`, {
        method: "DELETE",
      });
      if (res.ok) {
        setAvailableModels((prev) => prev.filter((m) => m.value !== modelValue));
        setSelectedModels((prev) => prev.filter((m) => m !== modelValue));
      } else {
        const err = await res.json().catch(() => ({}));
        setError(err.error || "Failed to delete model");
      }
    } catch (err) {
      setError(`Delete error: ${err.message}`);
    }
  };

  // Fetch available models
  useEffect(() => {
    const fetchModels = async () => {
      try {
        const res = await fetch(`${API_URL}/api/models/`);
        if (res.ok) {
          const data = await res.json();
          const models = [];
          staticModels.forEach((m) => {
            models.push({ name: m.name, value: m.path, type: "server" });
          });
          const uploaded = data.available_models || [];
          const serverNames = staticModels.map((m) => m.path);
          uploaded.forEach((m) => {
            const name = typeof m === "string" ? m : m.folder_name || m.name;
            if (name && !serverNames.includes(name)) {
              models.push({ name: `${name} (Uploaded)`, value: `uploaded:${name}`, type: "uploaded" });
            }
          });
          setAvailableModels(models);
        }
      } catch (err) {
        console.warn("Could not fetch models:", err);
      }
    };
    fetchModels();
  }, []);

  // Fetch available datasets
  useEffect(() => {
    const fetchDatasets = async () => {
      try {
        const res = await fetch(`${API_URL}/api/benchmark/datasets/`);
        if (res.ok) {
          const data = await res.json();
          setAvailableDatasets(data.datasets || {});
        }
      } catch (err) {
        console.warn("Could not fetch datasets:", err);
      }
    };
    fetchDatasets();
  }, []);

  // Fetch past benchmarks
  useEffect(() => {
    const fetchPast = async () => {
      try {
        const res = await fetch(`${API_URL}/api/benchmark/list/`);
        if (res.ok) {
          const data = await res.json();
          setPastBenchmarks(data.benchmarks || []);
        }
      } catch (err) {
        console.warn("Could not fetch past benchmarks:", err);
      }
    };
    fetchPast();
  }, [results]);

  const toggleModel = (modelValue) => {
    setSelectedModels((prev) =>
      prev.includes(modelValue)
        ? prev.filter((m) => m !== modelValue)
        : [...prev, modelValue]
    );
  };

  // Run benchmark
  const handleRunBenchmark = async () => {
    if (selectedModels.length < 1) {
      setError("Select at least 1 model to benchmark.");
      return;
    }
    setLoading(true);
    setError("");
    setResults(null);

    try {
      let res;
      if (benchmarkMode === "dataset") {
        if (!selectedDataset) {
          setError("Please select a dataset.");
          setLoading(false);
          return;
        }
        res = await fetch(`${API_URL}/api/benchmark/run-dataset/`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            models: selectedModels,
            datasetId: selectedDataset,
            selectedProteinRegion: selectedRegion || null,
            seed: Number(seed) || 42,
          }),
        });
      } else {
        if (!nodeId) {
          setError("Please enter a Node ID.");
          setLoading(false);
          return;
        }
        res = await fetch(`${API_URL}/api/benchmark/run/`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            models: selectedModels,
            nodeId,
            elapsedDay: Number(elapsedDay) || 60,
            selectedProteinRegion: selectedRegion || null,
            seed: Number(seed) || 42,
          }),
        });
      }

      if (res.ok) {
        const data = await res.json();
        setResults(data);
      } else {
        const errData = await res.json().catch(() => ({}));
        setError(errData.error || `Benchmark failed (HTTP ${res.status})`);
      }
    } catch (err) {
      setError(`Network error: ${err.message}`);
    } finally {
      setLoading(false);
    }
  };

  // Export
  const handleExport = async (format) => {
    if (!results?.benchmark_id) return;
    try {
      const res = await fetch(`${API_URL}/api/benchmark/export/`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ id: results.benchmark_id, format }),
      });
      if (!res.ok) {
        const errText = await res.text();
        setError(`Export failed: ${errText}`);
        return;
      }
      const blob = await res.blob();
      const ext = format === 'json' ? 'json' : format === 'csv' ? 'csv' : 'html';
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${results.benchmark_id}.${ext}`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      window.URL.revokeObjectURL(url);
    } catch (err) {
      setError(`Export error: ${err.message}`);
    }
  };

  // Load past benchmark
  const loadPastBenchmark = async (benchmarkId) => {
    try {
      const res = await fetch(`${API_URL}/api/benchmark/results/?id=${benchmarkId}`);
      if (res.ok) {
        const data = await res.json();
        setResults(data);
        setSelectedModels(Object.keys(data.models || {}));
      }
    } catch (err) {
      setError(`Failed to load benchmark: ${err.message}`);
    }
  };

  // Successful models
  const successfulModels = useMemo(() => {
    if (!results) return [];
    if (results.type === "dataset") {
      return Object.keys(results.aggregated?.models || {});
    }
    return Object.entries(results.models || {})
      .filter(([_, v]) => v.status === "success")
      .map(([k]) => k);
  }, [results]);

  const proteinRegionOptions = [
    "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"
  ];

  // Helper: get summary for best model highlighting
  const getBestIdx = (metricKey, higher) => {
    const values = successfulModels.map((m) => results.models[m]?.metrics?.[metricKey]);
    const valid = values.filter((v) => v !== null && v !== undefined);
    if (higher === null || valid.length === 0) return -1;
    const best = higher ? Math.max(...valid) : Math.min(...valid);
    return values.indexOf(best);
  };

  return (
    <div className="min-h-screen bg-gray-50 py-6 px-4">
      <div className="max-w-7xl mx-auto">
        {/* Header */}
        <div className="mb-6">
          <button
            onClick={() => navigate("/")}
            className="flex items-center gap-2 text-gray-600 hover:text-blue-600 mb-3 transition-colors"
          >
            <MdArrowBack /> Back to Home
          </button>
          <h1 className="text-2xl font-bold text-gray-800 flex items-center gap-3">
            <MdCompareArrows className="text-blue-600" />
            Model Benchmark & Comparison
          </h1>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
          {/* ============================================ */}
          {/* LEFT PANEL - Configuration */}
          {/* ============================================ */}
          <div className="lg:col-span-4 space-y-4">
            {/* Model Selection */}
            <div className="bg-white p-5 rounded-xl shadow-sm border border-gray-100">
              <h3 className="font-bold text-gray-800 mb-3 flex items-center gap-2 text-sm">
                <MdScience className="text-blue-500" /> Select Models
              </h3>
              <div className="space-y-1.5 max-h-52 overflow-y-auto">
                {availableModels.map((model, idx) => (
                  <label
                    key={idx}
                    className={`flex items-center gap-2.5 p-2 rounded-lg cursor-pointer transition-colors text-sm ${
                      selectedModels.includes(model.value)
                        ? "bg-blue-50 border border-blue-200"
                        : "hover:bg-gray-50 border border-transparent"
                    }`}
                  >
                    <input
                      type="checkbox"
                      checked={selectedModels.includes(model.value)}
                      onChange={() => toggleModel(model.value)}
                      className="rounded text-blue-600"
                    />
                    <span className="text-gray-700 flex-1 truncate">{model.name}</span>
                    <span className={`text-xs px-1.5 py-0.5 rounded-full ${
                      model.type === "uploaded" ? "bg-green-100 text-green-700" : "bg-gray-100 text-gray-500"
                    }`}>
                      {model.type}
                    </span>
                    {model.type === "uploaded" && (
                      <button
                        onClick={(e) => { e.preventDefault(); e.stopPropagation(); handleDeleteModel(model.value, model.name); }}
                        className="text-red-400 hover:text-red-600 transition-colors"
                        title="Delete model"
                      >
                        <MdError size={14} />
                      </button>
                    )}
                  </label>
                ))}
                {availableModels.length === 0 && (
                  <p className="text-sm text-gray-400 text-center py-4">Loading models...</p>
                )}
              </div>
              <p className="text-xs text-gray-400 mt-2">
                {selectedModels.length} model(s) selected
              </p>
            </div>

            {/* Parameters */}
            <div className="bg-white p-5 rounded-xl shadow-sm border border-gray-100">
              <h3 className="font-bold text-gray-800 mb-3 text-sm">Parameters</h3>

              {/* Benchmark Mode Toggle */}
              <div className="flex gap-2 mb-4">
                <button
                  onClick={() => setBenchmarkMode("single")}
                  className={`flex-1 py-2 px-3 rounded-lg text-xs font-semibold transition-colors ${
                    benchmarkMode === "single" ? "bg-blue-600 text-white" : "bg-gray-100 text-gray-600 hover:bg-gray-200"
                  }`}
                >
                  Single Variant
                </button>
                <button
                  onClick={() => setBenchmarkMode("dataset")}
                  className={`flex-1 py-2 px-3 rounded-lg text-xs font-semibold transition-colors ${
                    benchmarkMode === "dataset" ? "bg-blue-600 text-white" : "bg-gray-100 text-gray-600 hover:bg-gray-200"
                  }`}
                >
                  Dataset (Multi)
                </button>
              </div>

              <div className="space-y-3">
                {benchmarkMode === "single" ? (
                  <>
                    <div>
                      <label className="text-xs font-bold text-gray-500 uppercase mb-1 block">Node ID</label>
                      <input
                        type="text"
                        value={nodeId}
                        onChange={(e) => setNodeId(e.target.value)}
                        className="w-full border rounded-lg px-3 py-2 text-sm"
                        placeholder="Variant Node ID"
                      />
                    </div>
                    <div>
                      <label className="text-xs font-bold text-gray-500 uppercase mb-1 block">Elapsed Days</label>
                      <input
                        type="number"
                        value={elapsedDay}
                        onChange={(e) => setElapsedDay(e.target.value)}
                        className="w-full border rounded-lg px-3 py-2 text-sm"
                        min={1}
                      />
                    </div>
                  </>
                ) : (
                  <div>
                    <label className="text-xs font-bold text-gray-500 uppercase mb-1 block">Dataset</label>
                    <select
                      value={selectedDataset}
                      onChange={(e) => setSelectedDataset(e.target.value)}
                      className="w-full border rounded-lg px-3 py-2 text-sm"
                    >
                      <option value="">Select a dataset...</option>
                      {Object.entries(availableDatasets).map(([key, ds]) => (
                        <option key={key} value={key}>
                          {ds.name} ({ds.num_variants} variant{ds.num_variants !== 1 ? "s" : ""})
                        </option>
                      ))}
                    </select>
                  </div>
                )}

                <div className="grid grid-cols-2 gap-3">
                  <div>
                    <label className="text-xs font-bold text-gray-500 uppercase mb-1 block">Region</label>
                    <select
                      value={selectedRegion}
                      onChange={(e) => setSelectedRegion(e.target.value)}
                      className="w-full border rounded-lg px-3 py-2 text-sm"
                    >
                      <option value="">Whole Genome</option>
                      {proteinRegionOptions.map((pr) => (
                        <option key={pr} value={pr}>{pr}</option>
                      ))}
                    </select>
                  </div>
                  <div>
                    <label className="text-xs font-bold text-gray-500 uppercase mb-1 block">Seed</label>
                    <input
                      type="number"
                      value={seed}
                      onChange={(e) => setSeed(e.target.value)}
                      className="w-full border rounded-lg px-3 py-2 text-sm"
                      min={0}
                    />
                  </div>
                </div>
              </div>

              {/* Run Button */}
              <button
                onClick={handleRunBenchmark}
                disabled={loading || selectedModels.length < 1}
                className="w-full mt-4 py-3 bg-blue-600 hover:bg-blue-700 disabled:bg-blue-300 text-white font-semibold rounded-xl flex items-center justify-center gap-2 transition-colors text-sm"
              >
                {loading ? (
                  <>
                    <div className="animate-spin rounded-full h-4 w-4 border-2 border-white border-t-transparent" />
                    Running...
                  </>
                ) : (
                  <>
                    <MdPlayArrow className="text-lg" />
                    Run Benchmark ({selectedModels.length} model{selectedModels.length !== 1 ? "s" : ""})
                  </>
                )}
              </button>

              {/* Visual Compare Button */}
              {benchmarkMode === "single" && selectedModels.length >= 2 && nodeId && (
                <button
                  onClick={() => {
                    const params = new URLSearchParams({
                      models: selectedModels.join(","),
                      nodeId,
                      elapsedDay: elapsedDay || "60",
                      proteinRegion: selectedRegion || "",
                    });
                    navigate(`/compare?${params.toString()}`);
                  }}
                  className="w-full mt-2 py-2.5 bg-purple-600 hover:bg-purple-700 text-white font-semibold rounded-xl flex items-center justify-center gap-2 transition-colors text-sm"
                >
                  <MdCompareArrows className="text-lg" />
                  Visual Compare
                </button>
              )}

              {error && (
                <div className="mt-3 p-3 bg-red-50 border border-red-200 rounded-lg text-red-700 text-sm">
                  {error}
                </div>
              )}
            </div>

            {/* Past Benchmarks */}
            <div className="bg-white p-4 rounded-xl shadow-sm border border-gray-100">
              <button
                onClick={() => setShowPast(!showPast)}
                className="w-full flex justify-between items-center"
              >
                <h3 className="font-bold text-gray-800 text-sm">Past Benchmarks</h3>
                {showPast ? <MdExpandLess /> : <MdExpandMore />}
              </button>
              {showPast && (
                <div className="mt-3 space-y-1.5 max-h-36 overflow-y-auto">
                  {pastBenchmarks.map((b, idx) => (
                    <button
                      key={idx}
                      onClick={() => loadPastBenchmark(b.benchmark_id)}
                      className="w-full text-left p-2 hover:bg-blue-50 rounded-lg text-sm transition-colors"
                    >
                      <p className="font-medium text-gray-700 text-xs">{b.model_names?.join(" vs ")}</p>
                      <p className="text-xs text-gray-400">{b.timestamp}</p>
                    </button>
                  ))}
                  {pastBenchmarks.length === 0 && (
                    <p className="text-sm text-gray-400 text-center py-2">No benchmarks yet</p>
                  )}
                </div>
              )}
            </div>
          </div>

          {/* ============================================ */}
          {/* RIGHT PANEL - Results */}
          {/* ============================================ */}
          <div className="lg:col-span-8 space-y-4">
            {/* Empty State */}
            {!results && !loading && (
              <div className="bg-white rounded-2xl shadow-sm border border-gray-100 h-80 flex items-center justify-center">
                <div className="text-center opacity-60">
                  <MdBarChart size={48} className="mx-auto text-blue-300 mb-3" />
                  <h3 className="text-lg font-bold text-gray-700 mb-1">Ready to Compare</h3>
                  <p className="text-gray-500 text-sm max-w-sm">
                    Select models and run a benchmark to see results.
                  </p>
                </div>
              </div>
            )}

            {/* Loading */}
            {loading && (
              <div className="bg-white rounded-2xl shadow-sm border border-gray-100 h-80 flex items-center justify-center">
                <div className="text-center">
                  <div className="animate-spin rounded-full h-10 w-10 border-3 border-blue-500 border-t-transparent mx-auto mb-4" />
                  <p className="text-lg font-semibold text-blue-600">Running Benchmark...</p>
                  <p className="text-sm text-gray-400 mt-1">This may take a few minutes</p>
                </div>
              </div>
            )}

            {/* Results */}
            {results && !loading && (
              <>
                {/* Benchmark Header - Compact */}
                <div className="bg-white px-5 py-3 rounded-xl shadow-sm border border-gray-100 flex justify-between items-center flex-wrap gap-2">
                  <div>
                    <h3 className="font-bold text-gray-800 text-sm">
                      {results.type === "dataset"
                        ? `Dataset: ${results.dataset?.name}`
                        : `Benchmark: ${results.benchmark_id?.substring(0, 20)}`
                      }
                    </h3>
                    <p className="text-xs text-gray-400">
                      {results.timestamp}
                      {results.type === "dataset"
                        ? ` | ${results.dataset?.num_variants} variants`
                        : ` | ${results.parameters?.num_mutations} mutations`
                      }
                    </p>
                  </div>
                  <div className="flex gap-2 items-center">
                    <div className="flex gap-1">
                      {["json", "csv", "html"].map((fmt) => (
                        <button key={fmt} onClick={() => handleExport(fmt)}
                          className="px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded font-mono transition-colors">
                          {fmt.toUpperCase()}
                        </button>
                      ))}
                    </div>
                    <div className="flex gap-1 ml-2">
                      {successfulModels.map((name, idx) => (
                        <span key={name} className="px-2 py-0.5 rounded-full text-xs font-bold text-white"
                              style={{ backgroundColor: MODEL_COLORS[idx % MODEL_COLORS.length].border }}>
                          {name}
                        </span>
                      ))}
                    </div>
                  </div>
                </div>

                {/* Key Metrics - Compact Cards */}
                <Section title="Key Metrics" icon={<MdBarChart className="text-blue-500" />} defaultOpen={true}>
                  <div className="overflow-x-auto">
                    <table className="w-full text-sm">
                      <thead>
                        <tr className="border-b border-gray-200">
                          <th className="text-left py-2 px-3 text-gray-500 uppercase text-xs">Metric</th>
                          {successfulModels.map((name, idx) => (
                            <th key={name} className="text-center py-2 px-3">
                              <span className="px-2 py-0.5 rounded text-xs font-bold text-white"
                                    style={{ backgroundColor: MODEL_COLORS[idx % MODEL_COLORS.length].border }}>
                                {name}
                              </span>
                            </th>
                          ))}
                        </tr>
                      </thead>
                      <tbody>
                        {[
                          { key: "auroc", label: "AUROC", higher: true },
                          { key: "auprc", label: "AUPRC", higher: true },
                          { key: "brier_score", label: "Brier Score", higher: false },
                          { key: "ece", label: "ECE", higher: false },
                          { key: "runtime_seconds", label: "Runtime (s)", higher: false },
                        ].map((metric) => {
                          const values = successfulModels.map(
                            (m) => results.models[m]?.metrics?.[metric.key]
                          );
                          const bestIdx = getBestIdx(metric.key, metric.higher);

                          return (
                            <tr key={metric.key} className="border-b border-gray-100 hover:bg-gray-50">
                              <td className="py-2 px-3 font-medium text-gray-700 text-xs">
                                {metric.label}
                                <span className="ml-1 text-gray-400">
                                  {metric.higher ? "↑" : "↓"}
                                </span>
                              </td>
                              {values.map((val, idx) => (
                                <td key={idx} className={`text-center py-2 px-3 font-mono text-xs ${
                                  idx === bestIdx ? "font-bold text-green-600" : "text-gray-600"
                                }`}>
                                  {val !== null && val !== undefined ? val.toFixed(4) : "—"}
                                  {idx === bestIdx && " ★"}
                                </td>
                              ))}
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                  </div>
                </Section>

                {/* Prediction Chart */}
                {successfulModels.some((m) => results.models?.[m]?.prediction_curve) && (
                  <Section
                    title="Prediction Comparison"
                    icon={<MdCompareArrows className="text-purple-500" />}
                    defaultOpen={true}
                  >
                    <div className="flex justify-end mb-2">
                      <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                        <input
                          type="checkbox"
                          checked={useLogScale}
                          onChange={(e) => setUseLogScale(e.target.checked)}
                          className="rounded text-blue-600"
                        />
                        Log Scale
                      </label>
                    </div>
                    <div style={{ height: "320px" }}>
                      <PredictionOverlayChart models={successfulModels} results={results} useLogScale={useLogScale} />
                    </div>
                    <p className="text-xs text-gray-400 mt-2 text-center">
                      Lines show predicted mutation probability. Red triangles mark actual mutations.
                    </p>
                  </Section>
                )}

                {/* Charts Row - Calibration & Per-Protein */}
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <Section title="Calibration" icon={<MdInfo className="text-green-500" />} defaultOpen={true}>
                    <div style={{ height: "280px" }}>
                      <CalibrationChart models={successfulModels} results={results} />
                    </div>
                  </Section>

                  <Section title="Per-Protein" icon={<MdBarChart className="text-amber-500" />} defaultOpen={true}>
                    <div className="flex justify-end mb-2">
                      <select
                        value={perProteinMetric}
                        onChange={(e) => setPerProteinMetric(e.target.value)}
                        className="text-xs border rounded px-2 py-1"
                      >
                        <option value="mean_prediction">Mean Prediction</option>
                        <option value="auroc">AUROC</option>
                        <option value="auprc">AUPRC</option>
                        <option value="brier_score">Brier Score</option>
                      </select>
                    </div>
                    <div style={{ height: "260px" }}>
                      <PerProteinChart
                        models={successfulModels}
                        results={results}
                        metric={perProteinMetric}
                      />
                    </div>
                  </Section>
                </div>

                {/* Model Agreement - Collapsed by default */}
                {results.model_agreement && Object.keys(results.model_agreement).length > 0 && (
                  <Section title="Model Agreement" icon={<MdInfo className="text-blue-500" />} defaultOpen={false} badge="Correlation">
                    <div className="overflow-x-auto">
                      <table className="w-full text-sm">
                        <thead>
                          <tr className="border-b">
                            <th className="text-left py-2 px-3"></th>
                            {successfulModels.map((m) => (
                              <th key={m} className="text-center py-2 px-3 text-xs text-gray-600">{m}</th>
                            ))}
                          </tr>
                        </thead>
                        <tbody>
                          {successfulModels.map((rowModel) => (
                            <tr key={rowModel} className="border-b border-gray-100">
                              <td className="py-2 px-3 font-medium text-gray-700 text-xs">{rowModel}</td>
                              {successfulModels.map((colModel) => {
                                const val = results.model_agreement?.[rowModel]?.[colModel];
                                const isSelf = rowModel === colModel;
                                return (
                                  <td key={colModel} className={`text-center py-2 px-3 font-mono text-xs ${
                                    isSelf ? "bg-gray-100 text-gray-400" :
                                    val > 0.8 ? "bg-green-50 text-green-700" :
                                    val > 0.5 ? "bg-yellow-50 text-yellow-700" :
                                    "bg-red-50 text-red-700"
                                  }`}>
                                    {val !== null ? val.toFixed(3) : "—"}
                                  </td>
                                );
                              })}
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </Section>
                )}

                {/* Per-Protein Detail - Collapsed by default */}
                {successfulModels.length > 0 && results.models[successfulModels[0]]?.per_protein && (
                  <Section title="Per-Protein Details" icon={<MdScience className="text-gray-500" />} defaultOpen={false} badge="Advanced">
                    <div className="overflow-x-auto">
                      <table className="w-full text-xs">
                        <thead>
                          <tr className="border-b">
                            <th className="text-left py-2 px-2">Region</th>
                            <th className="text-center py-2 px-2">Pos</th>
                            <th className="text-center py-2 px-2">Mut</th>
                            {successfulModels.map((m, idx) => (
                              <th key={m} className="text-center py-2 px-2" colSpan="3">
                                <span className="px-1 py-0.5 rounded text-white"
                                      style={{ backgroundColor: MODEL_COLORS[idx % MODEL_COLORS.length].border }}>
                                  {m}
                                </span>
                                <div className="flex justify-center gap-1 mt-1 text-gray-400 font-normal">
                                  <span>Mean</span><span>AUROC</span><span>Brier</span>
                                </div>
                              </th>
                            ))}
                          </tr>
                        </thead>
                        <tbody>
                          {Object.keys(results.models[successfulModels[0]]?.per_protein || {}).map((protein) => (
                            <tr key={protein} className="border-b border-gray-100 hover:bg-gray-50">
                              <td className="py-1.5 px-2 font-medium">{protein}</td>
                              <td className="text-center py-1.5 px-2">
                                {results.models[successfulModels[0]].per_protein[protein]?.num_positions}
                              </td>
                              <td className="text-center py-1.5 px-2">
                                {results.models[successfulModels[0]].per_protein[protein]?.num_mutations}
                              </td>
                              {successfulModels.map((m) => {
                                const pp = results.models[m]?.per_protein?.[protein];
                                return (
                                  <React.Fragment key={m}>
                                    <td className="text-center py-1.5 px-1 font-mono">{pp?.mean_prediction?.toFixed(4) ?? "—"}</td>
                                    <td className="text-center py-1.5 px-1 font-mono">{pp?.auroc?.toFixed(3) ?? "—"}</td>
                                    <td className="text-center py-1.5 px-1 font-mono">{pp?.brier_score?.toFixed(4) ?? "—"}</td>
                                  </React.Fragment>
                                );
                              })}
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </Section>
                )}
              </>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default BenchmarkDashboard;
