import React, { useState, useEffect, useCallback, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import Chart from "chart.js/auto";
import { useRef } from "react";
import { modelList as staticModels } from "../data/modelList";
import CompareModels from "./CompareModels";
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
  MdUploadFile,
  MdClose,
  MdOpenInFull,
  MdCloseFullscreen,
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
                    metric === "auroc" ? "AUROC" :
                    metric === "auprc" ? "AUPRC" :
                    metric === "precision" ? "Precision" :
                    metric === "recall" ? "Recall" :
                    metric === "f1_score" ? "F1 Score" :
                    metric === "mcc" ? "MCC" :
                    metric === "brier_score" ? "Brier Score" : metric,
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
const PredictionOverlayChart = ({ models, results, useLogScale, gtMinThreshold = 0, colorOffset = 0 }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current || !results?.models) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const datasets = [];
    models.forEach((modelName, idx) => {
      const curve = results.models?.[modelName]?.prediction_curve;
      if (!curve) return;
      const color = MODEL_COLORS[(idx + colorOffset) % MODEL_COLORS.length];
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
    const mutValues = results.parameters?.mutation_values || [];
    const gtType = results.parameters?.ground_truth_type || "binary";
    if (mutPositions.length > 0) {
      const isProbabilityGt = gtType === "probability";

      // Filter by threshold (probability GT only) and scale point radius by proportion
      // so that important mutations (proportion near 1.0) pop visually while
      // low-proportion "noise" fades out.
      const filtered = [];
      const radii = [];
      for (let i = 0; i < mutPositions.length; i++) {
        const val = isProbabilityGt ? (mutValues[i] ?? 1.0) : 1.0;
        if (isProbabilityGt && val < gtMinThreshold) continue;
        filtered.push({ x: mutPositions[i], y: val });
        // Stars are drawn as strokes (not filled), so they read smaller than a
        // filled circle of the same radius — scale them up a bit. Range ~4 (low
        // proportion) to ~8 (proportion ≈ 1) so high-confidence GT clearly pops.
        radii.push(isProbabilityGt ? Math.max(4, Math.min(8, 4 + val * 4)) : 6);
      }

      // Distinct dark color for GT to avoid clashing with model colors (blue/red/green/amber/purple)
      const GT_COLOR = "#111827"; // gray-900
      datasets.push({
        label: isProbabilityGt ? "Actual Proportion (GT)" : "Actual Mutations",
        data: filtered,
        borderColor: GT_COLOR,
        backgroundColor: GT_COLOR,
        borderWidth: 1.5,
        pointRadius: radii,
        pointHoverRadius: radii.map((r) => r + 2),
        pointStyle: "star",
        showLine: false,
        order: 0, // Draw GT on top so models don't fully hide it
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
            // Extend slightly beyond [0,1] in linear mode so GT dots at y≈1 and y≈0 aren't clipped.
            ...(useLogScale ? { min: 0.0001 } : { min: -0.03, max: 1.05 }),
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
  }, [results, models, useLogScale, gtMinThreshold]);

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
  // Prediction Comparison view mode: "overlay" (all models on one chart) or
  // "stacked" (one chart per model). Both reuse the same prediction_curve data
  // already in the benchmark result — no re-run required.
  const [predViewMode, setPredViewMode] = useState("overlay");
  // Top-level results view: "benchmark" (metrics/charts) or "visualization"
  // (the full multi-model genome visualization, embedded inline so the user can
  // switch back and forth without opening a new tab).
  const [resultsView, setResultsView] = useState("benchmark");
  // Snapshot of the parameters that produced the current results, so the inline
  // visualization matches what was benchmarked even if the form is edited after.
  const [vizParams, setVizParams] = useState(null);
  // When true, the left configuration panel is hidden and the visualization
  // expands to the full page width (a near-fullscreen view inside the page).
  const [vizFullWidth, setVizFullWidth] = useState(false);
  const [benchmarkMode, setBenchmarkMode] = useState("single");
  const [availableDatasets, setAvailableDatasets] = useState({});
  const [selectedDataset, setSelectedDataset] = useState("");
  const [seed, setSeed] = useState("42");
  const [useLogScale, setUseLogScale] = useState(false);

  // Ground truth source: "api" (cov-spectrum.org LAPIS, default) or "csv" (upload)
  const [groundTruthMode, setGroundTruthMode] = useState("api");
  const [groundTruthCsv, setGroundTruthCsv] = useState(null);
  const [csvPreview, setCsvPreview] = useState(null); // { rowCount, firstMutation }

  // cov-spectrum.org API parameters — variant (Pango lineage) + optional date range
  const [variantQuery, setVariantQuery] = useState("");
  const [variantResults, setVariantResults] = useState([]);
  const [variantSearchLoading, setVariantSearchLoading] = useState(false);
  const [variantSearchError, setVariantSearchError] = useState("");
  const [selectedVariant, setSelectedVariant] = useState(null); // { lineage, count }
  const [apiDateFrom, setApiDateFrom] = useState("");
  const [apiDateTo, setApiDateTo] = useState("");

  // Which metric's info popover is currently open (null = none)
  const [openMetricInfo, setOpenMetricInfo] = useState(null);

  // Minimum proportion for GT dots to be visible in the overlay chart
  // (cov-spectrum.org CSVs often contain tens of thousands of low-proportion
  // mutations that create visual noise — default hides anything below 5%).
  const [gtMinThreshold, setGtMinThreshold] = useState(0.05);

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

  // Handle CSV file selection (cov-spectrum.org format)
  const handleCsvUpload = async (file) => {
    if (!file) {
      setGroundTruthCsv(null);
      setCsvPreview(null);
      return;
    }
    setGroundTruthCsv(file);
    setError("");
    try {
      const text = await file.text();
      const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
      if (lines.length < 2) {
        setError("CSV is empty or has only a header.");
        setCsvPreview(null);
        return;
      }
      const header = lines[0].toLowerCase();
      if (!header.includes("mutation") || !header.includes("proportion")) {
        setError("CSV must contain 'mutation' and 'proportion' columns.");
        setCsvPreview(null);
        return;
      }
      const firstRow = lines[1].split(",");
      setCsvPreview({
        rowCount: lines.length - 1,
        firstMutation: firstRow[0],
        firstProportion: firstRow[1],
      });
    } catch (err) {
      setError(`Failed to read CSV: ${err.message}`);
      setCsvPreview(null);
    }
  };

  const clearCsv = () => {
    setGroundTruthCsv(null);
    setCsvPreview(null);
  };

  // Load ALL Pango lineages (cached 1h server-side) once; client-side filter afterwards.
  const loadAllVariants = useCallback(async () => {
    if (variantResults.length > 0 || variantSearchLoading) return;
    setVariantSearchLoading(true);
    setVariantSearchError("");
    try {
      const res = await fetch(`${API_URL}/api/benchmark/cov-spectrum/search-variants/?q=&limit=10000`);
      if (!res.ok) {
        const err = await res.json().catch(() => ({}));
        throw new Error(err.error || `Failed to load variants (HTTP ${res.status})`);
      }
      const data = await res.json();
      setVariantResults(data.variants || []);
    } catch (err) {
      setVariantSearchError(err.message);
    } finally {
      setVariantSearchLoading(false);
    }
  }, [variantResults.length, variantSearchLoading]);

  // Auto-load variants when user switches to API mode for the first time.
  useEffect(() => {
    if (groundTruthMode === "api" && variantResults.length === 0 && !variantSearchLoading && !variantSearchError) {
      loadAllVariants();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [groundTruthMode]);

  // Client-side filter for the variant combobox.
  const filteredVariants = useMemo(() => {
    const q = variantQuery.trim().toUpperCase();
    if (!q) return variantResults;
    return variantResults.filter((v) => v.lineage.toUpperCase().includes(q));
  }, [variantResults, variantQuery]);

  const clearSelectedVariant = () => {
    setSelectedVariant(null);
    setVariantQuery("");
  };

  // Quick-preset helper for date range. preset = "30d" | "3mo" | "1y" | "all".
  const applyDatePreset = (preset) => {
    if (preset === "all") {
      setApiDateFrom("");
      setApiDateTo("");
      return;
    }
    const today = new Date();
    const fromDate = new Date(today);
    if (preset === "30d") fromDate.setDate(today.getDate() - 30);
    else if (preset === "3mo") fromDate.setMonth(today.getMonth() - 3);
    else if (preset === "1y") fromDate.setFullYear(today.getFullYear() - 1);
    const toIso = today.toISOString().slice(0, 10);
    const fromIso = fromDate.toISOString().slice(0, 10);
    setApiDateFrom(fromIso);
    setApiDateTo(toIso);
  };

  // Which preset (if any) is currently active — for button highlighting.
  const activeDatePreset = useMemo(() => {
    if (!apiDateFrom && !apiDateTo) return "all";
    if (!apiDateFrom || !apiDateTo) return null;
    const today = new Date();
    const todayIso = today.toISOString().slice(0, 10);
    if (apiDateTo !== todayIso) return null;
    const fromDate = new Date(apiDateFrom);
    const diffDays = Math.round((today - fromDate) / (1000 * 60 * 60 * 24));
    if (diffDays === 30) return "30d";
    if (diffDays >= 89 && diffDays <= 92) return "3mo";
    if (diffDays >= 364 && diffDays <= 366) return "1y";
    return null;
  }, [apiDateFrom, apiDateTo]);

  // Run benchmark
  const handleRunBenchmark = async () => {
    if (selectedModels.length < 1) {
      setError("Select at least 1 model to benchmark.");
      return;
    }
    setLoading(true);
    setError("");
    setResults(null);
    setResultsView("benchmark"); // always land on the metrics view first
    // Snapshot the params that define this run so the inline visualization
    // stays consistent even if the form is edited afterwards.
    setVizParams({
      models: [...selectedModels],
      nodeId,
      elapsedDay: elapsedDay || 60,
      proteinRegion: selectedRegion || "",
    });

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

        if (groundTruthMode === "csv") {
          if (!groundTruthCsv) {
            setError("Please upload a mutation CSV file (mutation, proportion, count, jaccard).");
            setLoading(false);
            return;
          }
          // Multipart for file upload
          const form = new FormData();
          form.append("models", JSON.stringify(selectedModels));
          form.append("nodeId", nodeId);
          form.append("elapsedDay", String(Number(elapsedDay) || 60));
          form.append("selectedProteinRegion", selectedRegion || "");
          form.append("seed", String(Number(seed) || 42));
          form.append("groundTruthSource", "cov_spectrum_csv");
          form.append("groundTruthCsv", groundTruthCsv);
          res = await fetch(`${API_URL}/api/benchmark/run/`, {
            method: "POST",
            body: form,
          });
        } else {
          // API mode — variant required, dates optional (empty = all time)
          if (!selectedVariant) {
            setError("Please search and select a variant first.");
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
              groundTruthSource: "cov_spectrum_api",
              covSpectrumLineage: selectedVariant.lineage,
              covSpectrumDateFrom: apiDateFrom || null,
              covSpectrumDateTo: apiDateTo || null,
            }),
          });
        }
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

  // Protein region options: use per_protein keys from results when available
  // (organism-aware), otherwise fall back to the built-in COVID list.
  const COVID_PROTEIN_REGIONS = ["ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"];
  const proteinRegionOptions = useMemo(() => {
    if (!results) return COVID_PROTEIN_REGIONS;
    const models = results.type === "dataset"
      ? Object.values(results.per_variant?.[0]?.models || {})
      : Object.values(results.models || {});
    const firstWithRegions = models.find((m) => m.per_protein && Object.keys(m.per_protein).length > 0);
    if (firstWithRegions) return Object.keys(firstWithRegions.per_protein);
    return COVID_PROTEIN_REGIONS;
  }, [results]);

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
      <div className={`${vizFullWidth ? "max-w-[1800px]" : "max-w-7xl"} mx-auto transition-[max-width] duration-200`}>
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
          <div className={`${vizFullWidth ? "hidden" : "lg:col-span-4"} space-y-4`}>
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

              {/* Ground Truth (single variant mode only) */}
              {benchmarkMode === "single" && (
                <div className="mt-4 pt-4 border-t border-gray-100">
                  <label className="text-xs font-bold text-gray-500 uppercase mb-2 block">
                    Ground Truth
                  </label>

                  {/* Mode toggle: API fetch (default) vs CSV upload */}
                  <div className="flex gap-2 mb-3">
                    <button
                      type="button"
                      onClick={() => setGroundTruthMode("api")}
                      className={`flex-1 py-2 px-2 rounded-lg text-xs font-semibold transition-colors ${
                        groundTruthMode === "api"
                          ? "bg-blue-600 text-white"
                          : "bg-gray-100 text-gray-600 hover:bg-gray-200"
                      }`}
                    >
                      Fetch from cov-spectrum.org
                    </button>
                    <button
                      type="button"
                      onClick={() => setGroundTruthMode("csv")}
                      className={`flex-1 py-2 px-2 rounded-lg text-xs font-semibold transition-colors ${
                        groundTruthMode === "csv"
                          ? "bg-blue-600 text-white"
                          : "bg-gray-100 text-gray-600 hover:bg-gray-200"
                      }`}
                    >
                      Upload Mutation CSV
                    </button>
                  </div>

                  {groundTruthMode === "csv" && (
                    <div className="space-y-2">
                      {/* Format guide — what file to upload */}
                      <div className="p-3 bg-blue-50 border border-blue-200 rounded-lg text-xs text-gray-700 leading-snug">
                        <p className="font-semibold text-blue-900 mb-1.5">
                          Required CSV format
                        </p>
                        <p className="text-gray-700 mb-2">
                          A comma-separated file with one mutation per row. The first row must be a header
                          with these column names:
                        </p>
                        <div className="bg-white border border-blue-100 rounded p-2 font-mono text-[11px] overflow-x-auto whitespace-pre">
{`mutation,proportion,count,jaccard
A23403G,0.987,1542,0.95
C14408T,0.959,1498,0.92
T3-,0.401,623,0.55`}
                        </div>
                        <ul className="mt-2 space-y-0.5 text-[11px] text-gray-600 list-disc list-inside">
                          <li><code className="bg-white px-1 rounded">mutation</code> — e.g. <code className="bg-white px-1 rounded">A23403G</code> (substitution) or <code className="bg-white px-1 rounded">T3-</code> (deletion)</li>
                          <li><code className="bg-white px-1 rounded">proportion</code> — value in [0, 1] used as probability ground truth</li>
                          <li><code className="bg-white px-1 rounded">count</code>, <code className="bg-white px-1 rounded">jaccard</code> — currently informational only (can be empty)</li>
                        </ul>
                        <p className="mt-2 text-[11px] text-gray-500">
                          Tip: cov-spectrum.org's "Mutations" tab exports this exact format directly.
                        </p>
                      </div>

                      {/* Upload dropzone */}
                      <label
                        className={`flex items-center gap-2 p-3 border-2 border-dashed rounded-lg cursor-pointer transition-colors text-xs ${
                          groundTruthCsv
                            ? "border-green-300 bg-green-50"
                            : "border-gray-300 hover:border-blue-400 hover:bg-blue-50"
                        }`}
                      >
                        <MdUploadFile className={groundTruthCsv ? "text-green-600" : "text-gray-400"} size={20} />
                        <div className="flex-1 min-w-0">
                          {groundTruthCsv ? (
                            <>
                              <p className="font-semibold text-green-700 truncate">{groundTruthCsv.name}</p>
                              {csvPreview && (
                                <p className="text-gray-500">
                                  {csvPreview.rowCount} rows · e.g. {csvPreview.firstMutation} ({csvPreview.firstProportion})
                                </p>
                              )}
                            </>
                          ) : (
                            <>
                              <p className="font-semibold text-gray-700">Click to choose a CSV file</p>
                              <p className="text-gray-400">Format: see above</p>
                            </>
                          )}
                        </div>
                        {groundTruthCsv && (
                          <button
                            type="button"
                            onClick={(e) => { e.preventDefault(); e.stopPropagation(); clearCsv(); }}
                            className="text-gray-400 hover:text-red-500"
                            title="Remove"
                          >
                            <MdClose size={16} />
                          </button>
                        )}
                        <input
                          type="file"
                          accept=".csv,text/csv"
                          className="hidden"
                          onChange={(e) => handleCsvUpload(e.target.files?.[0] || null)}
                        />
                      </label>
                    </div>
                  )}

                  {groundTruthMode === "api" && (
                    <div className="space-y-3">
                      {/* Step 1: Variant combobox (all variants loaded once, filtered client-side) */}
                      <div>
                        <label className="text-xs font-semibold text-gray-500 mb-1 block">
                          Variant (Pango Lineage)
                          {variantResults.length > 0 && (
                            <span className="text-gray-400 font-normal ml-1">
                              — {variantResults.length.toLocaleString()} available
                            </span>
                          )}
                        </label>

                        {selectedVariant ? (
                          // Selected state
                          <div className="flex items-center gap-2 p-2.5 border-2 border-green-300 bg-green-50 rounded-lg">
                            <div className="flex-1 min-w-0">
                              <p className="font-bold text-green-800 text-sm truncate">
                                {selectedVariant.lineage}
                              </p>
                              <p className="text-xs text-green-700">
                                {selectedVariant.count.toLocaleString()} sequences available
                              </p>
                            </div>
                            <button
                              type="button"
                              onClick={clearSelectedVariant}
                              className="text-gray-400 hover:text-red-500"
                              title="Change variant"
                            >
                              <MdClose size={16} />
                            </button>
                          </div>
                        ) : (
                          // Combobox: always-visible scrollable dropdown + live filter
                          <>
                            <input
                              type="text"
                              value={variantQuery}
                              onChange={(e) => setVariantQuery(e.target.value)}
                              className="w-full border rounded-lg px-3 py-2 text-sm"
                              placeholder={
                                variantSearchLoading
                                  ? "Loading variants…"
                                  : "Filter: BA.1, XBB, JN.1…"
                              }
                              disabled={variantSearchLoading}
                            />

                            {variantSearchError && (
                              <div className="mt-1 flex items-center gap-2">
                                <p className="text-xs text-red-600 flex-1">{variantSearchError}</p>
                                <button
                                  type="button"
                                  onClick={loadAllVariants}
                                  className="text-xs text-blue-600 hover:underline"
                                >
                                  Retry
                                </button>
                              </div>
                            )}

                            <div className="mt-2 max-h-56 overflow-y-auto border rounded-lg bg-white">
                              {variantSearchLoading && (
                                <div className="px-3 py-4 text-xs text-gray-400 text-center">
                                  Loading variants from cov-spectrum.org…
                                </div>
                              )}
                              {!variantSearchLoading && filteredVariants.length === 0 && variantResults.length > 0 && (
                                <div className="px-3 py-4 text-xs text-gray-400 text-center">
                                  No variants match "{variantQuery}"
                                </div>
                              )}
                              {!variantSearchLoading && filteredVariants.slice(0, 500).map((v) => (
                                <button
                                  key={v.lineage}
                                  type="button"
                                  onClick={() => { setSelectedVariant(v); setVariantQuery(""); }}
                                  className="w-full flex justify-between items-center px-3 py-1.5 text-xs hover:bg-blue-50 text-left border-b last:border-b-0 border-gray-100"
                                >
                                  <span className="font-semibold text-gray-700">{v.lineage}</span>
                                  <span className="text-gray-400 font-mono">{v.count.toLocaleString()}</span>
                                </button>
                              ))}
                              {!variantSearchLoading && filteredVariants.length > 500 && (
                                <div className="px-3 py-2 text-xs text-gray-400 text-center border-t bg-gray-50">
                                  +{(filteredVariants.length - 500).toLocaleString()} more — keep typing to narrow down
                                </div>
                              )}
                            </div>
                          </>
                        )}
                      </div>

                      {/* Step 2: Date range with quick presets (default = all time) */}
                      <div>
                        <label className="text-xs font-semibold text-gray-500 mb-1 block">
                          Date Range <span className="text-gray-400 font-normal">(default: all time)</span>
                        </label>
                        <div className="flex flex-wrap gap-1 mb-2">
                          {[
                            { key: "30d", label: "Last 30 days" },
                            { key: "3mo", label: "Last 3 months" },
                            { key: "1y", label: "Last year" },
                            { key: "all", label: "All time" },
                          ].map((preset) => (
                            <button
                              key={preset.key}
                              type="button"
                              onClick={() => applyDatePreset(preset.key)}
                              className={`px-2 py-1 text-xs rounded-md border transition ${
                                activeDatePreset === preset.key
                                  ? "bg-blue-600 border-blue-600 text-white font-semibold"
                                  : "bg-white border-gray-300 text-gray-600 hover:bg-gray-50"
                              }`}
                            >
                              {preset.label}
                            </button>
                          ))}
                        </div>
                        <div className="grid grid-cols-2 gap-2">
                          <input
                            type="date"
                            value={apiDateFrom}
                            onChange={(e) => setApiDateFrom(e.target.value)}
                            className="w-full border rounded-lg px-2 py-2 text-sm"
                            placeholder="From"
                          />
                          <input
                            type="date"
                            value={apiDateTo}
                            onChange={(e) => setApiDateTo(e.target.value)}
                            className="w-full border rounded-lg px-2 py-2 text-sm"
                            placeholder="To"
                          />
                        </div>
                      </div>

                      <p className="text-xs text-gray-400 leading-snug">
                        Fetches mutations live from <code className="bg-gray-100 px-1 rounded">cov-spectrum.org</code>.
                        The <code className="bg-gray-100 px-1 rounded">proportion</code> of each mutation is used as probability ground truth.
                      </p>
                    </div>
                  )}
                </div>
              )}

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

              {/* Visualization is now available inline via the "Model Visualization"
                  toggle that appears above the benchmark results — no separate page. */}

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
          <div className={`${vizFullWidth ? "lg:col-span-12" : "lg:col-span-8"} space-y-4`}>
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
                      {(results.parameters?.ground_truth_source === "cov_spectrum_csv" ||
                        results.parameters?.ground_truth_source === "cov_spectrum_api") && (
                        <span className="ml-2 px-1.5 py-0.5 rounded bg-purple-100 text-purple-700 font-semibold">
                          GT: {results.parameters.ground_truth_source === "cov_spectrum_api" ? "cov-spectrum API" : "cov-spectrum CSV"}
                          {results.parameters?.ground_truth_meta?.num_mutations_parsed !== undefined && (
                            <> ({results.parameters.ground_truth_meta.num_mutations_parsed} parsed)</>
                          )}
                        </span>
                      )}
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

                {/* Top-level view toggle: Benchmark Results vs full Model Visualization.
                    The visualization is embedded inline so the user can flip back and
                    forth without opening a new tab. Only offered when at least two
                    models were benchmarked (the visualization compares models). */}
                {benchmarkMode === "single" && vizParams && vizParams.models.length >= 2 && vizParams.nodeId && (
                  <div className="flex bg-gray-100 rounded-xl p-1 w-full sm:w-auto self-start">
                    {[
                      { key: "benchmark", label: "Benchmark Results" },
                      { key: "visualization", label: "Model Visualization" },
                    ].map((v) => (
                      <button
                        key={v.key}
                        type="button"
                        onClick={() => {
                          setResultsView(v.key);
                          // Going back to the metrics view always restores the
                          // configuration panel so the form is reachable again.
                          if (v.key === "benchmark") setVizFullWidth(false);
                        }}
                        className={`flex-1 sm:flex-none px-4 py-2 text-sm font-semibold rounded-lg transition-colors ${
                          resultsView === v.key
                            ? "bg-purple-600 text-white shadow-sm"
                            : "text-gray-600 hover:text-gray-800"
                        }`}
                      >
                        {v.label}
                      </button>
                    ))}
                  </div>
                )}

                {/* ---- Inline Model Visualization (embedded CompareModels) ---- */}
                {resultsView === "visualization" && vizParams && (
                  <div className="bg-white rounded-xl shadow-sm border border-gray-100 overflow-hidden">
                    {/* Expand / collapse the left config panel for a near-fullscreen view */}
                    <div className="flex justify-end px-3 py-2 border-b border-gray-100 bg-gray-50">
                      <button
                        type="button"
                        onClick={() => setVizFullWidth((v) => !v)}
                        className="flex items-center gap-1.5 text-xs font-semibold text-gray-600 hover:text-purple-600 transition-colors"
                        title={vizFullWidth ? "Show configuration panel" : "Expand to full width"}
                      >
                        {vizFullWidth ? <MdCloseFullscreen size={16} /> : <MdOpenInFull size={16} />}
                        {vizFullWidth ? "Collapse" : "Full width"}
                      </button>
                    </div>
                    <CompareModels
                      embedded
                      key={`${vizParams.models.join(",")}|${vizParams.nodeId}|${vizParams.elapsedDay}|${vizParams.proteinRegion}`}
                      models={vizParams.models}
                      nodeId={vizParams.nodeId}
                      elapsedDay={vizParams.elapsedDay}
                      proteinRegion={vizParams.proteinRegion}
                    />
                  </div>
                )}

                {/* ---- Benchmark results (metrics, charts, agreement) ---- */}
                {resultsView === "benchmark" && (
                <>
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
                          {
                            key: "auroc",
                            label: "AUROC",
                            higher: true,
                            info: "Area Under the ROC Curve. Measures the model's ability to rank mutated positions higher than non-mutated ones. 0.5 = random guessing, 1.0 = perfect separation. For probability ground truth, positions are binarized at a 0.5 threshold before computing.",
                          },
                          {
                            key: "auprc",
                            label: "AUPRC",
                            higher: true,
                            info: "Area Under the Precision-Recall Curve. More informative than AUROC for imbalanced classes (few mutation positions vs. tens of thousands of non-mutated positions). High = the model captures true mutations well with few false positives.",
                          },
                          {
                            key: "precision",
                            label: "Precision",
                            higher: true,
                            info: "Of all positions the model flagged as mutated (prediction ≥ 0.5), what fraction truly were mutated? TP / (TP + FP). High precision = few false alarms. 0 = every positive prediction was wrong, 1 = every positive prediction was right.",
                          },
                          {
                            key: "recall",
                            label: "Recall",
                            higher: true,
                            info: "Of all actually mutated positions, what fraction did the model successfully flag (prediction ≥ 0.5)? TP / (TP + FN). High recall = few missed mutations. 0 = missed every mutation, 1 = caught every mutation.",
                          },
                          {
                            key: "f1_score",
                            label: "F1 Score",
                            higher: true,
                            info: "Harmonic mean of Precision and Recall: 2·P·R / (P + R). Single number summarizing both. 0 = either precision or recall is zero, 1 = both are perfect. Useful when you need one balanced score.",
                          },
                          {
                            key: "mcc",
                            label: "MCC",
                            higher: true,
                            info: "Matthews Correlation Coefficient. A balanced score in [-1, 1] that stays reliable even when mutated positions are very rare (heavy class imbalance). +1 = perfect agreement, 0 = no better than random, -1 = total disagreement. Often the single most honest summary for imbalanced data.",
                          },
                          {
                            key: "brier_score",
                            label: "Brier Score",
                            higher: false,
                            info: "Mean squared error between predicted probabilities and ground-truth values: mean((pred − gt)²). 0 = perfect, higher = worse. Directly measures probability quality; works with both binary and probability ground truth.",
                          },
                          {
                            key: "ece",
                            label: "ECE",
                            higher: false,
                            info: "Expected Calibration Error. Measures how well predicted probabilities match the actual observed rates. E.g. among positions predicted at 80%, do roughly 80% actually show a mutation? 0 = perfectly calibrated, higher = over- or under-confident.",
                          },
                          {
                            key: "runtime_seconds",
                            label: "Runtime (s)",
                            higher: false,
                            info: "Total prediction time for the model (in seconds). Lower = faster. A practical performance metric when comparing multiple models.",
                          },
                        ].map((metric) => {
                          const values = successfulModels.map(
                            (m) => results.models[m]?.metrics?.[metric.key]
                          );
                          const bestIdx = getBestIdx(metric.key, metric.higher);
                          const isInfoOpen = openMetricInfo === metric.key;

                          return (
                            <tr key={metric.key} className="border-b border-gray-100 hover:bg-gray-50">
                              <td className="py-2 px-3 font-medium text-gray-700 text-xs relative">
                                <div className="flex items-center gap-1.5">
                                  <span>{metric.label}</span>
                                  <span className="text-gray-400">
                                    {metric.higher ? "↑" : "↓"}
                                  </span>
                                  <button
                                    type="button"
                                    onClick={() => setOpenMetricInfo(isInfoOpen ? null : metric.key)}
                                    className={`ml-0.5 text-gray-400 hover:text-blue-500 transition-colors ${isInfoOpen ? "text-blue-500" : ""}`}
                                    title="Metric description"
                                  >
                                    <MdInfo size={14} />
                                  </button>
                                </div>
                                {isInfoOpen && (
                                  <div className="absolute left-0 top-full mt-1 z-20 w-72 p-3 bg-gray-900 text-white text-xs rounded-lg shadow-lg normal-case font-normal leading-relaxed">
                                    <div className="flex justify-between items-start mb-1">
                                      <span className="font-bold text-blue-300">{metric.label}</span>
                                      <button
                                        type="button"
                                        onClick={() => setOpenMetricInfo(null)}
                                        className="text-gray-400 hover:text-white -mt-0.5 -mr-1"
                                      >
                                        <MdClose size={14} />
                                      </button>
                                    </div>
                                    {metric.info}
                                    <div className="mt-2 pt-2 border-t border-gray-700 text-gray-400">
                                      {metric.higher ? "Higher is better ↑" : "Lower is better ↓"}
                                    </div>
                                  </div>
                                )}
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
                    {(() => {
                      const gtType = results.parameters?.ground_truth_type;
                      const gtMeta = results.parameters?.ground_truth_meta || {};
                      const totalGtPoints = results.parameters?.mutation_positions?.length || 0;
                      const mutValues = results.parameters?.mutation_values || [];
                      const visibleGtCount = gtType === "probability"
                        ? mutValues.filter((v) => v >= gtMinThreshold).length
                        : totalGtPoints;

                      return (
                        <>
                          <div className="flex flex-wrap justify-between items-center gap-3 mb-2">
                            {gtType === "probability" ? (
                              <div className="flex items-center gap-2 text-xs text-gray-600 flex-1 min-w-[260px]">
                                <span className="font-semibold text-gray-700 whitespace-nowrap">
                                  GT ≥
                                </span>
                                <input
                                  type="range"
                                  min="0"
                                  max="1"
                                  step="0.01"
                                  value={gtMinThreshold}
                                  onChange={(e) => setGtMinThreshold(parseFloat(e.target.value))}
                                  className="flex-1 accent-gray-700"
                                />
                                <span className="font-mono text-gray-700 w-10 text-right">
                                  {gtMinThreshold.toFixed(2)}
                                </span>
                                <span className="text-gray-400 whitespace-nowrap">
                                  · {visibleGtCount.toLocaleString()} / {totalGtPoints.toLocaleString()} shown
                                </span>
                              </div>
                            ) : <div />}
                            <div className="flex items-center gap-3">
                              {/* View toggle: overlay (all models together) vs stacked (one chart each).
                                  Both reuse data already in the result — no benchmark re-run. */}
                              <div className="flex bg-gray-100 rounded-lg p-0.5">
                                {[
                                  { key: "overlay", label: "Overlay" },
                                  { key: "stacked", label: "Stacked" },
                                ].map((mode) => (
                                  <button
                                    key={mode.key}
                                    type="button"
                                    onClick={() => setPredViewMode(mode.key)}
                                    className={`px-3 py-1 text-xs font-semibold rounded-md transition-colors ${
                                      predViewMode === mode.key
                                        ? "bg-purple-600 text-white shadow-sm"
                                        : "text-gray-600 hover:text-gray-800"
                                    }`}
                                  >
                                    {mode.label}
                                  </button>
                                ))}
                              </div>
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
                          </div>
                          {predViewMode === "overlay" ? (
                            <div style={{ height: "500px" }}>
                              <PredictionOverlayChart
                                models={successfulModels}
                                results={results}
                                useLogScale={useLogScale}
                                gtMinThreshold={gtMinThreshold}
                              />
                            </div>
                          ) : (
                            <div className="space-y-3">
                              {successfulModels.map((modelName, idx) => {
                                const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
                                return (
                                  <div key={modelName}>
                                    <div className="flex items-center gap-2 mb-1 px-1">
                                      <span
                                        className="inline-block w-3 h-3 rounded-sm"
                                        style={{ backgroundColor: mc.border }}
                                      />
                                      <span className="text-xs font-bold" style={{ color: mc.label || mc.border }}>
                                        {modelName}
                                      </span>
                                    </div>
                                    <div style={{ height: "230px" }}>
                                      <PredictionOverlayChart
                                        models={[modelName]}
                                        results={results}
                                        useLogScale={useLogScale}
                                        gtMinThreshold={gtMinThreshold}
                                        colorOffset={idx}
                                      />
                                    </div>
                                  </div>
                                );
                              })}
                            </div>
                          )}
                          <p className="text-xs text-gray-400 mt-2 text-center">
                            {gtType === "probability" ? (
                              <>
                                Lines show predicted mutation probability. Black stars mark actual
                                mutation proportions from cov-spectrum.org (star size scales with proportion).
                                {gtMeta.num_mutations_parsed !== undefined && (
                                  <> Parsed {gtMeta.num_mutations_parsed.toLocaleString()} mutations from CSV.</>
                                )}
                              </>
                            ) : (
                              <>Lines show predicted mutation probability. Stars mark actual mutations.</>
                            )}
                          </p>
                        </>
                      );
                    })()}
                  </Section>
                )}

                {/* Charts Row - Calibration & Per-Protein */}
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <Section title="Calibration" icon={<MdInfo className="text-green-500" />} defaultOpen={true}>
                    <div className="flex justify-end mb-2 relative">
                      <button
                        type="button"
                        onClick={() => setOpenMetricInfo(openMetricInfo === "calibration_info" ? null : "calibration_info")}
                        className={`flex items-center gap-1 text-xs hover:text-blue-500 transition-colors ${openMetricInfo === "calibration_info" ? "text-blue-500" : "text-gray-400"}`}
                        title="What is this?"
                      >
                        <MdInfo size={14} />
                        <span>What is this?</span>
                      </button>
                      {openMetricInfo === "calibration_info" && (
                        <div className="absolute right-0 top-full mt-1 z-20 w-80 p-3 bg-gray-900 text-white text-xs rounded-lg shadow-lg normal-case font-normal leading-relaxed">
                          <div className="flex justify-between items-start mb-1">
                            <span className="font-bold text-blue-300">Calibration Plot</span>
                            <button type="button" onClick={() => setOpenMetricInfo(null)} className="text-gray-400 hover:text-white -mt-0.5 -mr-1">
                              <MdClose size={14} />
                            </button>
                          </div>
                          Compares predicted probabilities (x-axis) with the actual observed mutation fraction (y-axis). Predictions are grouped into 10 bins. The dashed diagonal is perfect calibration — e.g. among positions predicted at 70%, roughly 70% should truly mutate.
                          <div className="mt-2 pt-2 border-t border-gray-700 text-gray-400">
                            Below diagonal = overconfident · Above diagonal = underconfident
                          </div>
                        </div>
                      )}
                    </div>
                    <div style={{ height: "280px" }}>
                      <CalibrationChart models={successfulModels} results={results} />
                    </div>
                  </Section>

                  <Section title="Per-Protein" icon={<MdBarChart className="text-amber-500" />} defaultOpen={true}>
                    <div className="flex justify-between items-center mb-2 relative">
                      <button
                        type="button"
                        onClick={() => setOpenMetricInfo(openMetricInfo === "per_protein_info" ? null : "per_protein_info")}
                        className={`flex items-center gap-1 text-xs hover:text-blue-500 transition-colors ${openMetricInfo === "per_protein_info" ? "text-blue-500" : "text-gray-400"}`}
                        title="What is this?"
                      >
                        <MdInfo size={14} />
                        <span>What is this?</span>
                      </button>
                      <select
                        value={perProteinMetric}
                        onChange={(e) => setPerProteinMetric(e.target.value)}
                        className="text-xs border rounded px-2 py-1"
                      >
                        <option value="mean_prediction">Mean Prediction</option>
                        <option value="auroc">AUROC</option>
                        <option value="auprc">AUPRC</option>
                        <option value="precision">Precision</option>
                        <option value="recall">Recall</option>
                        <option value="f1_score">F1 Score</option>
                        <option value="mcc">MCC</option>
                        <option value="brier_score">Brier Score</option>
                      </select>
                      {openMetricInfo === "per_protein_info" && (
                        <div className="absolute left-0 top-full mt-1 z-20 w-80 p-3 bg-gray-900 text-white text-xs rounded-lg shadow-lg normal-case font-normal leading-relaxed">
                          <div className="flex justify-between items-start mb-1">
                            <span className="font-bold text-blue-300">Per-Protein Performance</span>
                            <button type="button" onClick={() => setOpenMetricInfo(null)} className="text-gray-400 hover:text-white -mt-0.5 -mr-1">
                              <MdClose size={14} />
                            </button>
                          </div>
                          Shows the selected metric broken down by SARS-CoV-2 protein region (Spike, ORF1ab, N, …). Useful to see where each model predicts well or struggles. Some regions (e.g. Spike) mutate more frequently and are typically easier to score.
                          <div className="mt-2 pt-2 border-t border-gray-700 text-gray-400">
                            Switch metric with the dropdown on the right.
                          </div>
                        </div>
                      )}
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
                    <div className="flex justify-end mb-2 relative">
                      <button
                        type="button"
                        onClick={() => setOpenMetricInfo(openMetricInfo === "model_agreement_info" ? null : "model_agreement_info")}
                        className={`flex items-center gap-1 text-xs hover:text-blue-500 transition-colors ${openMetricInfo === "model_agreement_info" ? "text-blue-500" : "text-gray-400"}`}
                        title="What is this?"
                      >
                        <MdInfo size={14} />
                        <span>What is this?</span>
                      </button>
                      {openMetricInfo === "model_agreement_info" && (
                        <div className="absolute right-0 top-full mt-1 z-20 w-80 p-3 bg-gray-900 text-white text-xs rounded-lg shadow-lg normal-case font-normal leading-relaxed">
                          <div className="flex justify-between items-start mb-1">
                            <span className="font-bold text-blue-300">Model Agreement</span>
                            <button type="button" onClick={() => setOpenMetricInfo(null)} className="text-gray-400 hover:text-white -mt-0.5 -mr-1">
                              <MdClose size={14} />
                            </button>
                          </div>
                          Pearson correlation between each pair of models' predicted-probability vectors across all genome positions. 1.0 = predictions are identical; 0.0 = uncorrelated.
                          <div className="mt-2 pt-2 border-t border-gray-700 text-gray-400">
                            Two models with very high agreement (&gt;0.95) are effectively redundant — low agreement often means the models capture complementary signals (good ensemble candidates).
                          </div>
                        </div>
                      )}
                    </div>
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
                </>
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
