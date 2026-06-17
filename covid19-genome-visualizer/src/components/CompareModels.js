import React, { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { useSearchParams, useNavigate } from "react-router-dom";
import GenomeChart from "./Recharts";
import DoughnutChart from "./DoughnutChart";
import Chart from "chart.js/auto";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";
const GENOME_LENGTH = 29903;

// ============================================
// SCALAR / BINARY LINE CHART
// Fallback for non-categorical payloads (binary_per_position, scalar_per_position).
// ============================================
const ScalarLineChart = ({ payload, color }) => {
  const chartRef = useRef(null);
  const chartInst = useRef(null);

  useEffect(() => {
    if (!chartRef.current || !payload) return;
    if (chartInst.current) chartInst.current.destroy();

    const values = payload.predictions?.values;
    if (!Array.isArray(values) || values.length === 0) return;

    const region = payload.domain?.region;
    const offset = region?.start ?? 0;
    const totalLength = payload.domain?.total_length ?? GENOME_LENGTH;

    // Flatten in case values are nested (e.g. [[0.1], [0.2]])
    const flat = values.map((v) => (Array.isArray(v) ? Math.max(...v) : v));

    // Downsample to ~500 points for rendering performance
    const step = Math.max(1, Math.floor(flat.length / 500));
    const data = [];
    for (let i = 0; i < flat.length; i += step) {
      data.push({ x: offset + i, y: flat[i] });
    }

    const taskKind = payload.task?.kind ?? "score";
    const label = taskKind === "binary_per_position" ? "Mutation probability"
                : taskKind === "scalar_per_position" ? "Score"
                : taskKind;

    chartInst.current = new Chart(chartRef.current, {
      type: "scatter",
      data: {
        datasets: [{
          label,
          data,
          borderColor: color,
          backgroundColor: "transparent",
          borderWidth: 1.5,
          pointRadius: 0,
          showLine: true,
          tension: 0.1,
        }],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false,
        scales: {
          x: {
            title: { display: true, text: "Genome Position", font: { size: 11 } },
            min: 0,
            max: totalLength,
            ticks: { font: { size: 10 } },
          },
          y: {
            title: { display: true, text: label, font: { size: 11 } },
            min: 0,
            max: 1,
            ticks: { font: { size: 10 } },
          },
        },
        plugins: {
          legend: { display: false },
          tooltip: {
            callbacks: {
              title: (items) => `Position: ${Math.round(items[0]?.parsed?.x ?? 0)}`,
              label: (ctx) => `${label}: ${ctx.parsed.y.toFixed(4)}`,
            },
          },
        },
      },
    });

    return () => { if (chartInst.current) chartInst.current.destroy(); };
  }, [payload, color]);

  return <canvas ref={chartRef} />;
};

const MODEL_COLORS = [
  { border: "#2563EB", bg: "rgba(37, 99, 235, 0.15)", label: "#1D4ED8", line: "rgba(37, 99, 235, 0.8)" },
  { border: "#DC2626", bg: "rgba(220, 38, 38, 0.15)", label: "#B91C1C", line: "rgba(220, 38, 38, 0.8)" },
  { border: "#059669", bg: "rgba(5, 150, 105, 0.15)", label: "#047857", line: "rgba(5, 150, 105, 0.8)" },
  { border: "#D97706", bg: "rgba(217, 119, 6, 0.15)", label: "#B45309", line: "rgba(217, 119, 6, 0.8)" },
];

// SARS-CoV-2 protein region coordinates (1-based genome positions)
const PROTEIN_REGION_RANGES = {
  ORF1ab: [266, 21555],
  S: [21563, 25384],
  ORF3a: [25393, 26220],
  E: [26245, 26472],
  M: [26523, 27191],
  ORF6: [27202, 27387],
  ORF7a: [27394, 27759],
  ORF7b: [27756, 27887],
  ORF8: [27894, 28259],
  N: [28274, 29533],
  ORF10: [29558, 29674],
};

// ============================================
// SHARED RANGE SLIDER — drives both charts via syncZoomRange
// ============================================
const RangeSlider = ({ value, onChange }) => {
  const [start, end] = value;
  const trackRef = useRef(null);
  const dragging = useRef(null);

  const handleMouseDown = useCallback((which, e) => {
    e.preventDefault();
    dragging.current = which;
    const onMove = (me) => {
      if (!trackRef.current) return;
      const rect = trackRef.current.getBoundingClientRect();
      let pct = (me.clientX - rect.left) / rect.width;
      pct = Math.max(0, Math.min(1, pct));
      onChange((prev) => {
        if (dragging.current === "start") return [Math.min(pct, prev[1] - 0.01), prev[1]];
        if (dragging.current === "end") return [prev[0], Math.max(pct, prev[0] + 0.01)];
        return prev;
      });
    };
    const onUp = () => {
      dragging.current = null;
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    };
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
  }, [onChange]);

  const leftPct = start * 100;
  const widthPct = (end - start) * 100;

  return (
    <div className="relative h-8 select-none">
      <div ref={trackRef} className="absolute inset-x-0 top-3 h-2 bg-slate-200 rounded-full">
        <div
          className="absolute h-full bg-indigo-400 rounded-full"
          style={{ left: `${leftPct}%`, width: `${widthPct}%` }}
        />
        <div
          className="absolute top-1/2 -translate-y-1/2 w-4 h-4 bg-white border-2 border-indigo-500 rounded-full cursor-ew-resize shadow-sm hover:scale-110 transition-transform"
          style={{ left: `calc(${leftPct}% - 8px)` }}
          onMouseDown={(e) => handleMouseDown("start", e)}
        />
        <div
          className="absolute top-1/2 -translate-y-1/2 w-4 h-4 bg-white border-2 border-indigo-500 rounded-full cursor-ew-resize shadow-sm hover:scale-110 transition-transform"
          style={{ left: `calc(${end * 100}% - 8px)` }}
          onMouseDown={(e) => handleMouseDown("end", e)}
        />
      </div>
    </div>
  );
};


// ============================================
// OVERLAY CHART - All models on one synchronized chart
// ============================================
const OverlayChart = ({ predictions, modelNames, useLogScale, nucBreakdown }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const datasets = [];
    // Per-nucleotide stroke styling — dashed lines so they stay readable when
    // overlapped with the same model color.
    const NUC_DASH = { 0: [], 1: [4, 3], 2: [2, 2], 3: [6, 3, 2, 3] }; // A, T, G, C

    modelNames.forEach((name, idx) => {
      const pred = predictions[name];
      if (!pred?.genomeDataRaw) return;

      const color = MODEL_COLORS[idx % MODEL_COLORS.length];
      const numPositions = pred.genomeDataRaw[0]?.length || 0;
      const step = Math.max(1, Math.floor(numPositions / 2000));

      if (nucBreakdown) {
        // One line per nucleotide × model (4 lines per model). Per-position value
        // is the predicted probability for that specific base.
        ["A", "T", "G", "C"].forEach((nuc, nucIdx) => {
          const series = [];
          for (let i = 0; i < numPositions; i += step) {
            series.push({ x: i, y: pred.genomeDataRaw[nucIdx]?.[i] || 0 });
          }
          datasets.push({
            label: `${name} · ${nuc}`,
            data: series,
            borderColor: color.line,
            backgroundColor: "transparent",
            borderWidth: 1.2,
            borderDash: NUC_DASH[nucIdx],
            pointRadius: 0,
            showLine: true,
            tension: 0.1,
            fill: false,
          });
        });
      } else {
        // MAX over A/T/G/C — the most informative single-line summary. Tells you
        // "what is the highest probability assigned to any mutation at this
        // position", which is the metric the model is actually scoring.
        const series = [];
        for (let i = 0; i < numPositions; i += step) {
          const a = pred.genomeDataRaw[0]?.[i] || 0;
          const t = pred.genomeDataRaw[1]?.[i] || 0;
          const g = pred.genomeDataRaw[2]?.[i] || 0;
          const c = pred.genomeDataRaw[3]?.[i] || 0;
          series.push({ x: i, y: Math.max(a, t, g, c) });
        }
        datasets.push({
          label: name,
          data: series,
          borderColor: color.line,
          backgroundColor: "transparent",
          borderWidth: 1.5,
          pointRadius: 0,
          showLine: true,
          tension: 0.1,
          fill: false,
        });
      }
    });

    if (datasets.length === 0) return;

    chartInstance.current = new Chart(chartRef.current, {
      type: "scatter",
      data: { datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false,
        interaction: {
          mode: "index",
          intersect: false,
        },
        scales: {
          x: {
            title: { display: true, text: "Genome Position", font: { size: 12 } },
            min: 0,
            ticks: { font: { size: 10 } },
          },
          y: {
            type: useLogScale ? "logarithmic" : "linear",
            title: {
              display: true,
              text: nucBreakdown
                ? (useLogScale ? "Per-Nucleotide Probability (log)" : "Per-Nucleotide Probability")
                : (useLogScale ? "Max Mutation Probability (log)" : "Max Mutation Probability"),
              font: { size: 12 },
            },
            ...(useLogScale ? { min: 0.0001, max: 1 } : { min: 0, max: 1 }),
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
          legend: {
            position: "top",
            labels: {
              usePointStyle: true,
              padding: 16,
              font: { size: 12, weight: "bold" },
            },
          },
          tooltip: {
            backgroundColor: "rgba(0,0,0,0.85)",
            titleFont: { size: 12 },
            bodyFont: { size: 11 },
            padding: 10,
            callbacks: {
              title: (items) => `Position: ${Math.round(items[0]?.parsed?.x || 0)}`,
              label: (ctx) => `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(6)}`,
            },
          },
        },
      },
    });

    return () => { if (chartInstance.current) chartInstance.current.destroy(); };
  }, [predictions, modelNames, useLogScale]);

  return <canvas ref={chartRef} />;
};

// ============================================
// DIFFERENCE CHART - Shows where models diverge
// ============================================
const DifferenceChart = ({ predictions, modelNames, useLogScale }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current || modelNames.length < 2) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const pred1 = predictions[modelNames[0]];
    const pred2 = predictions[modelNames[1]];
    if (!pred1?.genomeDataRaw || !pred2?.genomeDataRaw) return;

    const numPositions = Math.min(
      pred1.genomeDataRaw[0]?.length || 0,
      pred2.genomeDataRaw[0]?.length || 0
    );
    const step = Math.max(1, Math.floor(numPositions / 2000));

    // Difference is taken on the MAX-per-position metric (matching OverlayChart's
    // default summary line). Values lie in [-1, +1] so the chart stays compact.
    const diffData = [];
    for (let i = 0; i < numPositions; i += step) {
      const max1 = Math.max(
        pred1.genomeDataRaw[0]?.[i] || 0,
        pred1.genomeDataRaw[1]?.[i] || 0,
        pred1.genomeDataRaw[2]?.[i] || 0,
        pred1.genomeDataRaw[3]?.[i] || 0,
      );
      const max2 = Math.max(
        pred2.genomeDataRaw[0]?.[i] || 0,
        pred2.genomeDataRaw[1]?.[i] || 0,
        pred2.genomeDataRaw[2]?.[i] || 0,
        pred2.genomeDataRaw[3]?.[i] || 0,
      );
      diffData.push({ x: i, y: max1 - max2 });
    }

    chartInstance.current = new Chart(chartRef.current, {
      type: "scatter",
      data: {
        datasets: [{
          label: `${modelNames[0]} - ${modelNames[1]}`,
          data: diffData,
          borderColor: "rgba(139, 92, 246, 0.7)",
          backgroundColor: diffData.map(d => d.y >= 0 ? "rgba(37, 99, 235, 0.3)" : "rgba(220, 38, 38, 0.3)"),
          borderWidth: 1,
          pointRadius: 1.5,
          showLine: true,
          tension: 0.1,
          fill: true,
        }],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false,
        scales: {
          x: {
            title: { display: true, text: "Genome Position", font: { size: 12 } },
            min: 0,
            ticks: { font: { size: 10 } },
          },
          y: {
            // Sym range [-1, 1] keeps the diff chart compact and centered on 0
            // (max1 - max2 is bounded by [-1, +1] since both are probabilities).
            min: -1,
            max: 1,
            title: { display: true, text: "Δ Max Probability", font: { size: 11 } },
            ticks: { font: { size: 10 }, stepSize: 0.5 },
            grid: {
              color: (ctx) => ctx.tick.value === 0 ? "rgba(0,0,0,0.25)" : "rgba(0,0,0,0.05)",
              lineWidth: (ctx) => ctx.tick.value === 0 ? 1.5 : 1,
            },
          },
        },
        plugins: {
          legend: { display: false },
          tooltip: {
            callbacks: {
              title: (items) => `Position: ${Math.round(items[0]?.parsed?.x || 0)}`,
              label: (ctx) => {
                const val = ctx.parsed.y;
                const sign = val >= 0 ? "+" : "";
                return `Diff: ${sign}${val.toFixed(6)}`;
              },
            },
          },
        },
      },
    });

    return () => { if (chartInstance.current) chartInstance.current.destroy(); };
  }, [predictions, modelNames, useLogScale]);

  return <canvas ref={chartRef} />;
};


const CompareModels = ({
  embedded = false,
  models: modelsProp,
  nodeId: nodeIdProp,
  elapsedDay: elapsedDayProp,
  proteinRegion: proteinRegionProp,
} = {}) => {
  const [searchParams] = useSearchParams();
  const navigate = useNavigate();

  // When embedded inside the benchmark page we receive parameters as props;
  // as a standalone /compare route we read them from the URL query string.
  const models = modelsProp
    ? modelsProp.filter(Boolean)
    : (searchParams.get("models") || "").split(",").filter(Boolean);
  const nodeId = nodeIdProp ?? (searchParams.get("nodeId") || "");
  const elapsedDay = elapsedDayProp ?? (Number(searchParams.get("elapsedDay")) || 60);
  const proteinRegion = proteinRegionProp ?? (searchParams.get("proteinRegion") || "");

  const [predictions, setPredictions] = useState({});
  const [loading, setLoading] = useState({});
  const [errors, setErrors] = useState({});
  // View modes: "overlay" | "stacked" | "sidebyside"
  const [viewMode, setViewMode] = useState("overlay");
  const [useLogScale, setUseLogScale] = useState(false);
  const [showDoughnut, setShowDoughnut] = useState(true);
  // Overlay-mode toggle: show one summary line per model (MAX A/T/G/C) vs four
  // dashed lines per model (one per nucleotide).
  const [nucBreakdown, setNucBreakdown] = useState(false);
  // Sync zoom: when ON, zooming one chart updates the others (and the slider).
  const [zoomSyncEnabled, setZoomSyncEnabled] = useState(true);
  // Shared genome view range as fractional [start, end] in [0, 1] — drives the
  // syncZoomRange prop on every GenomeChart so they stay aligned.
  const [viewRange, setViewRange] = useState([0, 1]);
  // Tracks which chart initiated the most recent manual zoom so we don't echo
  // its own zoom event back to it.
  const zoomSourceRef = useRef(null);

  // Convert fractional viewRange → genome-position {min,max} that GenomeChart accepts.
  const sharedZoomRange = useMemo(() => {
    if (!zoomSyncEnabled) return null;
    return {
      min: Math.round(viewRange[0] * GENOME_LENGTH),
      max: Math.round(viewRange[1] * GENOME_LENGTH),
    };
  }, [viewRange, zoomSyncEnabled]);

  // When a chart manually zooms, update the slider so they all stay in sync.
  const handleChartZoom = useCallback((modelId) => (range) => {
    zoomSourceRef.current = modelId;
    setViewRange([
      Math.max(0, range.min / GENOME_LENGTH),
      Math.min(1, range.max / GENOME_LENGTH),
    ]);
  }, []);

  // Click a protein-region badge to zoom both charts to that region.
  const zoomToRegion = useCallback((region) => {
    const r = PROTEIN_REGION_RANGES[region];
    if (!r) return;
    zoomSourceRef.current = null;
    setViewRange([r[0] / GENOME_LENGTH, r[1] / GENOME_LENGTH]);
  }, []);

  useEffect(() => {
    if (!models.length || !nodeId) return;
    models.forEach(async (modelId) => {
      const displayName = modelId.startsWith("uploaded:") ? modelId.replace("uploaded:", "") : modelId;
      setLoading((p) => ({ ...p, [displayName]: true }));
      setErrors((p) => ({ ...p, [displayName]: null }));
      try {
        const fd = new FormData();
        fd.append("selectedModel", modelId);
        fd.append("nodeId", nodeId);
        fd.append("elapsedDay", elapsedDay);
        if (proteinRegion) fd.append("selectedProteinRegion", proteinRegion);

        const res = await fetch(`${API_URL}/api/predict/`, { method: "POST", body: fd });
        if (res.ok) {
          const data = await res.json();
          let genomeDataRaw = null;
          if (data.genomeData) {
            const gd = data.genomeData;
            if (Array.isArray(gd) && Array.isArray(gd[0]) && typeof gd[0][0] === "number") {
              genomeDataRaw = gd;
            } else if (Array.isArray(gd) && gd[0]?.mutationPoss) {
              genomeDataRaw = [
                gd.map((d) => d.mutationPoss?.A || 0),
                gd.map((d) => d.mutationPoss?.T || 0),
                gd.map((d) => d.mutationPoss?.G || 0),
                gd.map((d) => d.mutationPoss?.C || 0),
              ];
            } else if (Array.isArray(gd) && typeof gd[0] === "number") {
              genomeDataRaw = [gd, gd, gd, gd];
            }
          }
          setPredictions((p) => ({
            ...p,
            [displayName]: {
              genomeDataRaw,
              predictionPayload: data.predictionPayload || null,
              genomeSequence: data.genomeSequence || "",
              proteinMutationProbs: data.protein_mutation_probs || {},
            },
          }));
        } else {
          const err = await res.json().catch(() => ({}));
          setErrors((p) => ({ ...p, [displayName]: err.error || `HTTP ${res.status}` }));
        }
      } catch (e) {
        setErrors((p) => ({ ...p, [displayName]: e.message }));
      } finally {
        setLoading((p) => ({ ...p, [displayName]: false }));
      }
    });
    // eslint-disable-next-line
  }, []);

  const modelNames = models.map((m) => (m.startsWith("uploaded:") ? m.replace("uploaded:", "") : m));
  const allDone = modelNames.every((n) => !loading[n]);
  const anyLoading = modelNames.some((n) => loading[n]);
  const loadedModels = modelNames.filter((n) => predictions[n]?.genomeDataRaw || predictions[n]?.predictionPayload);

  return (
    <div className={embedded ? "bg-[#f6f7f9]" : "min-h-screen bg-[#f6f7f9]"}>
      {/* Header */}
      <div className={`bg-white border-b shadow-sm z-10 ${embedded ? "" : "sticky top-0"}`}>
        <div className="max-w-[1800px] mx-auto px-6 py-3 flex items-center justify-between flex-wrap gap-2">
          {!embedded && (
            <div className="flex items-center gap-4">
              <button onClick={() => navigate("/benchmark")} className="text-gray-400 hover:text-gray-700 text-sm">
                ← Back to Benchmark
              </button>
              <h1 className="font-bold text-gray-800 text-lg">Model Visualization</h1>
            </div>
          )}

          {/* Controls */}
          <div className="flex items-center gap-3 flex-wrap">
            {/* View Mode Toggle */}
            <div className="flex bg-gray-100 rounded-lg p-0.5">
              {[
                { key: "overlay", label: "Overlay" },
                { key: "stacked", label: "Stacked" },
                { key: "sidebyside", label: "Side by Side" },
              ].map((mode) => (
                <button
                  key={mode.key}
                  onClick={() => setViewMode(mode.key)}
                  className={`px-3 py-1.5 text-xs font-semibold rounded-md transition-colors ${
                    viewMode === mode.key
                      ? "bg-blue-600 text-white shadow-sm"
                      : "text-gray-600 hover:text-gray-800"
                  }`}
                >
                  {mode.label}
                </button>
              ))}
            </div>

            {/* Log Scale Toggle - only meaningful for OverlayChart */}
            {viewMode === "overlay" && (
              <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                <input
                  type="checkbox"
                  checked={useLogScale}
                  onChange={(e) => setUseLogScale(e.target.checked)}
                  className="rounded text-blue-600"
                />
                Log Scale
              </label>
            )}

            {/* Nucleotide breakdown - only in overlay mode */}
            {viewMode === "overlay" && (
              <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                <input
                  type="checkbox"
                  checked={nucBreakdown}
                  onChange={(e) => setNucBreakdown(e.target.checked)}
                  className="rounded text-indigo-600"
                />
                A/T/G/C breakdown
              </label>
            )}

            {/* Sync Zoom Toggle - charts zoom together in stacked/sidebyside */}
            {(viewMode === "stacked" || viewMode === "sidebyside") && (
              <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                <input
                  type="checkbox"
                  checked={zoomSyncEnabled}
                  onChange={(e) => setZoomSyncEnabled(e.target.checked)}
                  className="rounded text-purple-600"
                />
                Sync Zoom
              </label>
            )}

            {/* Doughnut Toggle */}
            {viewMode !== "overlay" && (
              <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                <input
                  type="checkbox"
                  checked={showDoughnut}
                  onChange={(e) => setShowDoughnut(e.target.checked)}
                  className="rounded text-blue-600"
                />
                Protein Chart
              </label>
            )}

            {/* Model Legend */}
            <div className="flex items-center gap-2 ml-2">
              {modelNames.map((name, idx) => {
                const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
                return (
                  <div key={name} className="flex items-center gap-1">
                    <div className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: mc.border }} />
                    <span className="text-xs font-semibold" style={{ color: mc.label }}>{name}</span>
                  </div>
                );
              })}
            </div>

            {/* Variant Info */}
            <div className="text-xs text-gray-400">
              {nodeId.substring(0, 30)}... | Day {elapsedDay} | {proteinRegion || "Full Genome"}
            </div>
          </div>
        </div>
      </div>

      <div className="max-w-[1800px] mx-auto px-4 py-4">
        {/* Loading */}
        {anyLoading && (
          <div className="text-center py-6 text-gray-400 animate-pulse flex items-center justify-center gap-2">
            <div className="animate-spin rounded-full h-5 w-5 border-2 border-blue-500 border-t-transparent" />
            Loading predictions... ({loadedModels.length}/{modelNames.length})
          </div>
        )}

        {/* ============================================ */}
        {/* OVERLAY MODE - All models on one chart. Overlay + Difference are */}
        {/* fused into a single card so both fit on screen without scrolling. */}
        {/* ============================================ */}
        {viewMode === "overlay" && allDone && loadedModels.length > 0 && (
          <div className="space-y-4">
            <div className="bg-white rounded-xl shadow-sm border overflow-hidden">
              <div className="px-4 py-2.5 border-b bg-gray-50 flex items-center justify-between">
                <div>
                  <h2 className="text-sm font-bold text-gray-700">
                    Prediction Overlay {nucBreakdown ? "— Per-Nucleotide Probability" : "— Max Mutation Probability"}
                  </h2>
                  <p className="text-[11px] text-gray-400 mt-0.5">
                    {nucBreakdown
                      ? "Four dashed lines per model (A solid, T dashed, G dotted, C dash-dot)."
                      : "Shows max(P(A), P(T), P(G), P(C)) per position — the most likely mutation at each site."}
                  </p>
                </div>
              </div>

              {/* Overlay chart — taller, primary view */}
              <div style={{ height: "42vh", minHeight: "320px" }} className="px-3 pt-3 pb-1">
                <OverlayChart
                  predictions={predictions}
                  modelNames={loadedModels}
                  useLogScale={useLogScale}
                  nucBreakdown={nucBreakdown}
                />
              </div>

              {/* Difference chart — compact, attached directly below overlay (no */}
              {/* card-in-card, sticks to the top so the user can see both at */}
              {/* once on a single screen). */}
              {loadedModels.length === 2 && (
                <>
                  <div className="px-4 pt-1 pb-0.5 text-[11px] text-gray-500 border-t border-gray-100 bg-gray-50/40 flex items-center justify-between">
                    <span className="font-semibold text-gray-600">
                      Δ Max Probability — {loadedModels[0]} − {loadedModels[1]}
                    </span>
                    <span className="text-[10px] text-gray-400">
                      Above 0 = first higher · Below 0 = second higher
                    </span>
                  </div>
                  <div style={{ height: "18vh", minHeight: "150px" }} className="px-3 pt-1 pb-3">
                    <DifferenceChart
                      predictions={predictions}
                      modelNames={loadedModels}
                      useLogScale={useLogScale}
                    />
                  </div>
                </>
              )}
            </div>

            {/* Doughnut Charts side by side in overlay mode */}
            {loadedModels.some(n => predictions[n]?.proteinMutationProbs && Object.keys(predictions[n].proteinMutationProbs).length > 0) && (
              <div className="bg-white rounded-xl shadow-sm border overflow-hidden">
                <div className="px-4 py-3 border-b bg-gray-50">
                  <h2 className="text-sm font-bold text-gray-700">Protein Region Distribution</h2>
                </div>
                <div className={`grid gap-4 p-4 ${loadedModels.length <= 2 ? "grid-cols-2" : loadedModels.length === 3 ? "grid-cols-3" : "grid-cols-4"}`}>
                  {loadedModels.map((name, idx) => {
                    const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
                    const pred = predictions[name];
                    if (!pred?.proteinMutationProbs || Object.keys(pred.proteinMutationProbs).length === 0) return null;
                    return (
                      <div key={name} className="text-center">
                        <div className="flex items-center justify-center gap-1.5 mb-2">
                          <div className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: mc.border }} />
                          <span className="text-sm font-bold" style={{ color: mc.label }}>{name}</span>
                        </div>
                        <DoughnutChart data={pred.proteinMutationProbs} />
                      </div>
                    );
                  })}
                </div>
              </div>
            )}
          </div>
        )}

        {/* ============================================ */}
        {/* SHARED CONTROLS for stacked/sidebyside modes */}
        {/*  - Protein Region badges (click to zoom both charts to that region) */}
        {/*  - Range slider (drag handles to zoom both charts to a custom range) */}
        {/* ============================================ */}
        {(viewMode === "stacked" || viewMode === "sidebyside") && allDone && loadedModels.length > 0 && (
          <div className="bg-white rounded-xl p-4 shadow-sm border border-gray-100 mb-4">
            <div className="flex items-center justify-between mb-3">
              <h2 className="text-sm font-bold text-gray-700">Protein Regions</h2>
              <span className="text-[11px] text-gray-400 font-mono">
                Position {Math.round(viewRange[0] * GENOME_LENGTH) + 1} — {Math.round(viewRange[1] * GENOME_LENGTH)}
                {" "}({Math.round((viewRange[1] - viewRange[0]) * GENOME_LENGTH).toLocaleString()} bp)
              </span>
            </div>
            <div className="flex flex-wrap gap-1.5 mb-4">
              {Object.entries(PROTEIN_REGION_RANGES).map(([region, [start, end]]) => (
                <button
                  key={region}
                  type="button"
                  onClick={() => zoomToRegion(region)}
                  className="text-xs font-medium px-2.5 py-1 rounded-full border border-gray-200 hover:border-gray-400 hover:shadow-sm transition-all"
                  style={{ backgroundColor: proteinRegionColorMap[region] || "#E5E7EB" }}
                  title={`${start}-${end}`}
                >
                  <span className="font-bold text-gray-800">{region}</span>
                  <span className="ml-1.5 text-[10px] text-gray-600 font-mono">{start}-{end}</span>
                </button>
              ))}
              <button
                type="button"
                onClick={() => { zoomSourceRef.current = null; setViewRange([0, 1]); }}
                className="text-xs font-medium px-2.5 py-1 rounded-full bg-gray-700 text-white hover:bg-gray-900 transition-colors ml-auto"
              >
                Reset View
              </button>
            </div>
            <div className="px-1">
              <RangeSlider
                value={viewRange}
                onChange={(updater) => {
                  zoomSourceRef.current = null;
                  setViewRange(updater);
                }}
              />
            </div>
            <div className="flex justify-between text-[10px] text-gray-400 font-mono mt-1 px-1">
              <span>1</span>
              <span>{GENOME_LENGTH.toLocaleString()}</span>
            </div>
          </div>
        )}

        {/* ============================================ */}
        {/* STACKED MODE - GenomeChart from main prediction page, no per-chart */}
        {/* sidebar (handled by shared widget above). Charts sync via the */}
        {/* shared range slider + Sync Zoom toggle. */}
        {/* ============================================ */}
        {viewMode === "stacked" && (
          <div className="space-y-4">
            {models.map((modelId, idx) => {
              const name = modelNames[idx];
              const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
              const pred = predictions[name];
              const isLoading = loading[name];
              const error = errors[name];

              return (
                <div key={modelId}>
                  <div
                    className="flex items-center gap-2.5 px-3 py-2 rounded-t-xl"
                    style={{ background: mc.bg, borderLeft: `4px solid ${mc.border}` }}
                  >
                    <div className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: mc.border }} />
                    <span className="font-bold text-sm" style={{ color: mc.label }}>{name}</span>
                    {isLoading && (
                      <div className="flex items-center gap-2 ml-auto">
                        <div className="animate-spin h-3 w-3 border-2 border-t-transparent rounded-full" style={{ borderColor: mc.border }} />
                        <span className="text-xs" style={{ color: mc.border }}>Computing...</span>
                      </div>
                    )}
                    {error && <span className="text-xs text-red-500 ml-auto">{error}</span>}
                  </div>

                  <div className="bg-white rounded-b-xl shadow-sm border overflow-hidden"
                       style={{ borderLeftColor: mc.border, borderLeftWidth: "4px" }}>
                    {isLoading ? (
                      <div className="flex items-center justify-center py-16">
                        <div className="animate-spin rounded-full h-8 w-8 border-2 border-blue-500 border-t-transparent" />
                      </div>
                    ) : pred?.genomeDataRaw ? (
                      <GenomeChart
                        key={`genome-chart-${modelId}`}
                        genomeData={pred.genomeDataRaw}
                        genomeSequence={pred.genomeSequence}
                        hideSidebar={true}
                        onZoomSync={zoomSyncEnabled ? handleChartZoom(modelId) : undefined}
                        syncZoomRange={zoomSyncEnabled && sharedZoomRange && zoomSourceRef.current !== modelId ? sharedZoomRange : undefined}
                      />
                    ) : pred?.predictionPayload ? (
                      <div style={{ height: 300, padding: "12px 4px" }}>
                        <ScalarLineChart payload={pred.predictionPayload} color={mc.border} />
                      </div>
                    ) : error ? (
                      <div className="py-12 text-center text-red-400 text-sm">{error}</div>
                    ) : null}
                  </div>
                </div>
              );
            })}
          </div>
        )}

        {/* ============================================ */}
        {/* SIDE BY SIDE MODE - same as stacked but in a horizontal grid */}
        {/* ============================================ */}
        {viewMode === "sidebyside" && (
          <div
            className={`grid gap-4 ${
              models.length <= 2 ? "grid-cols-2" : models.length === 3 ? "grid-cols-3" : "grid-cols-2"
            }`}
          >
            {models.map((modelId, idx) => {
              const name = modelNames[idx];
              const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
              const pred = predictions[name];
              const isLoading = loading[name];
              const error = errors[name];

              return (
                <div key={modelId} className="flex flex-col">
                  <div
                    className="flex items-center gap-2.5 px-3 py-2 rounded-t-xl"
                    style={{ background: mc.bg, borderTop: `3px solid ${mc.border}` }}
                  >
                    <div className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: mc.border }} />
                    <span className="font-bold text-sm" style={{ color: mc.label }}>{name}</span>
                    {isLoading && (
                      <div className="flex items-center gap-2 ml-auto">
                        <div className="animate-spin h-3 w-3 border-2 border-t-transparent rounded-full" style={{ borderColor: mc.border }} />
                        <span className="text-xs" style={{ color: mc.border }}>Computing...</span>
                      </div>
                    )}
                    {error && <span className="text-xs text-red-500 ml-auto">{error}</span>}
                  </div>

                  <div className="bg-white rounded-b-xl shadow-sm border overflow-hidden flex-1">
                    {isLoading ? (
                      <div className="flex items-center justify-center py-16">
                        <div className="animate-spin rounded-full h-8 w-8 border-2 border-blue-500 border-t-transparent" />
                      </div>
                    ) : pred?.genomeDataRaw ? (
                      <GenomeChart
                        key={`genome-chart-side-${modelId}`}
                        genomeData={pred.genomeDataRaw}
                        genomeSequence={pred.genomeSequence}
                        hideSidebar={true}
                        onZoomSync={zoomSyncEnabled ? handleChartZoom(modelId) : undefined}
                        syncZoomRange={zoomSyncEnabled && sharedZoomRange && zoomSourceRef.current !== modelId ? sharedZoomRange : undefined}
                      />
                    ) : pred?.predictionPayload ? (
                      <div style={{ height: 300, padding: "12px 4px" }}>
                        <ScalarLineChart payload={pred.predictionPayload} color={mc.border} />
                      </div>
                    ) : error ? (
                      <div className="py-12 text-center text-red-400 text-sm">{error}</div>
                    ) : null}
                  </div>
                </div>
              );
            })}
          </div>
        )}

        {/* ============================================ */}
        {/* PROTEIN MUTATION PROBABILITIES - same Doughnut grid as the main */}
        {/* prediction page, one card per model (replaces the comparison table). */}
        {/* ============================================ */}
        {(viewMode === "stacked" || viewMode === "sidebyside") && allDone && showDoughnut &&
          loadedModels.some((n) => predictions[n]?.proteinMutationProbs && Object.keys(predictions[n].proteinMutationProbs).length > 0) && (
            <div className="bg-white rounded-xl shadow-sm border overflow-hidden mt-4">
              <div className="px-4 py-3 border-b bg-gray-50">
                <h2 className="text-sm font-bold text-gray-700">Protein Mutation Probabilities</h2>
              </div>
              <div className={`grid gap-4 p-4 ${loadedModels.length <= 2 ? "grid-cols-2" : loadedModels.length === 3 ? "grid-cols-3" : "grid-cols-4"}`}>
                {loadedModels.map((name, idx) => {
                  const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
                  const pred = predictions[name];
                  if (!pred?.proteinMutationProbs || Object.keys(pred.proteinMutationProbs).length === 0) return null;
                  return (
                    <div key={name} className="text-center">
                      <div className="flex items-center justify-center gap-1.5 mb-2">
                        <div className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: mc.border }} />
                        <span className="text-sm font-bold" style={{ color: mc.label }}>{name}</span>
                      </div>
                      <DoughnutChart data={pred.proteinMutationProbs} />
                    </div>
                  );
                })}
              </div>
            </div>
          )}
      </div>
    </div>
  );
};

export default CompareModels;
