import React, { useState, useEffect, useRef, useMemo } from "react";
import { useSearchParams, useNavigate } from "react-router-dom";
import GenomeChart from "./Recharts";
import DoughnutChart from "./DoughnutChart";
import Chart from "chart.js/auto";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

const MODEL_COLORS = [
  { border: "#2563EB", bg: "rgba(37, 99, 235, 0.15)", label: "#1D4ED8", line: "rgba(37, 99, 235, 0.8)" },
  { border: "#DC2626", bg: "rgba(220, 38, 38, 0.15)", label: "#B91C1C", line: "rgba(220, 38, 38, 0.8)" },
  { border: "#059669", bg: "rgba(5, 150, 105, 0.15)", label: "#047857", line: "rgba(5, 150, 105, 0.8)" },
  { border: "#D97706", bg: "rgba(217, 119, 6, 0.15)", label: "#B45309", line: "rgba(217, 119, 6, 0.8)" },
];

// ============================================
// OVERLAY CHART - All models on one synchronized chart
// ============================================
const OverlayChart = ({ predictions, modelNames, useLogScale }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const datasets = [];

    modelNames.forEach((name, idx) => {
      const pred = predictions[name];
      if (!pred?.genomeDataRaw) return;

      const color = MODEL_COLORS[idx % MODEL_COLORS.length];
      // Sum all 4 nucleotide probabilities per position for total mutation probability
      const totalProbs = [];
      const numPositions = pred.genomeDataRaw[0]?.length || 0;
      const step = Math.max(1, Math.floor(numPositions / 2000)); // Downsample for performance

      for (let i = 0; i < numPositions; i += step) {
        const total = (pred.genomeDataRaw[0][i] || 0) +
                      (pred.genomeDataRaw[1][i] || 0) +
                      (pred.genomeDataRaw[2][i] || 0) +
                      (pred.genomeDataRaw[3][i] || 0);
        totalProbs.push({ x: i, y: total });
      }

      datasets.push({
        label: name,
        data: totalProbs,
        borderColor: color.line,
        backgroundColor: "transparent",
        borderWidth: 1.5,
        pointRadius: 0,
        showLine: true,
        tension: 0.1,
        fill: false,
      });
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
              text: useLogScale ? "Total Mutation Probability (log)" : "Total Mutation Probability",
              font: { size: 12 },
            },
            ...(useLogScale ? { min: 0.0001 } : { min: 0 }),
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

    const diffData = [];
    for (let i = 0; i < numPositions; i += step) {
      const total1 = (pred1.genomeDataRaw[0][i] || 0) + (pred1.genomeDataRaw[1][i] || 0) +
                     (pred1.genomeDataRaw[2][i] || 0) + (pred1.genomeDataRaw[3][i] || 0);
      const total2 = (pred2.genomeDataRaw[0][i] || 0) + (pred2.genomeDataRaw[1][i] || 0) +
                     (pred2.genomeDataRaw[2][i] || 0) + (pred2.genomeDataRaw[3][i] || 0);
      diffData.push({ x: i, y: total1 - total2 });
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
            title: { display: true, text: "Probability Difference", font: { size: 12 } },
            ticks: { font: { size: 10 } },
          },
        },
        plugins: {
          legend: { position: "top", labels: { font: { size: 12 } } },
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


const CompareModels = () => {
  const [searchParams] = useSearchParams();
  const navigate = useNavigate();

  const models = (searchParams.get("models") || "").split(",").filter(Boolean);
  const nodeId = searchParams.get("nodeId") || "";
  const elapsedDay = Number(searchParams.get("elapsedDay")) || 60;
  const proteinRegion = searchParams.get("proteinRegion") || "";

  const [predictions, setPredictions] = useState({});
  const [loading, setLoading] = useState({});
  const [errors, setErrors] = useState({});
  // View modes: "overlay" | "stacked" | "sidebyside"
  const [viewMode, setViewMode] = useState("overlay");
  const [useLogScale, setUseLogScale] = useState(false);
  const [showDoughnut, setShowDoughnut] = useState(true);
  const [zoomSyncEnabled, setZoomSyncEnabled] = useState(true);
  // Shared zoom range for synchronized charts - { min, max } genome positions
  const [sharedZoomRange, setSharedZoomRange] = useState(null);
  const zoomSourceRef = useRef(null); // tracks which chart initiated zoom

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
  const loadedModels = modelNames.filter((n) => predictions[n]?.genomeDataRaw);

  return (
    <div className="min-h-screen bg-[#f6f7f9]">
      {/* Header */}
      <div className="bg-white border-b shadow-sm sticky top-0 z-10">
        <div className="max-w-[1800px] mx-auto px-6 py-3 flex items-center justify-between flex-wrap gap-2">
          <div className="flex items-center gap-4">
            <button onClick={() => navigate("/benchmark")} className="text-gray-400 hover:text-gray-700 text-sm">
              ← Back to Benchmark
            </button>
            <h1 className="font-bold text-gray-800 text-lg">Visual Model Comparison</h1>
          </div>

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

            {/* Log Scale Toggle */}
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

            {/* Zoom Sync Toggle */}
            {(viewMode === "stacked" || viewMode === "sidebyside") && (
              <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                <input
                  type="checkbox"
                  checked={zoomSyncEnabled}
                  onChange={(e) => {
                    setZoomSyncEnabled(e.target.checked);
                    if (!e.target.checked) setSharedZoomRange(null);
                  }}
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
        {/* OVERLAY MODE - All models on one chart */}
        {/* ============================================ */}
        {viewMode === "overlay" && allDone && loadedModels.length > 0 && (
          <div className="space-y-4">
            {/* Main Overlay Chart */}
            <div className="bg-white rounded-xl shadow-sm border overflow-hidden">
              <div className="px-4 py-3 border-b bg-gray-50">
                <h2 className="text-sm font-bold text-gray-700">
                  Prediction Overlay - Total Mutation Probability
                </h2>
                <p className="text-xs text-gray-400 mt-0.5">
                  All models shown on the same chart. Hover to compare values at each position.
                </p>
              </div>
              <div style={{ height: "45vh", minHeight: "350px" }} className="p-3">
                <OverlayChart
                  predictions={predictions}
                  modelNames={loadedModels}
                  useLogScale={useLogScale}
                />
              </div>
            </div>

            {/* Difference Chart (only for 2 models) */}
            {loadedModels.length === 2 && (
              <div className="bg-white rounded-xl shadow-sm border overflow-hidden">
                <div className="px-4 py-3 border-b bg-gray-50">
                  <h2 className="text-sm font-bold text-gray-700">
                    Prediction Difference: {loadedModels[0]} vs {loadedModels[1]}
                  </h2>
                  <p className="text-xs text-gray-400 mt-0.5">
                    Blue = first model predicts higher, Red = second model predicts higher.
                  </p>
                </div>
                <div style={{ height: "25vh", minHeight: "200px" }} className="p-3">
                  <DifferenceChart
                    predictions={predictions}
                    modelNames={loadedModels}
                    useLogScale={useLogScale}
                  />
                </div>
              </div>
            )}

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
        {/* STACKED MODE - Models vertically stacked */}
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
                  <div className="flex items-center gap-2 mb-1.5 px-2">
                    <div className="w-3 h-3 rounded-full" style={{ backgroundColor: mc.border }} />
                    <span className="font-bold text-sm" style={{ color: mc.label }}>{name}</span>
                    {isLoading && <span className="text-xs text-gray-400 animate-pulse ml-2">Loading...</span>}
                    {error && <span className="text-xs text-red-500 ml-2">{error}</span>}
                  </div>

                  <div className="bg-white rounded-xl shadow-sm border overflow-hidden"
                       style={{ borderLeftColor: mc.border, borderLeftWidth: "4px" }}>
                    {isLoading ? (
                      <div className="flex items-center justify-center py-16">
                        <div className="animate-spin rounded-full h-8 w-8 border-2 border-blue-500 border-t-transparent" />
                      </div>
                    ) : pred?.genomeDataRaw ? (
                      <div className="flex">
                        <div className="flex-1" style={{ height: "40vh", maxHeight: "400px", overflow: "hidden" }}>
                          <GenomeChart
                            key={`genome-chart-${modelId}`}
                            genomeData={pred.genomeDataRaw}
                            genomeSequence={pred.genomeSequence}
                            onZoomSync={zoomSyncEnabled ? (range) => {
                              zoomSourceRef.current = modelId;
                              setSharedZoomRange({ ...range });
                            } : undefined}
                            syncZoomRange={zoomSyncEnabled && sharedZoomRange && zoomSourceRef.current !== modelId ? sharedZoomRange : undefined}
                          />
                        </div>
                        {showDoughnut && pred.proteinMutationProbs && Object.keys(pred.proteinMutationProbs).length > 0 && (
                          <div className="w-64 flex-shrink-0 flex justify-center items-start pt-4">
                            <DoughnutChart data={pred.proteinMutationProbs} />
                          </div>
                        )}
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
        {/* SIDE BY SIDE MODE - Models horizontally */}
        {/* ============================================ */}
        {viewMode === "sidebyside" && (
          <div className={`grid gap-4 ${models.length <= 2 ? "grid-cols-2" : models.length === 3 ? "grid-cols-3" : "grid-cols-2"}`}>
            {models.map((modelId, idx) => {
              const name = modelNames[idx];
              const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
              const pred = predictions[name];
              const isLoading = loading[name];
              const error = errors[name];

              return (
                <div key={modelId} className="flex flex-col">
                  <div className="flex items-center gap-2 mb-1.5 px-2">
                    <div className="w-3 h-3 rounded-full" style={{ backgroundColor: mc.border }} />
                    <span className="font-bold text-sm" style={{ color: mc.label }}>{name}</span>
                    {isLoading && <span className="text-xs text-gray-400 animate-pulse ml-2">Loading...</span>}
                    {error && <span className="text-xs text-red-500 ml-2">{error}</span>}
                  </div>

                  <div className="bg-white rounded-xl shadow-sm border overflow-hidden flex-1"
                       style={{ borderTopColor: mc.border, borderTopWidth: "3px" }}>
                    {isLoading ? (
                      <div className="flex items-center justify-center py-16">
                        <div className="animate-spin rounded-full h-8 w-8 border-2 border-blue-500 border-t-transparent" />
                      </div>
                    ) : pred?.genomeDataRaw ? (
                      <div>
                        <div style={{ height: "35vh", minHeight: "280px", overflow: "hidden" }}>
                          <GenomeChart
                            key={`genome-chart-side-${modelId}`}
                            genomeData={pred.genomeDataRaw}
                            genomeSequence={pred.genomeSequence}
                            onZoomSync={zoomSyncEnabled ? (range) => {
                              zoomSourceRef.current = modelId;
                              setSharedZoomRange({ ...range });
                            } : undefined}
                            syncZoomRange={zoomSyncEnabled && sharedZoomRange && zoomSourceRef.current !== modelId ? sharedZoomRange : undefined}
                          />
                        </div>
                        {showDoughnut && pred.proteinMutationProbs && Object.keys(pred.proteinMutationProbs).length > 0 && (
                          <div className="flex justify-center py-3 border-t">
                            <DoughnutChart data={pred.proteinMutationProbs} />
                          </div>
                        )}
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
      </div>
    </div>
  );
};

export default CompareModels;
