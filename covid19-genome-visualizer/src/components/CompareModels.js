import React, { useState, useEffect, useCallback } from "react";
import { useSearchParams, useNavigate } from "react-router-dom";
import GenomeChart from "./Recharts";
import DoughnutChart from "./DoughnutChart";
import { proteinRegions } from "../data/proteinRegions";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

const MODEL_COLORS = [
  { border: "#2563EB", bg: "rgba(37, 99, 235, 0.15)", label: "#1D4ED8" },
  { border: "#DC2626", bg: "rgba(220, 38, 38, 0.15)", label: "#B91C1C" },
  { border: "#059669", bg: "rgba(5, 150, 105, 0.15)", label: "#047857" },
  { border: "#D97706", bg: "rgba(217, 119, 6, 0.15)", label: "#B45309" },
];

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

  // Layout mode
  const [layoutMode, setLayoutMode] = useState("stacked");

  // Zoom sync
  const [syncZoomRange, setSyncZoomRange] = useState(null);
  const [zoomSyncEnabled, setZoomSyncEnabled] = useState(true);
  const [syncSource, setSyncSource] = useState(null);

  // Shared protein region focus - controls all charts simultaneously
  const [sharedFocusedProtein, setSharedFocusedProtein] = useState(null);

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

  const createZoomSyncHandler = useCallback((modelName) => {
    return (range) => {
      if (!zoomSyncEnabled) return;
      setSyncSource(modelName);
      setSyncZoomRange({ ...range, _ts: Date.now() });
    };
  }, [zoomSyncEnabled]);

  // Handle protein region click - toggle focus
  const handleProteinClick = (proteinName) => {
    setSharedFocusedProtein((prev) => (prev === proteinName ? null : proteinName));
  };

  const modelNames = models.map((m) => (m.startsWith("uploaded:") ? m.replace("uploaded:", "") : m));
  const anyLoading = modelNames.some((n) => loading[n]);
  const loadedModels = modelNames.filter((n) => predictions[n]?.genomeDataRaw);
  const hasDoughnuts = loadedModels.some(n =>
    predictions[n]?.proteinMutationProbs && Object.keys(predictions[n].proteinMutationProbs).length > 0
  );

  const getGridClass = () => {
    if (layoutMode === "sidebyside") {
      return "grid grid-cols-2 gap-4";
    }
    return "space-y-4";
  };

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

          <div className="flex items-center gap-4 flex-wrap">
            {/* Layout Toggle */}
            <div className="flex items-center bg-gray-100 rounded-lg p-0.5">
              <button
                onClick={() => setLayoutMode("stacked")}
                className={`px-3 py-1.5 text-xs font-medium rounded-md transition-all ${
                  layoutMode === "stacked" ? "bg-white text-gray-800 shadow-sm" : "text-gray-500 hover:text-gray-700"
                }`}
              >
                Stacked
              </button>
              <button
                onClick={() => setLayoutMode("sidebyside")}
                className={`px-3 py-1.5 text-xs font-medium rounded-md transition-all ${
                  layoutMode === "sidebyside" ? "bg-white text-gray-800 shadow-sm" : "text-gray-500 hover:text-gray-700"
                }`}
              >
                Side by Side
              </button>
            </div>

            {/* Zoom Sync Toggle */}
            <button
              onClick={() => setZoomSyncEnabled((v) => !v)}
              className={`flex items-center gap-1.5 px-3 py-1.5 text-xs font-medium rounded-lg border transition-all ${
                zoomSyncEnabled ? "bg-blue-50 border-blue-300 text-blue-700" : "bg-gray-50 border-gray-300 text-gray-500"
              }`}
            >
              <svg xmlns="http://www.w3.org/2000/svg" className="h-3.5 w-3.5" viewBox="0 0 20 20" fill="currentColor">
                <path d="M8 5a1 1 0 011-1h2a1 1 0 110 2H9a1 1 0 01-1-1zM8 15a1 1 0 011-1h2a1 1 0 110 2H9a1 1 0 01-1-1z" />
                <path fillRule="evenodd" d="M10 2a1 1 0 011 1v2.586l1.707-1.293a1 1 0 011.286 1.414L12 7.414V12.586l1.993 1.707a1 1 0 01-1.286 1.414L11 14.414V17a1 1 0 11-2 0v-2.586l-1.993 1.293a1 1 0 01-1.286-1.414L8 12.586V7.414L6.007 5.707a1 1 0 011.286-1.414L9 5.586V3a1 1 0 011-1z" clipRule="evenodd" />
              </svg>
              Zoom Sync {zoomSyncEnabled ? "ON" : "OFF"}
            </button>

            {/* Model Legend */}
            <div className="flex items-center gap-3">
              {modelNames.map((name, idx) => {
                const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
                return (
                  <div key={name} className="flex items-center gap-1.5">
                    <div className="w-3 h-3 rounded-full" style={{ backgroundColor: mc.border }} />
                    <span className="text-sm font-bold" style={{ color: mc.label }}>{name}</span>
                  </div>
                );
              })}
            </div>

            <div className="text-xs text-gray-400 ml-2">
              {nodeId.substring(0, 30)}... | Day {elapsedDay} | {proteinRegion || "Full Genome"}
            </div>
          </div>
        </div>
      </div>

      <div className="max-w-[1800px] mx-auto px-4 py-4">
        {/* Loading */}
        {anyLoading && (
          <div className="text-center py-8 text-gray-400 animate-pulse flex items-center justify-center gap-2">
            <div className="animate-spin rounded-full h-5 w-5 border-2 border-blue-500 border-t-transparent" />
            Loading predictions... ({loadedModels.length}/{modelNames.length})
          </div>
        )}

        {/* Protein Regions - Clickable Shared Panel (above charts) */}
        {loadedModels.length > 0 && (
          <div className="mb-4 bg-white rounded-xl shadow-sm border overflow-hidden">
            <div className="px-4 py-2.5 border-b bg-gray-50 flex items-center justify-between">
              <div>
                <h2 className="text-sm font-bold text-gray-700">Protein Regions</h2>
                <p className="text-xs text-gray-400">Click a region to zoom all charts to that protein</p>
              </div>
              {sharedFocusedProtein && (
                <button
                  onClick={() => setSharedFocusedProtein(null)}
                  className="text-xs bg-red-50 text-red-600 border border-red-200 px-3 py-1 rounded-lg hover:bg-red-100 transition-colors font-medium"
                >
                  Reset View
                </button>
              )}
            </div>
            <div className="p-3">
              <div className="flex flex-wrap gap-1.5">
                {Object.entries(proteinRegions).map(([name, range]) => {
                  const isActive = sharedFocusedProtein === name;
                  return (
                    <button
                      key={name}
                      onClick={() => handleProteinClick(name)}
                      className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg border text-sm font-medium transition-all cursor-pointer ${
                        isActive
                          ? "border-2 shadow-md scale-105"
                          : "border-gray-200 bg-gray-50 hover:bg-gray-100 hover:border-gray-300"
                      }`}
                      style={isActive ? {
                        backgroundColor: (proteinRegionColorMap[name] || "#ccc") + "30",
                        borderColor: proteinRegionColorMap[name] || "#ccc",
                      } : {}}
                    >
                      <div
                        className="w-2.5 h-2.5 rounded-full flex-shrink-0"
                        style={{ backgroundColor: proteinRegionColorMap[name] || "#ccc" }}
                      />
                      <span className={isActive ? "text-gray-900" : "text-gray-700"}>{name}</span>
                      <span className="text-xs text-gray-400 hidden sm:inline">{range}</span>
                    </button>
                  );
                })}
              </div>
            </div>
          </div>
        )}

        {/* Charts */}
        <div className={getGridClass()}>
          {models.map((modelId, idx) => {
            const name = modelNames[idx];
            const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
            const pred = predictions[name];
            const isLoading = loading[name];
            const error = errors[name];

            return (
              <div key={modelId}>
                {/* Model Label */}
                <div className="flex items-center gap-2 mb-1.5 px-2">
                  <div className="w-3 h-3 rounded-full" style={{ backgroundColor: mc.border }} />
                  <span className="font-bold text-sm" style={{ color: mc.label }}>{name}</span>
                  {isLoading && <span className="text-xs text-gray-400 animate-pulse ml-2">Loading...</span>}
                  {error && <span className="text-xs text-red-500 ml-2">{error}</span>}
                </div>

                {/* Chart Card */}
                <div className="bg-white rounded-xl shadow-sm border overflow-hidden"
                     style={{ borderLeftColor: mc.border, borderLeftWidth: "4px" }}>
                  {isLoading ? (
                    <div className="flex items-center justify-center py-20">
                      <div className="animate-spin rounded-full h-8 w-8 border-2 border-blue-500 border-t-transparent" />
                    </div>
                  ) : pred?.genomeDataRaw ? (
                    <GenomeChart
                      key={`genome-chart-${modelId}-${sharedFocusedProtein || 'full'}`}
                      genomeData={pred.genomeDataRaw}
                      genomeSequence={pred.genomeSequence}
                      compact={true}
                      onZoomSync={zoomSyncEnabled ? createZoomSyncHandler(name) : undefined}
                      syncZoomRange={zoomSyncEnabled && syncSource !== name ? syncZoomRange : undefined}
                      externalFocusedProtein={sharedFocusedProtein}
                    />
                  ) : error ? (
                    <div className="py-16 text-center text-red-400 text-sm">{error}</div>
                  ) : null}
                </div>
              </div>
            );
          })}
        </div>

        {/* Protein Region Mutation Distribution - Shared Section */}
        {hasDoughnuts && loadedModels.length > 0 && (
          <div className="mt-6 bg-white rounded-xl shadow-sm border overflow-hidden">
            <div className="px-5 py-3 border-b bg-gray-50">
              <h2 className="text-sm font-bold text-gray-700">Mutation Distribution by Protein Region</h2>
              <p className="text-xs text-gray-400 mt-0.5">
                Each doughnut shows mutation probability distribution. Use the Normalize switch for per-base density comparison.
              </p>
            </div>
            <div className={`grid gap-4 p-4 ${loadedModels.length <= 2 ? "grid-cols-1 md:grid-cols-2" : loadedModels.length === 3 ? "grid-cols-1 md:grid-cols-3" : "grid-cols-1 md:grid-cols-2 lg:grid-cols-4"}`}>
              {loadedModels.map((name, idx) => {
                const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
                const pred = predictions[name];
                if (!pred?.proteinMutationProbs || Object.keys(pred.proteinMutationProbs).length === 0) return null;
                return (
                  <div key={name}>
                    <div className="flex items-center justify-center gap-2 mb-2">
                      <div className="w-3 h-3 rounded-full" style={{ backgroundColor: mc.border }} />
                      <span className="text-sm font-bold" style={{ color: mc.label }}>{name}</span>
                    </div>
                    <div className="flex justify-center">
                      <DoughnutChart data={pred.proteinMutationProbs} />
                    </div>
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
