import React, { useState, useEffect, useRef } from "react";
import { useSearchParams, useNavigate } from "react-router-dom";
import GenomeChart from "./Recharts";
import DoughnutChart from "./DoughnutChart";

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
  const anyLoading = modelNames.some((n) => loading[n]);
  const loadedModels = modelNames.filter((n) => predictions[n]?.genomeDataRaw);
  const hasDoughnuts = loadedModels.some(n =>
    predictions[n]?.proteinMutationProbs && Object.keys(predictions[n].proteinMutationProbs).length > 0
  );

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

          <div className="flex items-center gap-3 flex-wrap">
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

        {/* Charts - Stacked vertically, full width, no side panels */}
        <div className="space-y-4">
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
                      key={`genome-chart-${modelId}`}
                      genomeData={pred.genomeDataRaw}
                      genomeSequence={pred.genomeSequence}
                      compact={true}
                    />
                  ) : error ? (
                    <div className="py-16 text-center text-red-400 text-sm">{error}</div>
                  ) : null}
                </div>
              </div>
            );
          })}
        </div>

        {/* Doughnut Charts - Shared section at the bottom */}
        {hasDoughnuts && loadedModels.length > 0 && (
          <div className="mt-6 bg-white rounded-xl shadow-sm border overflow-hidden">
            <div className="px-5 py-3 border-b bg-gray-50">
              <h2 className="text-sm font-bold text-gray-700">Protein Region Distribution Comparison</h2>
            </div>
            <div className={`grid gap-6 p-6 ${loadedModels.length <= 2 ? "grid-cols-2" : loadedModels.length === 3 ? "grid-cols-3" : "grid-cols-4"}`}>
              {loadedModels.map((name, idx) => {
                const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
                const pred = predictions[name];
                if (!pred?.proteinMutationProbs || Object.keys(pred.proteinMutationProbs).length === 0) return null;
                return (
                  <div key={name} className="text-center">
                    <div className="flex items-center justify-center gap-2 mb-3">
                      <div className="w-3 h-3 rounded-full" style={{ backgroundColor: mc.border }} />
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
