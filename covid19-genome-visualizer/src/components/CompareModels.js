import React, { useState, useEffect } from "react";
import { useSearchParams, useNavigate } from "react-router-dom";
import GenomeChart from "./Recharts";
import DoughnutChart from "./DoughnutChart";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

const MODEL_COLORS = [
  { border: "#2563EB", bg: "#EFF6FF", label: "#1D4ED8" },
  { border: "#DC2626", bg: "#FEF2F2", label: "#B91C1C" },
  { border: "#059669", bg: "#ECFDF5", label: "#047857" },
  { border: "#D97706", bg: "#FFFBEB", label: "#B45309" },
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
  const allDone = modelNames.every((n) => !loading[n]);

  return (
    <div className="min-h-screen bg-[#f6f7f9]">
      {/* Header */}
      <div className="bg-white border-b shadow-sm sticky top-0 z-10">
        <div className="max-w-[1600px] mx-auto px-6 py-3 flex items-center justify-between">
          <div className="flex items-center gap-4">
            <button onClick={() => navigate("/benchmark")} className="text-gray-400 hover:text-gray-700 text-sm">
              ← Back to Benchmark
            </button>
            <h1 className="font-bold text-gray-800 text-lg">Visual Model Comparison</h1>
          </div>
          <div className="text-xs text-gray-400">
            {nodeId.substring(0, 40)}... · Day {elapsedDay} · {proteinRegion || "Full Genome"}
          </div>
        </div>
      </div>

      {/* Models stacked */}
      <div className="max-w-[1600px] mx-auto px-4 py-4 space-y-6">
        {models.map((modelId, idx) => {
          const name = modelNames[idx];
          const mc = MODEL_COLORS[idx % MODEL_COLORS.length];
          const pred = predictions[name];
          const isLoading = loading[name];
          const error = errors[name];

          return (
            <div key={modelId}>
              {/* Model Label */}
              <div className="flex items-center gap-2 mb-2 px-2">
                <div className="w-3 h-3 rounded-full" style={{ backgroundColor: mc.border }} />
                <span className="font-bold text-lg" style={{ color: mc.label }}>{name}</span>
                {isLoading && (
                  <span className="text-sm text-gray-400 animate-pulse ml-2">Loading...</span>
                )}
                {error && <span className="text-sm text-red-500 ml-2">{error}</span>}
              </div>

              {/* Chart - SAME as main prediction page */}
              <div className="bg-white rounded-xl shadow-sm border overflow-hidden"
                   style={{ borderLeftColor: mc.border, borderLeftWidth: "4px" }}>
                {isLoading ? (
                  <div className="flex items-center justify-center py-24">
                    <div className="animate-spin rounded-full h-10 w-10 border-3 border-blue-500 border-t-transparent" />
                  </div>
                ) : pred?.genomeDataRaw ? (
                  <div className="flex">
                    {/* GenomeChart - same component as main page, isolated with unique key */}
                    <div className="flex-1" style={{ height: "70vh", maxHeight: "70vh", overflow: "hidden" }}>
                      <GenomeChart
                        key={`genome-chart-${modelId}`}
                        genomeData={pred.genomeDataRaw}
                        genomeSequence={pred.genomeSequence}
                      />
                    </div>

                    {/* DoughnutChart - same as main page */}
                    {pred.proteinMutationProbs && Object.keys(pred.proteinMutationProbs).length > 0 && (
                      <div className="w-80 flex-shrink-0 flex justify-center items-start pt-8">
                        <DoughnutChart data={pred.proteinMutationProbs} />
                      </div>
                    )}
                  </div>
                ) : error ? (
                  <div className="py-16 text-center text-red-400">{error}</div>
                ) : null}
              </div>
            </div>
          );
        })}
      </div>

      {/* Loading indicator */}
      {!allDone && (
        <div className="text-center py-8 text-gray-400 animate-pulse">
          Loading predictions...
        </div>
      )}
    </div>
  );
};

export default CompareModels;