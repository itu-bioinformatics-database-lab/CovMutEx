import React, { useState, useEffect, useRef, useCallback } from "react";
import { useSearchParams, useNavigate } from "react-router-dom";
import Chart from "chart.js/auto";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

const NUC_COLORS = {
  A: { solid: "#DC2626", light: "rgba(220,38,38,0.15)" },
  T: { solid: "#16A34A", light: "rgba(22,163,74,0.15)" },
  G: { solid: "#CA8A04", light: "rgba(202,138,4,0.15)" },
  C: { solid: "#2563EB", light: "rgba(37,99,235,0.15)" },
};
const NUCS = ["A", "T", "G", "C"];

const MODEL_STYLES = [
  { border: "#2563EB", bg: "#EFF6FF", label: "#1D4ED8" },
  { border: "#DC2626", bg: "#FEF2F2", label: "#B91C1C" },
  { border: "#059669", bg: "#ECFDF5", label: "#047857" },
  { border: "#D97706", bg: "#FFFBEB", label: "#B45309" },
];

// ============================================
// UTILS
// ============================================
function movingAvg(arr, windowSize) {
  const result = new Array(arr.length);
  const half = Math.floor(windowSize / 2);
  for (let i = 0; i < arr.length; i++) {
    let sum = 0, count = 0;
    for (let j = Math.max(0, i - half); j <= Math.min(arr.length - 1, i + half); j++) {
      sum += arr[j];
      count++;
    }
    result[i] = sum / count;
  }
  return result;
}

function downsample(arr, targetLen) {
  if (arr.length <= targetLen) return { data: arr, indices: arr.map((_, i) => i) };
  const step = arr.length / targetLen;
  const data = [], indices = [];
  for (let i = 0; i < targetLen; i++) {
    const idx = Math.round(i * step);
    data.push(arr[Math.min(idx, arr.length - 1)]);
    indices.push(idx);
  }
  return { data, indices };
}

// ============================================
// GENOME CHART — clean area chart, no zoom plugin
// ============================================
const GenomeChart = ({ genomeData, viewRange, showNucBreakdown, height = 280 }) => {
  const canvasRef = useRef(null);
  const chartRef = useRef(null);

  useEffect(() => {
    if (!canvasRef.current || !genomeData?.length) return;
    if (chartRef.current) chartRef.current.destroy();

    const numPos = genomeData[0]?.length || 0;
    const [rangeStart, rangeEnd] = viewRange;
    const start = Math.floor(rangeStart * numPos);
    const end = Math.ceil(rangeEnd * numPos);
    const sliceLen = end - start;

    // Compute total mutation probability per position (sum of all nucleotides)
    const totalProb = [];
    for (let i = start; i < end; i++) {
      let sum = 0;
      for (let n = 0; n < 4; n++) sum += genomeData[n]?.[i] || 0;
      totalProb.push(Math.min(sum, 1));
    }

    // Smoothing window based on visible range
    const windowSize = Math.max(3, Math.floor(sliceLen / 200));
    const smoothed = movingAvg(totalProb, windowSize);

    // Downsample to ~800 points
    const { data: dsTotal, indices } = downsample(smoothed, Math.min(800, sliceLen));
    const labels = indices.map((idx) => start + idx + 1);

    const datasets = [];

    if (showNucBreakdown) {
      // Show only total (sum) — less useful since A+T+G+C ≈ 1 always
      datasets.push({
        label: "Total Mutation Prob.",
        data: dsTotal,
        borderColor: "#6366F1",
        backgroundColor: "rgba(99,102,241,0.12)",
        borderWidth: 1.8,
        pointRadius: 0,
        fill: true,
        tension: 0.3,
      });
    } else {
      // Default: show individual nucleotides
      NUCS.forEach((nuc, nucIdx) => {
        const raw = [];
        for (let i = start; i < end; i++) raw.push(genomeData[nucIdx]?.[i] || 0);
        const sm = movingAvg(raw, windowSize);
        const { data: ds } = downsample(sm, Math.min(800, sliceLen));

        datasets.push({
          label: nuc,
          data: ds,
          borderColor: NUC_COLORS[nuc].solid,
          backgroundColor: NUC_COLORS[nuc].light,
          borderWidth: 1.5,
          pointRadius: 0,
          fill: false,
          tension: 0.3,
          order: nucIdx + 1,
        });
      });
    }

    chartRef.current = new Chart(canvasRef.current, {
      type: "line",
      data: { labels, datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false,
        interaction: { mode: "index", intersect: false },
        layout: { padding: { top: 6, right: 12, bottom: 2, left: 4 } },
        scales: {
          x: {
            grid: { display: false },
            ticks: {
              maxTicksLimit: 10,
              font: { size: 10 },
              color: "#94A3B8",
              callback: function (val) {
                const v = this.getLabelForValue(val);
                return Number(v) >= 1000 ? `${(Number(v) / 1000).toFixed(1)}k` : v;
              },
            },
            title: { display: true, text: "Genome Position", font: { size: 10 }, color: "#94A3B8" },
          },
          y: {
            type: "logarithmic",
            min: 0.001,
            max: 1,
            grid: { color: "rgba(0,0,0,0.04)" },
            ticks: {
              font: { size: 9 },
              color: "#94A3B8",
              callback: (v) => {
                if (v === 1) return "1.0";
                if (v === 0.1) return "0.1";
                if (v === 0.01) return "0.01";
                if (v === 0.001) return "0.001";
                return "";
              },
            },
            title: { display: true, text: "Mutation Prob. (log)", font: { size: 10 }, color: "#64748B" },
          },
        },
        plugins: {
          legend: { display: false },
          tooltip: {
            backgroundColor: "rgba(15,23,42,0.92)",
            titleFont: { size: 11 },
            bodyFont: { size: 10 },
            callbacks: {
              title: (items) => `Position ${items[0]?.label || ""}`,
              label: (ctx) => ` ${ctx.dataset.label}: ${ctx.parsed.y?.toFixed(5)}`,
            },
          },
        },
      },
    });

    return () => { if (chartRef.current) chartRef.current.destroy(); };
  }, [genomeData, viewRange, showNucBreakdown]);

  return <canvas ref={canvasRef} style={{ height: `${height}px` }} />;
};

// ============================================
// RANGE SLIDER
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
        if (dragging.current === "start") return [Math.min(pct, prev[1] - 0.02), prev[1]];
        if (dragging.current === "end") return [prev[0], Math.max(pct, prev[0] + 0.02)];
        return prev;
      });
    };
    const onUp = () => { dragging.current = null; window.removeEventListener("mousemove", onMove); window.removeEventListener("mouseup", onUp); };
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
  }, [onChange]);

  const leftPct = start * 100;
  const widthPct = (end - start) * 100;

  return (
    <div className="relative h-8 select-none">
      <div ref={trackRef} className="absolute inset-x-0 top-3 h-2 bg-slate-200 rounded-full">
        {/* Selected range */}
        <div
          className="absolute h-full bg-indigo-400 rounded-full"
          style={{ left: `${leftPct}%`, width: `${widthPct}%` }}
        />
        {/* Left handle */}
        <div
          className="absolute top-1/2 -translate-y-1/2 w-4 h-4 bg-white border-2 border-indigo-500 rounded-full cursor-ew-resize shadow-sm hover:scale-110 transition-transform"
          style={{ left: `calc(${leftPct}% - 8px)` }}
          onMouseDown={(e) => handleMouseDown("start", e)}
        />
        {/* Right handle */}
        <div
          className="absolute top-1/2 -translate-y-1/2 w-4 h-4 bg-white border-2 border-indigo-500 rounded-full cursor-ew-resize shadow-sm hover:scale-110 transition-transform"
          style={{ left: `calc(${(start + (end - start)) * 100}% - 8px)` }}
          onMouseDown={(e) => handleMouseDown("end", e)}
        />
      </div>
      {/* Labels */}
      <div className="absolute -bottom-3 left-0 text-[10px] text-slate-400">{Math.round(start * 29903) + 1}</div>
      <div className="absolute -bottom-3 right-0 text-[10px] text-slate-400">{Math.round(end * 29903)}</div>
    </div>
  );
};

// ============================================
// MAIN
// ============================================
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

  // Shared controls
  const [viewRange, setViewRange] = useState([0, 1]);
  const [showNucBreakdown, setShowNucBreakdown] = useState(false);

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
          let raw = null;
          if (data.genomeData) {
            const gd = data.genomeData;
            if (Array.isArray(gd) && Array.isArray(gd[0]) && typeof gd[0][0] === "number") raw = gd;
            else if (Array.isArray(gd) && gd[0]?.mutationPoss) {
              raw = [gd.map(d => d.mutationPoss?.A || 0), gd.map(d => d.mutationPoss?.T || 0),
                     gd.map(d => d.mutationPoss?.G || 0), gd.map(d => d.mutationPoss?.C || 0)];
            } else if (Array.isArray(gd) && typeof gd[0] === "number") raw = [gd, gd, gd, gd];
          }
          setPredictions(p => ({ ...p, [displayName]: { ...data, genomeDataRaw: raw, proteinMutationProbs: data.protein_mutation_probs } }));
        } else {
          const err = await res.json().catch(() => ({}));
          setErrors(p => ({ ...p, [displayName]: err.error || `HTTP ${res.status}` }));
        }
      } catch (e) {
        setErrors(p => ({ ...p, [displayName]: e.message }));
      } finally {
        setLoading(p => ({ ...p, [displayName]: false }));
      }
    });
  }, []);

  const modelNames = models.map(m => m.startsWith("uploaded:") ? m.replace("uploaded:", "") : m);
  const allDone = modelNames.every(n => !loading[n]);
  const regions = ["ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"];

  return (
    <div className="min-h-screen bg-slate-50">
      {/* Header */}
      <div className="bg-white border-b sticky top-0 z-10 shadow-sm">
        <div className="max-w-[1500px] mx-auto px-6 py-3 flex items-center justify-between">
          <div className="flex items-center gap-4">
            <button onClick={() => navigate("/benchmark")} className="text-slate-400 hover:text-slate-700 text-sm">← Back</button>
            <h1 className="font-bold text-slate-800 text-lg">Visual Comparison</h1>
          </div>
          <div className="flex items-center gap-4">
            {/* Nuc breakdown toggle */}
            <label className="flex items-center gap-2 cursor-pointer">
              <span className="text-xs text-slate-500">Show total only</span>
              <div className="relative" onClick={() => setShowNucBreakdown(v => !v)}>
                <div className={`w-9 h-5 rounded-full transition-colors ${showNucBreakdown ? "bg-indigo-500" : "bg-slate-300"}`} />
                <div className={`absolute top-0.5 w-4 h-4 bg-white rounded-full shadow transition-transform ${showNucBreakdown ? "translate-x-4" : "translate-x-0.5"}`} />
              </div>
            </label>
            {!showNucBreakdown && (
              <div className="flex items-center gap-2">
                {NUCS.map(n => (
                  <div key={n} className="flex items-center gap-1">
                    <div className="w-3 h-1.5 rounded-sm" style={{ backgroundColor: NUC_COLORS[n].solid }} />
                    <span className="text-[11px] font-bold text-slate-500">{n}</span>
                  </div>
                ))}
              </div>
            )}
            <button onClick={() => setViewRange([0, 1])} className="px-3 py-1.5 bg-slate-100 hover:bg-slate-200 rounded-lg text-xs font-medium">
              Reset View
            </button>
          </div>
        </div>
      </div>

      {/* Info + Range Slider */}
      <div className="max-w-[1500px] mx-auto px-6 pt-4 pb-2">
        <div className="flex items-center gap-2 text-xs text-slate-400 mb-4">
          <span className="bg-slate-100 px-2 py-0.5 rounded font-mono text-[11px]">
            {nodeId.length > 55 ? nodeId.substring(0, 55) + "..." : nodeId}
          </span>
          <span>· Day {elapsedDay} · {proteinRegion || "Full Genome"}</span>
        </div>

        {/* Shared Range Slider */}
        <div className="bg-white rounded-xl p-4 shadow-sm border border-slate-100 mb-4">
          <div className="flex items-center justify-between mb-2">
            <span className="text-xs font-semibold text-slate-600">Genome Range</span>
            <span className="text-xs text-slate-400 font-mono">
              {Math.round(viewRange[0] * 29903) + 1} — {Math.round(viewRange[1] * 29903)}
              {" "}({Math.round((viewRange[1] - viewRange[0]) * 29903).toLocaleString()} positions)
            </span>
          </div>
          <RangeSlider value={viewRange} onChange={setViewRange} />
        </div>
      </div>

      {/* Charts */}
      <div className="max-w-[1500px] mx-auto px-6 space-y-3 pb-4">
        {models.map((modelId, idx) => {
          const name = modelNames[idx];
          const mc = MODEL_STYLES[idx % MODEL_STYLES.length];
          const pred = predictions[name];
          const isLoading = loading[name];
          const error = errors[name];

          return (
            <div key={modelId} className="rounded-xl overflow-hidden bg-white" style={{ border: `1px solid ${mc.border}18` }}>
              <div className="flex items-center gap-2.5 px-5 py-2" style={{ background: mc.bg, borderBottom: `2px solid ${mc.border}20` }}>
                <div className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: mc.border }} />
                <span className="font-bold text-sm" style={{ color: mc.label }}>{name}</span>
                {isLoading && (
                  <div className="flex items-center gap-2 ml-auto">
                    <div className="animate-spin h-3 w-3 border-2 border-t-transparent rounded-full" style={{ borderColor: mc.border }} />
                    <span className="text-xs" style={{ color: mc.border }}>Computing...</span>
                  </div>
                )}
              </div>
              <div style={{ height: "300px", padding: "4px 8px" }}>
                {isLoading ? (
                  <div className="flex items-center justify-center h-full animate-pulse text-slate-300 text-sm">Running prediction...</div>
                ) : pred?.genomeDataRaw ? (
                  <GenomeChart genomeData={pred.genomeDataRaw} viewRange={viewRange} showNucBreakdown={showNucBreakdown} height={290} />
                ) : error ? (
                  <div className="flex items-center justify-center h-full text-red-300 text-sm">{error}</div>
                ) : null}
              </div>
            </div>
          );
        })}
      </div>

      {/* Protein Table */}
      {allDone && modelNames.some(n => predictions[n]?.proteinMutationProbs) && (
        <div className="max-w-[1500px] mx-auto px-6 pb-8">
          <div className="bg-white rounded-xl shadow-sm border border-slate-100 overflow-hidden">
            <div className="px-5 py-3 border-b bg-slate-50">
              <h2 className="font-bold text-slate-700 text-sm">Protein Region Comparison</h2>
            </div>
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-slate-100">
                  <th className="text-left px-5 py-2 text-xs font-semibold text-slate-500 uppercase w-28">Region</th>
                  {modelNames.map((name, idx) => (
                    <th key={name} className="text-right px-5 py-2">
                      <span className="text-xs font-semibold" style={{ color: MODEL_STYLES[idx % MODEL_STYLES.length].label }}>{name}</span>
                    </th>
                  ))}
                  {modelNames.length === 2 && <th className="text-right px-5 py-2 text-xs text-slate-400">Diff</th>}
                </tr>
              </thead>
              <tbody>
                {regions.map(region => {
                  const vals = modelNames.map(n => { const p = predictions[n]?.proteinMutationProbs?.[region]; return typeof p === "number" ? p : null; });
                  const valid = vals.filter(v => v !== null);
                  const mx = valid.length > 1 ? Math.max(...valid) : null;
                  const mn = valid.length > 1 ? Math.min(...valid) : null;
                  const diff = vals.length === 2 && vals[0] != null && vals[1] != null ? ((vals[0] - vals[1]) * 100).toFixed(2) : null;

                  return (
                    <tr key={region} className="border-b border-slate-50 hover:bg-slate-50/50">
                      <td className="px-5 py-2 font-semibold text-slate-700">{region}</td>
                      {vals.map((v, i) => (
                        <td key={i} className="text-right px-5 py-2 font-mono text-sm">
                          {v != null ? (
                            <span className={v === mx && valid.length > 1 ? "text-red-600 font-bold" : v === mn && valid.length > 1 ? "text-emerald-600 font-bold" : "text-slate-600"}>
                              {(v * 100).toFixed(2)}%
                            </span>
                          ) : "—"}
                        </td>
                      ))}
                      {modelNames.length === 2 && (
                        <td className="text-right px-5 py-2 font-mono text-xs">
                          {diff != null ? (
                            <span className={Number(diff) > 0 ? "text-blue-500" : Number(diff) < 0 ? "text-orange-500" : "text-slate-400"}>
                              {Number(diff) > 0 ? "+" : ""}{diff}pp
                            </span>
                          ) : "—"}
                        </td>
                      )}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {!allDone && <div className="text-center py-8 text-sm text-slate-400 animate-pulse">Loading...</div>}
    </div>
  );
};

export default CompareModels;