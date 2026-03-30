import React, { useState, useEffect, useRef } from "react";
import Chart from "chart.js/auto";

const API_URL = process.env.REACT_APP_API_URL || "http://127.0.0.1:8000";

// ============================================
// PROTEIN SIGNIFICANCE BAR
// ============================================
const SignificanceBar = ({ value, label }) => {
  const pct = Math.round((value || 0) * 100);
  const color = pct >= 70 ? "#EF4444" : pct >= 40 ? "#F59E0B" : "#10B981";
  return (
    <div className="flex items-center gap-2 text-xs">
      <span className="w-16 text-gray-600 truncate" title={label}>{label}</span>
      <div className="flex-1 bg-gray-200 rounded-full h-2">
        <div className="h-2 rounded-full transition-all" style={{ width: `${pct}%`, backgroundColor: color }} />
      </div>
      <span className="w-8 text-right font-mono text-gray-500">{pct}%</span>
    </div>
  );
};

// ============================================
// CASE TIMELINE CHART
// ============================================
const CaseTimelineChart = ({ monthlyData, title }) => {
  const chartRef = useRef(null);
  const chartInstance = useRef(null);

  useEffect(() => {
    if (!chartRef.current || !monthlyData) return;
    if (chartInstance.current) chartInstance.current.destroy();

    const labels = Object.keys(monthlyData).sort();
    const values = labels.map((k) => monthlyData[k]);

    chartInstance.current = new Chart(chartRef.current, {
      type: "line",
      data: {
        labels,
        datasets: [{
          label: title || "Cases",
          data: values,
          borderColor: "#3B82F6",
          backgroundColor: "rgba(59,130,246,0.1)",
          fill: true,
          tension: 0.3,
          pointRadius: 3,
        }],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          y: { beginAtZero: true, ticks: { callback: (v) => v >= 1000000 ? `${(v / 1000000).toFixed(1)}M` : v >= 1000 ? `${(v / 1000).toFixed(0)}K` : v } },
          x: { ticks: { maxRotation: 45 } },
        },
        plugins: { legend: { display: false }, title: { display: true, text: title || "Timeline", font: { size: 12 } } },
      },
    });
    return () => { if (chartInstance.current) chartInstance.current.destroy(); };
  }, [monthlyData, title]);

  return <canvas ref={chartRef} />;
};

// ============================================
// MAIN CONTEXT OVERLAY COMPONENT
// ============================================
const ContextOverlay = ({ nodeId, elapsedDay, selectedProteinRegion }) => {
  const [isOpen, setIsOpen] = useState(false);
  const [contextData, setContextData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

  // Fetch context when opened or params change
  useEffect(() => {
    if (!isOpen || !nodeId) return;

    const fetchContext = async () => {
      setLoading(true);
      setError("");
      try {
        const res = await fetch(
          `${API_URL}/api/context/variant/?nodeId=${encodeURIComponent(nodeId)}&elapsedDay=${elapsedDay || 60}`
        );
        if (res.ok) {
          const data = await res.json();
          setContextData(data.context);
        } else {
          setError("Failed to load context data");
        }
      } catch (err) {
        setError(`Network error: ${err.message}`);
      } finally {
        setLoading(false);
      }
    };
    fetchContext();
  }, [isOpen, nodeId, elapsedDay]);

  const lineage = contextData?.lineage_info;
  const timeCtx = contextData?.time_context;
  const proteinSig = contextData?.protein_significance;

  return (
    <>
      {/* Toggle Button */}
      <button
        onClick={() => setIsOpen(!isOpen)}
        className={`fixed right-4 top-1/2 -translate-y-1/2 z-50 px-2 py-4 rounded-l-lg shadow-lg transition-all ${
          isOpen ? "bg-blue-600 text-white" : "bg-white text-blue-600 border border-blue-200"
        }`}
        style={{ writingMode: "vertical-rl", textOrientation: "mixed" }}
        title="Toggle Context Overlay"
      >
        {isOpen ? "✕ Close" : "📊 Context"}
      </button>

      {/* Overlay Panel */}
      <div
        className={`fixed right-0 top-0 h-full bg-white shadow-2xl border-l border-gray-200 z-40 transition-transform duration-300 overflow-y-auto ${
          isOpen ? "translate-x-0" : "translate-x-full"
        }`}
        style={{ width: "360px" }}
      >
        <div className="p-4">
          <h2 className="text-lg font-bold text-gray-800 mb-1">Contextual Data</h2>
          <p className="text-xs text-gray-400 mb-4">Epidemiological & clinical context</p>

          {loading && (
            <div className="flex items-center justify-center py-8">
              <div className="animate-spin rounded-full h-8 w-8 border-2 border-blue-500 border-t-transparent" />
            </div>
          )}

          {error && <p className="text-red-500 text-sm p-2 bg-red-50 rounded">{error}</p>}

          {contextData && !loading && (
            <div className="space-y-4">

              {/* Lineage Info */}
              {lineage ? (
                <div className="bg-blue-50 rounded-lg p-3">
                  <h3 className="font-bold text-blue-800 text-sm mb-2">
                    {lineage.who_label || lineage.common_name} ({lineage.lineage})
                  </h3>
                  <div className="grid grid-cols-2 gap-2 text-xs">
                    <div>
                      <span className="text-gray-500">Transmissibility</span>
                      <p className="font-bold text-gray-800">
                        +{((lineage.transmissibility_increase || 0) * 100).toFixed(0)}%
                      </p>
                    </div>
                    <div>
                      <span className="text-gray-500">Severity</span>
                      <p className={`font-bold ${(lineage.severity_increase || 0) > 0 ? "text-red-600" : "text-green-600"}`}>
                        {(lineage.severity_increase || 0) > 0 ? "+" : ""}
                        {((lineage.severity_increase || 0) * 100).toFixed(0)}%
                      </p>
                    </div>
                    <div>
                      <span className="text-gray-500">Spike Mutations</span>
                      <p className="font-bold text-gray-800">{lineage.spike_mutations_count || "N/A"}</p>
                    </div>
                    <div>
                      <span className="text-gray-500">Key Mutations</span>
                      <p className="font-bold text-gray-800 text-[10px]">
                        {(lineage.key_mutations || []).join(", ")}
                      </p>
                    </div>
                  </div>
                </div>
              ) : (
                <div className="bg-gray-50 rounded-lg p-3 text-sm text-gray-500">
                  Lineage could not be identified from this variant. Context data is limited.
                </div>
              )}

              {/* Time Context */}
              {timeCtx && (
                <div className="bg-white border rounded-lg p-3">
                  <h3 className="font-bold text-gray-700 text-sm mb-2">Time Context</h3>
                  <p className="text-xs text-gray-500 mb-2">
                    Target month: <span className="font-mono font-bold">{timeCtx.target_month}</span>
                  </p>
                  {timeCtx.cases_at_time && (
                    <p className="text-xs">
                      Est. cases: <span className="font-bold">{(timeCtx.cases_at_time / 1000).toFixed(0)}K</span>
                    </p>
                  )}
                  {timeCtx.severity_at_time && (
                    <p className="text-xs">
                      Severity index: <span className="font-bold">{timeCtx.severity_at_time.toFixed(2)}</span>
                    </p>
                  )}

                  {/* Case Timeline Chart */}
                  {timeCtx.monthly_cases_timeline && (
                    <div className="mt-2" style={{ height: "140px" }}>
                      <CaseTimelineChart monthlyData={timeCtx.monthly_cases_timeline} title="Monthly Cases" />
                    </div>
                  )}

                  {/* Severity Timeline Chart */}
                  {timeCtx.severity_timeline && (
                    <div className="mt-2" style={{ height: "120px" }}>
                      <CaseTimelineChart monthlyData={timeCtx.severity_timeline} title="Severity Index" />
                    </div>
                  )}
                </div>
              )}

              {/* Protein Region Significance */}
              {proteinSig && (
                <div className="bg-white border rounded-lg p-3">
                  <h3 className="font-bold text-gray-700 text-sm mb-2">Protein Region Significance</h3>
                  <p className="text-[10px] text-gray-400 mb-2">Clinical significance based on literature review</p>
                  <div className="space-y-1.5">
                    {Object.entries(proteinSig)
                      .sort((a, b) => (b[1].clinical_significance || 0) - (a[1].clinical_significance || 0))
                      .map(([protein, info]) => (
                        <div key={protein}>
                          <SignificanceBar
                            value={info.clinical_significance}
                            label={protein}
                          />
                          {info.drug_target && (
                            <span className="text-[9px] ml-18 text-purple-600 font-bold">💊 Drug target</span>
                          )}
                        </div>
                      ))}
                  </div>
                </div>
              )}

              {/* Selected Protein Detail */}
              {selectedProteinRegion && proteinSig?.[selectedProteinRegion] && (
                <div className="bg-amber-50 border border-amber-200 rounded-lg p-3">
                  <h3 className="font-bold text-amber-800 text-sm mb-1">
                    Selected: {selectedProteinRegion}
                  </h3>
                  <p className="text-xs text-gray-600">
                    {proteinSig[selectedProteinRegion].description}
                  </p>
                  {proteinSig[selectedProteinRegion].known_drug_resistance_positions?.length > 0 && (
                    <p className="text-xs mt-1 text-red-600">
                      Known resistance positions: {proteinSig[selectedProteinRegion].known_drug_resistance_positions.join(", ")}
                    </p>
                  )}
                </div>
              )}
            </div>
          )}

          {/* Empty state */}
          {!contextData && !loading && !error && (
            <div className="text-center py-8 text-gray-400">
              <p className="text-sm">Run a prediction first to see contextual data.</p>
            </div>
          )}
        </div>
      </div>
    </>
  );
};

export default ContextOverlay;