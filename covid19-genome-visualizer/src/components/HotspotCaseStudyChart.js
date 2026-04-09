import React, { useRef, useMemo } from "react";
import {
  CartesianGrid,
  ComposedChart,
  Line,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

const formatScore = (value) =>
  value === null || value === undefined ? "—" : Number(value).toFixed(3);

const CaseStudyTooltip = ({
  active,
  payload,
  selectedModel,
  displayNormalization,
}) => {
  if (!active || !payload || payload.length === 0) return null;
  const row = payload[0].payload;
  if (!row) return null;

  const usesDisplayNormalization =
    displayNormalization ===
    "min_max_across_all_spike_sites_in_the_selected_context";
  const displayScoreLabel = usesDisplayNormalization
    ? "Hotspot score"
    : `${selectedModel === "PRIEST" ? "PRIEST" : "Site"} score`;
  const rankingScoreLabel = usesDisplayNormalization ? "Raw score" : "Ranking score";

  return (
    <div className="rounded-2xl border border-slate-200 bg-white px-4 py-3 shadow-lg">
      <p className="text-sm font-semibold text-slate-900">Spike site {row.aa_position}</p>
      <p className="mt-1 text-xs text-slate-600">Rank {row.rank}</p>
      <p className="mt-2 text-sm text-slate-700">
        {displayScoreLabel}:{" "}
        <span className="font-semibold text-slate-900">{formatScore(row.site_score)}</span>
      </p>
      <p className="mt-1 text-sm text-slate-700">
        {rankingScoreLabel}:{" "}
        <span className="font-semibold text-slate-900">{formatScore(row.raw_site_score)}</span>
      </p>
      <p className="mt-1 text-xs text-slate-500">
        Ref AA: {row.reference_aa || "—"} | Codon positions{" "}
        {Array.isArray(row.codon_genome_positions)
          ? row.codon_genome_positions.join(", ")
          : "—"}
      </p>
      <div className="mt-3 flex flex-wrap gap-2 text-[11px]">
        {row.is_top_k ? (
          <span className="rounded-full bg-sky-50 px-2 py-1 font-medium text-sky-700">Hotspot</span>
        ) : null}
        {row.is_comparison_site ?? row.is_omicron_site ? (
          <span className="rounded-full bg-amber-50 px-2 py-1 font-medium text-amber-700">Known variant site</span>
        ) : null}
        {row.is_overlap ? (
          <span className="rounded-full bg-rose-50 px-2 py-1 font-medium text-rose-700">Overlap hit</span>
        ) : null}
        {row.is_proximity_overlap && !row.is_overlap ? (
          <span className="rounded-full bg-teal-50 px-2 py-1 font-medium text-teal-700">Near-hit</span>
        ) : null}
      </div>
    </div>
  );
};

const buildScatterPoints = (scoreSeries, predicate) =>
  scoreSeries.filter(predicate).map((row) => ({ ...row, x: row.aa_position, y: row.site_score }));

const ScatterCircle = ({ cx, cy, fill }) => <circle cx={cx} cy={cy} r={4} fill={fill} />;

const ScatterDiamond = ({ cx, cy, fill }) => (
  <rect x={cx - 4} y={cy - 4} width={8} height={8} fill={fill} transform={`rotate(45 ${cx} ${cy})`} />
);

const buildStarPoints = (cx, cy, outerRadius = 6, innerRadius = 2.8) => {
  const points = [];
  for (let i = 0; i < 10; i++) {
    const angle = (Math.PI / 5) * i - Math.PI / 2;
    const r = i % 2 === 0 ? outerRadius : innerRadius;
    points.push(`${cx + Math.cos(angle) * r},${cy + Math.sin(angle) * r}`);
  }
  return points.join(" ");
};

const ScatterStar = ({ cx, cy, fill }) => <polygon points={buildStarPoints(cx, cy)} fill={fill} />;

const HotspotCaseStudyChart = ({
  scoreSeries,
  showTopK,
  showOmicron,
  showOverlap,
  showProximity,
  selectedModel,
  displayNormalization,
  variantLabel,
  topK,
}) => {
  const chartRef = useRef(null);

  const usesDisplayNormalization =
    displayNormalization === "min_max_across_all_spike_sites_in_the_selected_context";
  const chartTitle = selectedModel === "PRIEST"
    ? "PRIEST Spike Site Prevalence"
    : "Explorer Spike Hotspot Rankings";
  const chartDescription = usesDisplayNormalization
    ? "Hotspot scores across all Spike positions for the selected lineage. Top-K ranked positions are highlighted; known Spike mutation sites for this lineage are shown as an optional reference layer."
    : "PRIEST site prevalence scores across Spike positions for the selected lineage context. Top-K positions and known variant sites are shown as reference overlays.";
  const yAxisLabel = usesDisplayNormalization ? "Hotspot score" : "Site score";

  const chartData = useMemo(
    () =>
      (scoreSeries || []).map((row) => ({
        ...row,
        topKOnlyScore: row.is_top_k && !row.is_overlap ? Number(row.site_score) : null,
        omicronOnlyScore:
          (row.is_comparison_site ?? row.is_omicron_site) && !row.is_overlap
            ? Number(row.site_score)
            : null,
        overlapScore: row.is_overlap ? Number(row.site_score) : null,
      })),
    [scoreSeries]
  );

  const topKOnlyPoints = useMemo(
    () => buildScatterPoints(scoreSeries || [], (r) => r.is_top_k && !r.is_overlap && !r.is_proximity_overlap),
    [scoreSeries]
  );
  const omicronOnlyPoints = useMemo(
    () => buildScatterPoints(scoreSeries || [], (r) => (r.is_comparison_site ?? r.is_omicron_site) && !r.is_overlap),
    [scoreSeries]
  );
  const overlapPoints = useMemo(
    () => buildScatterPoints(scoreSeries || [], (r) => r.is_overlap),
    [scoreSeries]
  );
  const proximityOnlyPoints = useMemo(
    () => buildScatterPoints(scoreSeries || [], (r) => r.is_proximity_overlap && !r.is_overlap),
    [scoreSeries]
  );

  const handleDownloadPng = () => {
    const svgEl = chartRef.current?.querySelector("svg");
    if (!svgEl) return;

    const scale = 2;
    const width = svgEl.clientWidth || 800;
    const height = svgEl.clientHeight || 450;

    const svgClone = svgEl.cloneNode(true);
    svgClone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    svgClone.setAttribute("width", width);
    svgClone.setAttribute("height", height);

    const svgString = new XMLSerializer().serializeToString(svgClone);
    const blob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const img = new Image();

    img.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = width * scale;
      canvas.height = height * scale;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.scale(scale, scale);
      ctx.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);

      const safeVariant = (variantLabel || "variant").replace(/[^a-zA-Z0-9.]/g, "_").replace(/_+/g, "_");
      const safeModel = (selectedModel || "model").replace(/[^a-zA-Z0-9]/g, "_");
      const filename = `covmutex_case_study_${safeVariant}_${safeModel}_k${topK || "K"}.png`;

      const link = document.createElement("a");
      link.download = filename;
      link.href = canvas.toDataURL("image/png");
      link.click();
    };

    img.onerror = () => URL.revokeObjectURL(url);
    img.src = url;
  };

  if (!chartData.length) return null;

  const maxPosition = Math.max(...chartData.map((r) => r.aa_position));

  return (
    <div className="rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <div>
          <h2 className="text-lg font-semibold text-slate-900">{chartTitle}</h2>
          <p className="mt-1 text-sm text-slate-600">{chartDescription}</p>
        </div>
        <button
          type="button"
          onClick={handleDownloadPng}
          className="flex items-center gap-2 rounded-full border border-slate-300 bg-white px-4 py-2 text-sm font-medium text-slate-700 transition hover:bg-slate-50 hover:border-slate-400"
          title="Download chart as PNG"
        >
          <svg
            xmlns="http://www.w3.org/2000/svg"
            viewBox="0 0 20 20"
            fill="currentColor"
            className="h-4 w-4 text-slate-500"
          >
            <path d="M10.75 2.75a.75.75 0 0 0-1.5 0v8.614L6.295 8.235a.75.75 0 1 0-1.09 1.03l4.25 4.5a.75.75 0 0 0 1.09 0l4.25-4.5a.75.75 0 0 0-1.09-1.03l-2.955 3.129V2.75Z" />
            <path d="M3.5 12.75a.75.75 0 0 0-1.5 0v2.5A2.75 2.75 0 0 0 4.75 18h10.5A2.75 2.75 0 0 0 18 15.25v-2.5a.75.75 0 0 0-1.5 0v2.5c0 .69-.56 1.25-1.25 1.25H4.75c-.69 0-1.25-.56-1.25-1.25v-2.5Z" />
          </svg>
          Download PNG
        </button>
      </div>

      <div ref={chartRef} className="mt-5 h-[28rem] w-full" data-testid="retrospective-chart">
        <ResponsiveContainer width="100%" height="100%">
          <ComposedChart data={chartData} margin={{ top: 12, right: 18, bottom: 12, left: 0 }}>
            <CartesianGrid stroke="#e2e8f0" strokeDasharray="4 4" />
            <XAxis
              type="number"
              dataKey="aa_position"
              domain={[1, maxPosition]}
              tick={{ fill: "#64748b", fontSize: 12 }}
              tickLine={false}
              axisLine={{ stroke: "#cbd5e1" }}
              label={{
                value: "Spike amino-acid position",
                position: "insideBottom",
                offset: -6,
                style: { fill: "#475569", fontSize: 12 },
              }}
            />
            <YAxis
              type="number"
              dataKey="site_score"
              domain={[0, 1]}
              tick={{ fill: "#64748b", fontSize: 12 }}
              tickLine={false}
              axisLine={{ stroke: "#cbd5e1" }}
              label={{
                value: yAxisLabel,
                angle: -90,
                position: "insideLeft",
                style: { fill: "#475569", fontSize: 12 },
              }}
            />
            <Tooltip
              content={
                <CaseStudyTooltip selectedModel={selectedModel} displayNormalization={displayNormalization} />
              }
            />
            <Line
              type="linear"
              dataKey="site_score"
              stroke="#0f172a"
              strokeWidth={1.75}
              dot={false}
              isAnimationActive={false}
            />
            {showTopK && topKOnlyPoints.length ? (
              <Scatter data={topKOnlyPoints} fill="#0ea5e9" line={false} shape={<ScatterCircle />} isAnimationActive={false} />
            ) : null}
            {showOmicron && omicronOnlyPoints.length ? (
              <Scatter data={omicronOnlyPoints} fill="#f59e0b" line={false} shape={<ScatterDiamond />} isAnimationActive={false} />
            ) : null}
            {showOverlap && overlapPoints.length ? (
              <Scatter data={overlapPoints} fill="#e11d48" line={false} shape={<ScatterStar />} isAnimationActive={false} />
            ) : null}
            {showOverlap
              ? overlapPoints.map((point) => (
                  <ReferenceLine
                    key={`overlap-line-${point.x}`}
                    x={point.x}
                    stroke="#e11d48"
                    strokeDasharray="4 3"
                    strokeOpacity={0.45}
                    strokeWidth={1.5}
                  />
                ))
              : null}
            {showProximity && proximityOnlyPoints.length ? (
              <Scatter data={proximityOnlyPoints} fill="#14b8a6" line={false} shape={<ScatterCircle />} isAnimationActive={false} />
            ) : null}
            {showProximity
              ? proximityOnlyPoints.map((point) => (
                  <ReferenceLine
                    key={`proximity-line-${point.x}`}
                    x={point.x}
                    stroke="#14b8a6"
                    strokeDasharray="6 4"
                    strokeOpacity={0.3}
                    strokeWidth={1}
                  />
                ))
              : null}
          </ComposedChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
};

export default HotspotCaseStudyChart;
