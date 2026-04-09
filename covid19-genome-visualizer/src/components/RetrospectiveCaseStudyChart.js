import React, { useMemo } from "react";
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
  if (!active || !payload || payload.length === 0) {
    return null;
  }

  const row = payload[0].payload;
  if (!row) {
    return null;
  }

  const usesDisplayNormalization =
    displayNormalization ===
    "min_max_across_all_spike_sites_in_the_selected_context";
  const displayScoreLabel = usesDisplayNormalization
    ? "Hotspot score"
    : `${selectedModel === "PRIEST" ? "PRIEST" : "Site"} score`;
  const rankingScoreLabel = usesDisplayNormalization
    ? "Raw score"
    : "Ranking score";

  return (
    <div className="rounded-2xl border border-slate-200 bg-white px-4 py-3 shadow-lg">
      <p className="text-sm font-semibold text-slate-900">
        Spike site {row.aa_position}
      </p>
      <p className="mt-1 text-xs text-slate-600">Rank {row.rank}</p>
      <p className="mt-2 text-sm text-slate-700">
        {displayScoreLabel}:{" "}
        <span className="font-semibold text-slate-900">
          {formatScore(row.site_score)}
        </span>
      </p>
      <p className="mt-1 text-sm text-slate-700">
        {rankingScoreLabel}:{" "}
        <span className="font-semibold text-slate-900">
          {formatScore(row.raw_site_score)}
        </span>
      </p>
      <p className="mt-1 text-xs text-slate-500">
        Ref AA: {row.reference_aa || "—"} | Codon positions{" "}
        {Array.isArray(row.codon_genome_positions)
          ? row.codon_genome_positions.join(", ")
          : "—"}
      </p>
      <div className="mt-3 flex flex-wrap gap-2 text-[11px]">
        {row.is_top_k ? (
          <span className="rounded-full bg-sky-50 px-2 py-1 font-medium text-sky-700">
            Hotspot
          </span>
        ) : null}
        {row.is_comparison_site ?? row.is_omicron_site ? (
          <span className="rounded-full bg-amber-50 px-2 py-1 font-medium text-amber-700">
            Known variant site
          </span>
        ) : null}
        {row.is_overlap ? (
          <span className="rounded-full bg-rose-50 px-2 py-1 font-medium text-rose-700">
            Overlap hit
          </span>
        ) : null}
        {row.is_proximity_overlap && !row.is_overlap ? (
          <span className="rounded-full bg-teal-50 px-2 py-1 font-medium text-teal-700">
            Near-hit
          </span>
        ) : null}
      </div>
    </div>
  );
};

const buildScatterPoints = (scoreSeries, predicate) =>
  scoreSeries
    .filter(predicate)
    .map((row) => ({
      ...row,
      x: row.aa_position,
      y: row.site_score,
    }));

const ScatterCircle = ({ cx, cy, fill }) => (
  <circle cx={cx} cy={cy} r={4} fill={fill} />
);

const ScatterDiamond = ({ cx, cy, fill }) => (
  <rect
    x={cx - 4}
    y={cy - 4}
    width={8}
    height={8}
    fill={fill}
    transform={`rotate(45 ${cx} ${cy})`}
  />
);

const buildStarPoints = (cx, cy, outerRadius = 6, innerRadius = 2.8) => {
  const points = [];

  for (let index = 0; index < 10; index += 1) {
    const angle = ((Math.PI / 5) * index) - (Math.PI / 2);
    const radius = index % 2 === 0 ? outerRadius : innerRadius;
    points.push(`${cx + (Math.cos(angle) * radius)},${cy + (Math.sin(angle) * radius)}`);
  }

  return points.join(" ");
};

const ScatterStar = ({ cx, cy, fill }) => (
  <polygon points={buildStarPoints(cx, cy)} fill={fill} />
);

const RetrospectiveCaseStudyChart = ({
  scoreSeries,
  showTopK,
  showOmicron,
  showOverlap,
  showProximity,
  selectedModel,
  displayNormalization,
}) => {
  const usesDisplayNormalization =
    displayNormalization ===
    "min_max_across_all_spike_sites_in_the_selected_context";
  const chartTitle =
    selectedModel === "PRIEST"
      ? "PRIEST Spike Hotspot View"
      : "Explorer Spike Hotspot View";
  const chartDescription = usesDisplayNormalization
    ? "The explorer scores each Spike amino-acid position under the selected post-2022 lineage context, then highlights the top-K hotspots alongside the selected lineage's known Spike mutation sites and their overlap."
    : "PRIEST site prevalence scores are plotted directly on 1-based Spike amino-acid coordinates, with overlays for the selected top-K hotspots, known variant sites, and their overlap.";
  const yAxisLabel = usesDisplayNormalization
    ? "Hotspot score"
    : "Site score";
  const chartData = useMemo(
    () =>
      (scoreSeries || []).map((row) => ({
        ...row,
        topKOnlyScore:
          row.is_top_k && !row.is_overlap ? Number(row.site_score) : null,
        omicronOnlyScore:
          (row.is_comparison_site ?? row.is_omicron_site) && !row.is_overlap
            ? Number(row.site_score)
            : null,
        overlapScore: row.is_overlap ? Number(row.site_score) : null,
      })),
    [scoreSeries]
  );

  const topKOnlyPoints = useMemo(
    () => buildScatterPoints(scoreSeries || [], (row) => row.is_top_k && !row.is_overlap && !row.is_proximity_overlap),
    [scoreSeries]
  );
  const omicronOnlyPoints = useMemo(
    () =>
      buildScatterPoints(
        scoreSeries || [],
        (row) => (row.is_comparison_site ?? row.is_omicron_site) && !row.is_overlap
      ),
    [scoreSeries]
  );
  const overlapPoints = useMemo(
    () => buildScatterPoints(scoreSeries || [], (row) => row.is_overlap),
    [scoreSeries]
  );
  const proximityOnlyPoints = useMemo(
    () => buildScatterPoints(scoreSeries || [], (row) => row.is_proximity_overlap && !row.is_overlap),
    [scoreSeries]
  );

  if (!chartData.length) {
    return null;
  }

  const maxPosition = Math.max(...chartData.map((row) => row.aa_position));

  return (
    <div className="rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <div>
          <h2 className="text-lg font-semibold text-slate-900">
            {chartTitle}
          </h2>
          <p className="mt-1 text-sm text-slate-600">{chartDescription}</p>
        </div>
      </div>

      <div className="mt-5 h-[28rem] w-full" data-testid="retrospective-chart">
        <ResponsiveContainer width="100%" height="100%">
          <ComposedChart
            data={chartData}
            margin={{ top: 12, right: 18, bottom: 12, left: 0 }}
          >
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
                <CaseStudyTooltip
                  selectedModel={selectedModel}
                  displayNormalization={displayNormalization}
                />
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
              <Scatter
                data={topKOnlyPoints}
                fill="#0ea5e9"
                line={false}
                shape={<ScatterCircle />}
                isAnimationActive={false}
              />
            ) : null}
            {showOmicron && omicronOnlyPoints.length ? (
              <Scatter
                data={omicronOnlyPoints}
                fill="#f59e0b"
                line={false}
                shape={<ScatterDiamond />}
                isAnimationActive={false}
              />
            ) : null}
            {showOverlap && overlapPoints.length ? (
              <Scatter
                data={overlapPoints}
                fill="#e11d48"
                line={false}
                shape={<ScatterStar />}
                isAnimationActive={false}
              />
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
              <Scatter
                data={proximityOnlyPoints}
                fill="#14b8a6"
                line={false}
                shape={<ScatterCircle />}
                isAnimationActive={false}
              />
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

export default RetrospectiveCaseStudyChart;
