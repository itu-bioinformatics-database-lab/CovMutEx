import React from "react";
import { Card, Typography } from "@material-tailwind/react";
import {
  CartesianGrid,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

const formatScore = (score) =>
  score === null || score === undefined ? "—" : Number(score).toFixed(3);

const formatPercent = (value) =>
  value === null || value === undefined ? "—" : `${(value * 100).toFixed(1)}%`;

const formatSource = (source) => {
  if (source === "period") {
    return "Matched period";
  }
  if (source === "global") {
    return "Global fallback";
  }
  return "—";
};

const getScoreBand = (score, threshold) => {
  if (score === null || score === undefined) {
    return {
      label: "Unavailable",
      className: "bg-slate-100 text-slate-600",
    };
  }
  if (score >= threshold) {
    return {
      label: "High support",
      className: "bg-emerald-50 text-emerald-700",
    };
  }
  if (score >= threshold / 2) {
    return {
      label: "Moderate",
      className: "bg-amber-50 text-amber-700",
    };
  }
  return {
    label: "Low",
    className: "bg-slate-100 text-slate-700",
  };
};

const buildSummaryLine = ({
  annotation,
  totalMutations,
  nonSpikeMutations,
  spikeNtMutations,
  nonsynonymousCount,
  synonymousCount,
  annotatedCount,
}) => {
  if (!totalMutations) {
    return "No reconstructed mutations were detected for this node.";
  }

  if (!spikeNtMutations) {
    return `${totalMutations} reconstructed mutation${
      totalMutations === 1 ? "" : "s"
    } were identified for this node, and all ${nonSpikeMutations} occurred outside Spike, so no PRIEST site annotation is applicable for ${
      annotation.priest_period || "the mapped period"
    }.`;
  }

  return `${totalMutations} reconstructed mutation${
    totalMutations === 1 ? "" : "s"
  }, including ${spikeNtMutations} Spike nucleotide event${
    spikeNtMutations === 1 ? "" : "s"
  }, collapsed to ${nonsynonymousCount} nonsynonymous and ${synonymousCount} synonymous Spike site${
    synonymousCount === 1 ? "" : "s"
  }. In ${
    annotation.priest_period || "the mapped PRIEST interval"
  }, ${annotatedCount} of ${nonsynonymousCount} nonsynonymous Spike site${
    nonsynonymousCount === 1 ? "" : "s"
  } were matched to PRIEST site-level scores.`;
};

const MetadataItem = ({ label, value, tone = "default" }) => {
  const toneClass =
    tone === "strong"
      ? "text-slate-900"
      : tone === "muted"
      ? "text-slate-500"
      : "text-slate-700";

  return (
    <div className="rounded-2xl border border-slate-200 bg-white px-4 py-3">
      <Typography
        variant="small"
        className="font-semibold uppercase tracking-[0.14em] text-slate-400"
      >
        {label}
      </Typography>
      <Typography className={`mt-2 break-all text-sm ${toneClass}`}>
        {value || "—"}
      </Typography>
    </div>
  );
};

const MetricCard = ({ label, value }) => (
  <div className="rounded-2xl border border-slate-200 bg-white px-4 py-4 shadow-sm">
    <Typography
      variant="small"
      className="font-semibold uppercase tracking-[0.14em] text-slate-400"
    >
      {label}
    </Typography>
    <Typography className="mt-3 text-3xl font-semibold leading-none text-slate-900">
      {value}
    </Typography>
  </div>
);

const SectionHeading = ({ title, subtitle }) => (
  <div>
    <Typography variant="h6" className="text-slate-900">
      {title}
    </Typography>
    {subtitle ? (
      <Typography variant="small" className="mt-1 max-w-3xl text-slate-600">
        {subtitle}
      </Typography>
    ) : null}
  </div>
);

const RegionDistribution = ({ regions }) => {
  if (!regions || regions.length === 0) {
    return null;
  }

  return (
    <div className="rounded-2xl border border-slate-200 bg-white px-5 py-4 shadow-sm">
      <SectionHeading
        title="Observed Mutation Distribution"
        subtitle="Reconstructed mutations grouped by genomic region."
      />
      <div className="mt-3 flex flex-wrap gap-2">
        {regions.map((region) => (
          <div
            key={`${region.gene}-${region.count}`}
            className="rounded-full bg-slate-100 px-3 py-1 text-sm text-slate-700"
          >
            <span className="font-medium">{region.gene}</span>: {region.count}
          </div>
        ))}
      </div>
    </div>
  );
};

const PriestScoreChart = ({
  nonsynonymousMutations,
  synonymousMutations,
  threshold,
}) => {
  const nonsynonymousRows = (nonsynonymousMutations || [])
    .filter(
      (mutation) =>
        mutation.priest_score !== null && mutation.priest_score !== undefined
    )
    .map((mutation) => ({
      aa_mutation: mutation.aa_mutation,
      aa_position: Number(mutation.aa_position),
      priest_score: Number(mutation.priest_score),
      site_class: "Nonsynonymous",
      priest_score_source: mutation.priest_score_source,
    }));

  const synonymousRows = (synonymousMutations || [])
    .filter(
      (mutation) =>
        mutation.priest_score !== null && mutation.priest_score !== undefined
    )
    .map((mutation) => ({
      aa_mutation: mutation.aa_mutation,
      aa_position: Number(mutation.aa_position),
      priest_score: Number(mutation.priest_score),
      site_class: "Synonymous",
      priest_score_source: mutation.priest_score_source,
    }));

  const allRows = [...nonsynonymousRows, ...synonymousRows];
  if (allRows.length === 0) {
    return null;
  }

  const positions = allRows.map((row) => row.aa_position);
  const minPosition = Math.min(...positions);
  const maxPosition = Math.max(...positions);
  const positionPadding =
    minPosition === maxPosition
      ? 6
      : Math.max(6, Math.round((maxPosition - minPosition) * 0.08));

  return (
    <div className="rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
      <SectionHeading
        title="PRIEST Site View"
        subtitle="Observed Spike sites positioned along Spike and scored within the mapped temporal window."
      />

      <div className="mt-3 flex flex-wrap gap-2 text-xs text-slate-600">
        {nonsynonymousRows.length > 0 ? (
          <div className="inline-flex items-center gap-2 rounded-full bg-slate-100 px-3 py-1">
            <span className="h-2.5 w-2.5 rounded-full bg-blue-600" />
            Nonsynonymous sites
          </div>
        ) : null}
        {synonymousRows.length > 0 ? (
          <div className="inline-flex items-center gap-2 rounded-full bg-slate-100 px-3 py-1">
            <span className="h-2.5 w-2.5 rounded-full bg-slate-500" />
            Synonymous sites
          </div>
        ) : null}
        <div className="inline-flex items-center gap-2 rounded-full bg-amber-50 px-3 py-1 text-amber-800">
          <span className="h-px w-4 border-t border-dashed border-amber-500" />
          Threshold {formatScore(threshold)}
        </div>
      </div>

      <div className="mt-4 h-72 w-full">
        <ResponsiveContainer width="100%" height="100%">
          <ScatterChart margin={{ top: 12, right: 20, bottom: 12, left: 0 }}>
            <CartesianGrid stroke="#e2e8f0" strokeDasharray="3 3" />
            <XAxis
              type="number"
              dataKey="aa_position"
              name="Spike site"
              domain={[minPosition - positionPadding, maxPosition + positionPadding]}
              label={{
                value: "Spike amino-acid position",
                position: "insideBottom",
                offset: -6,
                style: { fill: "#475569", fontSize: 12 },
              }}
              tick={{ fill: "#64748b", fontSize: 12 }}
              tickLine={false}
              axisLine={{ stroke: "#cbd5e1" }}
            />
            <YAxis
              type="number"
              dataKey="priest_score"
              name="PRIEST score"
              domain={[0, 1]}
              label={{
                value: "PRIEST site score",
                angle: -90,
                position: "insideLeft",
                style: { fill: "#475569", fontSize: 12, textAnchor: "middle" },
              }}
              tick={{ fill: "#64748b", fontSize: 12 }}
              tickFormatter={(value) => Number(value).toFixed(1)}
              tickLine={false}
              axisLine={{ stroke: "#cbd5e1" }}
            />
            <ReferenceLine
              y={threshold}
              stroke="#f59e0b"
              strokeDasharray="4 4"
              ifOverflow="extendDomain"
            />
            <Tooltip
              cursor={{ stroke: "#cbd5e1", strokeDasharray: "4 4" }}
              content={({ active, payload }) => {
                if (!active || !payload || payload.length === 0) {
                  return null;
                }
                const row = payload[0].payload;
                return (
                  <div className="rounded-xl border border-slate-200 bg-white px-3 py-2 text-sm shadow-lg">
                    <div className="font-semibold text-slate-900">
                      {row.aa_mutation}
                    </div>
                    <div className="mt-1 text-slate-600">
                      Spike site {row.aa_position}
                    </div>
                    <div className="text-slate-600">
                      PRIEST score {formatScore(row.priest_score)}
                    </div>
                    <div className="text-slate-500">{row.site_class}</div>
                    <div className="text-slate-500">
                      {formatSource(row.priest_score_source)}
                    </div>
                  </div>
                );
              }}
            />
            <Scatter data={nonsynonymousRows} fill="#2563eb" />
            <Scatter data={synonymousRows} fill="#64748b" />
          </ScatterChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
};

const MutationTable = ({
  title,
  subtitle,
  mutations,
  threshold,
  emptyMessage,
}) => (
  <div className="rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
    <SectionHeading title={title} subtitle={subtitle} />
    {!mutations || mutations.length === 0 ? (
      <Typography variant="small" className="mt-4 text-slate-600">
        {emptyMessage}
      </Typography>
    ) : (
      <div className="mt-5 overflow-x-auto rounded-2xl border border-slate-200">
        <table className="min-w-full divide-y divide-slate-200 text-left text-sm">
          <thead className="bg-slate-50 text-slate-600">
            <tr>
              <th className="px-4 py-3 font-semibold">Mutation</th>
              <th className="px-4 py-3 font-semibold">Spike site</th>
              <th className="px-4 py-3 font-semibold">Codon transition</th>
              <th className="px-4 py-3 font-semibold">Genome loci</th>
              <th className="px-4 py-3 font-semibold">PRIEST score</th>
              <th className="px-4 py-3 font-semibold">Relative support</th>
              <th className="px-4 py-3 font-semibold">Source</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-200 bg-white text-slate-700">
            {mutations.map((mutation) => {
              const scoreBand = getScoreBand(mutation.priest_score, threshold);
              return (
                <tr
                  key={`${title}-${mutation.aa_mutation}-${mutation.aa_position}`}
                >
                  <td className="px-4 py-3 font-semibold text-slate-900">
                    {mutation.aa_mutation}
                  </td>
                  <td className="px-4 py-3">{mutation.aa_position}</td>
                  <td className="px-4 py-3">
                    {mutation.ref_codon && mutation.alt_codon
                      ? `${mutation.ref_codon} → ${mutation.alt_codon}`
                      : "—"}
                  </td>
                  <td className="px-4 py-3">
                    {mutation.genome_positions?.length
                      ? mutation.genome_positions.join(", ")
                      : mutation.genome_position || "—"}
                  </td>
                  <td className="px-4 py-3">
                    {formatScore(mutation.priest_score)}
                  </td>
                  <td className="px-4 py-3">
                    <span
                      className={`rounded-full px-3 py-1 text-xs font-semibold ${scoreBand.className}`}
                    >
                      {scoreBand.label}
                    </span>
                  </td>
                  <td className="px-4 py-3">
                    {formatSource(mutation.priest_score_source)}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    )}
  </div>
);

const SpikeMutationPanel = ({ annotation }) => {
  if (!annotation) {
    return null;
  }

  const spikeMutations = annotation.spike_mutations || [];
  const synonymousMutations = annotation.synonymous_spike_mutations || [];
  const summary = annotation.priest_summary || {};
  const totalMutations = summary.num_total_mutations ?? 0;
  const nonSpikeMutations = summary.num_non_spike_mutations ?? 0;
  const spikeNtMutations = summary.num_spike_nucleotide_mutations ?? 0;
  const nonsynonymousCount = summary.num_spike_mutations ?? 0;
  const synonymousCount =
    summary.num_synonymous_spike_mutations ?? synonymousMutations.length;
  const annotatedCount = summary.num_priest_annotated ?? 0;
  const threshold = summary.priest_support_threshold ?? 0.6;
  const regionCounts = summary.mutation_region_counts || [];
  const allScoredSpikeRows = [...spikeMutations, ...synonymousMutations].filter(
    (mutation) => mutation.priest_score !== null && mutation.priest_score !== undefined
  );
  const displayMeanPriestScore =
    summary.mean_priest_score !== null && summary.mean_priest_score !== undefined
      ? summary.mean_priest_score
      : allScoredSpikeRows.length > 0
      ? allScoredSpikeRows.reduce(
          (total, mutation) => total + Number(mutation.priest_score),
          0
        ) / allScoredSpikeRows.length
      : null;
  const meanScoreLabel =
    summary.mean_priest_score !== null && summary.mean_priest_score !== undefined
      ? "Mean PRIEST score"
      : allScoredSpikeRows.length > 0 && synonymousCount > 0 && nonsynonymousCount === 0
      ? "Mean PRIEST score (synonymous-only)"
      : "Mean PRIEST score";
  const summaryLine = buildSummaryLine({
    annotation,
    totalMutations,
    nonSpikeMutations,
    spikeNtMutations,
    nonsynonymousCount,
    synonymousCount,
    annotatedCount,
  });

  return (
    <Card className="mx-4 mb-6 mt-16 overflow-hidden border border-slate-200 bg-white shadow-[0_20px_60px_-40px_rgba(15,23,42,0.35)]">
      <div className="border-b border-slate-200 bg-gradient-to-r from-slate-50 via-white to-slate-50 px-6 py-6">
        <div className="flex flex-col gap-4 xl:flex-row xl:items-start xl:justify-between">
          <div className="max-w-4xl">
            <Typography variant="h4" className="text-slate-900">
              Spike Mutation Annotation
            </Typography>
            <Typography className="mt-2 max-w-3xl text-sm leading-6 text-slate-600">
              Reconstructed Spike amino-acid substitutions for the selected
              SARS-CoV-2 variant, cross-referenced against PRIEST site-level
              scores.
            </Typography>
          </div>
          <div className="grid gap-3 sm:grid-cols-2 xl:w-[320px] xl:grid-cols-1">
            <MetadataItem
              label="Collection date"
              value={annotation.node_date}
              tone="strong"
            />
            <MetadataItem
              label="PRIEST temporal window"
              value={annotation.priest_period}
              tone="strong"
            />
          </div>
        </div>
      </div>

      <div className="space-y-6 p-6">
        <div className="grid gap-4 xl:grid-cols-[1.6fr_1fr_1fr]">
          <MetadataItem
            label="Selected node"
            value={annotation.selected_node}
            tone="strong"
          />
          <MetadataItem
            label="Accession"
            value={annotation.selected_node_accession}
          />
          <MetadataItem
            label="Lookup mode"
            value={
              annotation.priest_score_method === "csv_lookup"
                ? "Local PRIEST lookup table"
                : "Raw PRIEST period prevalence fallback"
            }
          />
        </div>

        <div className="rounded-2xl border border-slate-200 bg-slate-50 px-5 py-4">
          <Typography className="text-sm leading-6 text-slate-700">
            {summaryLine}
          </Typography>
        </div>

        <div className="grid gap-4 md:grid-cols-2 xl:grid-cols-4">
          <MetricCard label="Total mutations" value={totalMutations} />
          <MetricCard label="Spike AA mutations" value={nonsynonymousCount} />
          <MetricCard label="PRIEST sites" value={annotatedCount} />
          <MetricCard
            label={meanScoreLabel}
            value={formatScore(displayMeanPriestScore)}
          />
        </div>

        <RegionDistribution regions={regionCounts} />

        <PriestScoreChart
          nonsynonymousMutations={spikeMutations}
          synonymousMutations={synonymousMutations}
          threshold={threshold}
        />

        <MutationTable
          title="Nonsynonymous Spike Sites"
          subtitle="Observed Spike amino-acid substitutions with PRIEST site support."
          mutations={spikeMutations}
          threshold={threshold}
          emptyMessage={
            nonSpikeMutations > 0
              ? `No nonsynonymous Spike amino-acid substitutions were reconstructed for this node. CovMutEx identified ${totalMutations} total mutation${
                  totalMutations === 1 ? "" : "s"
                }, with ${nonSpikeMutations} ${
                  nonSpikeMutations === 1 ? "event" : "events"
                } occurring outside Spike.`
              : "No nonsynonymous Spike amino-acid substitutions were reconstructed for this node."
          }
        />

        <MutationTable
          title="Synonymous Spike Sites"
          subtitle="Observed Spike nucleotide events that preserve amino-acid identity."
          mutations={synonymousMutations}
          threshold={threshold}
          emptyMessage="No synonymous Spike-site events were reconstructed for this node."
        />
      </div>
    </Card>
  );
};

export default SpikeMutationPanel;
