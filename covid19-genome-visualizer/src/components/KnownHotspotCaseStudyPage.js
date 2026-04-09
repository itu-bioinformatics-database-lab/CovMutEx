import React, { useEffect, useState } from "react";

import logo from "../CovMutexLogo-removebg-preview.png";
import { modelList } from "../data/modelList";
import HotspotCaseStudyChart from "./HotspotCaseStudyChart";

const formatPercent = (value) =>
  value === null || value === undefined
    ? "—"
    : `${(Number(value) * 100).toFixed(1)}%`;

const formatScore = (value) =>
  value === null || value === undefined ? "—" : Number(value).toFixed(3);

const formatMass = (value) =>
  value === null || value === undefined ? "—" : Number(value).toFixed(1);

const formatThreshold = (value) =>
  value === null || value === undefined
    ? "—"
    : `${(Number(value) * 100).toFixed(0)}%`;

const defaultControls = {
  selectedModel: "balanced_data_model",
  elapsedDay: "0",
  topK: "50",
  nodeId: "",
};

const MetricCard = ({ label, value, helper }) => (
  <div className="rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
    <p className="text-xs font-semibold uppercase tracking-[0.16em] text-slate-400">
      {label}
    </p>
    <p className="mt-3 text-3xl font-semibold leading-none text-slate-900">
      {value}
    </p>
    {helper ? <p className="mt-2 text-sm text-slate-500">{helper}</p> : null}
  </div>
);

const MarkerSymbol = ({ shape, className }) => {
  if (shape === "circle") {
    return (
      <span
        aria-hidden="true"
        className={`inline-block h-2.5 w-2.5 shrink-0 rounded-full ${className}`}
      />
    );
  }

  if (shape === "diamond") {
    return (
      <span
        aria-hidden="true"
        className={`inline-block h-2.5 w-2.5 shrink-0 rotate-45 rounded-[2px] ${className}`}
      />
    );
  }

  return (
    <svg
      aria-hidden="true"
      viewBox="0 0 24 24"
      className={`h-3.5 w-3.5 ${className}`}
      fill="currentColor"
    >
      <path d="M12 2.8l2.67 5.42 5.98.87-4.33 4.23 1.02 5.96L12 16.46 6.66 19.28l1.02-5.96-4.33-4.23 5.98-.87L12 2.8z" />
    </svg>
  );
};

const MarkerToggle = ({ label, toneClass, active, onClick, icon }) => (
  <button
    type="button"
    aria-pressed={active}
    onClick={onClick}
    className={`rounded-full border px-3 py-1.5 text-sm font-medium transition ${
      active
        ? `${toneClass} border-current`
        : "border-slate-300 bg-white text-slate-500 hover:border-slate-400"
    }`}
  >
    <span className="flex items-center gap-2">
      {icon}
      <span>{label}</span>
    </span>
  </button>
);

const MetadataCard = ({ label, value }) => (
  <div className="rounded-2xl border border-slate-200 bg-white px-4 py-4 shadow-sm">
    <p className="text-xs font-semibold uppercase tracking-[0.14em] text-slate-400">
      {label}
    </p>
    <p className="mt-2 break-all text-sm text-slate-800">{value || "—"}</p>
  </div>
);

const KnownHotspotCaseStudyPage = () => {
  const API_URL = process.env.REACT_APP_API_URL;
  const [controls, setControls] = useState(defaultControls);
  const [analysis, setAnalysis] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const [showTopK, setShowTopK] = useState(true);
  const [showOmicron, setShowOmicron] = useState(true);
  const [showOverlap, setShowOverlap] = useState(true);
  const [showProximity, setShowProximity] = useState(true);
  const [showAllRows, setShowAllRows] = useState(false);

  const loadAnalysis = async (overrideControls) => {
    const nextControls = overrideControls || controls;
    setLoading(true);
    setError("");

    try {
      const response = await fetch(
        `${API_URL}/api/case-studies/known-hotspot/`,
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            selectedModel: nextControls.selectedModel,
            nodeId: nextControls.nodeId || undefined,
            elapsedDay: Number(nextControls.elapsedDay || 0),
            topK: Number(nextControls.topK || 50),
          }),
        }
      );

      const payload = await response.json();
      if (!response.ok) {
        throw new Error(payload.error || "Unable to load case study.");
      }

      setAnalysis(payload);
      setControls({
        selectedModel:
          payload.metadata?.scoring_context?.selected_model ||
          nextControls.selectedModel,
        elapsedDay: String(
          payload.metadata?.scoring_context?.elapsed_day ??
            nextControls.elapsedDay
        ),
        topK: String(payload.metadata?.applied_top_k ?? nextControls.topK),
        nodeId:
          payload.metadata?.scoring_context?.node_id || nextControls.nodeId,
      });
    } catch (fetchError) {
      setError(fetchError.message);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadAnalysis(defaultControls);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const handleControlChange = (event) => {
    const { name, value } = event.target;
    setControls((previous) => ({
      ...previous,
      [name]: value,
    }));
  };

  const handleSubmit = (event) => {
    event.preventDefault();
    setShowAllRows(false);
    loadAnalysis(controls);
  };

  const metadata = analysis?.metadata;
  const metrics = analysis?.metrics;
  const scoringContext = metadata?.scoring_context || {};
  const comparisonSet = metadata?.comparison_set || {};
  const precomputedVariantSummary =
    scoringContext?.precomputed_variant_summary || {};
  const availableDeltaContextNodes =
    metadata?.available_delta_context_nodes ||
    metadata?.available_context_nodes ||
    {};
  const isPriestSelected = controls.selectedModel === "PRIEST";
  const selectedScoringModel =
    scoringContext?.selected_model || controls.selectedModel;
  const isPriestAnalysis = selectedScoringModel === "PRIEST";
  const selectedVariantLabel =
    scoringContext?.variant_display_name ||
    comparisonSet?.variant_label ||
    scoringContext?.pangolin_lineage ||
    "Selected variant";
  const comparisonSiteCount =
    metrics?.comparison_site_count ?? metrics?.omicron_site_count ?? 0;
  const usesDisplayNormalization =
    scoringContext?.display_normalization ===
    "min_max_across_all_spike_sites_in_the_selected_context";
  const deltaContextOptions =
    availableDeltaContextNodes.options ||
    (scoringContext.node_id
      ? [{ node_id: scoringContext.node_id, label: scoringContext.node_id }]
      : []);
  const overlaySeriesLabel = usesDisplayNormalization
    ? "display-normalized hotspot score series"
    : "site-score series";
  const rankedRowsDescription = isPriestAnalysis
    ? "Spike positions ranked by PRIEST site-level prevalence score for the period implied by the selected lineage. Scores are on a 0 to 1 scale and displayed as-is."
    : "Spike positions ranked by the explorer's hotspot score for the selected lineage. A representative consensus genome is constructed from lineage-specific nucleotide mutation frequencies and scored by CovMutEx across all Spike positions. Hotspot Score is min-max normalized for display; Raw Score is the underlying value used for ranking.";
  const displayScoreLabel = usesDisplayNormalization
    ? "Hotspot Score (normalised)"
    : "Site Score";
  const rankingScoreLabel = usesDisplayNormalization
    ? "Raw Score (un-normalised)"
    : "Ranking Score";
  const rankedRows = analysis?.ranked_rows ?? [];
  const visibleRows = showAllRows ? rankedRows : rankedRows.slice(0, 100);

  return (
    <div className="min-h-screen bg-slate-50 px-4 py-8 sm:px-6 lg:px-8">
      <div className="mx-auto max-w-7xl">
        <div className="mb-8 flex flex-col gap-6 lg:flex-row lg:items-end lg:justify-between">
          <div className="flex items-center gap-4">
            <img
              src={logo}
              alt="CovMutEx logo"
              className="h-16 w-16 rounded-2xl bg-white p-2 shadow-sm"
            />
            <div>
              <p className="text-xs font-semibold uppercase tracking-[0.2em] text-sky-700">
                Explorer Feature
              </p>
              <h1 className="mt-2 text-3xl font-semibold text-slate-900">
                Variant Hotspot Explorer
              </h1>
              <p className="mt-2 max-w-3xl text-sm text-slate-600">
                Select a post-2022 SARS-CoV-2 lineage to explore which Spike
                positions CovMutEx ranks as mutational hotspots. The explorer
                constructs a representative consensus genome from
                lineage-specific nucleotide mutation frequencies and scores
                every Spike position for mutational pressure. Known mutation
                sites for each lineage are overlaid as a biological reference
                layer, letting you visually inspect where the explorer's
                top-ranked positions fall relative to documented Spike biology.
              </p>
            </div>
          </div>
        </div>

        <div className="rounded-3xl border border-sky-200 bg-sky-50 px-5 py-4 text-sm text-sky-900 shadow-sm">
          <span className="font-semibold">How to use: </span>
          The explorer scores every Spike amino-acid position for mutational
          pressure using the selected lineage's genomic context and ranks
          positions by hotspot score. Known Spike mutation sites for that
          lineage, derived from population-level sequencing data, are displayed
          as a biological reference layer. The overlap and proximity metrics
          quantify how the top-ranked positions correspond to documented
          mutation sites, supporting exploratory interpretation of the Spike
          hotspot landscape rather than serving as a measure of predictive accuracy.
        </div>

        <div className="mt-6 rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
          <div className="flex flex-wrap items-center justify-between gap-4">
            <div>
              <h2 className="text-lg font-semibold text-slate-900">
                Explorer settings
              </h2>
              <p className="mt-1 text-sm text-slate-600">
                Choose a lineage and a scoring model. The explorer will score
                every Spike position under that lineage's genomic context and
                display the top-K ranked hotspot positions.
              </p>
            </div>
          </div>

          <form
            className="mt-5 grid gap-4 md:grid-cols-2 xl:grid-cols-5"
            onSubmit={handleSubmit}
          >
            <label className="block">
              <span className="text-sm font-medium text-slate-700">
                Scoring model
              </span>
              <select
                className="mt-2 w-full rounded-2xl border border-slate-300 bg-white px-3 py-2 text-sm text-slate-800"
                name="selectedModel"
                value={controls.selectedModel}
                onChange={handleControlChange}
              >
                {modelList.map((model) => (
                  <option key={model.path} value={model.path}>
                    {model.name}
                  </option>
                ))}
              </select>
            </label>

            <label className="block md:col-span-2 xl:col-span-2">
              <span className="text-sm font-medium text-slate-700">
                Lineage context
              </span>
              <select
                className="mt-2 w-full rounded-2xl border border-slate-300 bg-white px-3 py-2 text-sm text-slate-800"
                name="nodeId"
                value={controls.nodeId}
                onChange={handleControlChange}
              >
                {deltaContextOptions.length ? (
                  deltaContextOptions.map((option) => (
                    <option key={option.node_id} value={option.node_id}>
                      {option.label}
                    </option>
                  ))
                ) : (
                  <option value="">Loading lineage contexts…</option>
                )}
              </select>
              {availableDeltaContextNodes.total_count ? (
                <p className="mt-2 text-xs text-slate-500">
                  {availableDeltaContextNodes.returned_count} precomputed
                  lineage context
                  {availableDeltaContextNodes.returned_count === 1 ? "" : "s"}
                  {" "}available.
                </p>
              ) : null}
            </label>

            <label className="block">
              <span className="text-sm font-medium text-slate-700">
                Elapsed days
              </span>
              <input
                className={`mt-2 w-full rounded-2xl border px-3 py-2 text-sm ${
                  isPriestSelected
                    ? "border-slate-200 bg-slate-100 text-slate-400"
                    : "border-slate-300 bg-white text-slate-800"
                }`}
                name="elapsedDay"
                type="number"
                min="0"
                value={controls.elapsedDay}
                onChange={handleControlChange}
                disabled={isPriestSelected}
              />
              <p className="mt-2 text-xs text-slate-500">
                {isPriestSelected
                  ? "PRIEST uses site-level prevalence scores derived from sequencing data; elapsed days are not applicable."
                  : "Elapsed days since emergence are passed to the CovMutEx model as a temporal feature."}
              </p>
            </label>

            <label className="block">
              <span className="text-sm font-medium text-slate-700">Top-K</span>
              <input
                className="mt-2 w-full rounded-2xl border border-slate-300 bg-white px-3 py-2 text-sm text-slate-800"
                name="topK"
                type="number"
                min="1"
                value={controls.topK}
                onChange={handleControlChange}
              />
            </label>

            <div className="flex items-end">
              <button
                type="submit"
                className="w-full rounded-2xl bg-slate-900 px-4 py-2.5 text-sm font-semibold text-white transition hover:bg-slate-800"
              >
                Run analysis
              </button>
            </div>
          </form>
        </div>

        {loading ? (
          <div className="mt-6 rounded-3xl border border-slate-200 bg-white px-5 py-6 text-sm text-slate-600 shadow-sm">
            Loading case study…
          </div>
        ) : null}

        {error ? (
          <div className="mt-6 rounded-3xl border border-rose-200 bg-rose-50 px-5 py-4 text-sm text-rose-700 shadow-sm">
            {error}
          </div>
        ) : null}

        {!loading && analysis ? (
          <>
            <div className="mt-6 grid gap-4 md:grid-cols-3 xl:grid-cols-3">
              <MetadataCard
                label="Pango Lineage"
                value={scoringContext.pangolin_lineage || selectedVariantLabel}
              />
              <MetadataCard
                label="Common Name"
                value={scoringContext.variant_nickname || "—"}
              />
              <MetadataCard
                label="First Noted"
                value={scoringContext.emergence_label}
              />
            </div>

            {isPriestAnalysis ? (
              <div className="mt-4 rounded-2xl border border-sky-200 bg-sky-50 px-4 py-3 text-sm text-sky-800">
                <strong>Note on PRIEST:</strong> PRIEST was trained on
                sequencing data up to early 2022. The lineages shown here
                (2024 onward) are temporally out-of-sample. PRIEST provides
                quarterly population-level Spike-site prevalence scores rather
                than lineage-specific predictions, so its output is identical
                for any two lineages mapped to the same quarter.
              </div>
            ) : null}

            <div className="mt-6 grid gap-4 md:grid-cols-3">
              <MetricCard
                label="Exact Overlap"
                value={metrics?.overlap_count ?? 0}
                helper={`Top-${metrics?.top_k_count ?? "K"} hotspot positions that exactly match a known ${selectedVariantLabel} Spike site (${comparisonSiteCount} reference sites).`}
              />
              <MetricCard
                label="Precision @ K (Exact)"
                value={formatPercent(metrics?.precision_at_k)}
                helper="Share of top-K hotspot positions that fall on a reference Spike site."
              />
              <MetricCard
                label="Recall @ K (Exact)"
                value={
                  formatPercent(
                    metrics?.recall_against_comparison_sites ??
                      metrics?.recall_against_omicron_sites
                  )
                }
                helper="Share of reference Spike sites covered by the top-K hotspot positions."
              />
            </div>

            <div className="mt-3 grid gap-4 md:grid-cols-3">
              <MetricCard
                label={`Proximity Overlap (±${
                  metrics?.proximity_window ?? 3
                }aa)`}
                value={metrics?.proximity_overlap_count ?? 0}
                helper={`Top-K positions within ±${metrics?.proximity_window ?? 3}aa of a reference ${selectedVariantLabel} Spike site.`}
              />
              <MetricCard
                label="Proximity Precision @ K"
                value={formatPercent(metrics?.proximity_precision_at_k)}
                helper={`Share of top-K positions within ±${metrics?.proximity_window ?? 3}aa of a reference site.`}
              />
              <MetricCard
                label="Proximity Recall @ K"
                value={formatPercent(metrics?.proximity_recall)}
                helper={`Share of reference sites with an explorer hotspot within ±${metrics?.proximity_window ?? 3}aa.`}
              />
            </div>

            {/* Secondary support metrics — context strip */}
            {isPriestAnalysis ? (
              /* PRIEST: no synthetic genome, so skip injected-mutation stats.
                 Show comparison set size + overlap proportion instead.      */
              <div className="mt-4 rounded-2xl border border-slate-100 bg-slate-50 px-5 py-4">
                <p className="mb-3 text-xs font-semibold uppercase tracking-[0.14em] text-slate-400">
                  Context &amp; support data
                </p>
                <div className="flex flex-wrap gap-x-8 gap-y-3 text-sm text-slate-700">
                  <span>
                    <span className="font-medium text-slate-900">
                      {comparisonSiteCount}
                    </span>{" "}
                    known {selectedVariantLabel} Spike sites
                  </span>
                  <span>
                    Overlap mutations:{" "}
                    <span className="font-medium text-slate-900">
                      {metrics?.found_mutation_count ?? 0}
                    </span>
                  </span>
                  <span>
                    Mean overlap proportion:{" "}
                    <span className="font-medium text-slate-900">
                      {formatPercent(metrics?.found_mean_proportion)}
                    </span>
                  </span>
                  <span>
                    Median overlap proportion:{" "}
                    <span className="font-medium text-slate-900">
                      {formatPercent(metrics?.found_median_proportion)}
                    </span>
                  </span>
                </div>
              </div>
            ) : (
              /* CovMutEx ML models: show full synthetic genome injection stats */
              <div className="mt-4 rounded-2xl border border-slate-100 bg-slate-50 px-5 py-4">
                <p className="mb-3 text-xs font-semibold uppercase tracking-[0.14em] text-slate-400">
                  Context &amp; support data
                </p>
                <div className="flex flex-wrap gap-x-8 gap-y-3 text-sm text-slate-700">
                  <span>
                    <span className="font-medium text-slate-900">
                      {precomputedVariantSummary?.mutation_count ?? 0}
                    </span>{" "}
                    consensus mutations applied
                  </span>
                  <span>
                    <span className="font-medium text-slate-900">
                      {precomputedVariantSummary?.spike_site_count ?? 0}
                    </span>{" "}
                    Spike sites compared
                  </span>
                  <span>
                    Elapsed days:{" "}
                    <span className="font-medium text-slate-900">
                      {scoringContext?.elapsed_day ?? "—"}
                    </span>
                  </span>
                  <span>
                    Mean site proportion:{" "}
                    <span className="font-medium text-slate-900">
                      {formatPercent(
                        precomputedVariantSummary?.spike_support_summary
                          ?.mean_proportion
                      )}
                    </span>
                  </span>
                  <span>
                    Median site proportion:{" "}
                    <span className="font-medium text-slate-900">
                      {formatPercent(
                        precomputedVariantSummary?.spike_support_summary
                          ?.median_proportion
                      )}
                    </span>
                  </span>
                  <span>
                    Overlap mutations:{" "}
                    <span className="font-medium text-slate-900">
                      {metrics?.found_mutation_count ?? 0}
                    </span>
                  </span>
                  <span>
                    Mean overlap proportion:{" "}
                    <span className="font-medium text-slate-900">
                      {formatPercent(metrics?.found_mean_proportion)}
                    </span>
                  </span>
                  <span>
                    Median overlap proportion:{" "}
                    <span className="font-medium text-slate-900">
                      {formatPercent(metrics?.found_median_proportion)}
                    </span>
                  </span>
                </div>
              </div>
            )}

            <div className="mt-6 rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
              <div className="flex flex-wrap items-center justify-between gap-3">
                <div>
                  <h2 className="text-lg font-semibold text-slate-900">
                    Chart overlays
                  </h2>
                  <p className="mt-1 text-sm text-slate-600">
                    Toggle each marker layer to compare the explorer’s hotspot
                    scores (across all Spike positions) against known{" "}
                    {selectedVariantLabel} Spike mutation sites and their
                    overlap, plotted on the{" "}
                    {overlaySeriesLabel}.
                  </p>
                </div>
                <div className="flex flex-wrap gap-2">
                  <MarkerToggle
                    label="Explorer hotspots"
                    active={showTopK}
                    onClick={() => setShowTopK((previous) => !previous)}
                    toneClass="bg-sky-50 text-sky-700"
                    icon={
                      <MarkerSymbol shape="circle" className="bg-sky-500" />
                    }
                  />
                  <MarkerToggle
                    label="Known variant sites"
                    active={showOmicron}
                    onClick={() => setShowOmicron((previous) => !previous)}
                    toneClass="bg-amber-50 text-amber-700"
                    icon={
                      <MarkerSymbol shape="diamond" className="bg-amber-400" />
                    }
                  />
                  <MarkerToggle
                    label="Overlap hits"
                    active={showOverlap}
                    onClick={() => setShowOverlap((previous) => !previous)}
                    toneClass="bg-rose-50 text-rose-700"
                    icon={
                      <MarkerSymbol shape="star" className="text-rose-600" />
                    }
                  />
                  <MarkerToggle
                    label={`Near-hits (±${metrics?.proximity_window ?? 3}aa)`}
                    active={showProximity}
                    onClick={() => setShowProximity((previous) => !previous)}
                    toneClass="bg-teal-50 text-teal-700"
                    icon={
                      <MarkerSymbol shape="circle" className="bg-teal-400" />
                    }
                  />
                </div>
              </div>
            </div>

            <div className="mt-6">
              <HotspotCaseStudyChart
                scoreSeries={analysis.score_series}
                showTopK={showTopK}
                showOmicron={showOmicron}
                showOverlap={showOverlap}
                showProximity={showProximity}
                selectedModel={selectedScoringModel}
                displayNormalization={scoringContext?.display_normalization}
                variantLabel={scoringContext?.pangolin_lineage || "variant"}
                topK={controls.topK}
              />
            </div>

            <div className="mt-6 rounded-3xl border border-slate-200 bg-white px-5 py-5 shadow-sm">
              <div className="flex flex-wrap items-center justify-between gap-4">
                <div>
                  <h2 className="text-lg font-semibold text-slate-900">
                    Spike hotspot rankings
                  </h2>
                  <p className="mt-1 text-sm text-slate-600">
                    {rankedRowsDescription}
                  </p>
                </div>
                <button
                  type="button"
                  onClick={() => setShowAllRows((previous) => !previous)}
                  className="rounded-full border border-slate-300 px-3 py-1.5 text-sm font-medium text-slate-700 transition hover:bg-slate-50"
                >
                  {showAllRows
                    ? "Show first 100 rows"
                    : `Show all ${rankedRows.length} rows`}
                </button>
              </div>

              <div className="mt-5 overflow-x-auto">
                <table className="min-w-full border-separate border-spacing-0 text-left text-sm">
                  <thead>
                    <tr className="text-slate-500">
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        Rank
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        Spike Site
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        Ref AA
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        {displayScoreLabel}
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        {rankingScoreLabel}
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        Top-K
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        Variant Site
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        Overlap
                      </th>
                      <th className="border-b border-slate-200 px-3 py-3 font-semibold">
                        Codon Positions
                      </th>
                    </tr>
                  </thead>
                  <tbody>
                    {visibleRows.map((row) => (
                      <tr
                        key={`ranked-site-${row.aa_position}`}
                        className="text-slate-700"
                      >
                        <td className="border-b border-slate-100 px-3 py-3">
                          {row.rank}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3">
                          {row.aa_position}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3">
                          {row.reference_aa || "—"}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3 font-medium text-slate-900">
                          {formatScore(row.site_score)}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3 text-slate-700">
                          {formatScore(row.raw_site_score)}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3">
                          {row.is_top_k ? "Yes" : "No"}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3">
                          {(row.is_comparison_site ?? row.is_omicron_site)
                            ? "Yes"
                            : "No"}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3">
                          {row.is_overlap ? "Yes" : "No"}
                        </td>
                        <td className="border-b border-slate-100 px-3 py-3">
                          {Array.isArray(row.codon_genome_positions)
                            ? row.codon_genome_positions.join(", ")
                            : "—"}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </>
        ) : null}
      </div>
    </div>
  );
};

export default KnownHotspotCaseStudyPage;
