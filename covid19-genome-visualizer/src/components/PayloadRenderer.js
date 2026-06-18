import React from "react";
import { useSelector } from "react-redux";

import GenomeChart from "./Recharts";
import DoughnutChart from "./DoughnutChart";

/**
 * Dispatch the visualization based on the self-describing v2.0
 * PredictionPayload returned by the backend. The frontend stays organism-
 * and task-agnostic: it just renders what `task.kind` tells it to.
 *
 *   categorical_per_position -> existing ATGC GenomeChart (+ doughnut)
 *   binary_per_position      -> PerPositionTrack
 *   scalar_per_position      -> PerPositionTrack
 *
 * When `predictionPayload` is missing (legacy backend response) we fall back
 * to the legacy categorical view if the categorical arrays are populated, so
 * old responses still render.
 */
export default function PayloadRenderer({ scaleType }) {
  const {
    predictionPayload,
    genomeDataRaw,
    genome,
    proteinMutationProbs,
    selectedProteinRegion,
  } = useSelector((state) => state.genome);

  const kind = predictionPayload?.task?.kind ?? null;
  const hasLegacyCategorical = Array.isArray(genomeDataRaw) && genomeDataRaw.length > 0;

  if (!predictionPayload && !hasLegacyCategorical) {
    return (
      <div className="p-8 text-center text-gray-500">No prediction data yet.</div>
    );
  }

  if (kind === "binary_per_position" || kind === "scalar_per_position") {
    // Reuse GenomeChart's full chrome (SidePanel + canvas + pan/zoom) by
    // wrapping the scalar values as a single-channel "genomeData" shape.
    // GenomeChart's mode prop tells it to skip the 4-way ATGC dataset split
    // and draw a single bar per position instead.
    const values = predictionPayload?.predictions?.values ?? [];
    const region = predictionPayload?.domain?.region;
    const totalLength =
      predictionPayload?.domain?.total_length || (genome ? genome.length : values.length);
    const proteinRegionsAnn = predictionPayload?.annotations?.protein_regions ?? {};

    // Region-locked predictions only carry values for the picked window.
    // Pad with zeros so chart x-axis positions match the FULL genome and
    // the nucleotide labels (`1500-G`) stay aligned with the right base.
    let alignedValues = values;
    if (region && Number.isFinite(region.start) && Number.isFinite(region.end)) {
      const leadingZeros = Math.max(0, region.start - 1);
      const inRegionLength = region.end - region.start;
      const trailingZeros = Math.max(0, totalLength - (leadingZeros + inRegionLength));
      alignedValues = new Array(totalLength).fill(0);
      for (let i = 0; i < values.length && i < inRegionLength; i += 1) {
        alignedValues[leadingZeros + i] = values[i];
      }
      // trailingZeros stays implicit — the fill(0) covers it.
      void trailingZeros;
    }

    // Per-region aggregate: sum values within each protein region. For binary
    // this counts mutation positions; for scalar this sums the per-position
    // scores. DoughnutChart's "Normalize by region length" toggle divides by
    // region length to make regions of different sizes comparable. We skip
    // the doughnut entirely if the user picked a single region — only one
    // slice in the pie isn't informative.
    const doughnutData = !selectedProteinRegion
      ? aggregateValuesByRegion(alignedValues, proteinRegionsAnn)
      : {};
    const hasDoughnutData = Object.keys(doughnutData).length > 0;

    const chartMode = kind === "binary_per_position" ? "binary" : "scalar";
    return (
      <div className="p-4">
        <GenomeChart
          genomeData={[alignedValues]}
          genomeSequence={genome || "N".repeat(alignedValues.length)}
          scaleType={scaleType}
          mode={chartMode}
        />
        {hasDoughnutData && (
          <div className="mt-8 pb-12 flex justify-center">
            <DoughnutChart data={doughnutData} />
          </div>
        )}
      </div>
    );
  }

  if (kind && kind !== "categorical_per_position") {
    return (
      <div className="p-8 text-center text-gray-600">
        Unsupported task kind: <span className="font-mono">{kind}</span>
      </div>
    );
  }

  // categorical_per_position OR legacy response with no payload — render the
  // existing ATGC GenomeChart and the protein-doughnut below it.
  return (
    <LegacyCategoricalView
      genomeDataRaw={genomeDataRaw}
      genomeSequence={genome}
      proteinMutationProbs={proteinMutationProbs}
      selectedProteinRegion={selectedProteinRegion}
      scaleType={scaleType}
    />
  );
}

/**
 * Sum the per-position values that fall inside each declared protein
 * region. `annotations.protein_regions` uses 1-based inclusive coordinates
 * ({name: [start, end]}); we convert to 0-based half-open when slicing
 * `alignedValues` (which is itself indexed against the FULL genome).
 *
 * Returns {name: number} keyed by region name. Empty when no regions are
 * annotated. DoughnutChart converts each value to a percentage of the
 * total, so the absolute scale of the values doesn't matter.
 */
function aggregateValuesByRegion(alignedValues, proteinRegionsAnn) {
  if (!Array.isArray(alignedValues) || alignedValues.length === 0) return {};
  if (!proteinRegionsAnn || Object.keys(proteinRegionsAnn).length === 0) return {};
  const out = {};
  for (const [name, range] of Object.entries(proteinRegionsAnn)) {
    if (!Array.isArray(range) || range.length < 2) continue;
    const start = Math.max(0, range[0] - 1);                  // 1-based → 0-based
    const end = Math.min(alignedValues.length, range[1]);     // 1-based inclusive
    if (end <= start) continue;
    let sum = 0;
    for (let i = start; i < end; i += 1) {
      sum += Number(alignedValues[i]) || 0;
    }
    out[name] = sum;
  }
  return out;
}

function LegacyCategoricalView({
  genomeDataRaw,
  genomeSequence,
  proteinMutationProbs,
  selectedProteinRegion,
  scaleType,
}) {
  const hasGenomeData = Array.isArray(genomeDataRaw) && genomeDataRaw.length > 0;
  const hasDoughnutData =
    !selectedProteinRegion &&
    proteinMutationProbs &&
    Object.keys(proteinMutationProbs).length > 0;

  return (
    <div className="p-4">
      {hasGenomeData && (
        <GenomeChart
          genomeData={genomeDataRaw}
          genomeSequence={genomeSequence}
          scaleType={scaleType}
        />
      )}
      {hasDoughnutData && (
        <div className="mt-8 pb-12 flex justify-center">
          <DoughnutChart data={proteinMutationProbs} />
        </div>
      )}
    </div>
  );
}

