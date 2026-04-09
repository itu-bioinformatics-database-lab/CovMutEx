import React from "react";

export const knownHotspotCaseStudy = {
  selectedVariant: "Delta-like",
  spikeMutations: [
    { aa_position: 452, mutation: "L452R", priest_score: 0.696 },
    { aa_position: 478, mutation: "T478K", priest_score: null },
    { aa_position: 614, mutation: "D614G", priest_score: 0.95 },
  ],
  priestTopSites: [
    { aa_position: 413, score: 0.95 },
    { aa_position: 414, score: 0.778 },
    { aa_position: 446, score: 0.59 },
    { aa_position: 452, score: 0.696 },
  ],
  omicronMutations: [417, 440, 446, 478, 501],
};

const SPIKE_START = 21563;

const uniquePositions = (positions) => Array.from(new Set(positions));

export const spikeAaToGenomeWindow = (aaPosition) => {
  const start = SPIKE_START + (aaPosition - 1) * 3;
  return {
    start,
    end: start + 2,
  };
};

export const buildKnownHotspotCaseStudy = () => {
  const omicronSet = new Set(knownHotspotCaseStudy.omicronMutations);

  const sharedPositions = uniquePositions(
    knownHotspotCaseStudy.spikeMutations
      .map((mutation) => mutation.aa_position)
      .filter((position) => omicronSet.has(position))
  );

  const priestHotspotsSeenInOmicron = uniquePositions(
    knownHotspotCaseStudy.priestTopSites
      .map((site) => site.aa_position)
      .filter((position) => omicronSet.has(position))
  );

  const supportedMutations = knownHotspotCaseStudy.spikeMutations.filter(
    (mutation) =>
      mutation.priest_score !== null && mutation.priest_score >= 0.6
  );

  const strongMutations = knownHotspotCaseStudy.spikeMutations.filter(
    (mutation) =>
      mutation.priest_score !== null && mutation.priest_score >= 0.7
  );

  const meanScoreValues = knownHotspotCaseStudy.spikeMutations
    .map((mutation) => mutation.priest_score)
    .filter((score) => score !== null && score !== undefined);

  const meanScore = meanScoreValues.length
    ? (
        meanScoreValues.reduce((sum, score) => sum + Number(score), 0) /
        meanScoreValues.length
      ).toFixed(3)
    : "—";

  const highlightSites = [
    ...sharedPositions.map((aaPosition) => ({
      aaPosition,
      genomeWindow: spikeAaToGenomeWindow(aaPosition),
      category: "shared",
      label: `Shared with Omicron: ${aaPosition}`,
      tone: {
        border: "#b91c1c",
        background: "rgba(248, 113, 113, 0.18)",
        chip: "border-red-200 bg-red-50 text-red-700",
      },
    })),
    ...priestHotspotsSeenInOmicron
      .filter((aaPosition) => !sharedPositions.includes(aaPosition))
      .map((aaPosition) => ({
        aaPosition,
        genomeWindow: spikeAaToGenomeWindow(aaPosition),
        category: "priest",
        label: `PRIEST hotspot later in Omicron: ${aaPosition}`,
        tone: {
          border: "#c2410c",
          background: "rgba(251, 146, 60, 0.18)",
          chip: "border-orange-200 bg-orange-50 text-orange-700",
        },
      })),
  ];

  return {
    ...knownHotspotCaseStudy,
    sharedPositions,
    priestHotspotsSeenInOmicron,
    supportedMutations,
    strongMutations,
    meanScore,
    highlightSites,
  };
};

const PositionChip = ({ tone, children }) => (
  <span
    className={`inline-flex items-center rounded-full border px-2.5 py-1 text-[11px] font-semibold ${tone}`}
  >
    {children}
  </span>
);

export default function KnownHotspotCaseStudy({
  className = "",
  isVisible = true,
  onToggle,
  viewMode = "genome",
}) {
  const caseStudy = buildKnownHotspotCaseStudy();
  const outsideSpike = viewMode === "outside-spike";

  if (!isVisible) {
    return (
      <div
        className={`rounded-full border border-slate-200 bg-white/92 px-3 py-2 shadow-[0_16px_40px_-28px_rgba(15,23,42,0.45)] backdrop-blur ${className}`}
      >
        <div className="flex items-center justify-between gap-3">
          <p className="text-[11px] font-semibold uppercase tracking-[0.2em] text-slate-500">
            Known Hotspot Case Study
          </p>
          {onToggle ? (
            <button
              type="button"
              onClick={onToggle}
              className="rounded-full border border-slate-300 px-3 py-1 text-[11px] font-semibold uppercase tracking-[0.14em] text-slate-600 transition hover:bg-slate-50"
            >
              Show
            </button>
          ) : null}
        </div>
      </div>
    );
  }

  return (
    <div
      className={`rounded-[18px] border border-slate-200 bg-white/90 px-3 py-2.5 shadow-[0_16px_40px_-28px_rgba(15,23,42,0.45)] backdrop-blur ${className}`}
    >
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-[11px] font-semibold uppercase tracking-[0.2em] text-slate-500">
          Known Hotspot Case Study
        </span>
        <span className="hidden h-1 w-1 rounded-full bg-slate-300 sm:inline-block" />
        <span className="text-sm font-medium text-slate-700">
          {caseStudy.selectedVariant}
        </span>
        <PositionChip tone="border-red-200 bg-red-50 text-red-700">
          Shared {caseStudy.sharedPositions.join(", ") || "None"}
        </PositionChip>
        <PositionChip tone="border-orange-200 bg-orange-50 text-orange-700">
          PRIEST to Omicron {caseStudy.priestHotspotsSeenInOmicron.join(", ") || "None"}
        </PositionChip>
        {outsideSpike ? (
          <PositionChip tone="border-amber-200 bg-amber-50 text-amber-700">
            Focus S to view markers
          </PositionChip>
        ) : (
          <span className="text-[11px] text-slate-500">
            Two highlighted Spike sites
          </span>
        )}
        {onToggle ? (
          <button
            type="button"
            onClick={onToggle}
            className="ml-auto rounded-full border border-slate-300 px-3 py-1 text-[11px] font-semibold uppercase tracking-[0.14em] text-slate-600 transition hover:bg-slate-50"
          >
            Hide
          </button>
        ) : null}
      </div>
    </div>
  );
}
