import React from "react";
import { render, screen } from "@testing-library/react";

import RetrospectiveCaseStudyChart from "./RetrospectiveCaseStudyChart";

jest.mock("recharts", () => {
  const React = require("react");

  const wrap = (testId) =>
    function WrappedComponent({ children, ...props }) {
      return (
        <div
          data-testid={testId}
          data-has-data={Array.isArray(props.data) ? "true" : "false"}
        >
          {children}
        </div>
      );
    };

  return {
    CartesianGrid: wrap("cartesian-grid"),
    ComposedChart: wrap("composed-chart"),
    Line: wrap("line"),
    ReferenceLine: wrap("reference-line"),
    ResponsiveContainer: wrap("responsive-container"),
    Tooltip: wrap("tooltip"),
    XAxis: wrap("x-axis"),
    YAxis: wrap("y-axis"),
    Scatter: function Scatter(props) {
      return (
        <div
          data-testid="scatter"
          data-fill={props.fill}
          data-count={(props.data || []).length}
        />
      );
    },
  };
});

const baseScoreSeries = [
  {
    aa_position: 446,
    rank: 1,
    site_score: 0.9,
    raw_site_score: 1.4,
    reference_aa: "G",
    codon_genome_positions: [1, 2, 3],
    is_top_k: true,
    is_omicron_site: false,
    is_overlap: false,
  },
  {
    aa_position: 478,
    rank: 2,
    site_score: 0.7,
    raw_site_score: 1.1,
    reference_aa: "T",
    codon_genome_positions: [4, 5, 6],
    is_top_k: false,
    is_omicron_site: true,
    is_overlap: false,
  },
];

describe("RetrospectiveCaseStudyChart", () => {
  it("does not render an overlap scatter when there are no overlap points", () => {
    render(
      <RetrospectiveCaseStudyChart
        scoreSeries={baseScoreSeries}
        showTopK
        showOmicron
        showOverlap
        selectedModel="balanced_data_model"
        displayNormalization="min_max_across_all_spike_sites_in_the_selected_context"
      />
    );

    const scatters = screen.getAllByTestId("scatter");
    expect(scatters).toHaveLength(2);
    expect(scatters.map((node) => node.getAttribute("data-fill"))).toEqual(
      expect.arrayContaining(["#0ea5e9", "#f59e0b"])
    );
    expect(
      scatters.some((node) => node.getAttribute("data-fill") === "#e11d48")
    ).toBe(false);
  });

  it("renders the overlap scatter only for overlap rows", () => {
    render(
      <RetrospectiveCaseStudyChart
        scoreSeries={[
          ...baseScoreSeries,
          {
            aa_position: 501,
            rank: 3,
            site_score: 0.6,
            raw_site_score: 1.0,
            reference_aa: "N",
            codon_genome_positions: [7, 8, 9],
            is_top_k: true,
            is_omicron_site: true,
            is_overlap: true,
          },
        ]}
        showTopK
        showOmicron
        showOverlap
        selectedModel="balanced_data_model"
        displayNormalization="min_max_across_all_spike_sites_in_the_selected_context"
      />
    );

    const overlapScatter = screen
      .getAllByTestId("scatter")
      .find((node) => node.getAttribute("data-fill") === "#e11d48");

    expect(overlapScatter).toBeInTheDocument();
    expect(overlapScatter).toHaveAttribute("data-count", "1");
  });
});
