import React from "react";
import { act, fireEvent, render, screen, waitFor } from "@testing-library/react";

import DeltaOmicronRetrospectivePage from "./DeltaOmicronRetrospectivePage";

jest.mock("./RetrospectiveCaseStudyChart", () => (props) => (
  <div data-testid="mock-retrospective-chart">
    {JSON.stringify({
      showTopK: props.showTopK,
      showOmicron: props.showOmicron,
      showOverlap: props.showOverlap,
      pointCount: props.scoreSeries?.length || 0,
    })}
  </div>
));

const mockPayload = {
  metadata: {
    analysis_id: "known_hotspot_case_study_delta_omicron",
    analysis_label: "Known Hotspot Case Study",
    analysis_mode: "known hotspot case study",
    analysis_disclaimer:
      "This case study demonstrates how the CovMutEx explorer highlights mutational hotspots.",
    coordinate_system: {
      system: "Spike amino-acid positions",
      indexing: "1-based",
    },
    scoring_context: {
      node_id: "XBB.1.5|precomputed_consensus",
      node_date: "2022-12-01",
      emergence_label: "Dec 2022",
      nextstrain_clade: null,
      pangolin_lineage: "XBB.1.5",
      variant_display_name: "Kraken (XBB.1.5)",
      variant_nickname: "Kraken",
      source_file: "XBB1.5.nucleotide-mutations.csv",
      consensus_threshold: 0.5,
      selected_model: "balanced_data_model",
      elapsed_day: 0,
      display_normalization:
        "min_max_across_all_spike_sites_in_the_selected_context",
      precomputed_variant_summary: {
        mutation_count: 147,
        spike_site_count: 31,
        spike_support_summary: {
          mean_proportion: 0.94,
          median_proportion: 1.0,
        },
      },
    },
    comparison_set: {
      variant_label: "Kraken (XBB.1.5)",
      set_name: "Kraken (XBB.1.5) known Spike mutation sites",
      site_count: 4,
      site_positions: [417, 446, 478, 501],
    },
    available_delta_context_nodes: {
      total_count: 7,
      returned_count: 2,
      selected_node_id: "XBB.1.5|precomputed_consensus",
      options: [
        {
          node_id: "XBB.1.5|precomputed_consensus",
          label: "Kraken (XBB.1.5) | Dec 2022 | 31 Spike sites | 147 consensus mutations",
          is_default: true,
          is_selected: true,
        },
        {
          node_id: "XBB.1.16|precomputed_consensus",
          label: "Arcturus (XBB.1.16) | Jan 2023 | 33 Spike sites | 156 consensus mutations",
          is_default: false,
          is_selected: false,
        },
      ],
    },
    requested_top_k: 2,
    applied_top_k: 2,
  },
  metrics: {
    overlap_count: 1,
    precision_at_k: 0.5,
    recall_against_comparison_sites: 0.25,
    recall_against_omicron_sites: 0.25,
    top_k_count: 2,
    comparison_site_count: 4,
    omicron_site_count: 4,
    found_mutation_count: 2,
    found_total_proportion: 1.8,
    found_mean_proportion: 0.9,
    found_median_proportion: 0.9,
  },
  score_series: [
    {
      aa_position: 446,
      rank: 1,
      reference_aa: "G",
      site_score: 0.9123,
      raw_site_score: 1.487,
      codon_genome_positions: [22899, 22900, 22901],
      is_top_k: true,
      is_comparison_site: true,
      is_omicron_site: true,
      is_overlap: true,
    },
    {
      aa_position: 452,
      rank: 2,
      reference_aa: "L",
      site_score: 0.8123,
      raw_site_score: 1.221,
      codon_genome_positions: [22917, 22918, 22919],
      is_top_k: true,
      is_comparison_site: false,
      is_omicron_site: false,
      is_overlap: false,
    },
  ],
  top_k_positions: [446, 452],
  comparison_positions: [417, 446, 478, 501],
  omicron_positions: [417, 446, 478, 501],
  overlap_positions: [446],
  ranked_rows: [
    {
      aa_position: 446,
      rank: 1,
      reference_aa: "G",
      site_score: 0.9123,
      raw_site_score: 1.487,
      codon_genome_positions: [22899, 22900, 22901],
      is_top_k: true,
      is_comparison_site: true,
      is_omicron_site: true,
      is_overlap: true,
    },
    {
      aa_position: 452,
      rank: 2,
      reference_aa: "L",
      site_score: 0.8123,
      raw_site_score: 1.221,
      codon_genome_positions: [22917, 22918, 22919],
      is_top_k: true,
      is_comparison_site: false,
      is_omicron_site: false,
      is_overlap: false,
    },
  ],
};

describe("DeltaOmicronRetrospectivePage", () => {
  beforeEach(() => {
    global.fetch = jest.fn().mockResolvedValue({
      ok: true,
      json: async () => mockPayload,
    });
  });

  afterEach(() => {
    jest.resetAllMocks();
  });

  it("renders the case study description, metrics, toggles, and ranked rows", async () => {
    await act(async () => {
      render(<DeltaOmicronRetrospectivePage />);
    });

    expect(
      screen.getByText(/known hotspot case study/i)
    ).toBeInTheDocument();

    await screen.findByText("XBB.1.5|precomputed_consensus");
    await screen.findByText("Kraken (XBB.1.5)");

    expect(
      screen.getByText(/how to read this/i)
    ).toBeInTheDocument();
    expect(screen.getByText("50.0%")).toBeInTheDocument();
    expect(screen.getByText("25.0%")).toBeInTheDocument();
    expect(screen.getByText("446")).toBeInTheDocument();
    expect(screen.getByText("0.912")).toBeInTheDocument();
    expect(screen.getByText("1.487")).toBeInTheDocument();
    expect(
      screen.getByRole("columnheader", { name: /hotspot score/i })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("columnheader", { name: /raw score/i })
    ).toBeInTheDocument();
    expect(
      screen.getByRole("columnheader", { name: /variant site/i })
    ).toBeInTheDocument();
    expect(screen.getByText(/showing 2 selectable precomputed case-study nodes/i)).toBeInTheDocument();
    expect(screen.getByText("1.8")).toBeInTheDocument();
    expect(screen.getAllByText("90.0%").length).toBeGreaterThan(0);
    expect(
      screen.getByRole("option", {
        name: /Kraken \(XBB\.1\.5\).*31 Spike sites.*147 consensus mutations/i,
      })
    ).toBeInTheDocument();
    expect(screen.getByTestId("mock-retrospective-chart")).toHaveTextContent(
      '"pointCount":2'
    );
  });

  it("updates toggle state and forwards it to the chart", async () => {
    await act(async () => {
      render(<DeltaOmicronRetrospectivePage />);
    });

    await screen.findByText("XBB.1.5|precomputed_consensus");

    const topKButton = screen.getByRole("button", { name: /explorer hotspots/i });
    const omicronButton = screen.getByRole("button", {
      name: /known variant sites/i,
    });
    const overlapButton = screen.getByRole("button", {
      name: /overlap hits/i,
    });

    expect(topKButton).toHaveAttribute("aria-pressed", "true");
    expect(omicronButton).toHaveAttribute("aria-pressed", "true");
    expect(overlapButton).toHaveAttribute("aria-pressed", "true");

    fireEvent.click(topKButton);
    fireEvent.click(overlapButton);

    await waitFor(() => {
      expect(topKButton).toHaveAttribute("aria-pressed", "false");
      expect(overlapButton).toHaveAttribute("aria-pressed", "false");
    });

    expect(screen.getByTestId("mock-retrospective-chart")).toHaveTextContent(
      '"showTopK":false'
    );
    expect(screen.getByTestId("mock-retrospective-chart")).toHaveTextContent(
      '"showOverlap":false'
    );
    expect(screen.getByTestId("mock-retrospective-chart")).toHaveTextContent(
      '"showOmicron":true'
    );
  });

  it("submits PRIEST and the selected Delta context node on rerun", async () => {
    await act(async () => {
      render(<DeltaOmicronRetrospectivePage />);
    });

    await screen.findByText("XBB.1.5|precomputed_consensus");

    fireEvent.change(screen.getByLabelText(/scoring model/i), {
      target: { value: "PRIEST" },
    });
    fireEvent.change(screen.getByLabelText(/case-study lineage node id/i), {
      target: { value: "XBB.1.5|precomputed_consensus" },
    });

    expect(screen.getByLabelText(/elapsed days/i)).toBeDisabled();

    await act(async () => {
      fireEvent.click(screen.getByRole("button", { name: /re-run analysis/i }));
    });

    await waitFor(() => {
      expect(global.fetch).toHaveBeenCalledTimes(2);
    });

    expect(JSON.parse(global.fetch.mock.calls[1][1].body)).toMatchObject({
      selectedModel: "PRIEST",
      nodeId: "XBB.1.5|precomputed_consensus",
      elapsedDay: 0,
      topK: 2,
    });
  });
});
