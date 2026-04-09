# Known Hotspot Case Study — Delta → Omicron

## Purpose

This feature demonstrates the CovMutEx explorer's ability to highlight mutational hotspots.

When a user selects a Delta-era context, the explorer scores every Spike amino-acid position. This case study shows how many of the top-ranked hotspots correspond to sites that later acquired key mutations in Omicron BA.1.

It is framed as a feature of the explorer, not as a proof-of-prediction claim.

The case study asks:

- under a Delta-era, pre-Omicron scoring context
- which Spike amino-acid positions does the explorer flag as hotspots
- and how many of those hotspots overlap with known Omicron BA.1 Spike mutation sites

## Endpoint

- `POST /api/case-studies/delta-omicron-retrospective/`

### Request body

```json
{
  "selectedModel": "balanced_data_model",
  "elapsedDay": 0,
  "topK": 25,
  "nodeId": "optional Delta-era node id"
}
```

### Request notes

- `selectedModel` may be a CovMutEx mutation model or `PRIEST`.
- `nodeId` is optional. If omitted, the backend uses a built-in Delta-era default context node.
- The backend validates that the effective context is:
  - Delta-labeled
  - dated before `2021-11-01`
- The response also includes `metadata.available_delta_context_nodes`, which provides a curated dropdown-friendly subset of valid pre-Omicron Delta nodes plus the total eligible count.

## Coordinate conventions

- Canonical analysis coordinates are **Spike amino-acid positions**
- Indexing is **1-based**
- Spike genomic window is **21563-25384**
- Each amino-acid site corresponds to one Spike codon:
  - `codon_start = 21563 + (aa_position - 1) * 3`
  - `codon_genome_positions = [codon_start, codon_start + 1, codon_start + 2]`

## Score mapping

For CovMutEx mutation models, each Spike amino-acid site score is computed as:

- the raw A/T/G/C CovMutEx scores are kept in their original model output scale
- the non-reference score mass is taken at each of the three codon nucleotides
- those three nucleotide-level non-reference masses are averaged to produce `raw_site_score`
- `raw_site_score` is used for ranking
- a separate `site_score` field is then min-max normalized to `0-1` across all Spike amino-acid sites in the selected Delta-era context for plotting and display

This keeps the retrospective rank order tied to the underlying CovMutEx signal while still giving the frontend a consistent `0-1` plotting scale.

For `PRIEST`:

- site scores are taken directly from PRIEST Spike-site prevalence values for the quarter implied by the selected Delta context node date
- global PRIEST fallback scores are used when a quarter-specific value is unavailable
- `raw_site_score` and `site_score` are identical because PRIEST scores already lie on a `0-1` scale

## Response schema

```json
{
  "metadata": {
    "analysis_id": "known_hotspot_case_study_delta_omicron",
    "analysis_label": "Known Hotspot Case Study — Delta → Omicron",
    "analysis_mode": "known hotspot case study",
    "analysis_disclaimer": "Explorer feature framing text",
    "coordinate_system": {
      "system": "Spike amino-acid positions",
      "indexing": "1-based",
      "reference_sequence": "SARS-CoV-2 reference Spike coding sequence",
      "spike_gene": "S",
      "spike_genome_window": {
        "start": 21563,
        "end": 25384
      }
    },
    "scoring_context": {
      "node_id": "Delta-era node id",
      "node_date": "2021-06-22",
      "nextstrain_clade": "21A(Delta)",
      "pangolin_lineage": "AY.4",
      "country": "England",
      "pre_omicron_cutoff_date": "2021-11-01",
      "selected_model": "balanced_data_model",
      "elapsed_day": 0,
      "score_aggregation": "mean_non_reference_raw_nucleotide_probability_mass_across_the_three_nucleotides_of_each_spike_codon",
      "display_normalization": "min_max_across_all_spike_sites_in_the_selected_context"
    },
    "available_delta_context_nodes": {
      "total_count": 234284,
      "returned_count": 250,
      "selected_node_id": "Delta-era node id",
      "options": [
        {
          "node_id": "Delta-era node id",
          "label": "2021-06-22 | AY.4 | England | Delta-era node id",
          "is_default": true,
          "is_selected": true
        }
      ]
    },
    "comparison_set": {
      "variant_label": "Omicron BA.1",
      "set_name": "Curated Omicron BA.1 spike mutation site set",
      "site_count": 37,
      "site_positions": [67, 69, 70]
    },
    "requested_top_k": 25,
    "applied_top_k": 25,
    "ranked_row_count": 1273
  },
  "score_series": [
    {
      "aa_position": 446,
      "reference_aa": "G",
      "site_score": 0.912,
      "raw_site_score": 1.487,
      "rank": 1,
      "codon_start_genome_position": 22899,
      "codon_genome_positions": [22899, 22900, 22901],
      "nucleotide_rows": [],
      "is_top_k": true,
      "is_omicron_site": true,
      "is_overlap": true
    }
  ],
  "top_k_positions": [446, 452],
  "omicron_positions": [417, 446, 478, 501],
  "overlap_positions": [446],
  "metrics": {
    "overlap_count": 1,
    "precision_at_k": 0.5,
    "recall_against_omicron_sites": 0.25,
    "top_k_count": 2,
    "omicron_site_count": 4
  },
  "ranked_rows": [
    {
      "aa_position": 446,
      "rank": 1,
      "site_score": 0.912,
      "is_top_k": true,
      "is_omicron_site": true,
      "is_overlap": true
    }
  ]
}
```

## Response notes

- `score_series` is ordered by ascending Spike amino-acid position for plotting.
- `ranked_rows` is ordered by descending `raw_site_score`, then ascending amino-acid position.
- `site_score` is the display-ready `0-1` hotspot score used by the frontend chart and table.
- `raw_site_score` is the underlying explorer-derived site score used for ranking.
- When `selected_model` is `PRIEST`, `site_score` and `raw_site_score` are the same `0-1` PRIEST prevalence value.
- `precision_at_k` is `overlap_count / top_k_count` ("hotspot hit rate").
- `recall_against_omicron_sites` is `overlap_count / omicron_site_count` ("Omicron coverage").
