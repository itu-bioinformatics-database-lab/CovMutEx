import json
from unittest.mock import patch

import numpy as np
from django.test import Client, SimpleTestCase

from genome.delta_omicron_retrospective import (
    RetrospectiveCaseStudyValidationError,
    build_precomputed_variant_mutation_tuples,
    build_delta_omicron_retrospective_payload,
    compute_overlap_metrics,
    extract_spike_site_score_rows,
    load_delta_context_option_catalog,
    normalize_ranked_site_scores,
    parse_compact_nucleotide_mutation,
    rank_spike_site_rows,
    validate_delta_pre_omicron_context,
)


class DeltaOmicronRetrospectiveUnitTests(SimpleTestCase):
    def setUp(self):
        self.valid_context = {
            "node_id": "XBB.1.5|precomputed_consensus",
            "node_label": "Kraken (XBB.1.5)",
            "variant_label": "XBB.1.5",
            "variant_display_name": "Kraken (XBB.1.5)",
            "variant_nickname": "Kraken",
            "node_date": "2022-12-01",
            "emergence_label": "Dec 2022",
            "nextstrain_clade": None,
            "pangolin_lineage": "XBB.1.5",
            "country": None,
            "source_file": "XBB1.5.nucleotide-mutations.csv",
            "consensus_threshold": 0.5,
            "mutation_count": 147,
            "spike_site_count": 31,
        }

    def test_validate_delta_pre_omicron_context_rejects_unknown_contexts(self):
        with self.assertRaises(RetrospectiveCaseStudyValidationError):
            validate_delta_pre_omicron_context({"node_id": "unknown|node"})

        with self.assertRaises(RetrospectiveCaseStudyValidationError):
            validate_delta_pre_omicron_context({})

    def test_parse_compact_nucleotide_mutation_supports_substitutions_and_deletions(self):
        self.assertEqual(
            parse_compact_nucleotide_mutation("A56T"),
            {
                "reference_nucleotide": "A",
                "genome_position": 56,
                "alternate_nucleotide": "T",
                "mutation": "A56T",
            },
        )
        self.assertEqual(
            parse_compact_nucleotide_mutation("C21636-")["alternate_nucleotide"],
            "-",
        )

    def test_build_precomputed_variant_mutation_tuples_filters_reference_mismatches(self):
        with patch(
            "genome.delta_omicron_retrospective.load_precomputed_variant_mutation_catalog",
            return_value={
                "selected_mutations": (
                    {
                        "genome_position": 1,
                        "reference_nucleotide": "A",
                        "alternate_nucleotide": "T",
                        "mutation": "A1T",
                        "is_deletion": False,
                    },
                    {
                        "genome_position": 2,
                        "reference_nucleotide": "A",
                        "alternate_nucleotide": "-",
                        "mutation": "A2-",
                        "is_deletion": True,
                    },
                ),
                "consensus_threshold": 0.5,
                "deletion_count": 1,
                "total_rows": 10,
                "invalid_rows": 0,
                "below_threshold_count": 8,
                "selected_support_summary": {
                    "mutation_count": 2,
                    "total_proportion": 1.7,
                    "mean_proportion": 0.85,
                    "median_proportion": 0.85,
                    "min_proportion": 0.8,
                    "max_proportion": 0.9,
                },
            },
        ), patch(
            "genome.delta_omicron_retrospective.load_precomputed_variant_spike_site_set",
            return_value={
                "site_count": 1,
                "nucleotide_mutation_count": 1,
                "support_summary": {
                    "mutation_count": 1,
                    "total_proportion": 0.9,
                    "mean_proportion": 0.9,
                    "median_proportion": 0.9,
                    "min_proportion": 0.9,
                    "max_proportion": 0.9,
                },
            },
        ):
            result = build_precomputed_variant_mutation_tuples(
                reference_genome_sequence="ATG",
                node_context=self.valid_context,
            )

        self.assertEqual(result["mutation_count"], 1)
        self.assertEqual(result["reference_mismatch_count"], 1)
        self.assertEqual(result["deletion_count"], 1)
        self.assertEqual(result["mutation_tuples"][0][:3], (1, "A", "T"))

    def test_extract_spike_site_score_rows_maps_scores_to_one_based_spike_sites(self):
        predictions = np.array(
            [
                [0.1, 0.2, 0.3, 0.4],
                [0.5, 0.1, 0.2, 0.3],
                [0.3, 0.2, 0.1, 0.4],
                [0.2, 0.1, 0.5, 0.2],
                [0.3, 0.6, 0.1, 0.1],
                [0.2, 0.2, 0.4, 0.1],
            ]
        )

        rows = extract_spike_site_score_rows(
            predictions=predictions,
            reference_genome_sequence="ATGCCA",
            spike_start_genome_position=1,
            prediction_start_genome_position=1,
            reference_spike_amino_acids="MP",
        )

        self.assertEqual([row["aa_position"] for row in rows], [1, 2])
        self.assertEqual(rows[0]["codon_genome_positions"], [1, 2, 3])
        self.assertAlmostEqual(rows[0]["site_score"], 0.9333333333)
        self.assertAlmostEqual(rows[1]["site_score"], 0.8333333333)
        self.assertAlmostEqual(rows[0]["raw_site_score"], rows[0]["site_score"])

    def test_normalize_ranked_site_scores_scales_display_series_without_changing_order(self):
        normalized_rows = normalize_ranked_site_scores(
            [
                {"aa_position": 446, "rank": 1, "site_score": 1.4, "raw_site_score": 1.4},
                {"aa_position": 452, "rank": 2, "site_score": 1.1, "raw_site_score": 1.1},
                {"aa_position": 478, "rank": 3, "site_score": 0.8, "raw_site_score": 0.8},
            ]
        )

        self.assertEqual([row["aa_position"] for row in normalized_rows], [446, 452, 478])
        self.assertEqual([row["rank"] for row in normalized_rows], [1, 2, 3])
        self.assertAlmostEqual(normalized_rows[0]["site_score"], 1.0)
        self.assertAlmostEqual(normalized_rows[1]["site_score"], 0.5)
        self.assertAlmostEqual(normalized_rows[2]["site_score"], 0.0)
        self.assertEqual(normalized_rows[0]["raw_site_score"], 1.4)

    def test_rank_and_overlap_metrics_are_stable(self):
        ranked = rank_spike_site_rows(
            [
                {"aa_position": 452, "site_score": 0.9},
                {"aa_position": 446, "site_score": 0.9},
                {"aa_position": 478, "site_score": 0.2},
            ]
        )

        self.assertEqual([row["aa_position"] for row in ranked], [446, 452, 478])
        self.assertEqual([row["rank"] for row in ranked], [1, 2, 3])

        metrics = compute_overlap_metrics([446, 452], [417, 446, 478, 501])
        self.assertEqual(metrics["overlap_count"], 1)
        self.assertAlmostEqual(metrics["precision_at_k"], 0.5)
        self.assertAlmostEqual(metrics["recall_against_omicron_sites"], 0.25)
        self.assertEqual(metrics["overlap_positions"], [446])

    def test_context_catalog_surfaces_all_seven_precomputed_variants(self):
        catalog = load_delta_context_option_catalog()

        self.assertEqual(catalog["total_count"], 7)
        self.assertEqual(catalog["returned_count"], 7)
        self.assertEqual(len(catalog["options"]), 7)
        self.assertIn(
            "XBB.1.16|precomputed_consensus",
            {option["node_id"] for option in catalog["options"]},
        )
        self.assertIn(
            "BA.2.86|precomputed_consensus",
            {option["node_id"] for option in catalog["options"]},
        )

    def test_payload_contains_expected_schema_and_flags(self):
        predictions = np.array(
            [
                [0.1, 0.2, 0.3, 0.4],
                [0.5, 0.1, 0.2, 0.3],
                [0.3, 0.2, 0.1, 0.4],
                [0.2, 0.1, 0.5, 0.2],
                [0.3, 0.6, 0.1, 0.1],
                [0.2, 0.2, 0.4, 0.1],
            ]
        )

        with patch(
            "genome.delta_omicron_retrospective.load_precomputed_variant_spike_site_set",
            return_value={
                "set_name": "Test variant set",
                "variant_label": "Kraken (XBB.1.5)",
                "positions": (1, 4),
                "site_count": 2,
                "nucleotide_mutation_count": 3,
                "spike_mutations": (
                    {"aa_position": 1, "proportion": 0.9},
                    {"aa_position": 4, "proportion": 0.8},
                    {"aa_position": 4, "proportion": 0.7},
                ),
                "support_summary": {
                    "mutation_count": 3,
                    "total_proportion": 2.4,
                    "mean_proportion": 0.8,
                    "median_proportion": 0.8,
                    "min_proportion": 0.7,
                    "max_proportion": 0.9,
                },
                "sites": [
                    {"aa_position": 1, "mutation_count": 1},
                    {"aa_position": 4, "mutation_count": 2},
                ],
            },
        ), patch(
            "genome.delta_omicron_retrospective.validate_delta_pre_omicron_context",
            return_value=self.valid_context,
        ), patch(
            "genome.delta_omicron_retrospective.build_delta_context_option_metadata",
            return_value={
                "total_count": 1,
                "returned_count": 1,
                "selected_node_id": self.valid_context["node_id"],
                "options": [{"node_id": self.valid_context["node_id"]}],
            },
        ):
            payload = build_delta_omicron_retrospective_payload(
                predictions=predictions,
                reference_genome_sequence="ATGCCA",
                node_context=self.valid_context,
                selected_model="balanced_data_model",
                elapsed_day=0,
                top_k=1,
                spike_start_genome_position=1,
                prediction_start_genome_position=1,
                reference_spike_amino_acids="MP",
            )

        self.assertIn("metadata", payload)
        self.assertIn("score_series", payload)
        self.assertIn("top_k_positions", payload)
        self.assertIn("comparison_positions", payload)
        self.assertIn("omicron_positions", payload)
        self.assertIn("overlap_positions", payload)
        self.assertIn("metrics", payload)
        self.assertIn("ranked_rows", payload)
        self.assertEqual(payload["top_k_positions"], [1])
        self.assertEqual(payload["comparison_positions"], [1, 4])
        self.assertEqual(payload["overlap_positions"], [])
        self.assertEqual(payload["metrics"]["precision_at_k"], 0.0)
        self.assertEqual(payload["metrics"]["comparison_site_count"], 2)
        self.assertEqual(payload["metrics"]["found_mutation_count"], 0)
        self.assertTrue(
            all(
                {
                    "rank",
                    "is_top_k",
                    "is_comparison_site",
                    "is_omicron_site",
                    "is_overlap",
                }
                <= set(row.keys())
                for row in payload["ranked_rows"]
            )
        )

    def test_priest_payload_uses_identity_site_scores_and_selected_context_options(self):
        with patch(
            "genome.delta_omicron_retrospective.load_precomputed_variant_spike_site_set",
            return_value={
                "set_name": "Test variant set",
                "variant_label": "Kraken (XBB.1.5)",
                "positions": (2,),
                "site_count": 1,
                "nucleotide_mutation_count": 1,
                "spike_mutations": ({"aa_position": 2, "proportion": 0.7},),
                "support_summary": {
                    "mutation_count": 1,
                    "total_proportion": 0.7,
                    "mean_proportion": 0.7,
                    "median_proportion": 0.7,
                    "min_proportion": 0.7,
                    "max_proportion": 0.7,
                },
                "sites": [{"aa_position": 2, "mutation_count": 1}],
            },
        ), patch(
            "genome.delta_omicron_retrospective.load_priest_score_lookup",
            return_value={
                "period_scores": {"Q2-2021": {1: 0.7, 2: 0.2}},
                "global_scores": {},
                "available_periods": ("Q2-2021",),
            },
        ), patch(
            "genome.delta_omicron_retrospective.validate_delta_pre_omicron_context",
            return_value={**self.valid_context, "node_date": "2021-05-01"},
        ), patch(
            "genome.delta_omicron_retrospective.build_delta_context_option_metadata",
            return_value={
                "total_count": 1,
                "returned_count": 1,
                "selected_node_id": self.valid_context["node_id"],
                "options": [
                    {"node_id": self.valid_context["node_id"], "is_selected": True},
                ],
            },
        ):
            payload = build_delta_omicron_retrospective_payload(
                predictions=None,
                reference_genome_sequence=None,
                node_context={**self.valid_context, "node_date": "2021-05-01"},
                selected_model="PRIEST",
                elapsed_day=0,
                top_k=1,
                spike_start_genome_position=1,
                prediction_start_genome_position=1,
                reference_spike_amino_acids="MP",
            )

        self.assertEqual(payload["top_k_positions"], [1])
        self.assertEqual(payload["metrics"]["overlap_count"], 0)
        self.assertEqual(
            payload["metadata"]["scoring_context"]["display_normalization"],
            "identity_already_on_0_1_scale",
        )
        self.assertEqual(
            payload["metadata"]["scoring_context"]["score_source"],
            "PRIEST_spike_site_prevalence",
        )
        self.assertEqual(payload["ranked_rows"][0]["site_score"], 0.7)
        self.assertEqual(payload["ranked_rows"][0]["raw_site_score"], 0.7)
        self.assertIn("available_delta_context_nodes", payload["metadata"])


class DeltaOmicronRetrospectiveEndpointTests(SimpleTestCase):
    def setUp(self):
        self.client = Client()

    @patch("genome.views.build_delta_omicron_retrospective_payload")
    @patch("genome.views.predict_mutations", return_value=np.zeros((3, 4)))
    @patch("genome.views.get_model", return_value=object())
    @patch(
        "genome.views.resolve_model_name_and_path",
        return_value=("balanced_data_model", "/tmp/balanced_data_model.keras"),
    )
    @patch(
        "genome.views.resolve_delta_context",
        return_value={
            "node_id": "XBB.1.5|precomputed_consensus",
            "node_label": "Kraken (XBB.1.5)",
            "node_date": None,
            "nextstrain_clade": None,
            "pangolin_lineage": "XBB.1.5",
            "country": None,
            "source_file": "XBB1.5.nucleotide-mutations.csv",
            "consensus_threshold": 0.5,
        },
    )
    @patch(
        "genome.views.build_precomputed_variant_mutation_tuples",
        return_value={
            "mutation_tuples": [],
            "mutation_count": 0,
            "deletion_count": 0,
            "consensus_threshold": 0.5,
            "source_file": "XBB1.5.nucleotide-mutations.csv",
            "reference_mismatch_count": 0,
            "total_rows": 10,
            "invalid_rows": 0,
            "below_threshold_count": 10,
            "selected_support_summary": {
                "mutation_count": 0,
                "total_proportion": 0.0,
                "mean_proportion": None,
                "median_proportion": None,
                "min_proportion": None,
                "max_proportion": None,
            },
            "spike_site_count": 0,
            "spike_mutation_count": 0,
            "spike_support_summary": {
                "mutation_count": 0,
                "total_proportion": 0.0,
                "mean_proportion": None,
                "median_proportion": None,
                "min_proportion": None,
                "max_proportion": None,
            },
        },
    )
    @patch("genome.views.construct_variant_genome", return_value="ATG")
    @patch("genome.views.read_genome_sequence", return_value="ATG")
    def test_endpoint_returns_plot_ready_json_schema(
        self,
        _read_genome_sequence,
        _construct_variant_genome,
        _build_precomputed_variant_mutation_tuples,
        _resolve_delta_context,
        _resolve_model_name_and_path,
        _get_model,
        _predict_mutations,
        mock_build_payload,
    ):
        mock_build_payload.return_value = {
            "metadata": {
                "analysis_mode": "known hotspot case study",
                "scoring_context": {
                    "node_id": "XBB.1.5|precomputed_consensus"
                },
            },
            "score_series": [{"aa_position": 446, "site_score": 0.91}],
            "top_k_positions": [446],
            "omicron_positions": [417, 446, 478],
            "overlap_positions": [446],
            "metrics": {
                "overlap_count": 1,
                "precision_at_k": 1.0,
                "recall_against_omicron_sites": 0.3333333333,
                "top_k_count": 1,
                "omicron_site_count": 3,
            },
            "ranked_rows": [
                {
                    "aa_position": 446,
                    "rank": 1,
                    "site_score": 0.91,
                    "is_top_k": True,
                    "is_omicron_site": True,
                    "is_overlap": True,
                }
            ],
        }

        response = self.client.post(
            "/api/case-studies/delta-omicron-retrospective/",
            data=json.dumps(
                {
                    "selectedModel": "balanced_data_model",
                    "elapsedDay": 0,
                    "topK": 1,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload["metrics"]["overlap_count"], 1)
        self.assertEqual(payload["overlap_positions"], [446])
        self.assertEqual(
            payload["metadata"]["analysis_mode"],
            "known hotspot case study",
        )
