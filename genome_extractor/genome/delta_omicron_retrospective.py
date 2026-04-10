import csv
import json
from datetime import date
from functools import lru_cache
from pathlib import Path
import re

import numpy as np

from .priest_annotations import (
    SPIKE_END,
    SPIKE_GENE,
    SPIKE_START,
    _load_reference_spike_amino_acids,
    _lookup_priest_score,
    load_priest_score_lookup,
    node_date_to_priest_period,
)


ANALYSIS_ID = "known_hotspot_case_study_post_2022_variants"
ANALYSIS_LABEL = "Known Hotspot Case Study — Post-2022 Variants"
ANALYSIS_MODE = "known hotspot case study"
ANALYSIS_DISCLAIMER = (
    "This case study is presented as a feature demonstration of the CovMutEx explorer, "
    "not as a proof of predictive accuracy. Under each selected variant context, the "
    "explorer scores every Spike amino-acid position for mutational pressure using its "
    "trained model. The resulting hotspot rankings are then compared against the same "
    "lineage's known Spike mutation sites — derived from sequencing proportion data — to "
    "assess how well the tool identifies positions of genuine biological relevance. "
    "Agreement between explorer-ranked hotspots and known mutation sites illustrates "
    "the tool's capacity to surface biologically meaningful signals, independent of "
    "any temporal prediction claim."
)
DEFAULT_TOP_K = 50
DEFAULT_ELAPSED_DAY = 0
DEFAULT_MODEL_NAME = "balanced_data_model"
PRIEST_MODEL_NAME = "PRIEST"
DEFAULT_DELTA_CONTEXT_NODE_ID = "XFG|precomputed_consensus"
DELTA_CONTEXT_OPTION_LIMIT = 15
PRECOMPUTED_MUTATION_CONSENSUS_THRESHOLD = 0.5  # majority consensus: only mutations present in ≥50% of reported sequences define the lineage
NUCLEOTIDES = ("A", "T", "G", "C")

BASE_DIR = Path(__file__).resolve().parent
OMICRON_SPIKE_SITE_PATH = BASE_DIR / "omicron_spike_sites.json"
PRECOMPUTED_VARIANT_DIR = BASE_DIR.parent
COMPACT_NUCLEOTIDE_MUTATION_RE = re.compile(r"^(?P<ref>[ACGT])(?P<position>\d+)(?P<alt>[ACGT-])$")


def _build_precomputed_variant_context(
    *,
    pangolin_lineage,
    variant_nickname,
    emergence_label,
    node_date,
    source_file,
    sort_order,
    hidden_from_ui=False,
):
    node_id = f"{pangolin_lineage}|precomputed_consensus"
    variant_display_name = (
        f"{variant_nickname} ({pangolin_lineage})"
        if variant_nickname
        else pangolin_lineage
    )
    # Compute elapsed days from SARS-CoV-2 reference origin (2019-12-01)
    _SARS_COV2_ORIGIN = date(2019, 12, 1)
    _emergence_date = date.fromisoformat(node_date)
    elapsed_days = (_emergence_date - _SARS_COV2_ORIGIN).days
    return (
        node_id,
        {
            "node_id": node_id,
            "variant_label": pangolin_lineage,
            "variant_display_name": variant_display_name,
            "variant_nickname": variant_nickname,
            "node_date": node_date,
            "emergence_label": emergence_label,
            "nextstrain_clade": None,
            "pangolin_lineage": pangolin_lineage,
            "country": None,
            "source_file": source_file,
            "csv_path": PRECOMPUTED_VARIANT_DIR / source_file,
            "consensus_threshold": PRECOMPUTED_MUTATION_CONSENSUS_THRESHOLD,
            "sort_order": sort_order,
            "elapsed_days": elapsed_days,
            "hidden_from_ui": hidden_from_ui,
        },
    )


PRECOMPUTED_VARIANT_CONTEXTS = dict(
    [
        _build_precomputed_variant_context(
            pangolin_lineage="XBB.1.5",
            variant_nickname="Kraken",
            emergence_label="Oct 2022 · WHO VOI Jan 2023",
            node_date="2022-10-01",
            source_file="XBB1.5.nucleotide-mutations.csv",
            sort_order=1,
        ),
        _build_precomputed_variant_context(
            pangolin_lineage="XBB.1.16",
            variant_nickname="Arcturus",
            emergence_label="Jan 2023 · WHO VOI Apr 2023",
            node_date="2023-01-01",
            source_file="XBB.1.16nucleotide-mutations.csv",
            sort_order=2,
        ),
        _build_precomputed_variant_context(
            pangolin_lineage="BA.2.86",
            variant_nickname="Pirola",
            emergence_label="Jul 2023 · WHO VUM Aug 2023",
            node_date="2023-07-01",
            source_file="BA2.86.nucleotide-mutations.csv",
            sort_order=3,
        ),
        _build_precomputed_variant_context(
            pangolin_lineage="KP.2",
            variant_nickname="FLiRT",
            emergence_label="Apr 2024",
            node_date="2024-04-01",
            source_file="KP.2.nucleotide-mutations.csv",
            sort_order=4,
        ),
        _build_precomputed_variant_context(
            pangolin_lineage="KP.3",
            variant_nickname="FLiRT",
            emergence_label="May 2024",
            node_date="2024-05-01",
            source_file="KP.3nucleotide-mutations.csv",
            sort_order=5,
        ),
        _build_precomputed_variant_context(
            pangolin_lineage="XFG",
            variant_nickname="Stratus",
            emergence_label="Jan 2025 · WHO Jun 2025",
            node_date="2025-01-01",
            source_file="XFG.nucleotide-mutations.csv",
            sort_order=6,
        ),
        _build_precomputed_variant_context(
            pangolin_lineage="NB.1.8.1",
            variant_nickname="Nimbus",
            emergence_label="Jan 2025 · WHO May 2025",
            node_date="2025-01-01",
            source_file="NB.1.8.1.nucleotide-mutations.csv",
            sort_order=7,
        ),
    ]
)


class CaseStudyValidationError(ValueError):
    """Raised when the requested retrospective case-study context is invalid."""


def is_priest_model(selected_model):
    return (selected_model or "").upper() == PRIEST_MODEL_NAME


@lru_cache(maxsize=1)
def load_curated_omicron_spike_site_set():
    with OMICRON_SPIKE_SITE_PATH.open(encoding="utf-8") as handle:
        payload = json.load(handle)

    sites = [
        {
            "aa_position": int(site["aa_position"]),
            "mutation": site.get("mutation"),
            "event_type": site.get("event_type", "substitution"),
        }
        for site in payload.get("sites", [])
    ]
    positions = tuple(sorted({site["aa_position"] for site in sites}))

    return {
        **payload,
        "sites": sites,
        "positions": positions,
    }


def parse_compact_nucleotide_mutation(mutation_text):
    match = COMPACT_NUCLEOTIDE_MUTATION_RE.fullmatch((mutation_text or "").strip())
    if not match:
        raise CaseStudyValidationError(
            f"Unsupported compact nucleotide mutation: {mutation_text!r}"
        )

    return {
        "reference_nucleotide": match.group("ref"),
        "genome_position": int(match.group("position")),
        "alternate_nucleotide": match.group("alt"),
        "mutation": match.group(0),
    }


def _get_precomputed_variant_definition(node_id):
    definition = PRECOMPUTED_VARIANT_CONTEXTS.get(node_id)
    if definition is None:
        raise CaseStudyValidationError(
            "The requested case-study node was not found in the precomputed variant catalog."
        )
    return definition


def _precomputed_candidate_sort_key(candidate):
    return (
        float(candidate["proportion"]),
        int(candidate["count"]),
        candidate["mutation"],
    )


def _summarize_proportions(proportions):
    values = [float(value) for value in proportions]
    if not values:
        return {
            "mutation_count": 0,
            "total_proportion": 0.0,
            "mean_proportion": None,
            "median_proportion": None,
            "min_proportion": None,
            "max_proportion": None,
        }

    proportion_array = np.asarray(values, dtype=float)
    return {
        "mutation_count": len(values),
        "total_proportion": float(np.sum(proportion_array)),
        "mean_proportion": float(np.mean(proportion_array)),
        "median_proportion": float(np.median(proportion_array)),
        "min_proportion": float(np.min(proportion_array)),
        "max_proportion": float(np.max(proportion_array)),
    }


def _genome_position_to_spike_aa_position(genome_position):
    if genome_position < SPIKE_START or genome_position > SPIKE_END:
        return None
    return ((genome_position - SPIKE_START) // 3) + 1


@lru_cache(maxsize=8)
def load_precomputed_variant_mutation_catalog(node_id):
    definition = _get_precomputed_variant_definition(node_id)
    csv_path = definition["csv_path"]
    if not csv_path.exists():
        raise CaseStudyValidationError(
            f"Precomputed mutation CSV not found: {csv_path.name}"
        )

    best_rows_by_position = {}
    total_rows = 0
    invalid_rows = 0

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            total_rows += 1
            mutation_text = (row.get("mutation") or "").strip()
            if not mutation_text:
                invalid_rows += 1
                continue

            try:
                parsed_mutation = parse_compact_nucleotide_mutation(mutation_text)
            except CaseStudyValidationError:
                invalid_rows += 1
                continue

            try:
                proportion = float(row.get("proportion") or 0.0)
            except (TypeError, ValueError):
                proportion = 0.0
            try:
                count = int(float(row.get("count") or 0))
            except (TypeError, ValueError):
                count = 0

            candidate = {
                **parsed_mutation,
                "proportion": proportion,
                "count": count,
                "jaccard": row.get("jaccard"),
                "is_deletion": parsed_mutation["alternate_nucleotide"] == "-",
            }
            current_best = best_rows_by_position.get(candidate["genome_position"])
            if current_best is None or _precomputed_candidate_sort_key(candidate) > _precomputed_candidate_sort_key(current_best):
                best_rows_by_position[candidate["genome_position"]] = candidate

    consensus_threshold = float(definition["consensus_threshold"])
    selected_mutations = []
    below_threshold_count = 0
    for genome_position in sorted(best_rows_by_position):
        candidate = best_rows_by_position[genome_position]
        if candidate["proportion"] < consensus_threshold:
            below_threshold_count += 1
            continue
        selected_mutations.append(candidate)

    selected_support_summary = _summarize_proportions(
        mutation["proportion"] for mutation in selected_mutations
    )

    return {
        "node_id": node_id,
        "total_rows": total_rows,
        "invalid_rows": invalid_rows,
        "candidate_position_count": len(best_rows_by_position),
        "consensus_threshold": consensus_threshold,
        "selected_mutations": tuple(selected_mutations),
        "mutation_count": len(selected_mutations),
        "deletion_count": sum(1 for mutation in selected_mutations if mutation["is_deletion"]),
        "below_threshold_count": below_threshold_count,
        "selected_support_summary": selected_support_summary,
    }


@lru_cache(maxsize=8)
def load_precomputed_variant_spike_site_set(node_id):
    definition = _get_precomputed_variant_definition(node_id)
    mutation_catalog = load_precomputed_variant_mutation_catalog(node_id)
    spike_site_rows = {}
    spike_mutations = []

    for mutation in mutation_catalog["selected_mutations"]:
        aa_position = _genome_position_to_spike_aa_position(
            mutation["genome_position"]
        )
        if aa_position is None:
            continue

        enriched_mutation = {
            **mutation,
            "aa_position": aa_position,
        }
        spike_mutations.append(enriched_mutation)
        site_row = spike_site_rows.setdefault(
            aa_position,
            {
                "aa_position": aa_position,
                "mutation_labels": [],
                "proportions": [],
            },
        )
        site_row["mutation_labels"].append(mutation["mutation"])
        site_row["proportions"].append(mutation["proportion"])

    site_rows = []
    for aa_position in sorted(spike_site_rows):
        site_row = spike_site_rows[aa_position]
        support_summary = _summarize_proportions(site_row["proportions"])
        site_rows.append(
            {
                "aa_position": aa_position,
                "mutation_labels": tuple(site_row["mutation_labels"]),
                **support_summary,
            }
        )

    return {
        "variant_label": definition.get("variant_display_name")
        or definition.get("pangolin_lineage")
        or definition.get("variant_label"),
        "pangolin_lineage": definition.get("pangolin_lineage"),
        "variant_nickname": definition.get("variant_nickname"),
        "set_name": (
            f"{definition.get('variant_display_name') or definition.get('pangolin_lineage')} "
            "known Spike mutation sites"
        ),
        "positions": tuple(sorted(spike_site_rows)),
        "sites": tuple(site_rows),
        "spike_mutations": tuple(spike_mutations),
        "site_count": len(site_rows),
        "nucleotide_mutation_count": len(spike_mutations),
        "support_summary": _summarize_proportions(
            mutation["proportion"] for mutation in spike_mutations
        ),
    }


def build_precomputed_variant_mutation_tuples(reference_genome_sequence, node_context):
    mutation_catalog = load_precomputed_variant_mutation_catalog(node_context["node_id"])
    spike_site_set = load_precomputed_variant_spike_site_set(node_context["node_id"])
    mutation_tuples = []
    reference_mismatch_count = 0

    for mutation in mutation_catalog["selected_mutations"]:
        genome_position = mutation["genome_position"]
        reference_index = genome_position - 1
        if reference_index < 0 or reference_index >= len(reference_genome_sequence):
            reference_mismatch_count += 1
            continue

        reference_nucleotide = reference_genome_sequence[reference_index].upper()
        if reference_nucleotide != mutation["reference_nucleotide"]:
            reference_mismatch_count += 1
            continue

        mutation_tuples.append(
            (
                genome_position,
                mutation["reference_nucleotide"],
                mutation["alternate_nucleotide"],
                ((genome_position - 1) // 3) + 1,
                mutation["mutation"],
            )
        )

    return {
        "mutation_tuples": mutation_tuples,
        "mutation_count": len(mutation_tuples),
        "reference_mismatch_count": reference_mismatch_count,
        "consensus_threshold": mutation_catalog["consensus_threshold"],
        "deletion_count": mutation_catalog["deletion_count"],
        "source_file": node_context.get("source_file"),
        "total_rows": mutation_catalog["total_rows"],
        "invalid_rows": mutation_catalog["invalid_rows"],
        "below_threshold_count": mutation_catalog["below_threshold_count"],
        "selected_support_summary": mutation_catalog["selected_support_summary"],
        "spike_site_count": spike_site_set["site_count"],
        "spike_mutation_count": spike_site_set["nucleotide_mutation_count"],
        "spike_support_summary": spike_site_set["support_summary"],
    }


@lru_cache(maxsize=1)
def load_node_context_lookup():
    lookup = {}
    for node_id, definition in PRECOMPUTED_VARIANT_CONTEXTS.items():
        mutation_catalog = load_precomputed_variant_mutation_catalog(node_id)
        spike_site_set = load_precomputed_variant_spike_site_set(node_id)
        lookup[node_id] = {
            "node_id": node_id,
            "node_label": definition.get("variant_display_name")
            or definition["variant_label"],
            "variant_label": definition["variant_label"],
            "variant_display_name": definition.get("variant_display_name"),
            "variant_nickname": definition.get("variant_nickname"),
            "node_date": definition.get("node_date"),
            "emergence_label": definition.get("emergence_label"),
            "nextstrain_clade": definition.get("nextstrain_clade"),
            "pangolin_lineage": definition.get("pangolin_lineage"),
            "country": definition.get("country"),
            "source_type": "precomputed_variant_csv",
            "source_file": definition["source_file"],
            "hidden_from_ui": definition.get("hidden_from_ui", False),
            "consensus_threshold": mutation_catalog["consensus_threshold"],
            "mutation_count": mutation_catalog["mutation_count"],
            "deletion_count": mutation_catalog["deletion_count"],
            "below_threshold_count": mutation_catalog["below_threshold_count"],
            "selected_support_summary": mutation_catalog["selected_support_summary"],
            "spike_site_count": spike_site_set["site_count"],
            "spike_mutation_count": spike_site_set["nucleotide_mutation_count"],
            "spike_support_summary": spike_site_set["support_summary"],
            "sort_order": definition.get("sort_order", 999),
        }
    return lookup


def validate_delta_pre_omicron_context(node_context, cutoff=None):
    if not node_context:
        raise CaseStudyValidationError(
            "A case-study node context is required for this analysis."
        )

    node_id = (node_context.get("node_id") or "").strip()
    if not node_id:
        raise CaseStudyValidationError(
            "The selected case-study node is missing a node identifier."
        )

    catalog_context = load_node_context_lookup().get(node_id)
    if catalog_context is None:
        raise CaseStudyValidationError(
            "The selected case-study node is not available in the precomputed variant catalog."
        )

    return {
        **catalog_context,
        **{
            key: value
            for key, value in node_context.items()
            if key not in {"mutation_count", "deletion_count", "below_threshold_count"}
        },
    }


def resolve_delta_context(node_id=None):
    resolved_node_id = node_id or DEFAULT_DELTA_CONTEXT_NODE_ID
    return validate_delta_pre_omicron_context({"node_id": resolved_node_id})


def _build_delta_context_option(node_context, is_default=False):
    mutation_count = node_context.get("mutation_count", 0)
    label_parts = [
        node_context.get("variant_display_name")
        or node_context.get("node_label")
        or node_context.get("pangolin_lineage")
        or node_context.get("nextstrain_clade")
        or "Unknown case-study node",
    ]
    if node_context.get("emergence_label"):
        label_parts.append(node_context["emergence_label"])
    label_parts.extend(
        [
            f"{node_context.get('spike_site_count', 0)} Spike sites",
            f"{mutation_count} consensus mutations",
        ]
    )

    return {
        "node_id": node_context["node_id"],
        "label": " | ".join(label_parts),
        "node_date": node_context.get("node_date"),
        "emergence_label": node_context.get("emergence_label"),
        "variant_nickname": node_context.get("variant_nickname"),
        "variant_display_name": node_context.get("variant_display_name"),
        "nextstrain_clade": node_context.get("nextstrain_clade"),
        "pangolin_lineage": node_context.get("pangolin_lineage"),
        "country": node_context.get("country"),
        "source_file": node_context.get("source_file"),
        "consensus_threshold": node_context.get("consensus_threshold"),
        "mutation_count": mutation_count,
        "spike_site_count": node_context.get("spike_site_count", 0),
        "is_default": is_default,
    }


@lru_cache(maxsize=4)
def load_delta_context_option_catalog(limit=DELTA_CONTEXT_OPTION_LIMIT):
    valid_contexts = [
        ctx for ctx in load_node_context_lookup().values()
        if not ctx.get("hidden_from_ui", False)
    ]
    sorted_contexts = sorted(
        valid_contexts,
        key=lambda item: (
            item.get("sort_order", 999),
            item.get("pangolin_lineage") or "",
            item["node_id"],
        ),
    )
    default_context = resolve_delta_context(DEFAULT_DELTA_CONTEXT_NODE_ID)
    options = [_build_delta_context_option(default_context, is_default=True)]
    seen_node_ids = {default_context["node_id"]}

    for node_context in sorted_contexts:
        if len(options) >= limit:
            break
        if node_context["node_id"] in seen_node_ids:
            continue

        options.append(_build_delta_context_option(node_context))
        seen_node_ids.add(node_context["node_id"])

    return {
        "total_count": len(valid_contexts),
        "returned_count": len(options),
        "limit": limit,
        "options": options,
    }


def build_delta_context_option_metadata(node_context, limit=DELTA_CONTEXT_OPTION_LIMIT):
    catalog = load_delta_context_option_catalog(limit=limit)
    selected_node_id = node_context["node_id"]
    options = list(catalog["options"])

    if selected_node_id not in {option["node_id"] for option in options}:
        insert_at = 1 if options else 0
        options.insert(
            insert_at,
            _build_delta_context_option(
                node_context,
                is_default=selected_node_id == DEFAULT_DELTA_CONTEXT_NODE_ID,
            ),
        )

    return {
        **catalog,
        "returned_count": len(options),
        "selected_node_id": selected_node_id,
        "options": [
            {
                **option,
                "is_selected": option["node_id"] == selected_node_id,
            }
            for option in options
        ],
    }


def spike_aa_position_to_genome_positions(aa_position, spike_start_genome_position=SPIKE_START):
    codon_start = spike_start_genome_position + ((aa_position - 1) * 3)
    return [codon_start, codon_start + 1, codon_start + 2]


def extract_spike_site_score_rows(
    predictions,
    reference_genome_sequence,
    spike_start_genome_position=SPIKE_START,
    prediction_start_genome_position=SPIKE_START,
    reference_spike_amino_acids=None,
):
    prediction_array = np.asarray(predictions)
    if prediction_array.ndim != 2 or prediction_array.shape[1] != 4:
        raise CaseStudyValidationError(
            "Predictions must be a two-dimensional array with four nucleotide scores per position."
        )

    reference_spike = reference_spike_amino_acids or _load_reference_spike_amino_acids()
    site_rows = []

    for aa_position in range(1, len(reference_spike) + 1):
        genome_positions = spike_aa_position_to_genome_positions(
            aa_position,
            spike_start_genome_position=spike_start_genome_position,
        )

        nucleotide_rows = []
        for genome_position in genome_positions:
            prediction_index = genome_position - prediction_start_genome_position
            if prediction_index < 0 or prediction_index >= len(prediction_array):
                continue

            reference_index = genome_position - 1
            if reference_index < 0 or reference_index >= len(reference_genome_sequence):
                continue

            raw_scores = prediction_array[prediction_index]
            reference_nt = reference_genome_sequence[reference_index].upper()
            raw_candidate_scores = {
                nucleotide: float(raw_scores[index])
                for index, nucleotide in enumerate(NUCLEOTIDES)
            }
            raw_score_total = float(sum(raw_candidate_scores.values()))
            normalized_candidate_scores = {
                nucleotide: (
                    score / raw_score_total if raw_score_total > 0 else 0.0
                )
                for nucleotide, score in raw_candidate_scores.items()
            }

            raw_non_reference_scores = [
                score
                for nucleotide, score in raw_candidate_scores.items()
                if nucleotide != reference_nt
            ]
            normalized_non_reference_scores = [
                score
                for nucleotide, score in normalized_candidate_scores.items()
                if nucleotide != reference_nt
            ]
            raw_mutation_mass = float(sum(raw_non_reference_scores))
            normalized_mutation_mass = float(sum(normalized_non_reference_scores))

            nucleotide_rows.append(
                {
                    "genome_position": genome_position,
                    "reference_nucleotide": reference_nt,
                    "candidate_scores": raw_candidate_scores,
                    "raw_mutation_mass": raw_mutation_mass,
                    "mutation_mass": normalized_mutation_mass,
                    "max_non_reference_score": float(
                        max(raw_non_reference_scores, default=0.0)
                    ),
                }
            )

        if len(nucleotide_rows) != 3:
            continue

        raw_site_score = float(
            sum(row["raw_mutation_mass"] for row in nucleotide_rows) / len(nucleotide_rows)
        )
        site_rows.append(
            {
                "aa_position": aa_position,
                "reference_aa": reference_spike[aa_position - 1],
                "raw_site_score": raw_site_score,
                "site_score": raw_site_score,
                "codon_start_genome_position": genome_positions[0],
                "codon_genome_positions": genome_positions,
                "nucleotide_rows": nucleotide_rows,
            }
        )

    return site_rows


def extract_priest_spike_site_score_rows(
    node_context,
    spike_start_genome_position=SPIKE_START,
    reference_spike_amino_acids=None,
):
    reference_spike = reference_spike_amino_acids or _load_reference_spike_amino_acids()
    all_positions = tuple(range(1, len(reference_spike) + 1))
    lookup = load_priest_score_lookup(all_positions)
    priest_period = node_date_to_priest_period(node_context.get("node_date"))
    site_rows = []

    for aa_position in all_positions:
        genome_positions = spike_aa_position_to_genome_positions(
            aa_position,
            spike_start_genome_position=spike_start_genome_position,
        )
        priest_score, score_source = _lookup_priest_score(
            aa_position=aa_position,
            period=priest_period,
            period_scores=lookup["period_scores"],
            global_scores=lookup["global_scores"],
        )
        resolved_score = float(priest_score or 0.0)

        site_rows.append(
            {
                "aa_position": aa_position,
                "reference_aa": reference_spike[aa_position - 1],
                "raw_site_score": resolved_score,
                "site_score": resolved_score,
                "codon_start_genome_position": genome_positions[0],
                "codon_genome_positions": genome_positions,
                "nucleotide_rows": [],
                "score_source": score_source or "missing",
                "priest_period": priest_period,
            }
        )

    return site_rows, {
        "score_aggregation": (
            "priest_period_specific_spike_site_prevalence_with_global_fallback"
        ),
        "display_normalization": "identity_already_on_0_1_scale",
        "score_source": "PRIEST_spike_site_prevalence",
        "priest_period": priest_period,
        "priest_available_periods": list(lookup["available_periods"]),
    }


def rank_spike_site_rows(site_rows):
    ranked_rows = []

    for rank, row in enumerate(
        sorted(
            site_rows,
            key=lambda item: (-item.get("raw_site_score", item["site_score"]), item["aa_position"]),
        ),
        start=1,
    ):
        ranked_rows.append({**row, "rank": rank})

    return ranked_rows


def normalize_ranked_site_scores(site_rows):
    if not site_rows:
        return []

    raw_scores = [row.get("raw_site_score", row["site_score"]) for row in site_rows]
    min_score = min(raw_scores)
    max_score = max(raw_scores)
    score_range = max_score - min_score

    normalized_rows = []
    for row in site_rows:
        raw_score = row.get("raw_site_score", row["site_score"])
        normalized_score = (
            (raw_score - min_score) / score_range
            if score_range > 0
            else 0.5
        )
        normalized_rows.append(
            {
                **row,
                "site_score": float(normalized_score),
            }
        )

    return normalized_rows


def build_ranked_spike_site_rows(
    selected_model,
    node_context,
    predictions=None,
    reference_genome_sequence=None,
    spike_start_genome_position=SPIKE_START,
    prediction_start_genome_position=SPIKE_START,
    reference_spike_amino_acids=None,
):
    if is_priest_model(selected_model):
        site_rows, scoring_metadata = extract_priest_spike_site_score_rows(
            node_context=node_context,
            spike_start_genome_position=spike_start_genome_position,
            reference_spike_amino_acids=reference_spike_amino_acids,
        )
        return rank_spike_site_rows(site_rows), scoring_metadata

    if predictions is None or reference_genome_sequence is None:
        raise CaseStudyValidationError(
            "CovMutEx score extraction requires both predictions and a reference genome sequence."
        )

    site_rows = extract_spike_site_score_rows(
        predictions=predictions,
        reference_genome_sequence=reference_genome_sequence,
        spike_start_genome_position=spike_start_genome_position,
        prediction_start_genome_position=prediction_start_genome_position,
        reference_spike_amino_acids=reference_spike_amino_acids,
    )
    return normalize_ranked_site_scores(rank_spike_site_rows(site_rows)), {
        "score_aggregation": (
            "mean_non_reference_raw_nucleotide_probability_mass_"
            "across_the_three_nucleotides_of_each_spike_codon"
        ),
        "display_normalization": (
            "min_max_across_all_spike_sites_in_the_selected_context"
        ),
        "score_source": "CovMutEx_mutation_model",
    }


def summarize_found_mutation_support(node_context, aa_positions):
    spike_site_set = load_precomputed_variant_spike_site_set(node_context["node_id"])
    target_positions = set(aa_positions)
    matched_proportions = [
        mutation["proportion"]
        for mutation in spike_site_set["spike_mutations"]
        if mutation["aa_position"] in target_positions
    ]
    return _summarize_proportions(matched_proportions)


def compute_overlap_metrics(top_k_positions, omicron_positions, proximity_window=3):
    top_k_set = set(top_k_positions)
    omicron_position_list = sorted(set(omicron_positions))
    omicron_set = set(omicron_position_list)

    # --- Exact overlap ---
    exact_overlap_positions = sorted(top_k_set & omicron_set)
    exact_overlap_count = len(exact_overlap_positions)

    # --- Proximity overlap (within ±proximity_window amino acids) ---
    proximity_hits = []
    for hotspot_pos in sorted(top_k_set):
        for omicron_pos in omicron_position_list:
            if abs(hotspot_pos - omicron_pos) <= proximity_window:
                proximity_hits.append({
                    "hotspot_position": hotspot_pos,
                    "omicron_position": omicron_pos,
                    "distance": abs(hotspot_pos - omicron_pos),
                })
                break  # count each hotspot at most once

    proximity_overlap_count = len(proximity_hits)
    proximity_overlap_positions = sorted(
        {hit["hotspot_position"] for hit in proximity_hits}
    )

    precision_at_k = exact_overlap_count / len(top_k_positions) if top_k_positions else 0.0
    recall_against_omicron_sites = (
        exact_overlap_count / len(omicron_positions) if omicron_positions else 0.0
    )
    proximity_precision_at_k = (
        proximity_overlap_count / len(top_k_positions) if top_k_positions else 0.0
    )
    omicron_positions_covered_by_proximity = set()
    for hit in proximity_hits:
        omicron_positions_covered_by_proximity.add(hit["omicron_position"])
    proximity_recall = (
        len(omicron_positions_covered_by_proximity) / len(omicron_positions)
        if omicron_positions else 0.0
    )

    return {
        "overlap_count": exact_overlap_count,
        "precision_at_k": precision_at_k,
        "recall_against_omicron_sites": recall_against_omicron_sites,
        "top_k_count": len(top_k_positions),
        "omicron_site_count": len(omicron_positions),
        "overlap_positions": exact_overlap_positions,
        "proximity_window": proximity_window,
        "proximity_overlap_count": proximity_overlap_count,
        "proximity_precision_at_k": proximity_precision_at_k,
        "proximity_recall": proximity_recall,
        "proximity_overlap_positions": proximity_overlap_positions,
        "proximity_hits": proximity_hits,
    }


def build_top_k_metric_sweep(ranked_rows, comparison_positions, proximity_window=3):
    comparison_position_list = sorted(set(comparison_positions))
    comparison_position_set = set(comparison_position_list)
    comparison_site_count = len(comparison_position_list)
    exact_overlap_count = 0
    proximity_overlap_count = 0
    exact_overlap_positions_seen = set()
    covered_comparison_positions = set()
    sweep_rows = []

    for rank, row in enumerate(ranked_rows, start=1):
        aa_position = row["aa_position"]

        if (
            aa_position in comparison_position_set
            and aa_position not in exact_overlap_positions_seen
        ):
            exact_overlap_positions_seen.add(aa_position)
            exact_overlap_count += 1

        matched_comparison_position = next(
            (
                comparison_position
                for comparison_position in comparison_position_list
                if abs(aa_position - comparison_position) <= proximity_window
            ),
            None,
        )
        if matched_comparison_position is not None:
            proximity_overlap_count += 1
            covered_comparison_positions.add(matched_comparison_position)

        precision_at_k = exact_overlap_count / rank
        recall_against_comparison_sites = (
            exact_overlap_count / comparison_site_count if comparison_site_count else 0.0
        )
        proximity_precision_at_k = proximity_overlap_count / rank
        proximity_recall = (
            len(covered_comparison_positions) / comparison_site_count
            if comparison_site_count else 0.0
        )

        sweep_rows.append(
            {
                "top_k": rank,
                "overlap_count": exact_overlap_count,
                "precision_at_k": precision_at_k,
                "recall_against_comparison_sites": recall_against_comparison_sites,
                "proximity_overlap_count": proximity_overlap_count,
                "proximity_precision_at_k": proximity_precision_at_k,
                "proximity_recall": proximity_recall,
            }
        )

    return sweep_rows


def build_known_hotspot_case_study_payload(
    predictions,
    reference_genome_sequence,
    node_context,
    selected_model,
    elapsed_day=DEFAULT_ELAPSED_DAY,
    top_k=DEFAULT_TOP_K,
    spike_start_genome_position=SPIKE_START,
    prediction_start_genome_position=SPIKE_START,
    reference_spike_amino_acids=None,
):
    validated_context = validate_delta_pre_omicron_context(node_context)
    normalized_model_name = selected_model or DEFAULT_MODEL_NAME
    comparison_site_set = load_precomputed_variant_spike_site_set(
        validated_context["node_id"]
    )
    comparison_positions = list(comparison_site_set["positions"])

    ranked_rows, scoring_metadata = build_ranked_spike_site_rows(
        selected_model=normalized_model_name,
        node_context=validated_context,
        predictions=predictions,
        reference_genome_sequence=reference_genome_sequence,
        spike_start_genome_position=spike_start_genome_position,
        prediction_start_genome_position=prediction_start_genome_position,
        reference_spike_amino_acids=reference_spike_amino_acids,
    )
    if not ranked_rows:
        raise CaseStudyValidationError(
            "No Spike amino-acid score rows could be built from the provided prediction output."
        )

    requested_top_k = int(top_k or DEFAULT_TOP_K)
    if requested_top_k < 1:
        raise CaseStudyValidationError("topK must be at least 1.")

    applied_top_k = min(requested_top_k, len(ranked_rows))
    top_k_positions = [row["aa_position"] for row in ranked_rows[:applied_top_k]]
    metrics = compute_overlap_metrics(top_k_positions, comparison_positions)
    top_k_sweep = build_top_k_metric_sweep(
        ranked_rows,
        comparison_positions,
        proximity_window=metrics["proximity_window"],
    )
    found_support_summary = summarize_found_mutation_support(
        validated_context,
        metrics["overlap_positions"],
    )
    overlap_set = set(metrics["overlap_positions"])
    proximity_overlap_set = set(metrics["proximity_overlap_positions"])
    top_k_set = set(top_k_positions)
    comparison_position_set = set(comparison_positions)
    available_context_nodes = build_delta_context_option_metadata(validated_context)

    # --- Diagnostic logging ---
    top_k_with_scores = [
        (row["aa_position"], row.get("raw_site_score", row["site_score"]))
        for row in ranked_rows[:applied_top_k]
    ]
    print(
        f"[case-study-overlap] model={normalized_model_name} "
        f"top_k={applied_top_k} total_sites={len(ranked_rows)}"
    )
    print(f"  top-K positions (aa, raw_score): {top_k_with_scores}")
    print(f"  comparison positions: {sorted(comparison_positions)}")
    print(
        f"  exact overlap={metrics['overlap_positions']} "
        f"precision={metrics['precision_at_k']:.3f} "
        f"recall={metrics['recall_against_omicron_sites']:.3f}"
    )
    print(
        f"  proximity overlap (±{metrics['proximity_window']}aa): "
        f"{metrics['proximity_overlap_positions']} "
        f"count={metrics['proximity_overlap_count']} "
        f"precision={metrics['proximity_precision_at_k']:.3f} "
        f"recall={metrics['proximity_recall']:.3f}"
    )
    for hit in metrics.get("proximity_hits", []):
        print(
            f"    hotspot {hit['hotspot_position']} ↔ "
            f"comparison {hit['omicron_position']} "
            f"(distance={hit['distance']})"
        )
    if ranked_rows:
        scores = [r.get("raw_site_score", r["site_score"]) for r in ranked_rows]
        print(
            f"  score range: min={min(scores):.8f} max={max(scores):.8f} "
            f"top1/bottom1 ratio={max(scores)/max(min(scores), 1e-15):.1f}"
        )

    ranked_rows_with_flags = [
        {
            **row,
            "is_top_k": row["aa_position"] in top_k_set,
            "is_comparison_site": row["aa_position"] in comparison_position_set,
            "is_omicron_site": row["aa_position"] in comparison_position_set,
            "is_overlap": row["aa_position"] in overlap_set,
            "is_proximity_overlap": row["aa_position"] in proximity_overlap_set,
        }
        for row in ranked_rows
    ]

    row_by_position = {
        row["aa_position"]: row for row in ranked_rows_with_flags
    }
    score_series = [
        {
            **row_by_position[row["aa_position"]],
        }
        for row in sorted(ranked_rows_with_flags, key=lambda item: item["aa_position"])
    ]

    return {
        "metadata": {
            "analysis_id": ANALYSIS_ID,
            "analysis_label": ANALYSIS_LABEL,
            "analysis_mode": ANALYSIS_MODE,
            "analysis_disclaimer": ANALYSIS_DISCLAIMER,
            "coordinate_system": {
                "system": "Spike amino-acid positions",
                "indexing": "1-based",
                "reference_sequence": "SARS-CoV-2 reference Spike coding sequence",
                "spike_gene": SPIKE_GENE,
                "spike_genome_window": {
                    "start": spike_start_genome_position,
                    "end": SPIKE_END,
                },
            },
            "scoring_context": {
                **validated_context,
                "selected_model": normalized_model_name,
                "elapsed_day": int(elapsed_day or 0),
                **scoring_metadata,
            },
            "comparison_set": {
                "variant_label": comparison_site_set.get("variant_label"),
                "set_name": comparison_site_set.get("set_name"),
                "site_count": len(comparison_positions),
                "site_positions": comparison_positions,
                "nucleotide_mutation_count": comparison_site_set.get(
                    "nucleotide_mutation_count",
                    0,
                ),
                "support_summary": comparison_site_set.get("support_summary"),
            },
            "available_delta_context_nodes": available_context_nodes,
            "available_context_nodes": available_context_nodes,
            "requested_top_k": requested_top_k,
            "applied_top_k": applied_top_k,
            "ranked_row_count": len(ranked_rows_with_flags),
        },
        "score_series": score_series,
        "top_k_positions": top_k_positions,
        "comparison_positions": comparison_positions,
        "omicron_positions": comparison_positions,
        "overlap_positions": metrics["overlap_positions"],
        "proximity_overlap_positions": metrics["proximity_overlap_positions"],
        "top_k_sweep": top_k_sweep,
        "metrics": {
            "overlap_count": metrics["overlap_count"],
            "precision_at_k": metrics["precision_at_k"],
            "recall_against_omicron_sites": metrics["recall_against_omicron_sites"],
            "recall_against_comparison_sites": metrics["recall_against_omicron_sites"],
            "top_k_count": metrics["top_k_count"],
            "comparison_site_count": metrics["omicron_site_count"],
            "omicron_site_count": metrics["omicron_site_count"],
            "proximity_window": metrics["proximity_window"],
            "proximity_overlap_count": metrics["proximity_overlap_count"],
            "proximity_precision_at_k": metrics["proximity_precision_at_k"],
            "proximity_recall": metrics["proximity_recall"],
            "proximity_hits": metrics["proximity_hits"],
            "found_mutation_count": found_support_summary["mutation_count"],
            "found_total_proportion": found_support_summary["total_proportion"],
            "found_mean_proportion": found_support_summary["mean_proportion"],
            "found_median_proportion": found_support_summary["median_proportion"],
        },
        "ranked_rows": ranked_rows_with_flags,
    }
