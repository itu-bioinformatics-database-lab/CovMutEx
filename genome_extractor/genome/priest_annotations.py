import csv
import json
import re
from datetime import date
from functools import lru_cache
from pathlib import Path


SPIKE_START = 21563
SPIKE_END = 25384
SPIKE_GENE = "S"
PRIEST_SUPPORT_THRESHOLD = 0.6
UNKNOWN_AMINO_ACIDS = {"", "X", "_", "-"}

PROTEIN_REGIONS = {
    "ORF1ab": (266, 21555),
    "S": (21563, 25384),
    "ORF3a": (25393, 26220),
    "E": (26245, 26472),
    "M": (26523, 27191),
    "ORF6": (27202, 27387),
    "ORF7a": (27394, 27759),
    "ORF7b": (27756, 27887),
    "ORF8": (27894, 28259),
    "N": (28274, 29533),
    "ORF10": (29558, 29674),
}

BASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = BASE_DIR.parent.parent
RAW_DATA_DIR = REPO_ROOT / "PRIEST" / "src" / "PRIEST_data" / "Raw Data"
PRIEST_LOOKUP_OVERRIDE = BASE_DIR / "priest_site_scores.csv"
CODON_MAPPING_PATH = BASE_DIR / "codon_aa_mapping.json"
PRIEST_LOG_PREFIX = "[PRIEST]"
PRIEST_LOG_PREVIEW_LIMIT = 20


def _priest_log(event, **payload):
    """Print a structured PRIEST debug log line.

    Args:
        event: Short label describing the logging event.
        **payload: JSON-serializable metadata to include in the log.

    Returns:
        None.
    """
    try:
        message = json.dumps(payload, sort_keys=True, default=str)
    except TypeError:
        message = str(payload)
    print(f"{PRIEST_LOG_PREFIX} {event}: {message}")


def _mutation_log_snapshot(row):
    """Build a compact mutation snapshot for debug logging.

    Args:
        row: One reconstructed Spike amino-acid mutation row.

    Returns:
        A dictionary with the key fields used in PRIEST debug output.
    """
    return {
        "aa_mutation": row["aa_mutation"],
        "aa_position": row["aa_position"],
        "ref_aa": row["ref_aa"],
        "alt_aa": row["alt_aa"],
        "ref_codon": row["ref_codon"],
        "alt_codon": row["alt_codon"],
        "is_synonymous": row["is_synonymous"],
        "genome_positions": row.get("genome_positions", []),
    }


def _score_preview(scores, positions, limit=PRIEST_LOG_PREVIEW_LIMIT):
    """Create a small score preview for requested Spike positions.

    Args:
        scores: Mapping of amino-acid positions to PRIEST scores.
        positions: Ordered Spike amino-acid positions requested by the node.
        limit: Maximum number of positions to include in the preview.

    Returns:
        A list of preview dictionaries with position and rounded score.
    """
    preview = []
    for position in positions[:limit]:
        preview.append(
            {
                "aa_position": position,
                "score": _round_score(scores.get(position)),
            }
        )
    return preview


def parse_selected_node(node_id):
    """Extract node metadata from a CovMutEx node identifier string.

    Args:
        node_id: Raw selected node identifier from the frontend.

    Returns:
        A dictionary containing the original node id, node name, accession,
        and parsed collection date when present.
    """
    tokens = [token.strip() for token in (node_id or "").split("|") if token.strip()]
    node_date = None

    for token in reversed(tokens):
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", token):
            node_date = token
            break

    core_tokens = list(tokens)
    if node_date and core_tokens and core_tokens[-1] == node_date:
        core_tokens = core_tokens[:-1]

    accession = core_tokens[-1] if len(core_tokens) >= 2 else None
    node_name = "|".join(core_tokens[:-1]) if len(core_tokens) >= 2 else (core_tokens[0] if core_tokens else node_id)

    return {
        "selected_node": node_id,
        "selected_node_name": node_name,
        "selected_node_accession": accession,
        "node_date": node_date,
    }


def build_node_priest_annotation(node_id, mutation_tuples, reference_genome_sequence, variant_genome_sequence, threshold=PRIEST_SUPPORT_THRESHOLD):
    """Build the full PRIEST annotation payload for one selected node.

    Args:
        node_id: Selected node identifier.
        mutation_tuples: Reconstructed mutation tuples applied to the reference.
        reference_genome_sequence: Reference SARS-CoV-2 genome sequence.
        variant_genome_sequence: Reconstructed genome sequence for the node.
        threshold: PRIEST support cutoff used for summary statistics.

    Returns:
        A structured annotation dictionary for backend responses and frontend
        rendering, including node metadata, Spike-site rows, and summary stats.
    """
    node_info = parse_selected_node(node_id)
    all_mutations = _serialize_mutations(mutation_tuples)
    annotated_mutations = _annotate_mutation_regions(all_mutations)
    spike_nucleotide_mutations = [mutation for mutation in annotated_mutations if mutation["gene"] == SPIKE_GENE]
    non_spike_mutations = [mutation for mutation in annotated_mutations if mutation["gene"] != SPIKE_GENE]
    aa_rows = _collapse_spike_mutations_to_amino_acids(
        spike_nucleotide_mutations=spike_nucleotide_mutations,
        reference_genome_sequence=reference_genome_sequence,
        variant_genome_sequence=variant_genome_sequence,
    )

    requested_positions = tuple(sorted({row["aa_position"] for row in aa_rows}))
    period = node_date_to_priest_period(node_info["node_date"])
    _priest_log(
        "annotation_start",
        selected_node=node_info["selected_node"],
        selected_node_name=node_info["selected_node_name"],
        accession=node_info["selected_node_accession"],
        node_date=node_info["node_date"],
        priest_period=period,
        total_mutations=len(all_mutations),
        spike_nt_mutations=len(spike_nucleotide_mutations),
        spike_aa_rows=len(aa_rows),
        requested_positions=list(requested_positions[:PRIEST_LOG_PREVIEW_LIMIT]),
    )
    _priest_log(
        "spike_aa_rows",
        rows=[_mutation_log_snapshot(row) for row in aa_rows[:PRIEST_LOG_PREVIEW_LIMIT]],
        truncated=len(aa_rows) > PRIEST_LOG_PREVIEW_LIMIT,
    )

    lookup = load_priest_score_lookup(requested_positions)
    period_scores = lookup["period_scores"]
    global_scores = lookup["global_scores"]
    period_available = period in lookup["available_periods"] if period else False
    matched_period_scores = period_scores.get(period, {}) if period else {}

    _priest_log(
        "lookup_context",
        method=lookup["method"],
        override_file_exists=PRIEST_LOOKUP_OVERRIDE.exists(),
        requested_period=period,
        period_available=period_available,
        available_period_count=len(lookup["available_periods"]),
        available_periods_preview=list(lookup["available_periods"][:PRIEST_LOG_PREVIEW_LIMIT]),
        requested_period_score_preview=_score_preview(matched_period_scores, requested_positions),
        global_score_preview=_score_preview(global_scores, requested_positions),
    )
    if period and not period_available:
        _priest_log(
            "period_missing_using_global_fallback",
            requested_period=period,
            requested_positions=list(requested_positions[:PRIEST_LOG_PREVIEW_LIMIT]),
        )

    spike_mutations = []
    synonymous_spike_mutations = []
    for row in aa_rows:
        row = dict(row)
        period_score = matched_period_scores.get(row["aa_position"]) if period else None
        global_score = global_scores.get(row["aa_position"])
        score, score_source = _lookup_priest_score(
            aa_position=row["aa_position"],
            period=period,
            period_scores=period_scores,
            global_scores=global_scores,
        )
        row["priest_period"] = period
        row["priest_score"] = _round_score(score)
        row["priest_score_source"] = score_source
        _priest_log(
            "score_resolution",
            aa_mutation=row["aa_mutation"],
            aa_position=row["aa_position"],
            is_synonymous=row["is_synonymous"],
            genome_positions=row.get("genome_positions", []),
            period_score=_round_score(period_score),
            global_score=_round_score(global_score),
            assigned_score=row["priest_score"],
            assigned_source=score_source,
        )
        if row["is_synonymous"]:
            synonymous_spike_mutations.append(row)
        else:
            spike_mutations.append(row)

    summary = _build_priest_summary(
        spike_mutations=spike_mutations,
        all_mutations=all_mutations,
        annotated_mutations=annotated_mutations,
        non_spike_mutations=non_spike_mutations,
        spike_nucleotide_mutations=spike_nucleotide_mutations,
        synonymous_spike_mutations=synonymous_spike_mutations,
        threshold=threshold,
    )
    _priest_log(
        "summary",
        **summary,
    )
    if spike_mutations and all(row["priest_score"] is None for row in spike_mutations):
        _priest_log(
            "warning_all_scores_missing",
            aa_mutations=[row["aa_mutation"] for row in spike_mutations],
            requested_period=period,
        )
    scored_rows = [row for row in spike_mutations if row["priest_score"] is not None]
    if scored_rows and max(row["priest_score"] for row in scored_rows) == 0:
        _priest_log(
            "warning_all_scored_sites_zero",
            aa_mutations=[row["aa_mutation"] for row in scored_rows],
            requested_period=period,
        )

    return {
        **node_info,
        "priest_period": period,
        "priest_period_available": period_available,
        "priest_score_method": lookup["method"],
        "priest_available_periods": list(lookup["available_periods"]),
        "spike_mutations": spike_mutations,
        "synonymous_spike_mutations": synonymous_spike_mutations,
        "priest_summary": summary,
    }


def node_date_to_priest_period(node_date):
    """Map an ISO node date to the PRIEST quarter label.

    Args:
        node_date: Collection date in YYYY-MM-DD format.

    Returns:
        A quarter label such as ``Q2-2021`` or ``None`` if parsing fails.
    """
    if not node_date:
        return None

    try:
        parsed_date = date.fromisoformat(node_date)
    except ValueError:
        return None

    quarter = ((parsed_date.month - 1) // 3) + 1
    return f"Q{quarter}-{parsed_date.year}"


def _serialize_mutations(mutation_tuples):
    """Normalize reconstructed mutation tuples into dictionaries.

    Args:
        mutation_tuples: Iterable of mutation tuples from the reconstruction path.

    Returns:
        A list of mutation dictionaries with genome and recorded amino-acid fields.
    """
    serialized = []
    for position, ref_nt, alt_nt, aa_position, aa_change in mutation_tuples:
        serialized.append(
            {
                "genome_position": position,
                "ref_nt": ref_nt,
                "alt_nt": alt_nt,
                "recorded_aa_position": aa_position,
                "recorded_aa_change": aa_change,
            }
        )
    return serialized


def _annotate_mutation_regions(mutations):
    """Annotate each reconstructed mutation with its genomic region.

    Args:
        mutations: Serialized mutation dictionaries.

    Returns:
        A list of mutation dictionaries extended with a ``gene`` field.
    """
    annotated = []
    for mutation in mutations:
        annotated.append(
            {
                **mutation,
                "gene": _find_gene_for_position(mutation["genome_position"]),
            }
        )
    return annotated


def _find_gene_for_position(genome_position):
    """Assign a genomic position to a named SARS-CoV-2 region.

    Args:
        genome_position: One-based nucleotide position in the genome.

    Returns:
        The matching gene/region label, or ``Non-coding`` if none matches.
    """
    for gene, (start, end) in PROTEIN_REGIONS.items():
        if start <= genome_position <= end:
            return gene
    return "Non-coding"


def _collapse_spike_mutations_to_amino_acids(spike_nucleotide_mutations, reference_genome_sequence, variant_genome_sequence):
    """Collapse Spike nucleotide mutations into amino-acid site events.

    Args:
        spike_nucleotide_mutations: Reconstructed mutations that fall in Spike.
        reference_genome_sequence: Reference SARS-CoV-2 genome sequence.
        variant_genome_sequence: Reconstructed genome sequence for the node.

    Returns:
        A list of Spike amino-acid rows, one per affected amino-acid position.
    """
    codon_mapping = _load_codon_mapping()
    aa_rows = {}

    for mutation in spike_nucleotide_mutations:
        genome_position = mutation["genome_position"]
        spike_nt_offset = genome_position - SPIKE_START
        codon_index_zero_based = spike_nt_offset // 3
        aa_position = codon_index_zero_based + 1
        codon_start = SPIKE_START + (codon_index_zero_based * 3)

        reference_codon = reference_genome_sequence[codon_start - 1 : codon_start + 2].upper()
        mutated_codon = variant_genome_sequence[codon_start - 1 : codon_start + 2].upper()
        if len(reference_codon) != 3 or len(mutated_codon) != 3:
            continue

        ref_aa = _translate_codon(reference_codon, codon_mapping)
        alt_aa = _translate_codon(mutated_codon, codon_mapping)
        if not ref_aa or not alt_aa:
            continue

        row = aa_rows.setdefault(
            aa_position,
            {
                "gene": SPIKE_GENE,
                "aa_position": aa_position,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "aa_mutation": f"{ref_aa}{aa_position}{alt_aa}",
                "ref_codon": reference_codon,
                "alt_codon": mutated_codon,
                "is_synonymous": ref_aa == alt_aa,
                "supporting_nt_mutations": [],
            },
        )

        row["supporting_nt_mutations"].append(
            {
                "genome_position": genome_position,
                "ref_nt": mutation["ref_nt"],
                "alt_nt": mutation["alt_nt"],
                "spike_nt_offset": spike_nt_offset,
                "codon_nt_index": spike_nt_offset % 3,
            }
        )

        # The reconstructed genome is the single source of truth for the final codon state.
        row["ref_aa"] = ref_aa
        row["alt_aa"] = alt_aa
        row["aa_mutation"] = f"{ref_aa}{aa_position}{alt_aa}"
        row["ref_codon"] = reference_codon
        row["alt_codon"] = mutated_codon
        row["is_synonymous"] = ref_aa == alt_aa

    collapsed_rows = []
    for aa_position in sorted(aa_rows):
        row = dict(aa_rows[aa_position])
        supporting_positions = sorted(
            mutation["genome_position"] for mutation in row["supporting_nt_mutations"]
        )
        row["genome_positions"] = supporting_positions
        row["genome_position"] = supporting_positions[0] if supporting_positions else None
        collapsed_rows.append(row)

    return collapsed_rows


def _translate_codon(codon, codon_mapping):
    """Translate a codon string into a one-letter amino acid.

    Args:
        codon: Three-nucleotide codon string.
        codon_mapping: Dictionary mapping codons to amino acids.

    Returns:
        A one-letter amino-acid code, or ``None`` for invalid codons.
    """
    if len(codon) != 3:
        return None
    return codon_mapping.get(codon.upper())


@lru_cache(maxsize=1)
def _load_codon_mapping():
    """Load the codon-to-amino-acid translation table from disk.

    Args:
        None.

    Returns:
        A dictionary mapping uppercase codons to amino-acid symbols.
    """
    with CODON_MAPPING_PATH.open() as handle:
        return json.load(handle)


@lru_cache(maxsize=1)
def _load_reference_spike_amino_acids():
    """Translate the reference Spike coding sequence into amino acids.

    Args:
        None.

    Returns:
        The reference Spike amino-acid sequence without the terminal stop codon.
    """
    genome_path = BASE_DIR / "genome.txt"
    with genome_path.open() as handle:
        next(handle)
        genome_sequence = "".join(line.strip() for line in handle)

    spike_nt = genome_sequence[SPIKE_START - 1 : SPIKE_END]
    codon_mapping = _load_codon_mapping()
    amino_acids = []
    for index in range(0, len(spike_nt) - 2, 3):
        amino_acids.append(codon_mapping.get(spike_nt[index : index + 3].upper(), "X"))

    reference = "".join(amino_acids)
    return reference[:-1] if reference.endswith("*") else reference


def load_priest_score_lookup(aa_positions):
    """Load PRIEST scores for the requested Spike amino-acid positions.

    Args:
        aa_positions: Iterable of Spike amino-acid positions to score.

    Returns:
        A lookup dictionary containing method, period scores, global scores,
        and available periods.
    """
    target_positions = tuple(sorted({int(position) for position in aa_positions if position}))
    override = _load_override_scores()
    if override is not None:
        return _filter_lookup_to_positions(override, target_positions)
    return _build_raw_period_scores_for_positions(target_positions)


def _load_override_scores():
    """Load the precomputed local PRIEST lookup CSV when available.

    Args:
        None.

    Returns:
        A lookup dictionary with period/global scores, or ``None`` if the CSV
        does not exist.
    """
    if not PRIEST_LOOKUP_OVERRIDE.exists():
        return None

    period_scores = {}
    global_scores = {}

    with PRIEST_LOOKUP_OVERRIDE.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            aa_position = int(row["aa_position"])
            score_key = "priest_score" if "priest_score" in row else "score"
            score = float(row[score_key])
            period = (row.get("period") or "").strip()

            if period:
                period_scores.setdefault(period, {})[aa_position] = score
            else:
                global_scores[aa_position] = score

    available_periods = tuple(sorted(period_scores, key=_period_sort_key))
    return {
        "method": "csv_lookup",
        "period_scores": period_scores,
        "global_scores": global_scores,
        "available_periods": available_periods,
    }


@lru_cache(maxsize=64)
def _build_raw_period_scores_for_positions(target_positions):
    """Compute PRIEST site scores directly from the raw period sequence files.

    Args:
        target_positions: Tuple of Spike amino-acid positions to evaluate.

    Returns:
        A lookup dictionary with period-specific and global site prevalence scores.
    """
    reference_spike = _load_reference_spike_amino_acids()
    positions = tuple(
        position for position in target_positions if 1 <= position <= len(reference_spike)
    )
    period_scores = {}
    available_periods = list(_available_priest_periods())

    if not positions:
        return {
            "method": "priest_raw_period_mutation_prevalence",
            "period_scores": period_scores,
            "global_scores": {},
            "available_periods": tuple(available_periods),
        }

    global_mutated = {position: 0 for position in positions}
    global_observed = {position: 0 for position in positions}

    for period, files in _iter_priest_period_files():
        mutated = {position: 0 for position in positions}
        observed = {position: 0 for position in positions}
        sequence_count = 0

        for file_path in files:
            with file_path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    sequence = (row.get("Sequence") or "").strip().upper()
                    if not sequence:
                        continue

                    sequence_count += 1
                    for position in positions:
                        if position > len(sequence):
                            continue

                        alt_aa = sequence[position - 1]
                        if alt_aa in UNKNOWN_AMINO_ACIDS:
                            continue

                        observed[position] += 1
                        global_observed[position] += 1

                        if alt_aa != reference_spike[position - 1]:
                            mutated[position] += 1
                            global_mutated[position] += 1

        if sequence_count:
            period_scores[period] = _scores_from_counts(mutated, observed)

    global_scores = _scores_from_counts(global_mutated, global_observed)
    return {
        "method": "priest_raw_period_mutation_prevalence",
        "period_scores": period_scores,
        "global_scores": global_scores,
        "available_periods": tuple(available_periods),
    }


@lru_cache(maxsize=1)
def _available_priest_periods():
    """List all PRIEST periods available in the raw data directory.

    Args:
        None.

    Returns:
        A tuple of quarter labels sorted chronologically.
    """
    periods = [period for period, _ in _iter_priest_period_files()]
    return tuple(sorted(set(periods), key=_period_sort_key))


def _iter_priest_period_files():
    """Yield PRIEST raw-data files grouped by quarter.

    Args:
        None.

    Returns:
        An iterator of ``(period, [csv_paths])`` pairs for each available quarter.
    """
    quarter_raw_dir = RAW_DATA_DIR / "quarter_raw"
    if quarter_raw_dir.exists():
        for directory in sorted(quarter_raw_dir.iterdir()):
            if not directory.is_dir():
                continue
            match = re.fullmatch(r"year_(\d{4})_(\d)", directory.name)
            if not match:
                continue

            year = int(match.group(1))
            quarter = int(match.group(2)) + 1
            files = sorted(path for path in directory.glob("*.csv") if path.is_file())
            if files:
                yield f"Q{quarter}-{year}", files

    data_2023_dir = RAW_DATA_DIR / "2023 data"
    if data_2023_dir.exists():
        for file_path in sorted(data_2023_dir.glob("quarter_*_data.csv")):
            match = re.fullmatch(r"quarter_(\d+)_data\.csv", file_path.name)
            if match:
                yield f"Q{int(match.group(1))}-2023", [file_path]


def _scores_from_counts(mutated, observed):
    """Convert mutation and observation counts into site-level prevalence scores.

    Args:
        mutated: Mapping of positions to mutated-sequence counts.
        observed: Mapping of positions to observed-sequence counts.

    Returns:
        A dictionary mapping positions to mutation prevalence scores.
    """
    scores = {}
    for position in observed:
        if observed[position]:
            scores[position] = mutated[position] / observed[position]
    return scores


def _filter_lookup_to_positions(lookup, target_positions):
    """Restrict a PRIEST lookup object to the requested positions only.

    Args:
        lookup: Lookup dictionary with period and global PRIEST scores.
        target_positions: Tuple of Spike amino-acid positions to keep.

    Returns:
        A lookup dictionary containing only the requested positions.
    """
    if not target_positions:
        return {
            **lookup,
            "period_scores": {},
            "global_scores": {},
        }

    period_scores = {}
    for period, scores in lookup["period_scores"].items():
        filtered = {
            position: score
            for position, score in scores.items()
            if position in target_positions
        }
        if filtered:
            period_scores[period] = filtered

    global_scores = {
        position: score
        for position, score in lookup["global_scores"].items()
        if position in target_positions
    }

    return {
        **lookup,
        "period_scores": period_scores,
        "global_scores": global_scores,
    }


def _lookup_priest_score(aa_position, period, period_scores, global_scores):
    """Resolve the best PRIEST score for one Spike amino-acid position.

    Args:
        aa_position: Spike amino-acid position to score.
        period: Requested PRIEST quarter label.
        period_scores: Mapping of quarters to site-score dictionaries.
        global_scores: Mapping of site positions to global fallback scores.

    Returns:
        A ``(score, source)`` tuple using period-specific scores first and
        global fallback scores second.
    """
    if period and aa_position in period_scores.get(period, {}):
        return period_scores[period][aa_position], "period"
    if aa_position in global_scores:
        return global_scores[aa_position], "global"
    return None, None


def _build_mutation_region_counts(annotated_mutations):
    """Summarize reconstructed mutations by genomic region.

    Args:
        annotated_mutations: Reconstructed mutations that already include ``gene``.

    Returns:
        An ordered list of region/count dictionaries for frontend display.
    """
    counts = {}
    for mutation in annotated_mutations:
        gene = mutation["gene"]
        counts[gene] = counts.get(gene, 0) + 1

    ordered_genes = list(PROTEIN_REGIONS) + sorted(
        gene for gene in counts if gene not in PROTEIN_REGIONS
    )
    return [
        {"gene": gene, "count": counts[gene]}
        for gene in ordered_genes
        if gene in counts
    ]


def _build_priest_summary(
    spike_mutations,
    all_mutations,
    annotated_mutations,
    non_spike_mutations,
    spike_nucleotide_mutations,
    synonymous_spike_mutations,
    threshold,
):
    """Build summary statistics for the node-level PRIEST response.

    Args:
        spike_mutations: Nonsynonymous Spike amino-acid mutation rows.
        all_mutations: All reconstructed mutations for the node.
        annotated_mutations: All reconstructed mutations annotated by region.
        non_spike_mutations: Reconstructed mutations outside Spike.
        spike_nucleotide_mutations: All reconstructed Spike nucleotide events.
        synonymous_spike_mutations: Synonymous Spike amino-acid rows.
        threshold: Support threshold used for summary and ranking.

    Returns:
        A summary dictionary with counts, mean score, region burden, and top
        supported/unsupported mutations.
    """
    scored_rows = [row for row in spike_mutations if row["priest_score"] is not None]
    supported_rows = sorted(
        [row for row in scored_rows if row["priest_score"] >= threshold],
        key=lambda row: row["priest_score"],
        reverse=True,
    )
    unsupported_rows = [row for row in spike_mutations if row["priest_score"] is None]
    if not unsupported_rows:
        unsupported_rows = sorted(
            [row for row in scored_rows if row["priest_score"] < threshold],
            key=lambda row: (row["priest_score"], row["aa_position"]),
        )

    mean_score = None
    if scored_rows:
        mean_score = sum(row["priest_score"] for row in scored_rows) / len(scored_rows)

    return {
        "num_total_mutations": len(all_mutations),
        "num_non_spike_mutations": len(non_spike_mutations),
        "num_spike_nucleotide_mutations": len(spike_nucleotide_mutations),
        "num_spike_mutations": len(spike_mutations),
        "num_synonymous_spike_mutations": len(synonymous_spike_mutations),
        "mutation_region_counts": _build_mutation_region_counts(annotated_mutations),
        "num_priest_annotated": len(scored_rows),
        "num_priest_supported_ge_0_6": len(supported_rows),
        "priest_support_threshold": threshold,
        "mean_priest_score": _round_score(mean_score),
        "top_supported_mutations": [_summarize_mutation(row) for row in supported_rows[:5]],
        "top_unsupported_mutations": [_summarize_mutation(row) for row in unsupported_rows[:5]],
    }


def _summarize_mutation(row):
    """Reduce a full mutation row to the summary fields used in responses.

    Args:
        row: One annotated Spike mutation row.

    Returns:
        A compact mutation dictionary for summary sections.
    """
    return {
        "aa_mutation": row["aa_mutation"],
        "aa_position": row["aa_position"],
        "priest_score": row["priest_score"],
        "priest_score_source": row["priest_score_source"],
    }


def _period_sort_key(period):
    """Create a chronological sort key for PRIEST quarter labels.

    Args:
        period: Quarter label such as ``Q2-2021``.

    Returns:
        A tuple sortable by year then quarter, with invalid labels sorted last.
    """
    match = re.fullmatch(r"Q(\d)-(\d{4})", period or "")
    if not match:
        return (9999, 9)
    return (int(match.group(2)), int(match.group(1)))


def _round_score(score):
    """Round a PRIEST score to a stable response precision.

    Args:
        score: Floating-point PRIEST score or ``None``.

    Returns:
        The rounded score to six decimals, or ``None`` when absent.
    """
    return round(score, 6) if score is not None else None
