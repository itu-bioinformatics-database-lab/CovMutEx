# PRIEST Integration Rundown

## Purpose

This integration lets CovMutEx take a selected tree node, reconstruct its observed mutations, convert Spike nucleotide changes into Spike amino-acid substitutions, and annotate those Spike sites with PRIEST-derived site-level mutation propensity scores.

The bridge is intentionally:

```text
selected node
-> reconstructed mutations
-> Spike amino-acid positions
-> node date mapped to PRIEST period
-> PRIEST site score lookup
```

It is not based on raw variant ID matching.

## Files Involved

- `genome_extractor/genome/views.py`
  Django API entry point. Builds the reconstructed genome, runs the prediction model, and now injects PRIEST annotation into the response.
- `genome_extractor/genome/priest_annotations.py`
  Core integration logic for node parsing, Spike mutation collapsing, PRIEST period mapping, score lookup, and summary generation.
- `genome_extractor/genome/build_priest_lookup.py`
  Optional precompute utility that builds a fast local CSV lookup from PRIEST raw quarter data.
- `covid19-genome-visualizer/src/App.js`
  Stores the returned PRIEST annotation in React state and passes it into the UI.
- `covid19-genome-visualizer/src/components/SpikeMutationPanel.js`
  Renders the detailed mutation table and summary metrics.

## End-to-End Flow

```text
Frontend node selection
-> POST /api/predict/
-> CovMutEx parses node and reconstructs genome
-> PRIEST integration extracts Spike AA substitutions
-> node date is mapped to quarter
-> PRIEST score is looked up by (period, aa_position), then by aa_position
-> backend returns annotated Spike rows + summary
-> frontend renders a Spike Mutation Annotation panel
```

## Backend Request Lifecycle

### 1. Request enters `handle_prediction`

The frontend calls `POST /api/predict/` with:

- `nodeId`
- `elapsedDay`
- `selectedModel`
- `selectedProteinRegion`

In `genome_extractor/genome/views.py`, `handle_prediction()` does three things in parallel conceptually:

1. reconstructs the selected genome from the node mutations
2. runs the CovMutEx prediction model
3. builds the PRIEST annotation payload

The PRIEST annotation does not change the existing reconstruction pipeline. It uses the same mutation list and reconstructed genome that CovMutEx already produces.

### 2. Node parsing

`build_node_priest_annotation()` starts by calling `parse_selected_node(node_id)`.

Example node:

```text
EGY/CCHE57357_Wave_3_A029/2021|MZ380261.1|2021-05-11
```

Parsed output:

- `selected_node`: full original string
- `selected_node_name`: `EGY/CCHE57357_Wave_3_A029/2021`
- `selected_node_accession`: `MZ380261.1`
- `node_date`: `2021-05-11`

The date is later converted into the PRIEST time bucket.

### 3. Genome reconstruction stays unchanged

The existing CovMutEx path is preserved:

- `get_sample_depth(...)`
- `parse_mutations(nodeId)`
- `construct_variant_genome(genome_sequence, mutations)`
- `predict_mutations(...)`

The integration reuses:

- `mutations`
- `genome_sequence` as the reference genome
- `variant_genome_sequence` as the reconstructed selected genome

This is important because the amino-acid annotation is based on the final reconstructed codon state, not on a fragile one-mutation-at-a-time guess.

### 4. Mutation serialization

`parse_mutations(nodeId)` already returns mutation tuples. The PRIEST helper normalizes them into explicit dictionaries:

```python
{
    "genome_position": 22917,
    "ref_nt": "T",
    "alt_nt": "G",
    "recorded_aa_position": 452,
    "recorded_aa_change": "L452R"
}
```

The recorded amino-acid metadata is preserved, but the integration does not trust it as the final source of truth. It recomputes Spike amino-acid changes from the reconstructed genome.

### 5. Gene and region annotation

Each mutation is assigned to a region using `PROTEIN_REGIONS`:

- `ORF1ab`
- `S`
- `ORF3a`
- `E`
- `M`
- `ORF6`
- `ORF7a`
- `ORF7b`
- `ORF8`
- `N`
- `ORF10`

Anything outside those windows is labeled `Non-coding`.

Spike is defined as:

- start: `21563`
- end: `25384`

Only mutations whose genomic positions fall inside this interval are kept for PRIEST annotation.

### 6. Spike nucleotide mutations are collapsed into Spike amino-acid substitutions

This is the main robustness step.

For every Spike nucleotide mutation:

1. compute Spike offset:

```python
spike_nt_offset = genome_position - 21563
```

2. compute codon index and amino-acid position:

```python
codon_index_zero_based = spike_nt_offset // 3
aa_position = codon_index_zero_based + 1
```

3. compute the codon start inside the genome
4. extract the reference codon from `genome_sequence`
5. extract the mutated codon from `variant_genome_sequence`
6. translate both codons using `codon_aa_mapping.json`
7. build the amino-acid substitution

Example result:

```python
{
    "gene": "S",
    "aa_position": 452,
    "ref_aa": "L",
    "alt_aa": "R",
    "aa_mutation": "L452R",
    "ref_codon": "CTT",
    "alt_codon": "CGT",
    "is_synonymous": False,
    "supporting_nt_mutations": [
        {
            "genome_position": 22917,
            "ref_nt": "T",
            "alt_nt": "G",
            "spike_nt_offset": 1354,
            "codon_nt_index": 1
        }
    ]
}
```

### Why this codon-based collapse matters

Multiple nucleotide changes can land in the same codon. The code groups rows by `aa_position` and always uses the reconstructed genome as the final codon state. That means:

- combined codon changes are handled correctly
- the final amino-acid call reflects the selected node as reconstructed
- matching against PRIEST happens on the biologically cleaner unit: Spike amino-acid site

### 7. Synonymous vs nonsynonymous handling

After codon translation:

- synonymous Spike changes are tracked separately
- nonsynonymous Spike changes are returned in `spike_mutations`

So the main table only shows nonsynonymous Spike amino-acid substitutions, while the summary still keeps count of synonymous Spike events.

## Temporal Alignment

### 8. Node date to PRIEST period mapping

`node_date_to_priest_period(node_date)` maps ISO dates into quarter buckets:

- Jan-Mar -> `Q1-YYYY`
- Apr-Jun -> `Q2-YYYY`
- Jul-Sep -> `Q3-YYYY`
- Oct-Dec -> `Q4-YYYY`

Example:

```text
2021-05-11 -> Q2-2021
```

If the node has no valid ISO date, the PRIEST period becomes `None` and the lookup falls back to global site scores only.

## PRIEST Score Lookup

### 9. Lookup strategy

The score lookup is intentionally two-tiered:

1. preferred: local precomputed CSV lookup
2. fallback: derive scores from PRIEST raw quarter data for the requested positions

The lookup key is:

```python
(period, aa_position)
```

with fallback:

```python
aa_position
```

### 9a. Preferred path: `priest_site_scores.csv`

If `genome_extractor/genome/priest_site_scores.csv` exists, it is loaded once and cached.

Expected columns:

```text
period,aa_position,priest_score
```

Behavior:

- rows with a non-empty `period` become period-specific scores
- rows with an empty `period` become global fallback scores

This path is reported as:

```text
priest_score_method = "csv_lookup"
```

### 9b. Fallback path: raw PRIEST quarter scan

If the CSV is missing, `priest_annotations.py` scans the PRIEST raw data only for the amino-acid positions needed by the selected node.

It reads quarter files from:

- `PRIEST/src/PRIEST_data/Raw Data/quarter_raw/...`
- `PRIEST/src/PRIEST_data/Raw Data/2023 data/...`

For each requested position and period:

1. count how many valid PRIEST sequences have an observed amino acid at that site
2. count how many of those amino acids differ from the CovMutEx reference Spike amino acid
3. compute:

```python
priest_score = mutated_count / observed_count
```

This score is site-level mutation prevalence or propensity at that Spike position for that period. It is not a mutation-specific score for a particular substitution like `L452R` versus `L452Q`.

This path is reported as:

```text
priest_score_method = "priest_raw_period_mutation_prevalence"
```

### 9c. Preferred lookup and fallback logic

For each Spike amino-acid row:

1. try period-specific lookup:

```python
period_scores[period][aa_position]
```

2. if missing, try global fallback:

```python
global_scores[aa_position]
```

3. otherwise return `None`

The returned row gets:

- `priest_period`
- `priest_score`
- `priest_score_source`

Where `priest_score_source` is:

- `period`
- `global`
- `None`

Separately, the payload also reports `priest_period_available`, which tells you whether the selected node's mapped quarter actually exists in the loaded PRIEST period set. If it is `false`, the code can still return global fallback scores.

## Summary Generation

### 10. Summary metrics

The backend also returns a compact summary:

- `num_spike_nucleotide_mutations`
- `num_spike_mutations`
- `num_synonymous_spike_mutations`
- `num_priest_annotated`
- `num_priest_supported_ge_0_6`
- `priest_support_threshold`
- `mean_priest_score`
- `top_supported_mutations`
- `top_unsupported_mutations`

Current support threshold default:

```text
0.6
```

Interpretation:

- `num_spike_mutations` counts nonsynonymous Spike amino-acid substitutions
- `num_priest_annotated` counts how many of those got any PRIEST score
- `num_priest_supported_ge_0_6` counts scored rows at or above the threshold
- `mean_priest_score` is averaged over scored nonsynonymous Spike rows only
- `top_supported_mutations` are the highest-scoring rows
- `top_unsupported_mutations` are rows with no score, or if all rows have scores, the lowest-scoring rows below threshold

## API Response Shape

### 11. Returned payload

The Django response now includes the original fields plus the PRIEST block:

```json
{
  "selected_node": "EGY/CCHE57357_Wave_3_A029/2021|MZ380261.1|2021-05-11",
  "selected_node_name": "EGY/CCHE57357_Wave_3_A029/2021",
  "selected_node_accession": "MZ380261.1",
  "node_date": "2021-05-11",
  "priest_period": "Q2-2021",
  "priest_period_available": true,
  "priest_score_method": "csv_lookup",
  "priest_available_periods": ["Q1-2020", "Q2-2020", "Q3-2020"],
  "spike_mutations": [
    {
      "gene": "S",
      "aa_position": 452,
      "ref_aa": "L",
      "alt_aa": "R",
      "aa_mutation": "L452R",
      "genome_position": 22917,
      "genome_positions": [22917],
      "priest_period": "Q2-2021",
      "priest_score": 0.696,
      "priest_score_source": "period"
    }
  ],
  "priest_summary": {
    "num_spike_nucleotide_mutations": 2,
    "num_spike_mutations": 2,
    "num_synonymous_spike_mutations": 0,
    "num_priest_annotated": 1,
    "num_priest_supported_ge_0_6": 1,
    "priest_support_threshold": 0.6,
    "mean_priest_score": 0.696,
    "top_supported_mutations": [
      {
        "aa_mutation": "L452R",
        "aa_position": 452,
        "priest_score": 0.696,
        "priest_score_source": "period"
      }
    ],
    "top_unsupported_mutations": [
      {
        "aa_mutation": "T478K",
        "aa_position": 478,
        "priest_score": null,
        "priest_score_source": null
      }
    ]
  }
}
```

## Frontend Behavior

### 12. React state handling

In `covid19-genome-visualizer/src/App.js`:

1. the first `/api/predict/` call stores the PRIEST payload in `priestAnnotation`
2. if a protein region is selected, a second fetch retrieves the full genome for chart rendering
3. the PRIEST annotation is not recomputed on the client and is not tied to the second fetch

That means the annotation shown in the UI always comes from the selected node response itself.

### 13. `SpikeMutationPanel`

The panel renders:

- selected node metadata
- node date
- mapped PRIEST period
- summary cards
- detailed mutation table
- top supported mutations
- top unsupported mutations
- note explaining whether scores came from CSV lookup or raw-data-derived fallback

The detailed table shows:

- mutation
- Spike position
- PRIEST score
- source

Source labels:

- `period` -> `Matched period`
- `global` -> `Global fallback`

## Lookup Precompute Utility

### 14. `build_priest_lookup.py`

This utility exists to speed up runtime annotation.

It:

1. loads the CovMutEx reference Spike amino-acid sequence
2. iterates through PRIEST period files
3. compares each PRIEST Spike amino-acid sequence against the reference
4. computes site-level mutation prevalence for every Spike amino-acid position
5. writes a local CSV to:

```text
genome_extractor/genome/priest_site_scores.csv
```

The output contains:

- period-specific rows
- global fallback rows with empty `period`

Run it with:

```bash
python3 genome_extractor/genome/build_priest_lookup.py
```

## Current Runtime Characteristics

### 15. Fast path vs slow path

Fast path:

- `priest_site_scores.csv` exists
- lookup is loaded once and filtered in memory

Slow path:

- CSV is missing
- backend scans PRIEST raw quarter files for the requested positions

The fallback is scientifically consistent, but slower because the PRIEST raw corpus is large.

In practice, if `genome_extractor/genome/priest_site_scores.csv` has not been generated yet, the integration will still work, but requests may be noticeably slower until that precompute step is completed.

### 16. Current caveats

- PRIEST scoring is site-level, not substitution-specific. `L452R` and `L452Q` both map to site `452`.
- The annotation is Spike-only by design because PRIEST is a Spike model.
- If the node date is missing or malformed, temporal matching is skipped and only global site fallback is possible.
- Synonymous Spike changes are tracked but not shown in the main mutation table.
- Without the precomputed CSV, first-time requests can be noticeably slower.

## Paper-Safe Claim

Use wording like:

> For a selected variant, CovMutEx reconstructs observed mutations, identifies Spike amino-acid substitutions, and cross-references them with PRIEST-derived site-level mutation propensity scores, optionally within the corresponding temporal window.

Avoid claiming:

- per-variant PRIEST inference
- substitution-specific PRIEST probabilities
- exact lineage-level matching inside PRIEST

## Practical Mental Model

If you want to explain the integration in one sentence:

```text
CovMutEx reconstructs what changed in the selected node, converts those Spike changes into amino-acid sites, maps the node date into the matching PRIEST quarter, and asks whether PRIEST considers those Spike sites historically prone to mutation.
```
