# CovMutEx Plugin Contract

Contract version: `1.0`

This backend accepts convention-based upload bundles. A bundle is a directory in
`uploaded_models/<bundle_name>/` with these conventions:

- Required: `model.*`
- Optional: `feature_extractor.py`
- Optional: `custom_parameters.json`
- Optional: any helper files needed by the extractor or model

Helper files are not part of the contract surface. The runtime only guarantees
that they will be saved into the same bundle directory. Plugins should resolve
them relative to `base_dir`.

## Model contract

The runtime talks to wrapped models through these methods:

- `metadata()`
- `input_schema()`
- `preprocess(inputs)`
- `predict(batch)`
- `postprocess(raw, nucleotides_per_position=1)`

The runtime expects `postprocess()` to return predictions in position-oriented
form so they can be converted into frontend visualization data.

## Feature extractor contract

Uploaded `feature_extractor.py` should expose:

- `extract_features(...)`
- `get_feature_dimension()` (recommended)
- `get_feature_description()` (recommended)
- `get_metadata()` (recommended)

Canonical `extract_features(...)` signature:

```python
def extract_features(
    genome_seq,
    mutations,
    node_ids,
    elapsed_day=None,
    protein_regions=None,
    k=30,
    **kwargs,
):
    ...
```

Runtime guarantees:

- `node_ids` is always a list of strings
- `elapsed_day` is always `int | None`
- `protein_regions` is `dict | None`
- `base_dir` is available in `kwargs`
- merged custom parameters are passed through `kwargs`

## Naming conventions

- Model file is always stored as `model.<ext>`
- Extractor file is always stored as `feature_extractor.py`
- Helper files keep the requested names after sanitization
- `custom_parameters.json` stores the merged default parameters for the bundle

## Security and path rules

- Uploaded filenames are sanitized before writing to disk
- Helper files must use an allowed extension
- All uploaded files are written into one bundle directory
- Plugins must only access bundle-local helper files or server-provided inputs
