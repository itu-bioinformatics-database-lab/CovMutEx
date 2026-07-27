# CovMutEx Plugin Contract

Contract version: `2.0`

How to plug your own organism, feature extractor, and/or model into CovMutEx
without touching its source code.

The platform exposes **four plug points**. Each one is optional — every plug
point you don't supply, CovMutEx fills in with its built-in default.

```
                                 YOUR INPUTS
                                      |
                        ┌─────────────┴─────────────┐
                        ▼                           ▼
                 feature_extractor              model_adapter
                     (.py)                          (.py)
                        │                           │
[ORGANISM]    ┌─────────┴───────┐         ┌─────────┴──────────┐
genome.fasta ─→ extract_features │ ─────→  │ preprocess         │
protein_regions │                │         │ predict            │
      .csv    └─────────────────┘         │ postprocess ──→ PredictionPayload v2.0
                                          └────────────────────┘
```

---

## Bundle layout

When you upload a model, CovMutEx creates a directory at
`uploaded_models/<your_model_name>/` and stores your files there with
**canonical names**:

```
uploaded_models/my_model/
├── model.<ext>             (required) — your model file, stored as model.*
├── pytorch_model.bin       (auto-symlink) — only for .bin uploads; lets HuggingFace from_pretrained() find the file
├── feature_extractor.py    (optional) — your feature extractor
├── model_adapter.py        (optional) — your output transformer
├── bundle_metadata.json    (auto-generated) — { "organism": "...", ... }
├── custom_parameters.json  (auto-generated if you declared params)
├── genome.fasta            (optional, only for custom organisms)
├── protein_regions.csv     (optional, only for custom organisms)
└── <any helper files>      (optional, keep their original sanitized names)
```

**Model and code files use canonical names; helper files keep their original
(sanitized) names.** Inside any of your `.py` files you can reference
companions relative to `base_dir` — they are always neighbors:

### Supported file extensions

**Model file** (Model File slot):
```
.keras  .h5  .hdf5  .pkl  .pt  .pth  .bin  .onnx  .safetensors
```

**Helper files** (Helper Files slot):
```
.py  .json  .txt  .csv  .nwk  .tsv  .fasta  .fa  .yaml  .yml  .npy
```

### HuggingFace note

For `.bin` models, `from_pretrained(base_dir)` works out of the box — the platform
handles the necessary file setup automatically.

```python
# inside feature_extractor.py or model_adapter.py
import os
base_dir = kwargs["base_dir"]              # passed by CovMutEx
lookup_table = os.path.join(base_dir, "my_table.csv")
```

The `bundle_metadata.json` records organism dispatch and optional organism
helper file paths:

```jsonc
{
  "organism": "covid",               // "covid" | "influenza" | "influenza_h1n1" | "influenza_h3n2" | "influenza_h5n1" | "custom"
  "variant": "H1N1_A_California_07_2009",  // influenza only, optional at save time
  "genome_file": "genome.fasta",     // only for organism="custom"
  "protein_regions_file": "protein_regions.csv"  // optional, custom only
}
```

---

## Plug 1 — Organism

### Built-in organisms

The platform ships three organisms out of the box. Select them from the form
dropdown — no genome/protein-region upload needed:

| `bundle_metadata.json` value | Virus | Notes |
|---|---|---|
| `"covid"` | SARS-CoV-2 (29903 bp) | Default when `organism` is absent. UShER mutation overlay applied at predict time. |
| `"influenza"` | Influenza A (generic) | **Preferred** for new uploads. Variant (H1N1 / H3N2 / H5N1) resolved at predict time via the `variant` field. |
| `"influenza_h1n1"` | Influenza A / H1N1 | Subtype alias — normalises to `"influenza"`. Kept for backward compatibility. |
| `"influenza_h3n2"` | Influenza A / H3N2 | Subtype alias — normalises to `"influenza"`. |
| `"influenza_h5n1"` | Influenza A / H5N1 | Subtype alias — normalises to `"influenza"`. |
| `"custom"` | Your own organism | Requires `genome.fasta` in the bundle. See below. |

For `"influenza"` bundles, pass `variant="<strain_name>"` in every predict
request so the runtime can load the correct HA FASTA. If no variant is passed,
the predict endpoint returns a 400 error.

### Custom organism — what to upload

**When to use:** You are predicting on a virus other than the built-in ones
(SARS-CoV-2, Influenza A H1N1/H3N2/H5N1).

**What to upload:**
- Genome (FASTA) — your reference sequence → stored as `genome.fasta`
- Protein regions (CSV, optional) → stored as `protein_regions.csv`

**Form selection:** Target organism → **Other — upload my own genome**

**Protein regions CSV format** (`name,start,end`, 1-based inclusive):
```csv
name,start,end
S,1,1500
M,1501,2400
N,2401,3000
```

If you skip the protein regions file, predictions span the whole genome with
no protein-region breakdown. Coordinates are 1-based inclusive — the same
convention as standard biology. The `domain.region` in the payload uses
0-based half-open `[start, end)`, so subtract 1 from start when constructing
it programmatically.

---

## Plug 2 — Feature Extractor

**When to use:** Your model expects different features than the default (e.g.
different feature width, custom embeddings, no codon awareness).

**Filename stored as:** `feature_extractor.py`

Expose either a **class instance** or a **module-level function**:

```python
# feature_extractor.py — class pattern (recommended)
import numpy as np

class FeatureExtractor:
    def extract_features(self, genome_seq, **kwargs) -> np.ndarray:
        # genome_seq is a string of A/T/G/C/N.
        # kwargs swallows everything else — declare only what you use.
        base_dir = kwargs["base_dir"]   # your bundle directory
        one_hot = np.zeros((len(genome_seq), 4))
        for i, c in enumerate(genome_seq):
            idx = "ATGC".find(c.upper())
            if idx >= 0:
                one_hot[i, idx] = 1.0
        return one_hot

    def get_feature_dimension(self) -> int:   # recommended
        return 4

    def get_feature_description(self) -> dict:  # recommended
        return {"channels": "A, T, G, C one-hot"}

    def get_metadata(self) -> dict:             # recommended
        return {"name": "OneHotExtractor", "version": "1.0"}

# CovMutEx looks for a module-level instance named anything, OR a
# module-level extract_features() function.
ExtractorInstance = FeatureExtractor()
```

### Keyword arguments CovMutEx passes

Everything is passed as keyword args. Declare only the ones you use;
`**kwargs` swallows the rest:

| Kwarg | Type | What it carries |
|---|---|---|
| `genome_seq` | `str` | Resolved genome sequence (mutations/variant already applied) |
| `mutations` | `list` | COVID: parsed UShER mutations applied to reference. Empty for influenza / custom. |
| `node_ids` | `List[str]` | COVID: the node IDs selected in the UI. Empty for others. |
| `elapsed_day` | `int` | COVID-only temporal feature; 0 elsewhere. |
| `protein_regions` | `Dict[str, Tuple[int,int]]` | Active region dict (filtered to selection if any). |
| `k` | `int` | Window size; default 30. |
| `cache_path` | `str` | Optional cache dir for memoizing extraction. |
| `codon_mapper` | `str` | Path to codon → AA mapping JSON. |
| `config_file` | `Any` | Internal config object; safe to ignore. |
| `depth` | `float` | Sample depth (COVID-only). |
| `base_dir` | `str` | Your bundle directory. Use to open helper files. |
| `**custom_params` | `dict` | The `customParameters` JSON you declared at upload time. |

**Canonical `extract_features` signature** (module-level or instance method):

```python
def extract_features(
    genome_seq,
    mutations,
    node_ids,
    elapsed_day=None,
    protein_regions=None,
    k=30,
    **kwargs,          # catches base_dir, cache_path, custom_params, etc.
):
    ...
```

A module-level `nucleotides_per_position` attribute (default `1`) tells the
runtime how many channel rows the extractor emits per genome position; it is
forwarded to `postprocess()` via `context["nucleotides_per_position"]`.

---

## Plug 3 — Model Adapter

**When to use:** One of the following:
- Your model's input/output shape differs from the built-in Keras defaults
- Your model is not a Keras `.keras`/`.h5` file (PyTorch, sklearn, ONNX, …)
- You need custom preprocess/predict/postprocess logic

**Filename stored as:** `model_adapter.py`

Expose either `load_adapter(...)` or `ModelAdapter` class. **Only
`postprocess` is required.** Any method you omit falls back to the default
Keras wrapper.

### Minimal example — change only output format

```python
# model_adapter.py
class ModelAdapter:
    def __init__(self, **kwargs):
        pass  # model loading handled by default Keras wrapper

    def postprocess(self, raw, context=None):
        ctx = context or {}
        return {
            "contract_version": "2.0",
            "task": {"kind": "binary_per_position"},
            "domain": {
                "total_length": ctx.get("reference_length", len(raw)),
                "region": ctx.get("region"),
            },
            "predictions": {
                "indexing": "absolute",
                "values": raw.flatten().tolist(),
                "value_kind": "probability",
                "value_range": [0.0, 1.0],
            },
            "annotations": {},
        }
```

### Full example — PyTorch model with custom preprocess

```python
# model_adapter.py
import importlib.util, os, torch, numpy as np
from genome.covmutex_models import wrap_predictions_as_payload

class ModelAdapter:
    def __init__(self, model_path, base_dir=None, **kwargs):
        # Load model class from a helper file in the bundle
        spec = importlib.util.spec_from_file_location(
            "my_model", os.path.join(base_dir, "my_model.py")
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        self.model = mod.MyModel()
        self.model.load_state_dict(torch.load(model_path, map_location="cpu"))
        self.model.eval()

    def metadata(self):
        return {"name": "MyModel", "model_type": "pytorch", "version": "1.0"}

    def input_schema(self):
        return {"features": "(N, D) float32 array"}

    def preprocess(self, inputs):
        return torch.from_numpy(inputs["features"]).float().unsqueeze(0)

    def predict(self, batch):
        with torch.no_grad():
            return self.model(batch).squeeze(0).numpy()

    def postprocess(self, raw, context=None):
        return wrap_predictions_as_payload(
            np.asarray(raw),
            kind="scalar_per_position",
            total_length=(context or {}).get("reference_length", len(raw)),
            region=(context or {}).get("region"),
            value_kind="score",
            value_range=None,
        )

def load_adapter(model_path, base_dir=None, **kwargs):
    return ModelAdapter(model_path=model_path, base_dir=base_dir, **kwargs)
```

### Which methods are optional?

| Method | Required? | If omitted, default is… |
|---|---|---|
| `__init__` / `load_adapter` | If your class needs init | — |
| `metadata()` | No | Default Keras model metadata |
| `input_schema()` | No | Default Keras input schema |
| `preprocess(inputs)` | No | Default Keras reshape |
| `predict(batch)` | No | Default Keras `model.predict(batch)` |
| `postprocess(raw, context)` | **YES** | — |

### The `context` parameter

```python
context = {
    "nucleotides_per_position": 4,     # channel duplication hint from extractor
    "reference_length": 29903,         # full organism genome length
    "region": {"start": 100, "end": 500},  # only present if a region was selected
}
```

---

## Plug 4 — PredictionPayload v2.0

`postprocess()` must return a `PredictionPayload` dict. Build it with
`genome.covmutex_models.wrap_predictions_as_payload(...)` (recommended) or
by hand. It is validated by `validate_prediction_payload(payload)` before
it leaves the backend.

```jsonc
{
  "contract_version": "2.0",
  "task": {
    "kind": "categorical_per_position",   // "binary_per_position" | "scalar_per_position"
    "labels": ["A", "T", "G", "C"]        // REQUIRED for categorical; omit otherwise
  },
  "domain": {
    "total_length": 29903,                // positive int: full reference length
    "region": { "start": 28273, "end": 29533 }  // 0-based half-open, or null
  },
  "predictions": {
    "indexing": "absolute",               // "absolute" | "relative_to_region"
    "values": [ /* see shape rules */ ],
    "value_kind": "probability",          // "probability" | "logit" | "score"
    "value_range": [0.0, 1.0]            // or null if unbounded
  },
  "annotations": { /* free-form; e.g. protein_regions, region_of_interest */ }
}
```

### Task kinds, value shapes, and frontend rendering

Let `region_len = region.end − region.start` when `region` is set, otherwise
`region_len = total_length`.

| `task.kind` | `values` shape | Frontend renders |
|---|---|---|
| `categorical_per_position` | `(region_len, len(labels))` | GenomeChart (nucleotide bars) + DoughnutChart (protein breakdown) |
| `binary_per_position` | `(region_len,)` | PerPositionTrack — continuous probability lane (label: "Mutation probability") |
| `scalar_per_position` | `(region_len,)` | PerPositionTrack — continuous score lane (label: "Per-position score") |

### Validation invariants

- `task.kind` is one of the three kinds above.
- `categorical_per_position` requires a non-empty `task.labels` list.
- `domain.total_length` is a positive int.
- `domain.region`, when present, satisfies `0 <= start < end <= total_length`.
- `predictions.indexing` is `absolute` or `relative_to_region`.
- `predictions.value_kind` is `probability`, `logit`, or `score`.
- `predictions.values` is present and matches the shape rule for the task kind.

Shape mismatches surface as `400 An error occurred: <kind> values must have
shape (N, k), got (M, j)`. The most common cause is off-by-one in region
coordinates — `domain.region` uses **half-open `[start, end)`** while
`protein_regions.csv` uses **closed `[start, end]`**.

### `wrap_predictions_as_payload` examples

Whole-genome categorical (built-in COVID default):

```python
wrap_predictions_as_payload(
    per_position,              # shape (29903, 4)
    kind="categorical_per_position",
    labels=["A", "T", "G", "C"],
    total_length=29903,
    region=None,
    value_kind="probability",
    value_range=[0.0, 1.0],
)
```

Sub-region binary ("will this position mutate?" over spike protein):

```python
wrap_predictions_as_payload(
    scores,                    # shape (span_length,)
    kind="binary_per_position",
    total_length=29903,
    region={"start": 21562, "end": 25385},  # 0-based half-open
    value_kind="probability",
    value_range=[0.0, 1.0],
)
```

Whole-sequence scalar (influenza or custom organism):

```python
wrap_predictions_as_payload(
    scores,                    # shape (L,)
    kind="scalar_per_position",
    total_length=L,            # this organism's length, not 29903
    region=None,
    value_kind="score",
    value_range=None,
)
```

---

## Bring your own genome (custom organism)

To predict on an organism the platform doesn't ship, select
**Other — I'll upload my own genome** in the upload form:

| You upload | Stored in the bundle as | Required |
|---|---|---|
| Genome (FASTA) | `genome.fasta` | yes |
| Protein regions (CSV `name,start,end`, 1-based inclusive) | `protein_regions.csv` | no |

The runtime resolves the genome for you and passes it in as `genome_seq`, so
a custom-organism extractor usually does **not** need to read the file:

```python
def extract_features(genome_seq, mutations, node_ids, elapsed_day=None,
                     protein_regions=None, k=30, **kwargs):
    base_dir = kwargs["base_dir"]   # if you need sibling helper files
    # genome_seq is already the string content of your genome.fasta
    ...
```

Then declare the organism's real `total_length` in `postprocess()`. Nothing
in the contract assumes SARS-CoV-2 or length 29903.

---

## Custom parameters

If your extractor or model needs runtime knobs (threshold, window size, …),
declare them in the **Custom Parameters** section of the upload modal.
They land in two places:

1. Saved on disk as `custom_parameters.json` in your bundle
2. Passed as `**kwargs` to `feature_extractor.extract_features()`

Each parameter has a name, default value, type, and a `required` flag —
predict requests are rejected if a required parameter is missing.

---

## End-to-end prediction flow

```
1. Frontend → POST /api/predict/ with selectedModel="uploaded:my_model"
2. CovMutEx resolves your bundle on disk
3. Organism dispatch:
     covid     → built-in COVID genome + UShER mutation overlay
     influenza → chosen variant FASTA + reference protein_regions
     custom    → your genome.fasta + protein_regions.csv (if present)
4. feature_extractor.extract_features(genome_seq, mutations, node_ids, ...)
     OR default extractor if you didn't upload one
5. adapter.preprocess(inputs)
     OR default Keras preprocess
6. adapter.predict(batch)
     OR default Keras model.predict()
7. adapter.postprocess(raw, context)   ← required if adapter uploaded
     OR default Keras postprocess (produces categorical_per_position, A/T/G/C)
8. CovMutEx injects annotations.protein_regions, validates the payload,
   returns it as JSON.
9. Frontend PayloadRenderer dispatches on task.kind to the right chart.
```

---

## Common recipes

### "Model outputs one mutation-likelihood score per position"

Upload only `model.keras` plus a `model_adapter.py` with just `postprocess`:

```python
class ModelAdapter:
    def __init__(self, **kwargs): pass
    def postprocess(self, raw, context=None):
        ctx = context or {}
        return {
            "contract_version": "2.0",
            "task": {"kind": "scalar_per_position"},
            "domain": {"total_length": ctx.get("reference_length", len(raw)), "region": None},
            "predictions": {"indexing": "absolute", "values": raw.flatten().tolist(),
                            "value_kind": "score", "value_range": None},
            "annotations": {},
        }
```

### "Model wants k-mer counts, not the COVID feature pipeline"

Upload `feature_extractor.py`:

```python
from collections import Counter
import numpy as np

class FeatureExtractor:
    def extract_features(self, genome_seq, k=6, **kwargs):
        counts = Counter(genome_seq[i:i+k] for i in range(len(genome_seq)-k+1))
        return np.array(sorted(counts.values()), dtype=np.float32)
    def get_feature_dimension(self): return None
    def get_feature_description(self): return {"kind": "kmer-counts"}
    def get_metadata(self): return {"name": "KmerCounts", "version": "1.0"}

ExtractorInstance = FeatureExtractor()
```

### "My organism is HIV, not in your built-ins"

Select **Other** in the organism dropdown and upload:
- `genome.fasta` — HXB2 reference
- `protein_regions.csv` — HIV gene coordinates (gag, pol, env, …)

No extractor/adapter needed if your model accepts the default features and
outputs `(N, 4)` nucleotide probabilities. Otherwise add extractor/adapter
as described above.

### "Same model, multiple influenza strains"

Upload the model with `organism="influenza"` (the preferred generic value) and
no variant in the bundle. On every predict call, pass `variant="<strain_name>"`
in the request — CovMutEx looks up the strain across all influenza subtypes
(H1N1 / H3N2 / H5N1) and swaps in the correct FASTA before running your pipeline.

---

## What does NOT work

- **Renaming `feature_extractor.py` or `model_adapter.py`.** CovMutEx looks
  for exactly those filenames; alternatives are silently ignored.
- **Using a different output shape than what your declared `task.kind`
  promises.** The validator will return 400.
- **Mixing `organism="custom"` with a `variant` parameter.** Variants only
  apply to built-in influenza subtypes.
- **Reading files from outside your bundle directory at runtime.** Use only
  `base_dir`-relative paths.

---

## Naming conventions

- Model file → `model.<ext>` (always)
- Feature extractor → `feature_extractor.py` (always)
- Model adapter → `model_adapter.py` (always)
- Custom organism genome → `genome.fasta` (when uploaded through the form)
- Custom organism protein regions → `protein_regions.csv` (when uploaded)
- Helper files → keep original sanitized name
- `custom_parameters.json` → merged defaults, auto-written by the platform

## Security and path rules

- Uploaded filenames are sanitized before writing to disk
- Helper files must use an allowed extension
- All uploaded files are written into one bundle directory
- Plugins must only access bundle-local helper files or server-provided inputs

---

## Backward compatibility

- `MODEL_CONTRACT_VERSION = "2.0"` is reported in model metadata.
- Built-in COVID-19 models and the generic framework wrappers default to
  `categorical_per_position`, `labels = [A, T, G, C]`, `total_length = 29903`,
  `region = null`. Their payloads are unchanged in meaning.
- v1.0 bundles without a `model_adapter.py` keep working: the generic
  framework wrappers produce the payload, wrapping the model's raw array
  into a default `categorical_per_position` result.
- A bundle that ships its own `model_adapter.py` must return a
  `PredictionPayload` dict from `postprocess(raw, context=None)`.
