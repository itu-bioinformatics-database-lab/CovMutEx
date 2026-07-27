import json
import os
import re
import shutil
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .security_scanner import scan_folder_for_malware, validate_uploaded_files


PLUGIN_CONTRACT_VERSION = "2.0"


@dataclass
class HelperUpload:
    index: str
    file_obj: Any
    stored_name: str


@dataclass
class BundleResolution:
    source: str
    model_path: str
    extractor_path: Optional[str]
    adapter_path: Optional[str]
    bundle_dir: Optional[str]
    bundle_name: Optional[str]
    model_name: str
    selected_model: str
    custom_parameters: Dict[str, Any]
    saved_folder: Optional[str]
    # Organism dispatch — read from bundle_metadata.json when present.
    # Default "covid" keeps existing COVID bundles & built-in models backward compatible.
    organism: str = "covid"
    # When organism == "custom", these point at helper files inside bundle_dir.
    # For built-in organisms they are None and the runtime resolves genome/regions
    # from organisms/<organism>/ via genome.organism_registry.
    genome_file: Optional[str] = None
    protein_regions_file: Optional[str] = None
    # Optional built-in-organism variant override. When set, the runtime uses
    # the catalog variant's FASTA in place of the organism's reference genome;
    # protein_regions still come from the reference. Custom organisms can't
    # carry a variant — the bundle already ships its own genome.
    variant: Optional[str] = None


INVALID_PATH_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def sanitize_bundle_name(name: Optional[str], fallback: str = "upload_bundle") -> str:
    raw_name = (name or fallback).strip()
    sanitized = INVALID_PATH_CHARS.sub("_", raw_name).strip("._")
    return sanitized or fallback


def sanitize_uploaded_filename(name: Optional[str], fallback: str) -> str:
    base_name = os.path.basename((name or "").strip()) or fallback
    if "." not in base_name and "." in fallback:
        stem, ext = os.path.splitext(fallback)
        base_name = f"{sanitize_bundle_name(base_name, stem)}{ext}"
    stem, ext = os.path.splitext(base_name)
    safe_stem = sanitize_bundle_name(stem, os.path.splitext(fallback)[0])
    safe_ext = INVALID_PATH_CHARS.sub("", ext) if ext else os.path.splitext(fallback)[1]
    return f"{safe_stem}{safe_ext}"


def parse_custom_parameters(raw_value: Any) -> Dict[str, Any]:
    if raw_value in (None, "", {}):
        return {}
    if isinstance(raw_value, dict):
        return raw_value
    if isinstance(raw_value, str):
        try:
            parsed = json.loads(raw_value)
        except json.JSONDecodeError as exc:
            raise ValueError(f"customParameters must be valid JSON: {exc}") from exc
        if not isinstance(parsed, dict):
            raise ValueError("customParameters must decode to a JSON object")
        return parsed
    raise ValueError("customParameters must be a JSON string or object")


def normalize_node_ids(node_ids: Any) -> List[str]:
    if node_ids is None:
        return []
    if isinstance(node_ids, str):
        return [node_ids] if node_ids else []
    if isinstance(node_ids, (list, tuple, set)):
        return [str(item) for item in node_ids if item]
    return [str(node_ids)]


def normalize_elapsed_day(value: Any) -> Optional[int]:
    if value in (None, ""):
        return None
    return int(value)


def collect_helper_uploads(files: Any, data: Any) -> List[HelperUpload]:
    helper_uploads: List[HelperUpload] = []
    for key in files.keys():
        if not key.startswith("helperFile_"):
            continue
        index = key.replace("helperFile_", "")
        fallback_name = f"helper_{index}.bin"
        requested_name = data.get(f"helperFileName_{index}", fallback_name)
        helper_uploads.append(
            HelperUpload(
                index=index,
                file_obj=files[key],
                stored_name=sanitize_uploaded_filename(requested_name, fallback_name),
            )
        )
    return helper_uploads


def ensure_python_package(path: str) -> None:
    os.makedirs(path, exist_ok=True)
    init_path = os.path.join(path, "__init__.py")
    if not os.path.exists(init_path):
        with open(init_path, "w", encoding="utf-8") as handle:
            handle.write("")


def _write_uploaded_file(file_obj: Any, destination_path: str) -> None:
    with open(destination_path, "wb+") as destination:
        for chunk in file_obj.chunks():
            destination.write(chunk)


def _find_model_file(bundle_dir: str) -> Optional[str]:
    for filename in sorted(os.listdir(bundle_dir)):
        if filename.startswith("model."):
            candidate = os.path.join(bundle_dir, filename)
            if os.path.isfile(candidate):
                return candidate
    return None


def _find_optional_bundle_file(bundle_dir: str, filename: str) -> Optional[str]:
    candidate = os.path.join(bundle_dir, filename)
    if os.path.isfile(candidate):
        return candidate
    return None


ALLOWED_ORGANISM_VALUES = (
    "covid",
    # Generic influenza A: bundle declares "influenza" and the prediction-time
    # variant resolves to the right HA subtype (H1N1 / H3N2 / H5N1). This is
    # the preferred form for new uploads.
    "influenza",
    # Subtype-specific aliases kept for bundles that already declared one.
    # _normalize_organism_fields() folds them back to "influenza" when a
    # variant is provided, so dispatch can pick the right protein_regions.
    "influenza_h1n1",
    "influenza_h3n2",
    "influenza_h5n1",
    "custom",
)


def read_bundle_metadata(bundle_dir: Optional[str]) -> Dict[str, Any]:
    """Load ``bundle_metadata.json`` from a bundle directory, or return ``{}``.

    Bundles ship metadata describing things the runtime acts on (organism,
    genome_file, protein_regions_file, …). For backward compatibility a bundle
    without metadata is treated as ``organism = "covid"`` with no helper-file
    overrides.
    """
    if not bundle_dir:
        return {}
    path = os.path.join(bundle_dir, "bundle_metadata.json")
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as handle:
            loaded = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}
    return loaded if isinstance(loaded, dict) else {}


def _normalize_organism_fields(
    metadata: Dict[str, Any],
    bundle_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Pull the organism dispatch fields out of bundle metadata, defaulting
    safely. Unknown organism values fall back to ``"covid"``.

    For ``"custom"`` organisms the helper files (genome.fasta,
    protein_regions.csv) are detected on disk under canonical names — the
    manifest does not need to declare them. Legacy bundles that did declare
    ``genome_file`` / ``protein_regions_file`` are still respected.
    """
    organism_raw = (metadata.get("organism") or "covid").strip().lower()
    organism = organism_raw if organism_raw in ALLOWED_ORGANISM_VALUES else "covid"

    genome_file = metadata.get("genome_file") or None
    protein_regions_file = metadata.get("protein_regions_file") or None
    variant_raw = metadata.get("variant")
    variant = variant_raw.strip() if isinstance(variant_raw, str) and variant_raw.strip() else None

    # Canonical-name fallback: if metadata is silent but the file exists on
    # disk in the bundle dir under the documented canonical name, use it.
    if organism == "custom" and bundle_dir:
        if not genome_file:
            candidate = os.path.join(bundle_dir, "genome.fasta")
            if os.path.exists(candidate):
                genome_file = "genome.fasta"
        if not protein_regions_file:
            candidate = os.path.join(bundle_dir, "protein_regions.csv")
            if os.path.exists(candidate):
                protein_regions_file = "protein_regions.csv"

    if organism == "custom" and not genome_file:
        raise ValueError(
            'Bundle declares organism="custom" but no genome file was found '
            '(expected genome.fasta in the bundle directory).'
        )

    if organism == "custom" and variant:
        raise ValueError(
            'Custom organisms cannot declare a "variant"; the bundle already '
            'ships its own genome. Drop the variant field or switch to a '
            'built-in organism (covid / influenza_h1n1 / influenza_h5n1).'
        )

    return {
        "organism": organism,
        "genome_file": genome_file,
        "protein_regions_file": protein_regions_file,
        "variant": variant,
    }


RESERVED_BUNDLE_FILENAMES = frozenset({
    "bundle_metadata.json",     # auto-generated from the organism dropdown
    "custom_parameters.json",   # auto-generated from the customParameters form
    "feature_extractor.py",     # dedicated form slot
    "model_adapter.py",         # dedicated form slot
    "genome.fasta",             # dedicated form slot (custom organism)
    "protein_regions.csv",      # dedicated form slot (custom organism)
    "__init__.py",              # bundle Python-package init
})


def _assert_no_reserved_helper_names(helper_uploads: List[HelperUpload]) -> None:
    """Reject helper uploads that would silently overwrite a file the platform
    writes itself (or a file that has its own dedicated form slot). Names are
    matched case-insensitively, and ``model.*`` extensions are blocked too
    because the bundle's model file is auto-named ``model.<ext>`` from the
    modelFile slot.
    """
    for helper in helper_uploads:
        name = (helper.stored_name or "").lower()
        if not name:
            continue
        if name in RESERVED_BUNDLE_FILENAMES:
            raise ValueError(
                f'Cannot upload "{helper.stored_name}" as a helper file — that '
                f"name is reserved (the platform writes it from a dedicated "
                f"form slot or auto-generates it). Use the matching slot instead."
            )
        if name.startswith("model.") and name.split(".")[-1] in (
            "keras", "h5", "hdf5", "pb", "pt", "pth", "pkl", "onnx", "savedmodel", "bin", "safetensors"
        ):
            raise ValueError(
                f'Cannot upload "{helper.stored_name}" as a helper file — the '
                f'main model file is written from the "Model File" slot using '
                f'the canonical name "model.<ext>". Use that slot instead.'
            )


def save_uploaded_bundle(
    uploaded_models_dir: str,
    uploaded_model: Any,
    uploaded_extractor: Any,
    helper_uploads: List[HelperUpload],
    upload_folder_name: Optional[str],
    custom_parameters: Dict[str, Any],
    organism_metadata: Optional[Dict[str, Any]] = None,
) -> BundleResolution:
    helper_names = [item.stored_name for item in helper_uploads] or None
    is_valid, error_msg = validate_uploaded_files(
        model_file=uploaded_model,
        extractor_file=uploaded_extractor,
        helper_files=helper_names,
    )
    if not is_valid:
        raise ValueError(error_msg)

    # Block helper uploads from shadowing canonical filenames that the
    # platform writes itself (e.g. bundle_metadata.json, feature_extractor.py).
    # Helpers added via dedicated slots are tagged with stored_name set BY US
    # (e.g. CUSTOM_GENOME_FILENAME) — those are safe, only check user-supplied
    # names where stored_name came from the upload itself.
    _assert_no_reserved_helper_names([
        h for h in helper_uploads if h.index not in (
            "genome", "protein_regions", "adapter"
        )
    ])

    ensure_python_package(uploaded_models_dir)
    bundle_name = sanitize_bundle_name(upload_folder_name, "uploaded_model")
    bundle_dir = os.path.join(uploaded_models_dir, bundle_name)

    if os.path.exists(bundle_dir):
        shutil.rmtree(bundle_dir)
    ensure_python_package(bundle_dir)

    model_ext = os.path.splitext(uploaded_model.name)[1]
    canonical_model_name = f"model{model_ext}"
    model_path = os.path.join(bundle_dir, canonical_model_name)
    _write_uploaded_file(uploaded_model, model_path)

    # HuggingFace from_pretrained() looks for "pytorch_model.bin" by name.
    # Always create that symlink for .bin models, regardless of upload filename.
    # Also create a symlink under the original upload name in case the adapter
    # references it explicitly.
    original_model_name = uploaded_model.name
    for alias in {original_model_name, f"pytorch_model{model_ext}" if model_ext == ".bin" else None} - {None, canonical_model_name}:
        alias_path = os.path.join(bundle_dir, alias)
        if not os.path.exists(alias_path):
            os.symlink(canonical_model_name, alias_path)

    extractor_path = None
    if uploaded_extractor:
        extractor_path = os.path.join(bundle_dir, "feature_extractor.py")
        _write_uploaded_file(uploaded_extractor, extractor_path)

    for helper_upload in helper_uploads:
        helper_path = os.path.join(bundle_dir, helper_upload.stored_name)
        _write_uploaded_file(helper_upload.file_obj, helper_path)

    all_clean, scan_results = scan_folder_for_malware(bundle_dir)
    if not all_clean:
        shutil.rmtree(bundle_dir, ignore_errors=True)
        raise ValueError(
            json.dumps(
                {
                    "message": "One or more uploaded files contain malware",
                    "scan_results": scan_results,
                }
            )
        )

    if custom_parameters:
        params_path = os.path.join(bundle_dir, "custom_parameters.json")
        with open(params_path, "w", encoding="utf-8") as handle:
            json.dump(custom_parameters, handle, indent=2)

    # Organism dispatch metadata: write bundle_metadata.json with just the
    # organism marker. Genome / protein_regions / adapter / extractor presence
    # is detected from canonical filenames on disk, so the manifest stays
    # minimal.
    organism_metadata = dict(organism_metadata or {})
    organism_fields = _normalize_organism_fields(organism_metadata, bundle_dir=bundle_dir)
    if organism_metadata:
        meta_to_write = {"organism": organism_fields["organism"]}
        if organism_fields.get("variant"):
            meta_to_write["variant"] = organism_fields["variant"]
        meta_path = os.path.join(bundle_dir, "bundle_metadata.json")
        with open(meta_path, "w", encoding="utf-8") as handle:
            json.dump(meta_to_write, handle, indent=2)

    return BundleResolution(
        source="uploaded",
        model_path=model_path,
        extractor_path=extractor_path,
        adapter_path=_find_optional_bundle_file(bundle_dir, "model_adapter.py"),
        bundle_dir=bundle_dir,
        bundle_name=bundle_name,
        model_name=uploaded_model.name,
        selected_model=f"uploaded:{bundle_name}",
        custom_parameters=custom_parameters,
        saved_folder=bundle_name,
        organism=organism_fields["organism"],
        genome_file=organism_fields["genome_file"],
        protein_regions_file=organism_fields["protein_regions_file"],
        variant=organism_fields["variant"],
    )


def resolve_uploaded_bundle(
    uploaded_models_dir: str,
    selected_model: str,
    request_custom_parameters: Dict[str, Any],
) -> BundleResolution:
    bundle_name = sanitize_bundle_name(selected_model.replace("uploaded:", ""), "uploaded_model")
    bundle_dir = os.path.join(uploaded_models_dir, bundle_name)
    if not os.path.isdir(bundle_dir):
        raise FileNotFoundError(f"Uploaded model bundle not found: {bundle_name}")

    model_path = _find_model_file(bundle_dir)
    if not model_path:
        raise FileNotFoundError(f"No model.* file found in uploaded bundle: {bundle_name}")

    extractor_path = os.path.join(bundle_dir, "feature_extractor.py")
    if not os.path.exists(extractor_path):
        extractor_path = None

    adapter_path = _find_optional_bundle_file(bundle_dir, "model_adapter.py")

    saved_parameters: Dict[str, Any] = {}
    params_path = os.path.join(bundle_dir, "custom_parameters.json")
    if os.path.exists(params_path):
        with open(params_path, "r", encoding="utf-8") as handle:
            loaded = json.load(handle)
        if isinstance(loaded, dict):
            saved_parameters = loaded

    merged_parameters = saved_parameters.copy()
    merged_parameters.update(request_custom_parameters)

    organism_fields = _normalize_organism_fields(read_bundle_metadata(bundle_dir), bundle_dir=bundle_dir)

    return BundleResolution(
        source="uploaded",
        model_path=model_path,
        extractor_path=extractor_path,
        adapter_path=adapter_path,
        bundle_dir=bundle_dir,
        bundle_name=bundle_name,
        model_name=os.path.basename(model_path),
        selected_model=f"uploaded:{bundle_name}",
        custom_parameters=merged_parameters,
        saved_folder=bundle_name,
        organism=organism_fields["organism"],
        genome_file=organism_fields["genome_file"],
        protein_regions_file=organism_fields["protein_regions_file"],
        variant=organism_fields["variant"],
    )


def resolve_server_model(model_directory: str, selected_model: Optional[str]) -> BundleResolution:
    resolved_name = selected_model or "balanced_data_model"
    candidate_paths = [
        os.path.join(model_directory, f"{resolved_name}.keras"),
        os.path.join(model_directory, f"{resolved_name}.h5"),
        os.path.join(model_directory, f"{resolved_name}.hdf5"),
        os.path.join(model_directory, f"{resolved_name}.pb"),
        os.path.join(model_directory, f"{resolved_name}.pt"),
        os.path.join(model_directory, f"{resolved_name}.pth"),
        os.path.join(model_directory, f"{resolved_name}.pkl"),
    ]

    model_path = next((path for path in candidate_paths if os.path.exists(path)), None)
    if not model_path:
        raise FileNotFoundError(f"Model not found for selection: {resolved_name}")

    return BundleResolution(
        source="server",
        model_path=model_path,
        extractor_path=None,
        adapter_path=None,
        bundle_dir=None,
        bundle_name=None,
        model_name=resolved_name,
        selected_model=resolved_name,
        custom_parameters={},
        saved_folder=None,
        # Server-side built-in models target the COVID reference genome.
        organism="covid",
    )
