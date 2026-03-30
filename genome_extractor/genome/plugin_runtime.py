import json
import os
import re
import shutil
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .security_scanner import scan_folder_for_malware, validate_uploaded_files


PLUGIN_CONTRACT_VERSION = "1.0"


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
    bundle_dir: Optional[str]
    bundle_name: Optional[str]
    model_name: str
    selected_model: str
    custom_parameters: Dict[str, Any]
    saved_folder: Optional[str]


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


def save_uploaded_bundle(
    uploaded_models_dir: str,
    uploaded_model: Any,
    uploaded_extractor: Any,
    helper_uploads: List[HelperUpload],
    upload_folder_name: Optional[str],
    custom_parameters: Dict[str, Any],
) -> BundleResolution:
    helper_names = [item.stored_name for item in helper_uploads] or None
    is_valid, error_msg = validate_uploaded_files(
        model_file=uploaded_model,
        extractor_file=uploaded_extractor,
        helper_files=helper_names,
    )
    if not is_valid:
        raise ValueError(error_msg)

    ensure_python_package(uploaded_models_dir)
    bundle_name = sanitize_bundle_name(upload_folder_name, "uploaded_model")
    bundle_dir = os.path.join(uploaded_models_dir, bundle_name)

    if os.path.exists(bundle_dir):
        shutil.rmtree(bundle_dir)
    ensure_python_package(bundle_dir)

    model_ext = os.path.splitext(uploaded_model.name)[1]
    model_path = os.path.join(bundle_dir, f"model{model_ext}")
    _write_uploaded_file(uploaded_model, model_path)

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

    return BundleResolution(
        source="uploaded",
        model_path=model_path,
        extractor_path=extractor_path,
        bundle_dir=bundle_dir,
        bundle_name=bundle_name,
        model_name=uploaded_model.name,
        selected_model=f"uploaded:{bundle_name}",
        custom_parameters=custom_parameters,
        saved_folder=bundle_name,
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

    saved_parameters: Dict[str, Any] = {}
    params_path = os.path.join(bundle_dir, "custom_parameters.json")
    if os.path.exists(params_path):
        with open(params_path, "r", encoding="utf-8") as handle:
            loaded = json.load(handle)
        if isinstance(loaded, dict):
            saved_parameters = loaded

    merged_parameters = saved_parameters.copy()
    merged_parameters.update(request_custom_parameters)

    return BundleResolution(
        source="uploaded",
        model_path=model_path,
        extractor_path=extractor_path,
        bundle_dir=bundle_dir,
        bundle_name=bundle_name,
        model_name=os.path.basename(model_path),
        selected_model=f"uploaded:{bundle_name}",
        custom_parameters=merged_parameters,
        saved_folder=bundle_name,
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
        bundle_dir=None,
        bundle_name=None,
        model_name=resolved_name,
        selected_model=resolved_name,
        custom_parameters={},
        saved_folder=None,
    )
