import importlib
import os
import sys
from typing import Any, Optional

from .covmutex_models import MODEL_CONTRACT_VERSION


# All methods on the CovMutExModel protocol. The user adapter may implement
# any subset; the rest are delegated to a default Keras wrapper loaded from
# the bundle's model file.
_ADAPTER_METHODS = ("metadata", "input_schema", "preprocess", "predict", "postprocess")


def _load_uploaded_adapter_module(adapter_path: str) -> Any:
    adapter_dir = os.path.dirname(os.path.abspath(adapter_path))
    bundle_name = os.path.basename(adapter_dir)
    uploaded_models_dir = os.path.dirname(adapter_dir)
    project_root = os.path.dirname(uploaded_models_dir)

    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    module_name = f"uploaded_models.{bundle_name}.model_adapter"
    importlib.invalidate_caches()

    if module_name in sys.modules:
        return importlib.reload(sys.modules[module_name])
    return importlib.import_module(module_name)


class _HybridAdapter:
    """Wraps a partial user adapter. Any of the 5 protocol methods that the
    user did not implement get delegated to a default Keras model wrapper
    loaded lazily from the bundle's model file.

    This lets a user override just ``postprocess`` (the common case — change
    the output shape) without having to re-implement preprocess / predict /
    metadata / input_schema for a standard Keras model.
    """

    def __init__(self, user_adapter: Any, default_factory):
        self._user = user_adapter
        self._default_factory = default_factory
        self._default = None

    def _get_default(self):
        if self._default is None:
            self._default = self._default_factory()
        return self._default

    def _dispatch(self, method_name, *args, **kwargs):
        if hasattr(self._user, method_name):
            return getattr(self._user, method_name)(*args, **kwargs)
        return getattr(self._get_default(), method_name)(*args, **kwargs)

    def metadata(self):
        return self._dispatch("metadata")

    def input_schema(self):
        return self._dispatch("input_schema")

    def preprocess(self, inputs):
        return self._dispatch("preprocess", inputs)

    def predict(self, batch):
        return self._dispatch("predict", batch)

    def postprocess(self, raw, context=None):
        # postprocess is the only method we required from the user — but we
        # still check hasattr so future relaxations are forward-compatible.
        return self._dispatch("postprocess", raw, context=context)


def _instantiate_user_adapter(
    module: Any,
    model_path: str,
    model_name: Optional[str],
    description: Optional[str],
    source: str,
    base_dir: str,
) -> Any:
    """Try to instantiate the user's adapter, tolerating constructors that
    don't accept any of the standard kwargs (common when the user only
    overrides postprocess and doesn't need the model path)."""
    full_kwargs = dict(
        model_path=model_path,
        model_name=model_name,
        description=description,
        source=source,
        base_dir=base_dir,
    )
    if hasattr(module, "load_adapter"):
        try:
            return module.load_adapter(**full_kwargs)
        except TypeError:
            return module.load_adapter()
    if hasattr(module, "ModelAdapter"):
        try:
            return module.ModelAdapter(**full_kwargs)
        except TypeError:
            return module.ModelAdapter()
    raise ValueError(
        "Uploaded model adapter must expose either load_adapter(...) or ModelAdapter"
    )


def load_model_adapter(
    adapter_path: str,
    model_path: str,
    model_name: Optional[str] = None,
    description: Optional[str] = None,
    source: str = "uploaded",
):
    module = _load_uploaded_adapter_module(adapter_path)
    base_dir = os.path.dirname(os.path.abspath(adapter_path))

    user_instance = _instantiate_user_adapter(
        module, model_path, model_name, description, source, base_dir
    )

    if not hasattr(user_instance, "postprocess"):
        raise ValueError(
            "Uploaded model adapter must define postprocess(raw, context=None) → dict. "
            "Other methods (metadata, input_schema, preprocess, predict) are optional — "
            "they fall back to the default Keras model wrapper if omitted."
        )

    # Lazily load the default Keras wrapper as a fallback for any methods the
    # user didn't implement. Importing here avoids a circular import at the
    # top of the module.
    def _default_factory():
        from .covmutex_models import load_model as load_default_model

        return load_default_model(
            model_path=model_path,
            model_name=model_name,
            description=description,
            source=source,
        )

    hybrid = _HybridAdapter(user_instance, _default_factory)

    metadata = hybrid.metadata()
    if isinstance(metadata, dict):
        metadata.setdefault("contract_version", MODEL_CONTRACT_VERSION)
        metadata.setdefault("adapter_module_path", adapter_path)
        metadata.setdefault("base_directory", base_dir)

    return hybrid
