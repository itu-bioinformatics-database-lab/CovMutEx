# /genome_extractor/genome/runtime.py

from .legacy_adapter import LegacyModelAdapter
from .sdk_protocol import CovMutExModel # Protokolü ayrı bir dosyadan aldığımızı varsayalım

MODELS = {
    "legacy_v1": LegacyModelAdapter()
    # "new_pytorch_model_v1": PyTorchModelAdapter()
}

def get_model(model_name: str) -> CovMutExModel | None:
    return MODELS.get(model_name)