"""Unified Semantic Layer public API.

The frozen-reference backend is imported lazily so contracts/types/validation
modules remain independent of the historical runtime.
"""

from .types import SemanticRecord, SemanticResult, SemanticValidation, ValidationStatus
from .version import (
    RUNTIME_VERSION,
    SEMANTIC_ENGINE_VERSION,
    SEMANTIC_RULESET_VERSION,
    FGI_CONTRACT_VERSION,
)


def classify_reaction(*args, **kwargs):
    from .pipeline import classify_reaction as _classify_reaction
    return _classify_reaction(*args, **kwargs)


__all__ = [
    "classify_reaction",
    "SemanticRecord",
    "SemanticResult",
    "SemanticValidation",
    "ValidationStatus",
    "RUNTIME_VERSION",
    "SEMANTIC_ENGINE_VERSION",
    "SEMANTIC_RULESET_VERSION",
    "FGI_CONTRACT_VERSION",
]
