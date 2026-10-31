"""Self-contained Runtime-v6.6 semantic pipeline — M6C3.

Both supported backends are now Runtime-v6.6 native:

``validated`` (default)
    migrated non-FGI rules
      -> v6.6 validation/enforcement
      -> canonical redox/profile

``reference``
    migrated frozen-v6.5.6 non-FGI rules represented through the v6.6
    SemanticResult contract, without importing Runtime-v6.5.

No Runtime-v6.5 module is imported or executed by this file.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import is_dataclass, replace
import inspect
from typing import Any

from .validated_pipeline import *  # noqa: F401,F403
from .validated_pipeline import (
    classify_reaction as _native_stagee_classify_reaction,
)
from .diagnostic_provenance import normalize_diagnostic_provenance
from .objective import native_formed_bond_classes


_NATIVE_SIGNATURE = inspect.signature(
    _native_stagee_classify_reaction
)


def _backend_for_call(args, kwargs) -> str:
    try:
        bound = _NATIVE_SIGNATURE.bind_partial(*args, **kwargs)
        value = bound.arguments.get("backend", None)
        if value is None:
            param = _NATIVE_SIGNATURE.parameters.get("backend")
            if (
                param is not None
                and param.default is not inspect._empty
            ):
                value = param.default
    except Exception:
        value = kwargs.get("backend", None)

    if value is None:
        return "validated"
    return str(value)


def _reaction_for_call(args, kwargs) -> str:
    try:
        bound = _NATIVE_SIGNATURE.bind_partial(*args, **kwargs)
        for name in (
            "mapped_reaction",
            "reaction",
            "reaction_smiles",
            "mapped_reaction_smiles",
        ):
            if name in bound.arguments:
                return str(bound.arguments[name])
    except Exception:
        pass

    for name in (
        "mapped_reaction",
        "reaction",
        "reaction_smiles",
        "mapped_reaction_smiles",
    ):
        if name in kwargs:
            return str(kwargs[name])

    if args:
        return str(args[0])

    raise TypeError("Could not identify mapped reaction argument")


def _coerce_like(reference_value: Any, native_values: list[Any]) -> Any:
    if isinstance(reference_value, tuple):
        return tuple(native_values)
    if isinstance(reference_value, list):
        return list(native_values)
    if isinstance(reference_value, set):
        return set(native_values)
    return native_values


def _replace_field(result: Any, field: str, value: Any) -> Any:
    if isinstance(result, dict):
        out = deepcopy(result)
        out[field] = value
        return out

    if is_dataclass(result):
        return replace(result, **{field: value})

    if hasattr(result, "_replace"):
        try:
            return result._replace(**{field: value})
        except Exception:
            pass

    out = deepcopy(result)
    setattr(out, field, value)
    return out


def classify_reaction(*args, **kwargs):
    backend = _backend_for_call(args, kwargs)

    if backend not in {"validated", "reference"}:
        # Preserve the native Stage-E API as authority for unsupported backend
        # behavior/error reporting.
        return _native_stagee_classify_reaction(*args, **kwargs)

    result = _native_stagee_classify_reaction(*args, **kwargs)

    # Preserve M6C1B objective ownership explicitly for both public backends.
    reaction = _reaction_for_call(args, kwargs)
    native_fb = native_formed_bond_classes(reaction)
    current_fb = (
        result.get("formed_bond_classes")
        if isinstance(result, dict)
        else getattr(result, "formed_bond_classes")
    )
    native_fb = _coerce_like(current_fb, native_fb)
    result = _replace_field(
        result,
        "formed_bond_classes",
        native_fb,
    )

    current_diag = (
        result.get("diagnostics")
        if isinstance(result, dict)
        else getattr(result, "diagnostics")
    )
    normalized_diag = normalize_diagnostic_provenance(
        current_diag,
        backend=backend,
    )
    result = _replace_field(
        result,
        "diagnostics",
        normalized_diag,
    )

    return result
