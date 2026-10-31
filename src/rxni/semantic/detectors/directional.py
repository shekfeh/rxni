"""Unified directional-hypothesis orchestration.

The return value is raw detector evidence. Nothing returned here is semantic
truth until it passes the validation layer.
"""

from __future__ import annotations

from typing import Any

from .acyl import detect_acyl_transitions
from .extensions import detect_structural_semantic_extensions
from .functional_group import (
    detect_core_functional_group_transitions,
    detect_directional_functional_group_extensions,
    detect_nitrile_transitions,
    detect_nitro_hydroxylamine_transition,
)


DETECTORS = (
    ("functional_group", detect_core_functional_group_transitions),
    ("acyl", detect_acyl_transitions),
    ("nitrile", detect_nitrile_transitions),
    ("structural_extensions", detect_structural_semantic_extensions),
    ("directional_extensions", detect_directional_functional_group_extensions),
    ("nitro_hydroxylamine", detect_nitro_hydroxylamine_transition),
)


def run_directional_hypothesis_detectors(
    facts: dict[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    return {name: fn(facts) for name, fn in DETECTORS}


def directional_hypotheses(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for records in run_directional_hypothesis_detectors(facts).values():
        out.extend(
            r
            for r in records
            if r.get("role") == "directional_conversion"
        )
    return out
