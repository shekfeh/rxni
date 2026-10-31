"""Acyl-center semantic hypotheses only."""

from __future__ import annotations

from typing import Any


def detect_acyl_transitions(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    from . import _generated_fgi_rules as r
    return r.detect_acyl_conversions(facts)
