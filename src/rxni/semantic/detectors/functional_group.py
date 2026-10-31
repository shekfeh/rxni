"""Directional functional-group hypothesis detectors.

These functions propose structural hypotheses only. They do not validate them,
select primary ownership, derive canonical redox, or build retrieval profiles.
"""

from __future__ import annotations

from typing import Any


def _rules():
    from . import _generated_fgi_rules as r
    return r


def detect_core_functional_group_transitions(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    return _rules().detect_functional_group_conversions(facts)


def detect_nitrile_transitions(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    return _rules().detect_nitrile_thioamide_and_ester(facts)


def detect_directional_functional_group_extensions(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    return _rules().detect_directional_functional_group_extensions(facts)


def detect_nitro_hydroxylamine_transition(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    return _rules().detect_nitro_to_hydroxylamine(facts)
