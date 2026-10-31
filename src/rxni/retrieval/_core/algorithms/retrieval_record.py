"""Canonical System-I retrieval-record construction.

semantic_profile() is behavior-preserved from Tool 92.
"""

from __future__ import annotations

from typing import Any

from .semantic_adapter import classify_for_retrieval
from .retrieval_representation import (
    extract_reaction_centers,
    generate_localized_fingerprints,
    generate_multichannel_fingerprint,
)

def semantic_profile(sem: dict[str, Any]) -> dict[str, float]:
    """Canonical retrieval semantics, deliberately excluding generic bookkeeping.

    Template truth is never included. We use primary events + directional conversions
    + formed-bond classes. This also avoids letting generic fallback specific-events
    (e.g. carbonyl_substitution bookkeeping) dominate the first experiment.
    """
    tok: dict[str, float] = {}

    def add(key: str, weight: float):
        if key:
            tok[key] = max(tok.get(key, 0.0), float(weight))

    for rec in sem.get('primary_events', []) or []:
        ev = str(rec.get('event') or '').strip()
        st = str(rec.get('subtype') or '').strip()
        if ev:
            add(f'P:{ev}', 3.0)
            if st and st != ev:
                add(f'PS:{ev}|{st}', 4.0)

    for rec in sem.get('directional_conversions', []) or []:
        ev = str(rec.get('event') or '').strip()
        direction = str(rec.get('redox_direction') or '').strip()
        if ev:
            add(f'D:{ev}', 4.0)
        if direction in {'oxidation', 'reduction', 'neutral'}:
            add(f'R:{direction}', 0.75)

    for bond_class in sem.get('formed_bond_classes', []) or []:
        add(f'B:{bond_class}', 1.0)

    # Center multiplicity is a weak structural-semantic descriptor, not a family label.
    return dict(sorted(tok.items()))

def build_retrieval_record(
    mapped_reaction: str,
    *,
    classification: dict[str, Any] | None = None,
    radius: int = 2,
    fp_size: int = 2048,
    record_id: str = "",
) -> dict[str, Any]:
    """Build the representation contract consumed by Tool-93 scoring."""

    sem = (
        classification
        if classification is not None
        else classify_for_retrieval(mapped_reaction)
    )

    centers = extract_reaction_centers(
        mapped_reaction,
        classification=sem,
    )

    localized = generate_localized_fingerprints(
        mapped_reaction,
        classification=sem,
        radius=radius,
        fp_size=fp_size,
    )

    multi = generate_multichannel_fingerprint(
        mapped_reaction,
        classification=sem,
        radius=radius,
        fp_size=fp_size,
    )

    return {
        "record_id": record_id,
        "reaction": mapped_reaction,
        "semantic_profile": semantic_profile(sem),
        "semantics": sem,
        "reaction_centers": centers,
        "localized_fingerprints": localized,
        "multichannel_fingerprint": multi,
    }
