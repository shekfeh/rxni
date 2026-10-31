"""Canonical construction of retrieval semantic tokens.

Stage E invariant:
- P tokens come from canonical primary events after conversion filtering.
- D tokens come only from ACCEPTED directional conversions.
- R tokens come only from the canonical redox axis derived from those accepted
  directional conversions.
"""

from __future__ import annotations

from ..types import SemanticResult
from ..version import PROFILE_POLICY_CANONICAL, PROFILE_POLICY_PARITY


def _add(tok: dict[str, float], key: str, weight: float) -> None:
    if key:
        tok[key] = max(tok.get(key, 0.0), float(weight))


def build_semantic_profile(
    result: SemanticResult,
    *,
    policy: str = PROFILE_POLICY_CANONICAL,
) -> dict[str, float]:
    tok: dict[str, float] = {}

    for rec in result.primary_events:
        ev = rec.event.strip()
        st = (rec.subtype or "").strip()
        if ev:
            _add(tok, f"P:{ev}", 3.0)
            if st and st != ev:
                _add(tok, f"PS:{ev}|{st}", 4.0)

    for rec in result.directional_conversions:
        ev = rec.event.strip()
        if ev:
            _add(tok, f"D:{ev}", 4.0)

    if policy == PROFILE_POLICY_CANONICAL:
        active = {
            str(rec.event).strip()
            for rec in result.redox_axis
            if str(rec.event).strip() in {"oxidation", "reduction"}
        }
        if active:
            for direction in sorted(active):
                _add(tok, f"R:{direction}", 0.75)
        elif result.directional_conversions:
            _add(tok, "R:neutral", 0.75)

    elif policy == PROFILE_POLICY_PARITY:
        # Explicit compatibility mode only.  This is retained for audit
        # comparisons and is no longer the canonical Stage-E default.
        for rec in result.directional_conversions:
            direction = (rec.redox_direction or "").strip()
            if direction in {"oxidation", "reduction", "neutral"}:
                _add(tok, f"R:{direction}", 0.75)
    else:
        raise ValueError(f"Unknown semantic profile policy: {policy}")

    for bond_class in result.formed_bond_classes:
        _add(tok, f"B:{bond_class}", 1.0)

    return dict(sorted(tok.items()))
