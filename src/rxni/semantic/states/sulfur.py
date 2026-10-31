"""Sulfur local-state classification."""

from __future__ import annotations

from typing import Any

from .common import atom, attachments


def sulfur_state(
    facts: dict[str, Any],
    side: str,
    map_id: int,
) -> str:
    center = atom(facts, side, map_id)
    if center.get("element") != "S":
        return ""

    n_oxo = sum(
        1
        for x in attachments(facts, side, map_id)
        if x.get("element") == "O"
        and float(x.get("bond_order", 0.0)) >= 1.5
    )

    if n_oxo >= 2:
        return "sulfone"
    if n_oxo == 1:
        return "sulfoxide"
    return "sulfide"
