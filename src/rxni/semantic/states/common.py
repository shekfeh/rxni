"""Local structural-state utilities for Runtime v6.6.

These helpers inspect already-built structural facts. They do not assign
reaction labels, primary events, redox axes, or retrieval tokens.
"""

from __future__ import annotations

from typing import Any

HALOGENS = {"F", "Cl", "Br", "I"}


def attachment_map(att: dict[str, Any]) -> int | None:
    v = att.get("map")
    return int(v) if v is not None else None


def atom(facts: dict[str, Any], side: str, map_id: int) -> dict[str, Any]:
    return facts[side]["atoms"].get(int(map_id), {})


def attachments(
    facts: dict[str, Any],
    side: str,
    map_id: int,
) -> list[dict[str, Any]]:
    return facts[side]["attachments"].get(int(map_id), [])


def bond_order(
    facts: dict[str, Any],
    side: str,
    a: int,
    b: int,
) -> float:
    key = tuple(sorted((int(a), int(b))))
    return float(facts[side]["bonds"].get(key, 0.0) or 0.0)


def element_attachments(
    facts: dict[str, Any],
    side: str,
    center: int,
    element: str,
    min_order: float = 0.0,
    max_order: float = 99.0,
) -> list[dict[str, Any]]:
    return [
        x
        for x in attachments(facts, side, center)
        if x.get("element") == element
        and min_order <= float(x.get("bond_order", 0.0)) <= max_order
    ]


def support_maps(
    facts: dict[str, Any],
    center: int | None,
) -> list[int]:
    """Mapped local support around one conserved center."""
    if center is None:
        return []
    maps = {int(center)}
    for side in ("reactant", "product"):
        for x in attachments(facts, side, int(center)):
            if x.get("map") is not None:
                maps.add(int(x["map"]))
    return sorted(maps)


def hydrogen_delta(facts: dict[str, Any], map_id: int) -> int:
    r = atom(facts, "reactant", map_id)
    p = atom(facts, "product", map_id)
    if not r or not p:
        return 0
    return int(p.get("total_h", 0)) - int(r.get("total_h", 0))
