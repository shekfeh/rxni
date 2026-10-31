"""Nitrogen-state classification.

`nitrogen_state()` preserves v6.5.6 behavior exactly for migration parity.

`nitrogen_state_detailed()` is an additional non-authoritative descriptor for
the future native validator. It does not alter semantic labels in Milestone 2.
"""

from __future__ import annotations

from typing import Any

from .common import atom, attachment_map, attachments


def is_azide_n(
    facts: dict[str, Any],
    side: str,
    map_id: int,
) -> bool:
    if atom(facts, side, map_id).get("element") != "N":
        return False
    for att in attachments(facts, side, map_id):
        if att.get("element") != "N":
            continue
        if any(str(x[0]) == "N" for x in (att.get("neighbor_shell") or [])):
            return True
    return False


def nitrogen_state(
    facts: dict[str, Any],
    side: str,
    map_id: int,
) -> str:
    """Parity state used to compare with v6.5.6."""
    center = atom(facts, side, map_id)
    if center.get("element") != "N":
        return ""

    atts = attachments(facts, side, map_id)

    if is_azide_n(facts, side, map_id):
        return "azide_n"

    if any(
        x.get("element") == "C"
        and float(x.get("bond_order", 0.0)) >= 2.8
        for x in atts
    ):
        return "nitrile_n"

    oxy = [x for x in atts if x.get("element") == "O"]
    if len(oxy) >= 2 and sum(float(x.get("bond_order", 0.0)) for x in oxy) >= 2.5:
        return "nitro"

    if len(oxy) == 1 and int(center.get("charge", 0)) > 0:
        return "aza_n_oxide" if center.get("aromatic") else "amine_n_oxide"

    if not oxy and int(center.get("charge", 0)) <= 0:
        return "aza_nitrogen" if center.get("aromatic") else "amine_nitrogen"

    return ""


def _carbonyl_carbon_neighbor(
    facts: dict[str, Any],
    side: str,
    n_map: int,
) -> bool:
    for att in attachments(facts, side, n_map):
        if att.get("element") != "C":
            continue
        cm = attachment_map(att)
        if cm is None:
            continue
        if any(
            x.get("element") == "O"
            and float(x.get("bond_order", 0.0)) >= 1.8
            for x in attachments(facts, side, cm)
        ):
            return True
    return False


def nitrogen_state_detailed(
    facts: dict[str, Any],
    side: str,
    map_id: int,
) -> str:
    """Higher-resolution nitrogen descriptor for future validation.

    This is intentionally *not* used by reaction semantics in Milestone 2.
    """
    center = atom(facts, side, map_id)
    if center.get("element") != "N":
        return ""

    coarse = nitrogen_state(facts, side, map_id)
    if coarse in {"azide_n", "nitrile_n", "nitro", "aza_n_oxide", "amine_n_oxide"}:
        return coarse

    atts = attachments(facts, side, map_id)

    if center.get("aromatic"):
        # Aromatic [nH] versus pyridine-like aromatic N.
        if int(center.get("total_h", 0)) > 0:
            return "pyrrolic_n"
        return "aromatic_n"

    if _carbonyl_carbon_neighbor(facts, side, map_id):
        return "amide_n"

    if any(
        x.get("element") == "C"
        and 1.8 <= float(x.get("bond_order", 0.0)) <= 2.2
        for x in atts
    ):
        return "imine_n"

    carbon_neighbors = sum(
        1 for x in atts
        if x.get("element") == "C"
        and float(x.get("bond_order", 0.0)) <= 1.2
    )
    h = int(center.get("total_h", 0))

    if carbon_neighbors <= 1 and h >= 1:
        return "primary_amine"
    if carbon_neighbors == 2:
        return "secondary_amine"
    if carbon_neighbors >= 3:
        return "tertiary_amine"

    return coarse
