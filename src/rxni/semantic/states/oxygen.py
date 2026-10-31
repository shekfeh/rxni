"""Oxygen local descriptors.

There was no single authoritative v6.5.6 oxygen-state classifier corresponding
to the carbon/nitrogen/sulfur helpers. This module therefore provides only
descriptive states in Milestone 2 and is not used for parity-sensitive labels.
"""

from __future__ import annotations

from typing import Any

from .common import atom, attachments


def oxygen_state(
    facts: dict[str, Any],
    side: str,
    map_id: int,
) -> str:
    center = atom(facts, side, map_id)
    if center.get("element") != "O":
        return ""

    atts = attachments(facts, side, map_id)
    heavy = [x for x in atts if x.get("element") != "H"]

    if any(float(x.get("bond_order", 0.0)) >= 1.8 for x in heavy):
        return "oxo_oxygen"

    if int(center.get("charge", 0)) < 0:
        return "anionic_oxygen"

    degree = int(center.get("degree", len(heavy)))
    h = int(center.get("total_h", 0))

    if h > 0 or degree <= 1:
        return "hydroxyl_like_oxygen"

    if degree >= 2:
        return "ether_like_oxygen"

    return "oxygen"
