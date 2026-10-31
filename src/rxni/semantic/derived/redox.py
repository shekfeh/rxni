"""Canonical redox-axis derivation.

Stage E invariant:
redox may be derived only from canonical ACCEPTED directional conversions.
Rejected/ambiguous/raw hypotheses are never consulted.
"""

from __future__ import annotations

from collections import OrderedDict

from ..types import SemanticRecord


def derive_redox_axis(
    directional: list[SemanticRecord],
) -> list[SemanticRecord]:
    """Derive the active redox axis from accepted canonical directionals only."""
    rows: list[SemanticRecord] = []

    for rec in directional:
        direction = (rec.redox_direction or "").strip()
        if direction not in {"oxidation", "reduction"}:
            continue
        rows.append(
            SemanticRecord(
                event=direction,
                role="axis_event",
                subtype=rec.event,
                center_map=rec.center_map,
                support_maps=rec.support_maps,
                mapping_scope="derived_from_accepted_directional_conversion",
            )
        )

    uniq = OrderedDict()
    for rec in rows:
        key = (rec.event, rec.subtype, rec.center_map, rec.support_maps)
        uniq.setdefault(key, rec)
    return list(uniq.values())
