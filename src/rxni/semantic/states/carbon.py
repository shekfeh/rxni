"""Carbon-center state classification.

`carbon_state()` is the parity-preserving v6.6 port of the trusted v6.5.6
carbon functional-state classifier. It identifies local state only; it never
infers a reaction.
"""

from __future__ import annotations

from typing import Any

from .common import (
    HALOGENS,
    atom,
    attachment_map,
    attachments,
    element_attachments,
)


def anhydride_partner_centers(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> list[int]:
    """Mapped acyl centers joined through C(=O)-O-C(=O)."""
    out: list[int] = []
    for oatt in element_attachments(facts, side, center, "O", 0.8, 1.2):
        om = attachment_map(oatt)
        if om is None:
            continue
        for catt in attachments(facts, side, om):
            cm = attachment_map(catt)
            if cm is None or cm == center or catt.get("element") != "C":
                continue
            if not (0.8 <= float(catt.get("bond_order", 0.0)) <= 1.2):
                continue
            if any(
                x.get("element") == "O"
                and float(x.get("bond_order", 0.0)) >= 1.8
                for x in attachments(facts, side, cm)
            ):
                out.append(int(cm))
    return sorted(set(out))


def carbon_state(
    facts: dict[str, Any],
    side: str,
    map_id: int,
) -> str:
    center = atom(facts, side, map_id)
    if center.get("element") != "C":
        return ""

    atts = attachments(facts, side, map_id)

    # C#N
    if any(
        x.get("element") == "N"
        and float(x.get("bond_order", 0.0)) >= 2.8
        for x in atts
    ):
        return "nitrile"

    carbonyl_o = [
        x for x in atts
        if x.get("element") == "O"
        and float(x.get("bond_order", 0.0)) >= 1.8
    ]
    if carbonyl_o:
        single_o = [
            x for x in atts
            if x.get("element") == "O"
            and 0.8 <= float(x.get("bond_order", 0.0)) <= 1.2
        ]
        single_n = [
            x for x in atts
            if x.get("element") == "N"
            and float(x.get("bond_order", 0.0)) <= 1.2
        ]
        single_x = [
            x for x in atts
            if x.get("element") in HALOGENS
            and float(x.get("bond_order", 0.0)) <= 1.2
        ]

        if single_n:
            return "amide"

        if single_o and single_x:
            return "haloformate"

        if single_o:
            if anhydride_partner_centers(facts, side, map_id):
                return "anhydride"
            if any(int(x.get("degree", 0)) >= 2 for x in single_o):
                return "ester"
            return "carboxylic_acid"

        if single_x:
            return "acyl_halide"

        h = int(center.get("total_h", 0))
        return "aldehyde" if h >= 1 else "ketone"

    if any(
        x.get("element") in HALOGENS
        and float(x.get("bond_order", 0.0)) <= 1.2
        for x in atts
    ):
        return "aryl_halide" if center.get("aromatic") else "alkyl_halide"

    if center.get("aromatic"):
        return ""

    oh = [
        x for x in atts
        if x.get("element") == "O"
        and 0.8 <= float(x.get("bond_order", 0.0)) <= 1.2
        and int(x.get("degree", 0)) <= 1
    ]
    if oh:
        h = int(center.get("total_h", 0))
        if h >= 2:
            return "primary_alcohol"
        if h == 1:
            return "secondary_alcohol"
        return "tertiary_alcohol"

    single_n = [
        x for x in atts
        if x.get("element") == "N"
        and float(x.get("bond_order", 0.0)) <= 1.2
    ]
    if single_n and int(center.get("total_h", 0)) >= 1:
        return "amine_carbon"

    return ""
