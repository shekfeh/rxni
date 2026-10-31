"""Priority-1 directional-FGI validators for Runtime v6.6.

Audit-only milestone 3B1 implementation.

These validators consume already-built native structural facts plus one raw
SemanticRecord hypothesis. They do not mutate the hypothesis, resolve semantic
ownership, derive redox, or alter the semantic profile.
"""

from __future__ import annotations

from typing import Any, Iterable

from ..states.carbon import carbon_state
from ..states.common import HALOGENS, atom, attachment_map, attachments, hydrogen_delta
from ..states.nitrogen import nitrogen_state, nitrogen_state_detailed
from ..types import SemanticRecord, SemanticValidation, ValidationCheck, ValidationStatus
from ..version import FGI_CONTRACT_VERSION

VALIDATION_IMPLEMENTATION_VERSION = "priority1-audit-v1"

TARGET_EVENTS = {
    "azide_to_amine",
    "nitro_to_amine",
    "alkyl_halide_to_hydrogen",
    # Same VC5 contract; no independent 509 positive/negative set yet.
    "aryl_halide_to_hydrogen",
}


def _check(name: str, passed: bool | None, detail: str = "") -> ValidationCheck:
    return ValidationCheck(name=name, passed=passed, detail=detail)


def _finish(
    record: SemanticRecord,
    validator_class: str,
    checks: Iterable[ValidationCheck],
) -> SemanticValidation:
    checks = tuple(checks)
    failed = tuple(x.name for x in checks if x.passed is False)
    unknown = tuple(x.name for x in checks if x.passed is None)
    if failed:
        status = ValidationStatus.REJECTED
        veto = failed
    elif unknown:
        status = ValidationStatus.AMBIGUOUS
        veto = unknown
    else:
        status = ValidationStatus.ACCEPTED
        veto = ()
    return SemanticValidation(
        label=record.event,
        status=status,
        validator_class=validator_class,
        checks=checks,
        veto_reasons=veto,
        contract_version=FGI_CONTRACT_VERSION,
    )


def _mapped_neighbors(
    facts: dict[str, Any],
    side: str,
    center: int,
    *,
    elements: set[str] | None = None,
    exclude_elements: set[str] | None = None,
) -> tuple[set[tuple[int, str, float]], bool]:
    """Mapped heavy-neighbor signature and whether heavy mapping is complete."""
    out: set[tuple[int, str, float]] = set()
    complete = True
    for att in attachments(facts, side, center):
        el = str(att.get("element") or "")
        if not el or el == "H":
            continue
        if elements is not None and el not in elements:
            continue
        if exclude_elements is not None and el in exclude_elements:
            continue
        mn = attachment_map(att)
        if mn is None:
            complete = False
            continue
        out.add((int(mn), el, round(float(att.get("bond_order", 0.0) or 0.0), 3)))
    return out, complete


def _mapped_carbon_neighbors(
    facts: dict[str, Any], side: str, center: int
) -> tuple[set[int], bool]:
    sig, complete = _mapped_neighbors(facts, side, center, elements={"C"})
    return {mn for mn, _, _ in sig}, complete


def _element_count(facts: dict[str, Any], side: str, center: int, element: str) -> int:
    return sum(1 for x in attachments(facts, side, center) if x.get("element") == element)


def _local_edit_conflict(
    facts: dict[str, Any],
    center: int,
    local_maps: set[int],
) -> tuple[bool | None, str]:
    """Detect mapped heavy-atom rewiring near a hydrodehalogenation anchor.

    A clean hydrodehalogenation may break C-X at the anchor. Any other mapped
    formed/broken/bond-order edit touching the anchor or its retained mapped
    heavy neighbors is incompatible with the strict VC5 contract.
    """
    saw_schema = False
    conflicts: list[str] = []

    for key in ("formed_bonds", "broken_bonds", "bond_order_changes"):
        rows = facts.get(key)
        if rows is None:
            continue
        saw_schema = True
        for row in rows:
            try:
                a = int(row["a"])
                b = int(row["b"])
            except Exception:
                continue
            if a not in local_maps and b not in local_maps:
                continue

            # Allow only explicit mapped C-X loss at the center.
            if key == "broken_bonds" and center in {a, b}:
                other = b if a == center else a
                if atom(facts, "reactant", other).get("element") in HALOGENS:
                    continue

            conflicts.append(f"{key}:{a}-{b}")

    if not saw_schema:
        return None, "mapped edit collections unavailable"
    return (not conflicts), (";".join(conflicts) if conflicts else "no incompatible local mapped edits")


def _validate_azide_to_amine(
    facts: dict[str, Any], record: SemanticRecord
) -> SemanticValidation:
    center = record.center_map
    if center is None:
        return _finish(
            record,
            "VC8_HETEROATOM_OXSTATE",
            [_check("mapped_center_present", None, "hypothesis has no center_map")],
        )

    center = int(center)
    r_atom = atom(facts, "reactant", center)
    p_atom = atom(facts, "product", center)
    r_c, r_c_complete = _mapped_carbon_neighbors(facts, "reactant", center)
    p_c, p_c_complete = _mapped_carbon_neighbors(facts, "product", center)
    product_detail = nitrogen_state_detailed(facts, "product", center)

    checks = [
        _check(
            "mapped_center_present",
            bool(r_atom) and bool(p_atom),
            f"reactant={bool(r_atom)} product={bool(p_atom)}",
        ),
        _check(
            "source_organic_azide",
            nitrogen_state(facts, "reactant", center) == "azide_n" and len(r_c) == 1,
            f"source_state={nitrogen_state(facts, 'reactant', center)} carbon_neighbors={sorted(r_c)}",
        ),
        _check(
            "carbon_anchor_mapping_complete",
            True if (r_c_complete and p_c_complete) else None,
            f"reactant_complete={r_c_complete} product_complete={p_c_complete}",
        ),
        _check(
            "original_cn_anchor_retained",
            len(r_c) == 1 and p_c == r_c,
            f"reactant={sorted(r_c)} product={sorted(p_c)}",
        ),
        _check(
            "n3_topology_destroyed",
            _element_count(facts, "reactant", center, "N") >= 1
            and _element_count(facts, "product", center, "N") == 0,
            (
                f"reactant_N_neighbors={_element_count(facts, 'reactant', center, 'N')} "
                f"product_N_neighbors={_element_count(facts, 'product', center, 'N')}"
            ),
        ),
        _check(
            "target_primary_amine",
            product_detail == "primary_amine",
            f"product_detailed_state={product_detail}",
        ),
        _check(
            "no_new_second_carbon_attachment",
            len(p_c) == 1,
            f"product_carbon_neighbors={sorted(p_c)}",
        ),
        _check(
            "hydrogen_gain",
            hydrogen_delta(facts, center) > 0,
            f"delta_H={hydrogen_delta(facts, center)}",
        ),
    ]
    return _finish(record, "VC8_HETEROATOM_OXSTATE", checks)


def _validate_nitro_to_amine(
    facts: dict[str, Any], record: SemanticRecord
) -> SemanticValidation:
    center = record.center_map
    if center is None:
        return _finish(
            record,
            "VC8_HETEROATOM_OXSTATE",
            [_check("mapped_center_present", None, "hypothesis has no center_map")],
        )

    center = int(center)
    r_atom = atom(facts, "reactant", center)
    p_atom = atom(facts, "product", center)
    r_c, r_c_complete = _mapped_carbon_neighbors(facts, "reactant", center)
    p_c, p_c_complete = _mapped_carbon_neighbors(facts, "product", center)
    product_detail = nitrogen_state_detailed(facts, "product", center)
    r_oxygen = _element_count(facts, "reactant", center, "O")
    p_oxygen = _element_count(facts, "product", center, "O")

    checks = [
        _check(
            "mapped_center_present",
            bool(r_atom) and bool(p_atom),
            f"reactant={bool(r_atom)} product={bool(p_atom)}",
        ),
        _check(
            "source_nitro",
            nitrogen_state(facts, "reactant", center) == "nitro" and len(r_c) == 1,
            f"source_state={nitrogen_state(facts, 'reactant', center)} carbon_neighbors={sorted(r_c)}",
        ),
        _check(
            "carbon_anchor_mapping_complete",
            True if (r_c_complete and p_c_complete) else None,
            f"reactant_complete={r_c_complete} product_complete={p_c_complete}",
        ),
        _check(
            "original_cn_anchor_retained",
            len(r_c) == 1 and p_c == r_c,
            f"reactant={sorted(r_c)} product={sorted(p_c)}",
        ),
        _check(
            "nitro_oxygen_removed",
            r_oxygen >= 2 and p_oxygen == 0,
            f"reactant_O_neighbors={r_oxygen} product_O_neighbors={p_oxygen}",
        ),
        _check(
            "target_primary_amine",
            product_detail == "primary_amine",
            f"product_detailed_state={product_detail}",
        ),
        _check(
            "no_new_carbon_attachment",
            len(p_c) == 1,
            f"product_carbon_neighbors={sorted(p_c)}",
        ),
        _check(
            "hydrogen_gain",
            hydrogen_delta(facts, center) > 0,
            f"delta_H={hydrogen_delta(facts, center)}",
        ),
    ]
    return _finish(record, "VC8_HETEROATOM_OXSTATE", checks)


def _validate_halide_to_hydrogen(
    facts: dict[str, Any], record: SemanticRecord
) -> SemanticValidation:
    center = record.center_map
    if center is None:
        return _finish(
            record,
            "VC5_HYDROGEN_REPLACEMENT",
            [_check("mapped_center_present", None, "hypothesis has no center_map")],
        )

    center = int(center)
    r_atom = atom(facts, "reactant", center)
    p_atom = atom(facts, "product", center)
    r_state = carbon_state(facts, "reactant", center)
    p_state = carbon_state(facts, "product", center)
    expected_source = "aryl_halide" if record.event.startswith("aryl_") else "alkyl_halide"

    r_sig, r_complete = _mapped_neighbors(
        facts, "reactant", center, exclude_elements=set(HALOGENS)
    )
    p_sig, p_complete = _mapped_neighbors(
        facts, "product", center, exclude_elements=set(HALOGENS)
    )
    local_maps = {center} | {x[0] for x in r_sig} | {x[0] for x in p_sig}
    edit_ok, edit_detail = _local_edit_conflict(facts, center, local_maps)

    r_hal = sum(1 for x in attachments(facts, "reactant", center) if x.get("element") in HALOGENS)
    p_hal = sum(1 for x in attachments(facts, "product", center) if x.get("element") in HALOGENS)

    aromatic_ok = True
    if expected_source == "aryl_halide":
        aromatic_ok = bool(r_atom.get("aromatic")) and bool(p_atom.get("aromatic"))
    else:
        aromatic_ok = not bool(r_atom.get("aromatic")) and not bool(p_atom.get("aromatic"))

    checks = [
        _check(
            "mapped_center_present",
            bool(r_atom) and bool(p_atom),
            f"reactant={bool(r_atom)} product={bool(p_atom)}",
        ),
        _check(
            "source_halide_state",
            r_state == expected_source and r_hal >= 1,
            f"source_state={r_state} expected={expected_source} reactant_halogen_neighbors={r_hal}",
        ),
        _check(
            "halide_destroyed",
            p_hal == 0,
            f"product_halogen_neighbors={p_hal}",
        ),
        _check(
            "target_hydrocarbon_state",
            p_state == "",
            f"product_carbon_state={p_state!r}",
        ),
        _check(
            "heavy_neighbor_mapping_complete",
            True if (r_complete and p_complete) else None,
            f"reactant_complete={r_complete} product_complete={p_complete}",
        ),
        _check(
            "mapped_heavy_neighborhood_preserved",
            r_sig == p_sig,
            f"reactant={sorted(r_sig)} product={sorted(p_sig)}",
        ),
        _check(
            "aromaticity_class_preserved",
            aromatic_ok,
            f"reactant_aromatic={bool(r_atom.get('aromatic'))} product_aromatic={bool(p_atom.get('aromatic'))}",
        ),
        _check(
            "hydrogen_gain",
            hydrogen_delta(facts, center) > 0,
            f"delta_H={hydrogen_delta(facts, center)}",
        ),
        _check(
            "no_incompatible_local_mapped_edits",
            edit_ok,
            edit_detail,
        ),
    ]
    return _finish(record, "VC5_HYDROGEN_REPLACEMENT", checks)


def validate_priority1_directional(
    facts: dict[str, Any], record: SemanticRecord
) -> SemanticValidation | None:
    """Validate one supported priority-1 hypothesis, else return None."""
    if record.event == "azide_to_amine":
        return _validate_azide_to_amine(facts, record)
    if record.event == "nitro_to_amine":
        return _validate_nitro_to_amine(facts, record)
    if record.event in {"alkyl_halide_to_hydrogen", "aryl_halide_to_hydrogen"}:
        return _validate_halide_to_hydrogen(facts, record)
    return None
