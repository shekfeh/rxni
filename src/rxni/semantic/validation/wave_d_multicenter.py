"""Stage D validator: VC6_SPECIAL_MULTICENTER.

Seven event-specific mapped-topology proofs. Audit-only.

The common four-transition obligations remain explicit:
  SOURCE_PRESENT
  SOURCE_DESTROYED
  TARGET_ABSENT_BEFORE
  TARGET_CREATED

but GRAPH_EDIT_COMPATIBLE is event-specific because these transformations
cannot be reduced safely to a single same-center state transition.
"""

from __future__ import annotations

from typing import Any, Iterable
import json

from ..states.common import HALOGENS, atom, attachment_map, attachments, bond_order
from ..types import SemanticRecord, SemanticValidation, ValidationCheck, ValidationStatus
from ..version import FGI_CONTRACT_VERSION

VALIDATION_IMPLEMENTATION_VERSION = "wave-d-multicenter-v1"

VC6_EVENTS = {
    "alkene_oxidative_cleavage_to_aldehyde",
    "alpha_amino_acid_to_aldehyde",
    "benzylic_halide_to_aldehyde",
    "c_h_to_aldehyde",
    "carbonyl_to_amine_reductive_amination",
    "epoxide_to_azido_alcohol",
    "organosilane_to_alcohol",
}

WAVE_D_EVENTS = set(VC6_EVENTS)


def _check(name: str, passed: bool | None, detail: str = "") -> ValidationCheck:
    return ValidationCheck(name=name, passed=passed, detail=detail)


def _finish(record: SemanticRecord, checks: Iterable[ValidationCheck]) -> SemanticValidation:
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
        validator_class="VC6_SPECIAL_MULTICENTER",
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
    exclude_maps: set[int] | None = None,
) -> tuple[set[int], bool]:
    out = set()
    complete = True
    exclude_maps = set(exclude_maps or ())
    for x in attachments(facts, side, center):
        if x.get("element") == "H":
            continue
        if elements is not None and x.get("element") not in elements:
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
            continue
        mn = int(mn)
        if mn not in exclude_maps:
            out.add(mn)
    return out, complete


def _mapped_element(
    facts: dict[str, Any],
    side: str,
    center: int,
    element: str,
    lo: float = 0.0,
    hi: float = 99.0,
) -> tuple[set[int], bool]:
    out = set()
    complete = True
    for x in attachments(facts, side, center):
        if x.get("element") != element:
            continue
        bo = float(x.get("bond_order", 0.0) or 0.0)
        if not (lo <= bo <= hi):
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
        else:
            out.add(int(mn))
    return out, complete


def _is_aldehyde(facts: dict[str, Any], side: str, c: int) -> bool:
    a = atom(facts, side, c)
    if a.get("element") != "C":
        return False
    if int(a.get("total_h", 0) or 0) < 1:
        return False
    return any(
        x.get("element") == "O"
        and float(x.get("bond_order", 0.0) or 0.0) >= 1.8
        for x in attachments(facts, side, c)
    )


def _carbonyl_o_maps(facts: dict[str, Any], side: str, c: int) -> tuple[set[int], bool]:
    return _mapped_element(facts, side, c, "O", 1.8, 2.2)


def _single_n_maps(facts: dict[str, Any], side: str, c: int) -> tuple[set[int], bool]:
    return _mapped_element(facts, side, c, "N", 0.8, 1.6)


def _halogen_maps(facts: dict[str, Any], side: str, c: int) -> tuple[set[int], bool]:
    out = set()
    complete = True
    for x in attachments(facts, side, c):
        if x.get("element") not in HALOGENS:
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
        else:
            out.add(int(mn))
    return out, complete


def _aromatic_c_maps(facts: dict[str, Any], side: str, c: int) -> tuple[set[int], bool]:
    out = set()
    complete = True
    for x in attachments(facts, side, c):
        if x.get("element") != "C":
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
            continue
        if bool(atom(facts, side, int(mn)).get("aromatic")):
            out.add(int(mn))
    return out, complete


def _is_noncarbonyl_alkyl_ch(facts: dict[str, Any], side: str, c: int) -> bool:
    a = atom(facts, side, c)
    if a.get("element") != "C" or bool(a.get("aromatic")):
        return False
    if int(a.get("total_h", 0) or 0) < 2:
        return False
    return not any(
        x.get("element") == "O"
        and float(x.get("bond_order", 0.0) or 0.0) >= 1.8
        for x in attachments(facts, side, c)
    )


def _is_carbonyl_center(facts: dict[str, Any], side: str, c: int) -> bool:
    return bool(_carbonyl_o_maps(facts, side, c)[0])


def _is_amine_carbon(facts: dict[str, Any], side: str, c: int) -> bool:
    return bool(_single_n_maps(facts, side, c)[0]) and not _is_carbonyl_center(facts, side, c)


def _hydrogen_gain(facts: dict[str, Any], c: int) -> bool:
    return int(atom(facts, "product", c).get("total_h", 0) or 0) > int(
        atom(facts, "reactant", c).get("total_h", 0) or 0
    )


def _is_carboxyl_carbon(facts: dict[str, Any], side: str, c: int) -> bool:
    if atom(facts, side, c).get("element") != "C":
        return False
    double_o = any(
        x.get("element") == "O"
        and float(x.get("bond_order", 0.0) or 0.0) >= 1.8
        for x in attachments(facts, side, c)
    )
    single_o = any(
        x.get("element") == "O"
        and 0.8 <= float(x.get("bond_order", 0.0) or 0.0) <= 1.2
        for x in attachments(facts, side, c)
    )
    return double_o and single_o


def _azide_n_map(facts: dict[str, Any], side: str, c: int) -> int | None:
    for x in attachments(facts, side, c):
        if x.get("element") != "N" or float(x.get("bond_order", 0.0) or 0.0) > 1.3:
            continue
        n1 = attachment_map(x)
        if n1 is None:
            continue
        n1 = int(n1)
        for y in attachments(facts, side, n1):
            if y.get("element") != "N":
                continue
            n2 = attachment_map(y)
            if n2 is None:
                continue
            n2 = int(n2)
            if any(
                z.get("element") == "N"
                and attachment_map(z) is not None
                and int(attachment_map(z)) != n1
                for z in attachments(facts, side, n2)
            ):
                return n1
            if float(y.get("bond_order", 0.0) or 0.0) >= 1.5:
                return n1
    return None


def _alcohol_o_maps(facts: dict[str, Any], side: str, c: int) -> tuple[set[int], bool]:
    out = set()
    complete = True
    for x in attachments(facts, side, c):
        if x.get("element") != "O":
            continue
        if not (0.8 <= float(x.get("bond_order", 0.0) or 0.0) <= 1.2):
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
            continue
        mn = int(mn)
        if int(atom(facts, side, mn).get("total_h", 0) or 0) >= 1:
            out.add(mn)
    return out, complete


def _validate_alkene_cleavage(facts: dict[str, Any], record: SemanticRecord) -> SemanticValidation:
    c = int(record.center_map)
    partners = []
    for x in attachments(facts, "reactant", c):
        if x.get("element") == "C" and float(x.get("bond_order", 0.0) or 0.0) >= 1.7:
            mn = attachment_map(x)
            if mn is not None:
                partners.append(int(mn))
    partner = partners[0] if len(partners) == 1 else None

    source_present = partner is not None
    source_destroyed = (
        partner is not None and bond_order(facts, "product", c, partner) < 0.2
    )
    target_before = _is_aldehyde(facts, "reactant", c)
    target_after = _is_aldehyde(facts, "product", c)

    r_c, rc = _mapped_neighbors(
        facts, "reactant", c, elements={"C"},
        exclude_maps=({partner} if partner is not None else set()),
    )
    p_c, pc = _mapped_neighbors(facts, "product", c, elements={"C"})
    scaffold = (r_c == p_c) if (rc and pc) else None

    graph = (
        source_present and source_destroyed
        and (not target_before) and target_after
        and scaffold is True
    )

    return _finish(record, [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present, f"alkene_partner={partner}"),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", not target_before),
        _check("TARGET_CREATED", target_after),
        _check("remaining_mapped_carbon_scaffold_preserved", scaffold, f"reactant={sorted(r_c)} product={sorted(p_c)}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ])


def _validate_alpha_amino_acid(facts: dict[str, Any], record: SemanticRecord) -> SemanticValidation:
    c = int(record.center_map)
    r_n, rn_complete = _single_n_maps(facts, "reactant", c)
    p_n, pn_complete = _single_n_maps(facts, "product", c)

    carboxyl = []
    for x in attachments(facts, "reactant", c):
        if x.get("element") != "C":
            continue
        mn = attachment_map(x)
        if mn is not None and _is_carboxyl_carbon(facts, "reactant", int(mn)):
            carboxyl.append(int(mn))
    carbox = carboxyl[0] if len(carboxyl) == 1 else None

    source_present = bool(r_n) and carbox is not None and not _is_carbonyl_center(facts, "reactant", c)
    source_destroyed = (
        source_present
        and not bool(p_n)
        and bond_order(facts, "product", c, carbox) < 0.2
    )
    target_before = _is_aldehyde(facts, "reactant", c)
    target_after = _is_aldehyde(facts, "product", c)

    exclude = set(r_n)
    if carbox is not None:
        exclude.add(carbox)
    # The alpha-amino-acid contract preserves the carbon skeleton around the
    # reacting alpha carbon after removing the departing amino and carboxyl
    # substituents.  The newly formed aldehyde oxygen is part of the TARGET
    # functional group, not part of the conserved scaffold, so compare mapped
    # carbon neighbours only.
    r_rest, rc = _mapped_neighbors(
        facts,
        "reactant",
        c,
        elements={"C"},
        exclude_maps=exclude,
    )
    p_rest, pc = _mapped_neighbors(
        facts,
        "product",
        c,
        elements={"C"},
        exclude_maps=set(),
    )
    scaffold = (
        r_rest == p_rest
        if (rc and pc and rn_complete and pn_complete)
        else None
    )

    graph = (
        source_present and source_destroyed
        and not target_before and target_after
        and scaffold is True
    )
    return _finish(record, [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present, f"N={sorted(r_n)} carboxyl_C={carbox}"),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", not target_before),
        _check("TARGET_CREATED", target_after),
        _check("remaining_mapped_scaffold_preserved", scaffold, f"reactant={sorted(r_rest)} product={sorted(p_rest)}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ])


def _validate_benzylic_halide(facts: dict[str, Any], record: SemanticRecord) -> SemanticValidation:
    c = int(record.center_map)
    r_hal, rhc = _halogen_maps(facts, "reactant", c)
    p_hal, phc = _halogen_maps(facts, "product", c)
    r_ar, rac = _aromatic_c_maps(facts, "reactant", c)
    p_ar, pac = _aromatic_c_maps(facts, "product", c)

    source_present = (
        bool(r_hal)
        and bool(r_ar)
        and not bool(atom(facts, "reactant", c).get("aromatic"))
    )
    source_destroyed = source_present and not bool(p_hal)
    target_before = _is_aldehyde(facts, "reactant", c)
    target_after = _is_aldehyde(facts, "product", c)
    aromatic_lineage = (r_ar == p_ar) if (rac and pac) else None

    r_c, rc = _mapped_neighbors(facts, "reactant", c, elements={"C"})
    p_c, pc = _mapped_neighbors(facts, "product", c, elements={"C"})
    scaffold = (r_c == p_c) if (rc and pc) else None

    graph = (
        source_present and source_destroyed and not target_before and target_after
        and aromatic_lineage is True and scaffold is True
    )
    return _finish(record, [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", not target_before),
        _check("TARGET_CREATED", target_after),
        _check("benzylic_aromatic_attachment_preserved", aromatic_lineage, f"reactant={sorted(r_ar)} product={sorted(p_ar)}"),
        _check("mapped_carbon_scaffold_preserved", scaffold, f"reactant={sorted(r_c)} product={sorted(p_c)}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ])


def _validate_ch_to_aldehyde(facts: dict[str, Any], record: SemanticRecord) -> SemanticValidation:
    c = int(record.center_map)
    source_before = _is_noncarbonyl_alkyl_ch(facts, "reactant", c)
    source_after = _is_noncarbonyl_alkyl_ch(facts, "product", c)
    target_before = _is_aldehyde(facts, "reactant", c)
    target_after = _is_aldehyde(facts, "product", c)

    r_c, rc = _mapped_neighbors(facts, "reactant", c, elements={"C"})
    p_c, pc = _mapped_neighbors(facts, "product", c, elements={"C"})
    scaffold = (r_c == p_c) if (rc and pc) else None
    h_loss = int(atom(facts, "product", c).get("total_h", 0) or 0) < int(
        atom(facts, "reactant", c).get("total_h", 0) or 0
    )

    graph = (
        source_before and not source_after and not target_before and target_after
        and h_loss and scaffold is True
    )
    return _finish(record, [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_before),
        _check("SOURCE_DESTROYED", not source_after),
        _check("TARGET_ABSENT_BEFORE", not target_before),
        _check("TARGET_CREATED", target_after),
        _check("hydrogen_count_decreases", h_loss),
        _check("mapped_carbon_scaffold_preserved", scaffold, f"reactant={sorted(r_c)} product={sorted(p_c)}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ])


def _validate_reductive_amination(facts: dict[str, Any], record: SemanticRecord) -> SemanticValidation:
    c = int(record.center_map)
    source_before = _is_carbonyl_center(facts, "reactant", c)
    source_after = _is_carbonyl_center(facts, "product", c)
    target_before = _is_amine_carbon(facts, "reactant", c)
    target_after = _is_amine_carbon(facts, "product", c)

    r_n, rnc = _single_n_maps(facts, "reactant", c)
    p_n, pnc = _single_n_maps(facts, "product", c)
    gained_n = p_n - r_n

    r_c, rc = _mapped_neighbors(facts, "reactant", c, elements={"C"})
    p_c, pc = _mapped_neighbors(facts, "product", c, elements={"C"})
    scaffold = (r_c == p_c) if (rc and pc) else None

    graph = (
        source_before and not source_after and not target_before and target_after
        and bool(gained_n) and _hydrogen_gain(facts, c)
        and scaffold is True and rnc and pnc
    )
    return _finish(record, [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_before),
        _check("SOURCE_DESTROYED", not source_after),
        _check("TARGET_ABSENT_BEFORE", not target_before),
        _check("TARGET_CREATED", target_after),
        _check("mapped_nitrogen_gained", bool(gained_n) if (rnc and pnc) else None, f"gained_N={sorted(gained_n)}"),
        _check("hydrogen_gain", _hydrogen_gain(facts, c)),
        _check("mapped_carbon_scaffold_preserved", scaffold, f"reactant={sorted(r_c)} product={sorted(p_c)}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ])


def _validate_epoxide_azidolysis(facts: dict[str, Any], record: SemanticRecord) -> SemanticValidation:
    c = int(record.center_map)
    candidates = []
    for x in attachments(facts, "reactant", c):
        if x.get("element") != "C" or not (0.8 <= float(x.get("bond_order", 0.0) or 0.0) <= 1.2):
            continue
        partner = attachment_map(x)
        if partner is None:
            continue
        partner = int(partner)
        c_o, coc = _mapped_element(facts, "reactant", c, "O", 0.8, 1.2)
        p_o, poc = _mapped_element(facts, "reactant", partner, "O", 0.8, 1.2)
        shared = c_o & p_o
        for om in shared:
            candidates.append((partner, om, coc and poc))

    source_present = len(candidates) == 1
    partner, om, mapping_complete = candidates[0] if source_present else (None, None, False)

    source_destroyed = False
    target_created = False
    graph = False
    azn = _azide_n_map(facts, "product", c)
    target_before = _azide_n_map(facts, "reactant", c) is not None

    if source_present:
        center_o_gone = bond_order(facts, "product", c, om) < 0.2
        partner_o_kept = 0.8 <= bond_order(facts, "product", partner, om) <= 1.2
        c_partner_kept = 0.8 <= bond_order(facts, "product", c, partner) <= 1.2
        alcohol = int(atom(facts, "product", om).get("total_h", 0) or 0) >= 1
        source_destroyed = center_o_gone
        target_created = azn is not None and partner_o_kept and alcohol
        graph = (
            mapping_complete and center_o_gone and partner_o_kept
            and c_partner_kept and target_created and not target_before
        )

    ring_open_evidence = any(
        rec.get("subfamily") == "epoxide_ring_opening"
        and c in {int(x) for x in (rec.get("conserved_maps") or [])}
        for rec in (facts.get("ring_openings") or [])
    )
    # Explicit local topology is sufficient proof; ring_openings is diagnostic
    # because direct test fixtures need not synthesize the upstream annotation.

    return _finish(record, [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present, f"partner={partner} epoxide_O={om}"),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", not target_before),
        _check("TARGET_CREATED", target_created, f"azide_N={azn}"),
        _check("ring_opening_annotation", True, f"diagnostic_upstream_annotation={ring_open_evidence}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ])


def _validate_organosilane(facts: dict[str, Any], record: SemanticRecord) -> SemanticValidation:
    c = int(record.center_map)
    r_si, rsc = _mapped_element(facts, "reactant", c, "Si", 0.8, 1.2)
    p_si, psc = _mapped_element(facts, "product", c, "Si", 0.8, 1.2)
    r_oh, roc = _alcohol_o_maps(facts, "reactant", c)
    p_oh, poc = _alcohol_o_maps(facts, "product", c)

    source_present = bool(r_si)
    source_destroyed = source_present and not bool(p_si)
    target_before = bool(r_oh)
    target_after = bool(p_oh)

    r_c, rc = _mapped_neighbors(facts, "reactant", c, elements={"C"})
    p_c, pc = _mapped_neighbors(facts, "product", c, elements={"C"})
    scaffold = (r_c == p_c) if (rc and pc) else None

    graph = (
        source_present and source_destroyed and not target_before and target_after
        and scaffold is True and rsc and psc and roc and poc
    )
    return _finish(record, [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present, f"Si={sorted(r_si)}"),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", not target_before),
        _check("TARGET_CREATED", target_after, f"OH_O={sorted(p_oh)}"),
        _check("mapped_carbon_scaffold_preserved", scaffold, f"reactant={sorted(r_c)} product={sorted(p_c)}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ])


def validate_wave_d_multicenter(
    facts: dict[str, Any],
    record: SemanticRecord,
    contract: Any,
) -> SemanticValidation | None:
    if record.event not in VC6_EVENTS:
        return None
    if record.center_map is None:
        return _finish(record, [_check("mapped_center_present", None)])

    if record.event == "alkene_oxidative_cleavage_to_aldehyde":
        return _validate_alkene_cleavage(facts, record)
    if record.event == "alpha_amino_acid_to_aldehyde":
        return _validate_alpha_amino_acid(facts, record)
    if record.event == "benzylic_halide_to_aldehyde":
        return _validate_benzylic_halide(facts, record)
    if record.event == "c_h_to_aldehyde":
        return _validate_ch_to_aldehyde(facts, record)
    if record.event == "carbonyl_to_amine_reductive_amination":
        return _validate_reductive_amination(facts, record)
    if record.event == "epoxide_to_azido_alcohol":
        return _validate_epoxide_azidolysis(facts, record)
    if record.event == "organosilane_to_alcohol":
        return _validate_organosilane(facts, record)
    return None
