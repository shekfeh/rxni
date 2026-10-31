"""Stage C validators: VC3 + VC4 + VC5 + VC7.

Audit-only.

VC3: mapped C-C pair bond-order / oxygenation proof.
VC4: mapped attachment replacement at the same carbon center.
VC5: mapped substituent loss with hydrogen gain and scaffold conservation.
VC7: same-class acyl center retained while partner identity changes.

No canonical semantic behavior is changed here.
"""

from __future__ import annotations

from typing import Any, Iterable
import json

from ..states.common import HALOGENS, atom, attachment_map, attachments, bond_order
from ..types import SemanticRecord, SemanticValidation, ValidationCheck, ValidationStatus
from ..version import FGI_CONTRACT_VERSION

VALIDATION_IMPLEMENTATION_VERSION = "wave-c-bond-attachment-v1"

VC3_EVENTS = {
    "alkane_to_alkene",
    "alkene_to_alcohol",
    "alkene_to_alkane",
    "alkene_to_diol",
    "alkene_to_epoxide",
    "alkyne_to_alkane",
    "alkyne_to_alkene",
}

VC4_EVENTS = {
    "alcohol_to_azide",
    "alcohol_to_halide",
    "alkyl_halide_to_azide",
    "aryl_halide_to_azide",
    "halide_to_alcohol",
    "heteroaryl_halide_to_oxo",
    "sulfonate_to_azide",
}

VC5_EVENTS = {
    "alkyl_halide_to_hydrogen",
    "aryl_halide_to_hydrogen",
    "aryl_amine_to_hydrogen",
    "carbon_oxygen_to_hydrogen_deoxygenation",
}

VC7_EVENTS = {
    "ethyl_ester_to_methyl_ester",
    "transamidation",
    "transesterification",
}

WAVE_C_EVENTS = VC3_EVENTS | VC4_EVENTS | VC5_EVENTS | VC7_EVENTS


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


def _mapped_c_neighbors(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], bool]:
    out = set()
    complete = True
    for x in attachments(facts, side, center):
        if x.get("element") != "C":
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
        else:
            out.add(int(mn))
    return out, complete


def _scaffold_preserved(
    facts: dict[str, Any],
    center: int,
) -> tuple[bool | None, str]:
    r, rc = _mapped_c_neighbors(facts, "reactant", center)
    p, pc = _mapped_c_neighbors(facts, "product", center)
    if not (rc and pc):
        return None, f"reactant={sorted(r)} product={sorted(p)} mapping_incomplete"
    return r == p, f"reactant={sorted(r)} product={sorted(p)}"


def _mapped_element_neighbors(
    facts: dict[str, Any],
    side: str,
    center: int,
    element: str,
    *,
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


def _halogen_maps(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], bool]:
    out = set()
    complete = True
    for x in attachments(facts, side, center):
        if x.get("element") not in HALOGENS:
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
        else:
            out.add(int(mn))
    return out, complete


def _single_o_maps(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], bool]:
    return _mapped_element_neighbors(
        facts, side, center, "O", lo=0.8, hi=1.2
    )


def _single_n_maps(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], bool]:
    return _mapped_element_neighbors(
        facts, side, center, "N", lo=0.8, hi=1.6
    )


def _double_o_maps(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], bool]:
    return _mapped_element_neighbors(
        facts, side, center, "O", lo=1.8, hi=2.2
    )


def _is_alcohol_center(facts: dict[str, Any], side: str, center: int) -> bool:
    for x in attachments(facts, side, center):
        if x.get("element") != "O":
            continue
        if not (0.8 <= float(x.get("bond_order", 0.0) or 0.0) <= 1.2):
            continue
        mn = attachment_map(x)
        if mn is None:
            # One-sided mapping may be valid chemically but is insufficient
            # for strict local proof.
            continue
        oa = atom(facts, side, int(mn))
        if int(oa.get("total_h", 0) or 0) >= 1:
            return True
    return False


def _is_halide_center(facts: dict[str, Any], side: str, center: int) -> bool:
    return bool(_halogen_maps(facts, side, center)[0])


def _azide_n_map(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> int | None:
    for x in attachments(facts, side, center):
        if x.get("element") != "N":
            continue
        if float(x.get("bond_order", 0.0) or 0.0) > 1.3:
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
            # Minimal mapped organic azide representation C-N-N is accepted
            # if the second N has a multiple bond to N or an additional N
            # may be represented through the n1->n2 bond order.
            if float(y.get("bond_order", 0.0) or 0.0) >= 1.5:
                return n1
    return None


def _is_sulfonate_center(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:
    for ox in attachments(facts, side, center):
        if ox.get("element") != "O":
            continue
        if not (0.8 <= float(ox.get("bond_order", 0.0) or 0.0) <= 1.2):
            continue
        om = attachment_map(ox)
        if om is None:
            continue
        om = int(om)
        for sx in attachments(facts, side, om):
            if sx.get("element") != "S":
                continue
            sm = attachment_map(sx)
            if sm is None:
                continue
            sm = int(sm)
            oxo = sum(
                1 for q in attachments(facts, side, sm)
                if q.get("element") == "O"
                and float(q.get("bond_order", 0.0) or 0.0) >= 1.8
            )
            if oxo >= 2:
                return True
    return False


def _is_oxo_center(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:
    return bool(_double_o_maps(facts, side, center)[0])


def _pair_candidates(
    facts: dict[str, Any],
    side: str,
    center: int,
    lo: float,
    hi: float,
) -> list[int]:
    out = []
    for x in attachments(facts, side, center):
        if x.get("element") != "C":
            continue
        mn = attachment_map(x)
        if mn is None:
            continue
        bo = float(x.get("bond_order", 0.0) or 0.0)
        if lo <= bo <= hi:
            out.append(int(mn))
    return sorted(set(out))


def _pair_transition(
    facts: dict[str, Any],
    center: int,
    rlo: float,
    rhi: float,
    plo: float,
    phi: float,
) -> tuple[int | None, bool]:
    for partner in _pair_candidates(facts, "reactant", center, rlo, rhi):
        po = bond_order(facts, "product", center, partner)
        if plo <= po <= phi:
            return partner, True
    return None, False


def _gained_oh(
    facts: dict[str, Any],
    center: int,
) -> list[int]:
    out = []
    r_o, _ = _single_o_maps(facts, "reactant", center)
    p_o, _ = _single_o_maps(facts, "product", center)
    for om in sorted(p_o - r_o):
        if int(atom(facts, "product", om).get("total_h", 0) or 0) >= 1:
            out.append(om)
    return out


def _shared_epoxide_o(
    facts: dict[str, Any],
    a: int,
    b: int,
    side: str,
) -> set[int]:
    oa, _ = _single_o_maps(facts, side, a)
    ob, _ = _single_o_maps(facts, side, b)
    return oa & ob


def _validate_vc3(
    facts: dict[str, Any],
    record: SemanticRecord,
) -> SemanticValidation:
    c = record.center_map
    if c is None:
        return _finish(record, "VC3_PI_PAIR", [_check("mapped_center_present", None)])
    c = int(c)

    event = record.event
    source_present = target_absent = source_destroyed = target_created = graph = False
    detail = {}

    if event == "alkane_to_alkene":
        partner, ok = _pair_transition(facts, c, 0.8, 1.2, 1.7, 2.2)
        source_present = partner is not None or bool(_pair_candidates(facts, "reactant", c, 0.8, 1.2))
        target_absent = all(bond_order(facts, "reactant", c, p) < 1.7 for p in _pair_candidates(facts, "reactant", c, 0.8, 1.2))
        source_destroyed = ok
        target_created = ok
        graph = ok
        detail = {"partner": partner}
    elif event == "alkene_to_alkane":
        partner, ok = _pair_transition(facts, c, 1.7, 2.2, 0.8, 1.2)
        source_present = bool(_pair_candidates(facts, "reactant", c, 1.7, 2.2))
        target_absent = True
        source_destroyed = ok
        target_created = ok
        graph = ok
        detail = {"partner": partner}
    elif event == "alkyne_to_alkane":
        partner, ok = _pair_transition(facts, c, 2.8, 3.2, 0.8, 1.2)
        source_present = bool(_pair_candidates(facts, "reactant", c, 2.8, 3.2))
        target_absent = True
        source_destroyed = ok
        target_created = ok
        graph = ok
        detail = {"partner": partner}
    elif event == "alkyne_to_alkene":
        partner, ok = _pair_transition(facts, c, 2.8, 3.2, 1.7, 2.2)
        source_present = bool(_pair_candidates(facts, "reactant", c, 2.8, 3.2))
        target_absent = True
        source_destroyed = ok
        target_created = ok
        graph = ok
        detail = {"partner": partner}
    elif event in {"alkene_to_alcohol", "alkene_to_diol", "alkene_to_epoxide"}:
        partners = _pair_candidates(facts, "reactant", c, 1.7, 2.2)
        source_present = bool(partners)
        partner = None
        for p in partners:
            if 0.8 <= bond_order(facts, "product", c, p) <= 1.2:
                partner = p
                break
        source_destroyed = partner is not None

        if event == "alkene_to_alcohol":
            target_absent = not _is_alcohol_center(facts, "reactant", c)
            target_created = bool(_gained_oh(facts, c))
            graph = partner is not None and target_created
        elif event == "alkene_to_diol":
            target_absent = True
            target_created = (
                partner is not None
                and bool(_gained_oh(facts, c))
                and bool(_gained_oh(facts, partner))
            )
            graph = target_created
        else:
            before = _shared_epoxide_o(facts, c, partner, "reactant") if partner is not None else set()
            after = _shared_epoxide_o(facts, c, partner, "product") if partner is not None else set()
            target_absent = not before
            target_created = bool(after - before)
            graph = partner is not None and target_created
        detail = {"partner": partner}

    checks = [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", target_absent),
        _check("TARGET_CREATED", target_created),
        _check("GRAPH_EDIT_COMPATIBLE", graph, json.dumps(detail, sort_keys=True)),
    ]
    return _finish(record, "VC3_PI_PAIR", checks)


def _validate_vc4(
    facts: dict[str, Any],
    record: SemanticRecord,
) -> SemanticValidation:
    c = record.center_map
    if c is None:
        return _finish(record, "VC4_ATTACHMENT_REPLACEMENT", [_check("mapped_center_present", None)])
    c = int(c)
    event = record.event

    r_alc = _is_alcohol_center(facts, "reactant", c)
    p_alc = _is_alcohol_center(facts, "product", c)
    r_hal = _is_halide_center(facts, "reactant", c)
    p_hal = _is_halide_center(facts, "product", c)
    r_az = _azide_n_map(facts, "reactant", c) is not None
    p_az = _azide_n_map(facts, "product", c) is not None
    r_sulf = _is_sulfonate_center(facts, "reactant", c)
    p_sulf = _is_sulfonate_center(facts, "product", c)
    r_oxo = _is_oxo_center(facts, "reactant", c)
    p_oxo = _is_oxo_center(facts, "product", c)
    scaffold, scaffold_detail = _scaffold_preserved(facts, c)

    source_present = source_destroyed = target_absent = target_created = False
    guard = True
    if event == "alcohol_to_azide":
        source_present, source_destroyed = r_alc, not p_alc
        target_absent, target_created = not r_az, p_az
    elif event == "alcohol_to_halide":
        source_present, source_destroyed = r_alc, not p_alc
        target_absent, target_created = not r_hal, p_hal
    elif event == "halide_to_alcohol":
        source_present, source_destroyed = r_hal, not p_hal
        target_absent, target_created = not r_alc, p_alc
    elif event == "alkyl_halide_to_azide":
        source_present, source_destroyed = r_hal, not p_hal
        target_absent, target_created = not r_az, p_az
        guard = not bool(atom(facts, "reactant", c).get("aromatic"))
    elif event == "aryl_halide_to_azide":
        source_present, source_destroyed = r_hal, not p_hal
        target_absent, target_created = not r_az, p_az
        guard = bool(atom(facts, "reactant", c).get("aromatic"))
    elif event == "sulfonate_to_azide":
        source_present, source_destroyed = r_sulf, not p_sulf
        target_absent, target_created = not r_az, p_az
    elif event == "heteroaryl_halide_to_oxo":
        source_present, source_destroyed = r_hal, not p_hal
        target_absent, target_created = not r_oxo, p_oxo
        guard = bool(atom(facts, "reactant", c).get("aromatic"))

    graph = (
        source_present
        and source_destroyed
        and target_absent
        and target_created
        and guard
        and scaffold is True
    )
    checks = [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", target_absent),
        _check("TARGET_CREATED", target_created),
        _check("event_specific_guard", guard),
        _check("mapped_carbon_scaffold_preserved", scaffold, scaffold_detail),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ]
    return _finish(record, "VC4_ATTACHMENT_REPLACEMENT", checks)


def _hydrogen_gain(facts: dict[str, Any], center: int) -> bool:
    rh = int(atom(facts, "reactant", center).get("total_h", 0) or 0)
    ph = int(atom(facts, "product", center).get("total_h", 0) or 0)
    return ph > rh


def _mapped_heavy_non_source(
    facts: dict[str, Any],
    side: str,
    center: int,
    excluded_elements: set[str],
) -> tuple[set[int], bool]:
    out = set()
    complete = True
    for x in attachments(facts, side, center):
        el = str(x.get("element") or "")
        if el == "H" or el in excluded_elements:
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
        else:
            out.add(int(mn))
    return out, complete


def _validate_vc5(
    facts: dict[str, Any],
    record: SemanticRecord,
) -> SemanticValidation:
    c = record.center_map
    if c is None:
        return _finish(record, "VC5_HYDROGEN_REPLACEMENT", [_check("mapped_center_present", None)])
    c = int(c)
    event = record.event

    r_hal = _is_halide_center(facts, "reactant", c)
    p_hal = _is_halide_center(facts, "product", c)
    r_n, _ = _single_n_maps(facts, "reactant", c)
    p_n, _ = _single_n_maps(facts, "product", c)
    r_o, _ = _single_o_maps(facts, "reactant", c)
    p_o, _ = _single_o_maps(facts, "product", c)

    aromatic_required = event in {"aryl_halide_to_hydrogen", "aryl_amine_to_hydrogen"}
    aromatic_guard = (
        bool(atom(facts, "reactant", c).get("aromatic"))
        if aromatic_required
        else True
    )
    if event == "alkyl_halide_to_hydrogen":
        aromatic_guard = not bool(atom(facts, "reactant", c).get("aromatic"))

    if event in {"alkyl_halide_to_hydrogen", "aryl_halide_to_hydrogen"}:
        source_present = r_hal
        source_destroyed = not p_hal
        excluded = set(HALOGENS)
    elif event == "aryl_amine_to_hydrogen":
        source_present = bool(r_n)
        source_destroyed = not bool(p_n)
        excluded = {"N"}
    else:
        source_present = bool(r_o)
        source_destroyed = not bool(p_o)
        excluded = {"O"}

    r_skel, rc = _mapped_heavy_non_source(facts, "reactant", c, excluded)
    p_skel, pc = _mapped_heavy_non_source(facts, "product", c, excluded)
    scaffold = (r_skel == p_skel) if (rc and pc) else None

    target_absent = source_present
    target_created = source_destroyed and _hydrogen_gain(facts, c)
    graph = (
        source_present
        and source_destroyed
        and target_created
        and aromatic_guard
        and scaffold is True
    )

    checks = [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("SOURCE_PRESENT", source_present),
        _check("SOURCE_DESTROYED", source_destroyed),
        _check("TARGET_ABSENT_BEFORE", target_absent),
        _check("TARGET_CREATED", target_created, "implicit hydrogen gain"),
        _check("event_specific_guard", aromatic_guard),
        _check("mapped_scaffold_preserved", scaffold, f"reactant={sorted(r_skel)} product={sorted(p_skel)}"),
        _check("GRAPH_EDIT_COMPATIBLE", graph),
    ]
    return _finish(record, "VC5_HYDROGEN_REPLACEMENT", checks)


def _is_ester_center(facts: dict[str, Any], side: str, c: int) -> bool:
    return bool(_double_o_maps(facts, side, c)[0]) and bool(_single_o_maps(facts, side, c)[0])


def _is_amide_center(facts: dict[str, Any], side: str, c: int) -> bool:
    return bool(_double_o_maps(facts, side, c)[0]) and bool(_single_n_maps(facts, side, c)[0])


def _carbonyl_lineage(facts: dict[str, Any], c: int) -> bool | None:
    r, rc = _double_o_maps(facts, "reactant", c)
    p, pc = _double_o_maps(facts, "product", c)
    if not (rc and pc):
        return None
    return bool(r & p)


def _validate_vc7(
    facts: dict[str, Any],
    record: SemanticRecord,
) -> SemanticValidation:
    c = record.center_map
    if c is None:
        return _finish(record, "VC7_SAME_CLASS_EXCHANGE", [_check("mapped_center_present", None)])
    c = int(c)
    event = record.event
    scaffold, scaffold_detail = _scaffold_preserved(facts, c)
    lineage = _carbonyl_lineage(facts, c)

    checks = [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
    ]

    if event == "transesterification":
        r_class = _is_ester_center(facts, "reactant", c)
        p_class = _is_ester_center(facts, "product", c)
        r_partner, rc = _single_o_maps(facts, "reactant", c)
        p_partner, pc = _single_o_maps(facts, "product", c)
        exchanged = (r_partner != p_partner) if (rc and pc) else None
        checks += [
            _check("SAME_CLASS_CENTER_RETAINED", r_class and p_class),
            _check("SOURCE_PARTNER_PRESENT", bool(r_partner)),
            _check("SOURCE_PARTNER_DESTROYED", None if exchanged is None else exchanged),
            _check("NEW_PARTNER_CREATED", None if exchanged is None else exchanged),
            _check("PARTNER_IDENTITY_EXCHANGED", exchanged, f"reactant_O={sorted(r_partner)} product_O={sorted(p_partner)}"),
        ]
    elif event == "transamidation":
        r_class = _is_amide_center(facts, "reactant", c)
        p_class = _is_amide_center(facts, "product", c)
        r_partner, rc = _single_n_maps(facts, "reactant", c)
        p_partner, pc = _single_n_maps(facts, "product", c)
        exchanged = (r_partner != p_partner) if (rc and pc) else None
        checks += [
            _check("SAME_CLASS_CENTER_RETAINED", r_class and p_class),
            _check("SOURCE_PARTNER_PRESENT", bool(r_partner)),
            _check("SOURCE_PARTNER_DESTROYED", None if exchanged is None else exchanged),
            _check("NEW_PARTNER_CREATED", None if exchanged is None else exchanged),
            _check("PARTNER_IDENTITY_EXCHANGED", exchanged, f"reactant_N={sorted(r_partner)} product_N={sorted(p_partner)}"),
        ]
    else:
        # Same ester O and alpha carbon retained; alpha carbon loses its distal
        # carbon and gains H: -O-CH2-CH3 -> -O-CH3.
        r_class = _is_ester_center(facts, "reactant", c)
        p_class = _is_ester_center(facts, "product", c)
        r_o, ro_complete = _single_o_maps(facts, "reactant", c)
        p_o, po_complete = _single_o_maps(facts, "product", c)
        common_o = r_o & p_o
        ok = False
        detail = {}
        for om in sorted(common_o):
            r_c, rc = _mapped_element_neighbors(facts, "reactant", om, "C", lo=0.8, hi=1.2)
            p_c, pc = _mapped_element_neighbors(facts, "product", om, "C", lo=0.8, hi=1.2)
            alpha = sorted((r_c & p_c) - {c})
            for ac in alpha:
                r_dist, rd = _mapped_c_neighbors(facts, "reactant", ac)
                p_dist, pd = _mapped_c_neighbors(facts, "product", ac)
                r_dist.discard(c)
                p_dist.discard(c)
                hgain = _hydrogen_gain(facts, ac)
                if rd and pd and r_dist and not p_dist and hgain:
                    ok = True
                    detail = {"ester_O": om, "alpha_C": ac, "lost_distal_C": sorted(r_dist)}
                    break
            if ok:
                break
        checks += [
            _check("SAME_CLASS_CENTER_RETAINED", r_class and p_class),
            _check("SAME_ESTER_O_RETAINED", bool(common_o) if (ro_complete and po_complete) else None),
            _check("ETHYL_TO_METHYL_PARTNER_EDIT", ok, json.dumps(detail, sort_keys=True)),
        ]

    checks += [
        _check("carbonyl_oxygen_lineage", lineage),
        _check("mapped_carbon_scaffold_preserved", scaffold, scaffold_detail),
    ]
    graph_ok = all(c.passed is True for c in checks[1:] if c.name != "mapped_center_present")
    checks.append(_check("GRAPH_EDIT_COMPATIBLE", graph_ok))
    return _finish(record, "VC7_SAME_CLASS_EXCHANGE", checks)


def validate_wave_c_bond_attachment(
    facts: dict[str, Any],
    record: SemanticRecord,
    contract: Any,
) -> SemanticValidation | None:
    if record.event in VC3_EVENTS:
        return _validate_vc3(facts, record)
    if record.event in VC4_EVENTS:
        return _validate_vc4(facts, record)
    if record.event in VC5_EVENTS:
        return _validate_vc5(facts, record)
    if record.event in VC7_EVENTS:
        return _validate_vc7(facts, record)
    return None
