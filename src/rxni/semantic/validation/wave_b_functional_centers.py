"""Stage B validators: VC2 acyl, VC9 nitrile, VC8 heteroatom oxidation state.

Audit-only. These validators prove local mapped transitions; they do not alter
canonical semantics, redox, profiles, or retrieval.

Core rule:
  mapped anchor
  + source present and destroyed
  + target absent and created
  + event-compatible local graph edit
  + no incompatible mapped scaffold rewiring
"""

from __future__ import annotations

from typing import Any, Iterable
import json

from ..states.carbon import anhydride_partner_centers, carbon_state
from ..states.common import HALOGENS, atom, attachment_map, attachments, bond_order
from ..states.nitrogen import nitrogen_state, nitrogen_state_detailed
from ..states.sulfur import sulfur_state
from ..types import SemanticRecord, SemanticValidation, ValidationCheck, ValidationStatus
from ..version import FGI_CONTRACT_VERSION

VALIDATION_IMPLEMENTATION_VERSION = "wave-b-functional-centers-v1"

VC2_EVENTS = {
    "acyl_halide_to_amide",
    "acyl_halide_to_aryl_ketone",
    "acyl_halide_to_ketene",
    "carboxylic_acid_to_acyl_halide",
    "carboxylic_acid_to_amide",
    "carboxylic_acid_to_ester",
    "carboxylic_acids_to_anhydride",
    "ester_to_amide",
    "ester_to_carboxylic_acid",
    "haloformate_to_carbamate",
    "intramolecular_carboxylic_acid_to_ester",
    "weinreb_amide_to_ketone",
}

VC9_EVENTS = {
    "carboxylic_acid_to_nitrile",
    "nitrile_to_aldehyde",
    "nitrile_to_amide",
    "nitrile_to_amidine",
    "nitrile_to_amine",
    "nitrile_to_carboxylic_acid",
    "nitrile_to_ester",
    "nitrile_to_imidate",
    "nitrile_to_orthoester",
    "nitrile_to_thioamide",
}

VC8_EVENTS = {
    "aza_n_oxide_deoxygenation",
    "aza_n_oxide_to_aza_nitrogen",
    "aza_nitrogen_to_n_oxide",
    "azide_to_amine",
    "nitro_to_amine",
    "nitro_to_hydroxylamine",
    "sulfide_to_sulfone",
    "sulfide_to_sulfoxide",
    "sulfone_to_sulfide",
    "sulfone_to_sulfoxide",
    "sulfoxide_to_sulfide",
    "sulfoxide_to_sulfone",
    "tertiary_amine_n_oxide_to_amine",
    "tertiary_amine_to_n_oxide",
}

WAVE_B_EVENTS = VC2_EVENTS | VC9_EVENTS | VC8_EVENTS


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


def _mapped_neighbor_ids(
    facts: dict[str, Any],
    side: str,
    center: int,
    *,
    elements: set[str] | None = None,
    exclude_elements: set[str] | None = None,
) -> tuple[set[int], bool]:
    out: set[int] = set()
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
        out.add(int(mn))
    return out, complete


def _element_atts(
    facts: dict[str, Any],
    side: str,
    center: int,
    element: str,
) -> list[dict[str, Any]]:
    return [x for x in attachments(facts, side, center) if x.get("element") == element]


def _carbonyl_oxygen_maps(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], bool]:
    out: set[int] = set()
    complete = True
    for x in attachments(facts, side, center):
        if x.get("element") != "O" or float(x.get("bond_order", 0.0) or 0.0) < 1.8:
            continue
        mn = attachment_map(x)
        if mn is None:
            complete = False
        else:
            out.add(int(mn))
    return out, complete


def _is_ketene(facts: dict[str, Any], side: str, center: int) -> bool:
    atts = attachments(facts, side, center)
    return (
        any(x.get("element") == "O" and float(x.get("bond_order", 0.0) or 0.0) >= 1.8 for x in atts)
        and any(x.get("element") == "C" and float(x.get("bond_order", 0.0) or 0.0) >= 1.8 for x in atts)
    )


def _is_carbamate(facts: dict[str, Any], side: str, center: int) -> bool:
    if atom(facts, side, center).get("element") != "C":
        return False
    atts = attachments(facts, side, center)
    return (
        any(x.get("element") == "O" and float(x.get("bond_order", 0.0) or 0.0) >= 1.8 for x in atts)
        and any(x.get("element") == "O" and float(x.get("bond_order", 0.0) or 0.0) <= 1.2 for x in atts)
        and any(x.get("element") == "N" and float(x.get("bond_order", 0.0) or 0.0) <= 1.6 for x in atts)
    )


def _is_cyclic_ester(facts: dict[str, Any], side: str, center: int) -> bool:
    if carbon_state(facts, side, center) != "ester":
        return False
    if side != "product":
        return False
    if int(facts.get("common_cycle_rank_delta", 0) or 0) <= 0:
        return False
    for rec in facts.get("ring_formations", []) or []:
        maps = set(rec.get("conserved_maps") or rec.get("ring_maps") or [])
        if int(center) in {int(x) for x in maps}:
            return True
    return False


def _is_weinreb_amide(facts: dict[str, Any], side: str, center: int) -> bool:
    if carbon_state(facts, side, center) != "amide":
        return False
    for natt in _element_atts(facts, side, center, "N"):
        if float(natt.get("bond_order", 0.0) or 0.0) > 1.6:
            continue
        shell = list(natt.get("neighbor_shell") or [])
        elems = {str(x[0]) for x in shell}
        if "O" in elems and "C" in elems:
            return True
    return False


def _is_amidine(facts: dict[str, Any], side: str, center: int) -> bool:
    n = _element_atts(facts, side, center, "N")
    if len(n) != 2:
        return False
    orders = sorted(float(x.get("bond_order", 0.0) or 0.0) for x in n)
    return orders[0] <= 1.2 and orders[1] >= 1.5


def _is_imidate(facts: dict[str, Any], side: str, center: int) -> bool:
    n = _element_atts(facts, side, center, "N")
    o = _element_atts(facts, side, center, "O")
    return bool(n and o and any(float(x.get("bond_order", 0.0) or 0.0) >= 1.5 for x in n))


def _is_orthoester(facts: dict[str, Any], side: str, center: int) -> bool:
    n = _element_atts(facts, side, center, "N")
    o = _element_atts(facts, side, center, "O")
    return len(o) >= 3 and not n


def _is_thioamide(facts: dict[str, Any], side: str, center: int) -> bool:
    return (
        any(x.get("element") == "S" and float(x.get("bond_order", 0.0) or 0.0) >= 1.8 for x in attachments(facts, side, center))
        and any(x.get("element") == "N" and float(x.get("bond_order", 0.0) or 0.0) <= 1.2 for x in attachments(facts, side, center))
    )


def _is_two_carboxylic_acids(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:
    if carbon_state(facts, side, center) != "carboxylic_acid":
        return False
    maps = facts.get("mapping", {}).get("common_maps", []) or []
    return any(
        int(m) != int(center)
        and carbon_state(facts, side, int(m)) == "carboxylic_acid"
        for m in maps
    )


def _is_hydroxylamine(facts: dict[str, Any], side: str, center: int) -> bool:
    a = atom(facts, side, center)
    if a.get("element") != "N" or int(a.get("charge", 0) or 0) != 0:
        return False
    if int(a.get("total_h", 0) or 0) < 1:
        return False
    oxy = [
        x for x in attachments(facts, side, center)
        if x.get("element") == "O" and float(x.get("bond_order", 0.0) or 0.0) <= 1.2
    ]
    if len(oxy) != 1:
        return False
    om = attachment_map(oxy[0])
    if om is None:
        return True
    return int(atom(facts, side, int(om)).get("total_h", 0) or 0) >= 1


def _state_match(
    facts: dict[str, Any],
    side: str,
    center: int,
    state: str,
) -> bool | None:
    native_carbon = {
        "acyl_halide", "amide", "ketone", "carboxylic_acid", "ester",
        "anhydride", "haloformate", "nitrile", "aldehyde", "amine_carbon",
    }
    if state in native_carbon:
        return carbon_state(facts, side, center) == state
    if state == "ketene":
        return _is_ketene(facts, side, center)
    if state == "carbamate":
        return _is_carbamate(facts, side, center)
    if state == "cyclic_ester":
        return _is_cyclic_ester(facts, side, center)
    if state == "weinreb_amide":
        return _is_weinreb_amide(facts, side, center)
    if state == "amidine":
        return _is_amidine(facts, side, center)
    if state == "imidate":
        return _is_imidate(facts, side, center)
    if state == "orthoester":
        return _is_orthoester(facts, side, center)
    if state == "thioamide":
        return _is_thioamide(facts, side, center)
    if state == "two_carboxylic_acids":
        return _is_two_carboxylic_acids(facts, side, center)

    if state in {"aza_n_oxide", "aza_nitrogen", "azide_n", "nitro", "amine_n_oxide", "amine_nitrogen"}:
        return nitrogen_state(facts, side, center) == state
    if state == "tertiary_amine":
        return nitrogen_state_detailed(facts, side, center) == "tertiary_amine"
    if state == "hydroxylamine":
        return _is_hydroxylamine(facts, side, center)

    if state in {"sulfide", "sulfoxide", "sulfone"}:
        return sulfur_state(facts, side, center) == state

    return None


def _states_match(
    facts: dict[str, Any],
    side: str,
    center: int,
    states: tuple[str, ...],
) -> bool | None:
    if not states:
        return None
    saw_unknown = False
    for state in states:
        v = _state_match(facts, side, center, state)
        if v is True:
            return True
        if v is None:
            saw_unknown = True
    return None if saw_unknown else False


def _local_signature(facts: dict[str, Any], side: str, center: int) -> tuple:
    rows = []
    for x in attachments(facts, side, center):
        if x.get("element") == "H":
            continue
        rows.append((
            str(x.get("element") or ""),
            attachment_map(x),
            round(float(x.get("bond_order", 0.0) or 0.0), 3),
        ))
    return tuple(sorted(rows, key=lambda z: (z[0], -1 if z[1] is None else z[1], z[2])))


def _ring_touches(facts: dict[str, Any], center: int) -> bool:
    for rec in facts.get("ring_formations", []) or []:
        maps = set(rec.get("conserved_maps") or rec.get("ring_maps") or [])
        if int(center) in {int(x) for x in maps}:
            return True
    return False


def _vc2_graph_edit(
    facts: dict[str, Any],
    event: str,
    center: int,
) -> tuple[bool | None, str]:
    ra = attachments(facts, "reactant", center)
    pa = attachments(facts, "product", center)

    def count(rows, element, lo=0.0, hi=99.0):
        return sum(
            1 for x in rows
            if x.get("element") == element
            and lo <= float(x.get("bond_order", 0.0) or 0.0) <= hi
        )

    r_c, rc = _mapped_neighbor_ids(facts, "reactant", center, elements={"C"})
    p_c, pc = _mapped_neighbor_ids(facts, "product", center, elements={"C"})
    r_co, rco = _carbonyl_oxygen_maps(facts, "reactant", center)
    p_co, pco = _carbonyl_oxygen_maps(facts, "product", center)

    carbonyl_lineage = (
        bool(r_co & p_co)
        if (r_co and p_co)
        else None
    )

    graph = None
    detail = {}

    if event == "acyl_halide_to_amide":
        graph = count(ra, "Cl", 0.0, 1.2) + count(ra, "Br", 0.0, 1.2) + count(ra, "I", 0.0, 1.2) + count(ra, "F", 0.0, 1.2) >= 1 and count(pa, "N", 0.0, 1.6) >= 1 and r_c == p_c
    elif event == "acyl_halide_to_aryl_ketone":
        gained = p_c - r_c
        graph = len(gained) >= 1 and r_c.issubset(p_c)
    elif event == "acyl_halide_to_ketene":
        graph = r_c == p_c and any(x.get("element") == "C" and float(x.get("bond_order", 0.0) or 0.0) >= 1.8 for x in pa)
    elif event == "carboxylic_acid_to_acyl_halide":
        graph = count(pa, "Cl", 0.0, 1.2) + count(pa, "Br", 0.0, 1.2) + count(pa, "I", 0.0, 1.2) + count(pa, "F", 0.0, 1.2) >= 1 and r_c == p_c
    elif event == "carboxylic_acid_to_amide":
        graph = count(pa, "N", 0.0, 1.6) >= 1 and r_c == p_c
    elif event == "carboxylic_acid_to_ester":
        graph = carbon_state(facts, "product", center) == "ester" and r_c == p_c
    elif event == "carboxylic_acids_to_anhydride":
        partners = anhydride_partner_centers(facts, "product", center)
        graph = bool(partners) and any(carbon_state(facts, "reactant", int(x)) == "carboxylic_acid" for x in partners)
    elif event == "ester_to_amide":
        graph = count(pa, "N", 0.0, 1.6) >= 1 and r_c == p_c
    elif event == "ester_to_carboxylic_acid":
        graph = carbon_state(facts, "product", center) == "carboxylic_acid" and r_c == p_c
    elif event == "haloformate_to_carbamate":
        graph = (
            _is_carbamate(facts, "product", center)
            and count(pa, "N", 0.0, 1.6) >= 1
            and r_c == p_c
        )
    elif event == "intramolecular_carboxylic_acid_to_ester":
        graph = _is_cyclic_ester(facts, "product", center) and r_c == p_c
    elif event == "weinreb_amide_to_ketone":
        gained = p_c - r_c
        graph = len(gained) >= 1 and r_c.issubset(p_c)

    detail = {
        "reactant_C_neighbors": sorted(r_c),
        "product_C_neighbors": sorted(p_c),
        "reactant_carbonyl_O": sorted(r_co),
        "product_carbonyl_O": sorted(p_co),
        "carbonyl_lineage": carbonyl_lineage,
    }
    if graph is True and carbonyl_lineage is False:
        graph = False
    elif graph is True and carbonyl_lineage is None:
        graph = None

    return graph, json.dumps(detail, sort_keys=True)


def _validate_vc2(
    facts: dict[str, Any],
    record: SemanticRecord,
    contract: Any,
) -> SemanticValidation:
    c = record.center_map
    if c is None:
        return _finish(record, "VC2_ACYL_CENTER", [_check("mapped_center_present", None)])
    c = int(c)
    rs = _states_match(facts, "reactant", c, tuple(contract.source_states))
    ps = _states_match(facts, "product", c, tuple(contract.source_states))
    rt = _states_match(facts, "reactant", c, tuple(contract.target_states))
    pt = _states_match(facts, "product", c, tuple(contract.target_states))
    graph, graph_detail = _vc2_graph_edit(facts, record.event, c)
    checks = [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("mapped_center_is_carbon", atom(facts, "reactant", c).get("element") == "C" and atom(facts, "product", c).get("element") == "C"),
        _check("SOURCE_PRESENT", rs),
        _check("SOURCE_DESTROYED", None if ps is None else not ps),
        _check("TARGET_ABSENT_BEFORE", None if rt is None else not rt),
        _check("TARGET_CREATED", pt),
        _check("GRAPH_EDIT_COMPATIBLE", graph, graph_detail),
    ]
    return _finish(record, "VC2_ACYL_CENTER", checks)


def _vc9_graph_edit(
    facts: dict[str, Any],
    event: str,
    center: int,
) -> tuple[bool | None, str]:
    r_c, rc = _mapped_neighbor_ids(facts, "reactant", center, elements={"C"})
    p_c, pc = _mapped_neighbor_ids(facts, "product", center, elements={"C"})
    if not (rc and pc):
        return None, "mapped carbon-neighbor identity incomplete"
    if r_c != p_c:
        return False, f"carbon scaffold changed reactant={sorted(r_c)} product={sorted(p_c)}"
    if _ring_touches(facts, center) or atom(facts, "product", center).get("aromatic"):
        return False, "nitrile carbon incorporated into ring/aromatic topology"

    r_triple_n = [
        x for x in attachments(facts, "reactant", center)
        if x.get("element") == "N" and float(x.get("bond_order", 0.0) or 0.0) >= 2.8
    ]
    p_triple_n = [
        x for x in attachments(facts, "product", center)
        if x.get("element") == "N" and float(x.get("bond_order", 0.0) or 0.0) >= 2.8
    ]

    if event == "carboxylic_acid_to_nitrile":
        ok = not r_triple_n and bool(p_triple_n)
    else:
        ok = bool(r_triple_n) and not p_triple_n

    return ok, f"reactant_C={sorted(r_c)} product_C={sorted(p_c)} r_C#N={len(r_triple_n)} p_C#N={len(p_triple_n)}"


def _validate_vc9(
    facts: dict[str, Any],
    record: SemanticRecord,
    contract: Any,
) -> SemanticValidation:
    c = record.center_map
    if c is None:
        return _finish(record, "VC9_NITRILE_CENTER", [_check("mapped_center_present", None)])
    c = int(c)
    rs = _states_match(facts, "reactant", c, tuple(contract.source_states))
    ps = _states_match(facts, "product", c, tuple(contract.source_states))
    rt = _states_match(facts, "reactant", c, tuple(contract.target_states))
    pt = _states_match(facts, "product", c, tuple(contract.target_states))
    graph, detail = _vc9_graph_edit(facts, record.event, c)
    checks = [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("mapped_center_is_carbon", atom(facts, "reactant", c).get("element") == "C" and atom(facts, "product", c).get("element") == "C"),
        _check("SOURCE_PRESENT", rs),
        _check("SOURCE_DESTROYED", None if ps is None else not ps),
        _check("TARGET_ABSENT_BEFORE", None if rt is None else not rt),
        _check("TARGET_CREATED", pt),
        _check("GRAPH_EDIT_COMPATIBLE", graph, detail),
    ]
    return _finish(record, "VC9_NITRILE_CENTER", checks)


def _non_oxygen_heavy_ids(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], bool]:
    return _mapped_neighbor_ids(
        facts, side, center, exclude_elements={"O"}
    )


def _oxo_count(facts: dict[str, Any], side: str, center: int) -> int:
    return sum(
        1 for x in attachments(facts, side, center)
        if x.get("element") == "O"
        and float(x.get("bond_order", 0.0) or 0.0) >= 1.5
    )


def _oxygen_count(facts: dict[str, Any], side: str, center: int) -> int:
    return sum(1 for x in attachments(facts, side, center) if x.get("element") == "O")


def _vc8_graph_edit(
    facts: dict[str, Any],
    event: str,
    center: int,
) -> tuple[bool | None, str]:
    r_heavy, rc = _non_oxygen_heavy_ids(facts, "reactant", center)
    p_heavy, pc = _non_oxygen_heavy_ids(facts, "product", center)
    if not (rc and pc):
        return None, "non-oxygen heavy-atom mapping incomplete"
    if r_heavy != p_heavy:
        return False, f"non-O heavy identity changed reactant={sorted(r_heavy)} product={sorted(p_heavy)}"

    rO = _oxygen_count(facts, "reactant", center)
    pO = _oxygen_count(facts, "product", center)
    rOx = _oxo_count(facts, "reactant", center)
    pOx = _oxo_count(facts, "product", center)

    if event in {"aza_n_oxide_deoxygenation", "aza_n_oxide_to_aza_nitrogen", "tertiary_amine_n_oxide_to_amine"}:
        ok = rO >= 1 and pO == 0
    elif event in {"aza_nitrogen_to_n_oxide", "tertiary_amine_to_n_oxide"}:
        ok = rO == 0 and pO >= 1
    elif event == "nitro_to_hydroxylamine":
        # Exact mapped N=O -> N-O reduction where possible.
        common_o = []
        for x in attachments(facts, "reactant", center):
            if x.get("element") != "O" or float(x.get("bond_order", 0.0) or 0.0) < 1.8:
                continue
            om = attachment_map(x)
            if om is None:
                continue
            if 0.8 <= bond_order(facts, "product", center, int(om)) <= 1.2:
                common_o.append(int(om))
        ok = rO >= 2 and pO == 1 and bool(common_o)
    elif event == "sulfide_to_sulfoxide":
        ok = rOx == 0 and pOx == 1
    elif event == "sulfide_to_sulfone":
        ok = rOx == 0 and pOx >= 2
    elif event == "sulfoxide_to_sulfide":
        ok = rOx == 1 and pOx == 0
    elif event == "sulfoxide_to_sulfone":
        ok = rOx == 1 and pOx >= 2
    elif event == "sulfone_to_sulfide":
        ok = rOx >= 2 and pOx == 0
    elif event == "sulfone_to_sulfoxide":
        ok = rOx >= 2 and pOx == 1
    else:
        # azide_to_amine and nitro_to_amine are handled by priority1 before
        # reaching this validator.
        ok = True

    return ok, f"nonO={sorted(r_heavy)} O:{rO}->{pO} oxo:{rOx}->{pOx}"


def _validate_vc8(
    facts: dict[str, Any],
    record: SemanticRecord,
    contract: Any,
) -> SemanticValidation:
    c = record.center_map
    if c is None:
        return _finish(record, "VC8_HETEROATOM_OXSTATE", [_check("mapped_center_present", None)])
    c = int(c)
    rs = _states_match(facts, "reactant", c, tuple(contract.source_states))
    ps = _states_match(facts, "product", c, tuple(contract.source_states))
    rt = _states_match(facts, "reactant", c, tuple(contract.target_states))
    pt = _states_match(facts, "product", c, tuple(contract.target_states))
    graph, detail = _vc8_graph_edit(facts, record.event, c)
    checks = [
        _check("mapped_center_present", bool(atom(facts, "reactant", c)) and bool(atom(facts, "product", c))),
        _check("same_element_anchor", atom(facts, "reactant", c).get("element") == atom(facts, "product", c).get("element")),
        _check("SOURCE_PRESENT", rs),
        _check("SOURCE_DESTROYED", None if ps is None else not ps),
        _check("TARGET_ABSENT_BEFORE", None if rt is None else not rt),
        _check("TARGET_CREATED", pt),
        _check("GRAPH_EDIT_COMPATIBLE", graph, detail),
    ]
    return _finish(record, "VC8_HETEROATOM_OXSTATE", checks)


def validate_wave_b_functional_centers(
    facts: dict[str, Any],
    record: SemanticRecord,
    contract: Any,
) -> SemanticValidation | None:
    if record.event in VC2_EVENTS:
        return _validate_vc2(facts, record, contract)
    if record.event in VC9_EVENTS:
        return _validate_vc9(facts, record, contract)
    if record.event in VC8_EVENTS:
        return _validate_vc8(facts, record, contract)
    return None
