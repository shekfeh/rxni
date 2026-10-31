"""VC1 same-center directional FGI validation.

Milestone M3B2B.

This module validates directional conversions whose chemistry is representable
as a change in local functional state at one conserved mapped carbon center.

Validation is intentionally stricter than historical detection:

    SOURCE_PRESENT
    SOURCE_DESTROYED
    TARGET_ABSENT_BEFORE
    TARGET_CREATED
    mapped carbon scaffold preserved
    local structural change present

The module does not mutate semantic records and is not wired into canonical
resolution, redox, profiles, or retrieval.
"""

from __future__ import annotations

from typing import Any, Iterable
import json

from ..states.carbon import carbon_state
from ..states.common import atom, attachment_map, attachments
from ..types import (
    SemanticRecord,
    SemanticValidation,
    ValidationCheck,
    ValidationStatus,
)
from ..version import FGI_CONTRACT_VERSION


VALIDATION_IMPLEMENTATION_VERSION = "vc1-center-state-audit-v2"

# All 18 VC1 events now have executable audit validators. Most use the native
# carbon-state vocabulary; oxime_to_nitrile and enol_ether_to_carbonyl use
# explicit local topology predicates because forcing them into carbon_state()
# would blur structural-state responsibilities.
SUPPORTED_EVENTS = {
    "aldehyde_to_carboxylic_acid",
    "aldehyde_to_ester",
    "aldehyde_to_nitrile",
    "aldehyde_to_primary_alcohol",
    "amide_to_aldehyde",
    "amide_to_amine",
    "amide_to_nitrile",
    "carboxylic_acid_to_aldehyde",
    "carboxylic_acid_to_primary_alcohol",
    "ester_to_aldehyde",
    "ester_to_primary_alcohol",
    "ketone_to_secondary_alcohol",
    "methylene_to_ketone",
    "primary_alcohol_to_aldehyde",
    "primary_alcohol_to_carboxylic_acid",
    "secondary_alcohol_to_ketone",
    "enol_ether_to_carbonyl",
    "oxime_to_nitrile",
}

DEFERRED_EVENTS = set()


def _check(name: str, passed: bool | None, detail: str = "") -> ValidationCheck:
    return ValidationCheck(name=name, passed=passed, detail=detail)


NON_VETO_DIAGNOSTICS = {
    "mapped_carbon_bond_orders_preserved",
}


def _finish(
    record: SemanticRecord,
    checks: Iterable[ValidationCheck],
) -> SemanticValidation:
    checks = tuple(checks)
    failed = tuple(
        x.name
        for x in checks
        if x.passed is False and x.name not in NON_VETO_DIAGNOSTICS
    )
    unknown = tuple(
        x.name
        for x in checks
        if x.passed is None and x.name not in NON_VETO_DIAGNOSTICS
    )

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
        validator_class="VC1_CENTER_STATE",
        checks=checks,
        veto_reasons=veto,
        contract_version=FGI_CONTRACT_VERSION,
    )


def _mapped_carbon_neighbors(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[set[int], dict[int, float], bool]:
    """Mapped carbon-neighbor identity plus bond-order diagnostics.

    VC1 scaffold preservation is defined by *neighbor identity*, not by
    preservation of C-C bond order. Carbonyl reduction/oxidation embedded in
    conjugated or aromatic systems can legitimately change adjacent C-C bond
    orders while retaining the same carbon skeleton.

    Bond orders are still returned for diagnostics and later event-specific
    graph-edit compatibility checks.
    """
    ids: set[int] = set()
    orders: dict[int, float] = {}
    complete = True

    for att in attachments(facts, side, center):
        if att.get("element") != "C":
            continue
        mn = attachment_map(att)
        if mn is None:
            complete = False
            continue
        mn = int(mn)
        ids.add(mn)
        orders[mn] = round(
            float(att.get("bond_order", 0.0) or 0.0),
            3,
        )

    return ids, orders, complete


def _local_heavy_signature(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> tuple[tuple[str, int | None, float], ...]:
    rows = []
    for att in attachments(facts, side, center):
        el = str(att.get("element") or "")
        if not el or el == "H":
            continue
        mn = attachment_map(att)
        rows.append(
            (
                el,
                None if mn is None else int(mn),
                round(float(att.get("bond_order", 0.0) or 0.0), 3),
            )
        )
    return tuple(
        sorted(
            rows,
            key=lambda x: (
                x[0],
                -1 if x[1] is None else x[1],
                x[2],
            ),
        )
    )


def _plain_methylene(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:
    """Strict source predicate for `methylene_to_ketone`.

    A "methylene" source is a non-aromatic carbon with >=2 hydrogens whose
    heavy-atom attachments are carbon only. Ether/acetal, amine, sulfur,
    halide, silicon, etc. substituted CH2 centers are intentionally excluded.

    This definition is deliberately narrower than the historical detector,
    which emitted `methylene_to_ketone` whenever carbon_state()=="".
    """
    a = atom(facts, side, center)
    if a.get("element") != "C":
        return False
    if bool(a.get("aromatic")):
        return False
    if int(a.get("total_h", 0)) < 2:
        return False
    if carbon_state(facts, side, center) != "":
        return False

    heavy = [
        x for x in attachments(facts, side, center)
        if x.get("element") != "H"
    ]
    return bool(heavy) and all(x.get("element") == "C" for x in heavy)


def _bond_order(
    facts: dict[str, Any],
    side: str,
    a: int,
    b: int,
) -> float:
    return float(
        facts[side]["bonds"].get(
            tuple(sorted((int(a), int(b)))),
            0.0,
        )
        or 0.0
    )


def _heavy_attachment_maps(
    facts: dict[str, Any],
    side: str,
    center: int,
    *,
    exclude: set[int] | None = None,
) -> tuple[set[int], bool]:
    exclude = set(exclude or ())
    out: set[int] = set()
    complete = True
    for att in attachments(facts, side, center):
        if att.get("element") == "H":
            continue
        mn = attachment_map(att)
        if mn is None:
            complete = False
            continue
        mn = int(mn)
        if mn not in exclude:
            out.add(mn)
    return out, complete


def _aldoxime_n(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> int | None:
    """Return mapped oxime N for an aldoxime-like C=N-O center."""
    ca = atom(facts, side, center)
    if ca.get("element") != "C":
        return None
    if int(ca.get("total_h", 0) or 0) < 1:
        return None

    candidates: list[int] = []
    for att in attachments(facts, side, center):
        if att.get("element") != "N":
            continue
        mn = attachment_map(att)
        if mn is None:
            continue
        bo = float(att.get("bond_order", 0.0) or 0.0)
        if not (1.7 <= bo <= 2.2):
            continue
        n = int(mn)
        has_o = any(
            x.get("element") == "O"
            for x in attachments(facts, side, n)
        )
        if has_o:
            candidates.append(n)

    return candidates[0] if len(candidates) == 1 else None


def _is_aldoxime(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:
    return _aldoxime_n(facts, side, center) is not None


def _enol_ether_candidates(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> list[tuple[int, int]]:
    """Return (alkene_partner_C, ether_O) for O-alkyl enol ether centers.

    Structural definition:
      * mapped carbon center,
      * C=C to mapped carbon partner,
      * C-O single bond to mapped oxygen,
      * that oxygen has another carbon neighbor besides the enol carbon.

    This intentionally excludes silyl enol ethers and simple enols.
    """
    ca = atom(facts, side, center)
    if ca.get("element") != "C":
        return []

    partners: list[int] = []
    oxygens: list[int] = []

    for att in attachments(facts, side, center):
        mn = attachment_map(att)
        if mn is None:
            continue
        mn = int(mn)
        bo = float(att.get("bond_order", 0.0) or 0.0)

        if att.get("element") == "C" and 1.7 <= bo <= 2.2:
            partners.append(mn)

        if att.get("element") == "O" and 0.8 <= bo <= 1.2:
            other_carbon = any(
                x.get("element") == "C"
                and (
                    attachment_map(x) is None
                    or int(attachment_map(x)) != center
                )
                for x in attachments(facts, side, mn)
            )
            if other_carbon:
                oxygens.append(mn)

    return [
        (partner, oxygen)
        for partner in partners
        for oxygen in oxygens
    ]


def _is_o_alkyl_enol_ether(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:
    return bool(_enol_ether_candidates(facts, side, center))


def _is_carbonyl_target(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:
    return carbon_state(facts, side, center) in {"aldehyde", "ketone"}


def _validate_oxime_to_nitrile(
    facts: dict[str, Any],
    record: SemanticRecord,
) -> SemanticValidation:
    center = record.center_map
    if center is None:
        return _finish(
            record,
            [_check("mapped_center_present", None, "hypothesis has no center_map")],
        )

    center = int(center)
    ra = atom(facts, "reactant", center)
    pa = atom(facts, "product", center)

    source_before = _is_aldoxime(facts, "reactant", center)
    source_after = _is_aldoxime(facts, "product", center)
    target_before = carbon_state(facts, "reactant", center) == "nitrile"
    target_after = carbon_state(facts, "product", center) == "nitrile"

    n = _aldoxime_n(facts, "reactant", center)

    same_cn_pair = None
    cn_double_to_triple = None
    n_o_lost = None
    skeleton_preserved = None
    skeleton_complete = None

    if n is not None:
        same_cn_pair = bool(atom(facts, "product", n))
        r_cn = _bond_order(facts, "reactant", center, n)
        p_cn = _bond_order(facts, "product", center, n)
        cn_double_to_triple = (
            1.7 <= r_cn <= 2.2 and p_cn >= 2.8
        )

        r_no = [
            x for x in attachments(facts, "reactant", n)
            if x.get("element") == "O"
        ]
        p_no = [
            x for x in attachments(facts, "product", n)
            if x.get("element") == "O"
        ]
        n_o_lost = bool(r_no) and not p_no

        r_skel, r_complete = _heavy_attachment_maps(
            facts, "reactant", center, exclude={n}
        )
        p_skel, p_complete = _heavy_attachment_maps(
            facts, "product", center, exclude={n}
        )
        skeleton_complete = True if (r_complete and p_complete) else None
        skeleton_preserved = (
            r_skel == p_skel
            if (r_complete and p_complete)
            else None
        )

    checks = [
        _check(
            "mapped_center_present",
            bool(ra) and bool(pa),
            f"reactant={bool(ra)} product={bool(pa)}",
        ),
        _check(
            "mapped_center_is_carbon",
            ra.get("element") == "C" and pa.get("element") == "C",
            f"reactant={ra.get('element')!r} product={pa.get('element')!r}",
        ),
        _check("SOURCE_PRESENT", source_before, "aldoxime-like C=N-O at mapped carbon"),
        _check("SOURCE_DESTROYED", not source_after, "aldoxime topology absent after"),
        _check("TARGET_ABSENT_BEFORE", not target_before, "nitrile absent before"),
        _check("TARGET_CREATED", target_after, "nitrile present after"),
        _check("same_mapped_cn_pair", same_cn_pair, f"nitrogen_map={n}"),
        _check(
            "cn_double_to_triple",
            cn_double_to_triple,
            (
                ""
                if n is None
                else (
                    f"reactant={_bond_order(facts, 'reactant', center, n)} "
                    f"product={_bond_order(facts, 'product', center, n)}"
                )
            ),
        ),
        _check("nitrogen_oxygen_attachment_lost", n_o_lost, f"nitrogen_map={n}"),
        _check("mapped_scaffold_complete", skeleton_complete, ""),
        _check("mapped_scaffold_preserved", skeleton_preserved, ""),
    ]
    return _finish(record, checks)


def _validate_enol_ether_to_carbonyl(
    facts: dict[str, Any],
    record: SemanticRecord,
) -> SemanticValidation:
    center = record.center_map
    if center is None:
        return _finish(
            record,
            [_check("mapped_center_present", None, "hypothesis has no center_map")],
        )

    center = int(center)
    ra = atom(facts, "reactant", center)
    pa = atom(facts, "product", center)

    source_before = _is_o_alkyl_enol_ether(facts, "reactant", center)
    source_after = _is_o_alkyl_enol_ether(facts, "product", center)
    target_before = _is_carbonyl_target(facts, "reactant", center)
    target_after = _is_carbonyl_target(facts, "product", center)

    candidates = _enol_ether_candidates(
        facts, "reactant", center
    )

    graph_ok = False
    graph_detail = []
    for partner, oxygen in candidates:
        r_cc = _bond_order(facts, "reactant", center, partner)
        p_cc = _bond_order(facts, "product", center, partner)
        r_co = _bond_order(facts, "reactant", center, oxygen)
        p_co = _bond_order(facts, "product", center, oxygen)

        p_o_heavy = [
            x for x in attachments(facts, "product", oxygen)
            if x.get("element") != "H"
        ]
        o_only_center = (
            len(p_o_heavy) == 1
            and attachment_map(p_o_heavy[0]) is not None
            and int(attachment_map(p_o_heavy[0])) == center
        )

        ok = (
            1.7 <= r_cc <= 2.2
            and 0.8 <= p_cc <= 1.2
            and 0.8 <= r_co <= 1.2
            and p_co >= 1.8
            and o_only_center
        )
        graph_detail.append(
            {
                "partner": partner,
                "oxygen": oxygen,
                "r_cc": r_cc,
                "p_cc": p_cc,
                "r_co": r_co,
                "p_co": p_co,
                "product_oxygen_only_center": o_only_center,
                "compatible": ok,
            }
        )
        if ok:
            graph_ok = True

    r_c, _, r_complete = _mapped_carbon_neighbors(
        facts, "reactant", center
    )
    p_c, _, p_complete = _mapped_carbon_neighbors(
        facts, "product", center
    )

    checks = [
        _check(
            "mapped_center_present",
            bool(ra) and bool(pa),
            f"reactant={bool(ra)} product={bool(pa)}",
        ),
        _check(
            "mapped_center_is_carbon",
            ra.get("element") == "C" and pa.get("element") == "C",
            f"reactant={ra.get('element')!r} product={pa.get('element')!r}",
        ),
        _check("SOURCE_PRESENT", source_before, "O-alkyl enol ether at mapped carbon"),
        _check("SOURCE_DESTROYED", not source_after, "enol ether topology absent after"),
        _check("TARGET_ABSENT_BEFORE", not target_before, "aldehyde/ketone absent before"),
        _check("TARGET_CREATED", target_after, "aldehyde/ketone present after"),
        _check(
            "enol_bond_to_carbonyl_graph_edit",
            graph_ok,
            json.dumps(graph_detail, sort_keys=True),
        ),
        _check(
            "mapped_carbon_scaffold_complete",
            True if (r_complete and p_complete) else None,
            f"reactant_complete={r_complete} product_complete={p_complete}",
        ),
        _check(
            "mapped_carbon_scaffold_preserved",
            r_c == p_c if (r_complete and p_complete) else None,
            f"reactant_neighbors={sorted(r_c)} product_neighbors={sorted(p_c)}",
        ),
    ]
    return _finish(record, checks)



def _state_match(
    facts: dict[str, Any],
    side: str,
    center: int,
    states: tuple[str, ...],
) -> bool | None:
    if not states:
        return None

    # VC1 registry entries are single-source/single-target in the current
    # contract inventory, but tuple handling is kept explicit.
    observed = carbon_state(facts, side, center)

    for state in states:
        if state == "methylene":
            if _plain_methylene(facts, side, center):
                return True
        elif state in {
            "aldehyde",
            "ketone",
            "ester",
            "carboxylic_acid",
            "amide",
            "nitrile",
            "primary_alcohol",
            "secondary_alcohol",
            "tertiary_alcohol",
            "amine_carbon",
        }:
            if observed == state:
                return True
        else:
            # Unknown native state vocabulary must not become a False or True.
            return None

    return False


def validate_vc1_center_state(
    facts: dict[str, Any],
    record: SemanticRecord,
    contract: Any,
) -> SemanticValidation | None:
    """Validate one supported VC1 hypothesis, otherwise return None."""
    if getattr(contract, "validator_class", None) != "VC1_CENTER_STATE":
        return None

    if record.event not in SUPPORTED_EVENTS:
        return None

    if record.event == "oxime_to_nitrile":
        return _validate_oxime_to_nitrile(facts, record)

    if record.event == "enol_ether_to_carbonyl":
        return _validate_enol_ether_to_carbonyl(facts, record)

    center = record.center_map
    if center is None:
        return _finish(
            record,
            [_check("mapped_center_present", None, "hypothesis has no center_map")],
        )

    center = int(center)
    ra = atom(facts, "reactant", center)
    pa = atom(facts, "product", center)

    source_states = tuple(getattr(contract, "source_states", ()) or ())
    target_states = tuple(getattr(contract, "target_states", ()) or ())

    source_before = _state_match(
        facts, "reactant", center, source_states
    )
    source_after = _state_match(
        facts, "product", center, source_states
    )
    target_before = _state_match(
        facts, "reactant", center, target_states
    )
    target_after = _state_match(
        facts, "product", center, target_states
    )

    r_c, r_c_orders, r_complete = _mapped_carbon_neighbors(
        facts, "reactant", center
    )
    p_c, p_c_orders, p_complete = _mapped_carbon_neighbors(
        facts, "product", center
    )

    r_local = _local_heavy_signature(
        facts, "reactant", center
    )
    p_local = _local_heavy_signature(
        facts, "product", center
    )

    checks = [
        _check(
            "mapped_center_present",
            bool(ra) and bool(pa),
            f"reactant={bool(ra)} product={bool(pa)}",
        ),
        _check(
            "mapped_center_is_carbon",
            (
                ra.get("element") == "C"
                and pa.get("element") == "C"
            ),
            (
                f"reactant_element={ra.get('element')!r} "
                f"product_element={pa.get('element')!r}"
            ),
        ),
        _check(
            "SOURCE_PRESENT",
            source_before,
            (
                f"contract={source_states} "
                f"observed={carbon_state(facts, 'reactant', center)!r}"
            ),
        ),
        _check(
            "SOURCE_DESTROYED",
            (
                None if source_after is None
                else not source_after
            ),
            (
                f"contract={source_states} "
                f"product_observed={carbon_state(facts, 'product', center)!r}"
            ),
        ),
        _check(
            "TARGET_ABSENT_BEFORE",
            (
                None if target_before is None
                else not target_before
            ),
            (
                f"contract={target_states} "
                f"reactant_observed={carbon_state(facts, 'reactant', center)!r}"
            ),
        ),
        _check(
            "TARGET_CREATED",
            target_after,
            (
                f"contract={target_states} "
                f"observed={carbon_state(facts, 'product', center)!r}"
            ),
        ),
        _check(
            "mapped_carbon_scaffold_complete",
            True if (r_complete and p_complete) else None,
            f"reactant_complete={r_complete} product_complete={p_complete}",
        ),
        _check(
            "mapped_carbon_scaffold_preserved",
            (
                r_c == p_c
                if (r_complete and p_complete)
                else None
            ),
            (
                f"reactant_neighbors={sorted(r_c)} "
                f"product_neighbors={sorted(p_c)}"
            ),
        ),
        _check(
            "mapped_carbon_bond_orders_preserved",
            (
                r_c_orders == p_c_orders
                if (r_complete and p_complete and r_c == p_c)
                else None
            ),
            (
                f"reactant_orders={dict(sorted(r_c_orders.items()))} "
                f"product_orders={dict(sorted(p_c_orders.items()))}"
            ),
        ),
        _check(
            "local_heavy_environment_changed",
            r_local != p_local,
            f"reactant={r_local} product={p_local}",
        ),
    ]

    return _finish(record, checks)
