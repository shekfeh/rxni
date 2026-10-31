from __future__ import annotations

from typing import Any

from .mapped_graph import canonical_bond_class
from .rearrangement_extensions import detect_v67_rearrangement_extensions

DONOR_HANDLE_ELEMENTS = {"B", "Sn", "Zn", "Mg", "Si", "Cu", "Al", "Li"}
LEAVING_GROUP_ELEMENTS = {"F", "Cl", "Br", "I"}
HETERO_NU = {"N", "O", "S"}
HALOGENS = {"F", "Cl", "Br", "I"}


def _att_element(att: dict[str, Any]) -> str:
    return str(att.get("element", ""))


def _att_map(att: dict[str, Any]) -> int | None:
    v = att.get("map")
    return int(v) if v is not None else None


def _atom(facts: dict[str, Any], side: str, mn: int) -> dict[str, Any]:
    return facts[side]["atoms"].get(mn, {})


def _attachments(facts: dict[str, Any], side: str, mn: int) -> list[dict[str, Any]]:
    return facts[side]["attachments"].get(mn, [])


def _has_neighbor_element(
    facts: dict[str, Any], side: str, mn: int, elements: set[str]
) -> bool:
    return any(a.get("element") in elements for a in _attachments(facts, side, mn))


def _has_conserved_carbonyl(facts: dict[str, Any], mn: int) -> bool:
    if (
        _atom(facts, "reactant", mn).get("element") != "C"
        or _atom(facts, "product", mn).get("element") != "C"
    ):
        return False

    def has_co(side: str) -> bool:
        return any(
            a.get("element") == "O" and float(a.get("bond_order", 0)) >= 1.9
            for a in _attachments(facts, side, mn)
        )

    return has_co("reactant") and has_co("product")


def _has_conserved_sulfonyl(facts: dict[str, Any], mn: int) -> bool:
    if (
        _atom(facts, "reactant", mn).get("element") != "S"
        or _atom(facts, "product", mn).get("element") != "S"
    ):
        return False

    def n_so(side: str) -> int:
        return sum(
            1
            for a in _attachments(facts, side, mn)
            if a.get("element") == "O" and float(a.get("bond_order", 0)) >= 1.9
        )

    return n_so("reactant") >= 2 and n_so("product") >= 2


def _retained_mapped_oxygen_boron_bond(facts: dict[str, Any], oxygen: int) -> bool:
    """A conserved O-B connection on both sides identifies boronate O-alkyl loss.

    This is not an ether O-dealkylation event; keep its underlying O-C
    cleavage/fragment-loss structural evidence.
    """
    if any(_atom(facts, side, oxygen).get("element") != "O" for side in ("reactant", "product")):
        return False
    boron_before = {
        _att_map(a) for a in _attachments(facts, "reactant", oxygen)
        if a.get("element") == "B" and _att_map(a) is not None
    }
    boron_after = {
        _att_map(a) for a in _attachments(facts, "product", oxygen)
        if a.get("element") == "B" and _att_map(a) is not None
    }
    return bool(boron_before & boron_after)


def _external_gain(facts: dict[str, Any], center: int, att: dict[str, Any]) -> bool:
    gm = _att_map(att)
    if gm is None:
        return False
    return facts["reactant"]["map_to_component"].get(gm) != facts["reactant"][
        "map_to_component"
    ].get(center)


def _support(*maps: int | None) -> list[int]:
    return sorted({int(x) for x in maps if x is not None})


def _is_sp3_carbon(facts: dict[str, Any], side: str, mn: int) -> bool:
    a = _atom(facts, side, mn)
    return (
        a.get("element") == "C"
        and not a.get("aromatic")
        and a.get("hybridization") == "SP3"
    )


def _is_alpha_to_carbonyl(facts: dict[str, Any], side: str, mn: int) -> bool:
    if _atom(facts, side, mn).get("element") != "C":
        return False
    return any(
        x.get("element") == "C" and bool(x.get("carbonyl_like"))
        for x in _attachments(facts, side, mn)
    )


def _is_benzylic_carbon(facts: dict[str, Any], side: str, mn: int) -> bool:
    return any(
        x.get("element") == "C"
        and _att_map(x) is not None
        and bool(_atom(facts, side, int(_att_map(x))).get("aromatic"))
        for x in _attachments(facts, side, mn)
    )


def _same_reactant_origin(facts: dict[str, Any], *maps: int | None) -> bool:
    vals = [m for m in maps if m is not None]
    if not vals:
        return False
    comps = {facts["reactant"]["map_to_component"].get(int(m)) for m in vals}
    return None not in comps and len(comps) == 1


def _internal_same_element_rewire(
    facts: dict[str, Any], center: int, lost: dict[str, Any], gained: dict[str, Any]
) -> bool:
    """Detect mapped internal attachment permutations that must not become substitution.

    If a conserved center exchanges one mapped neighbor for another mapped neighbor of
    the same element, and all three atoms originate from the same reactant component,
    this is internal skeleton rewiring. It may represent a true rearrangement or an atom-
    mapping permutation, but it is not evidence for synthetic substitution by itself.
    """
    lm, gm = _att_map(lost), _att_map(gained)
    if lm is None or gm is None:
        return False
    if _att_element(lost) != _att_element(gained):
        return False
    if not _same_reactant_origin(facts, center, lm, gm):
        return False
    # Strongest artifact-prone context: conserved aromatic/ring skeleton. The general
    # same-origin/same-element exchange is still non-substitution even outside rings.
    return True


def _oxygen_has_carbonyl_neighbor(
    facts: dict[str, Any], side: str, center: int
) -> bool:
    if _atom(facts, side, center).get("element") != "O":
        return False
    return any(
        att.get("element") == "C" and bool(att.get("carbonyl_like"))
        for att in _attachments(facts, side, center)
    )


def _hydrogen_delta(facts: dict[str, Any], mn: int) -> int:
    r = _atom(facts, "reactant", mn)
    p = _atom(facts, "product", mn)
    if not r or not p:
        return 0
    return int(p.get("total_h", 0)) - int(r.get("total_h", 0))


def _is_nitrile_carbon(facts: dict[str, Any], mn: int) -> bool:
    if _atom(facts, "reactant", mn).get("element") != "C":
        return False
    return any(
        att.get("element") == "N" and float(att.get("bond_order", 0.0)) >= 2.8
        for att in _attachments(facts, "reactant", mn)
    )


def _consumed_pi_bond_at(facts: dict[str, Any], mn: int) -> list[int]:
    partners = []
    for rec in facts.get("bond_order_changes", []):
        if mn not in (rec.get("a"), rec.get("b")):
            continue
        if (
            float(rec.get("reactant_order", 0.0)) >= 1.8
            and float(rec.get("product_order", 0.0)) <= 1.2
        ):
            partners.append(int(rec["b"] if rec["a"] == mn else rec["a"]))
    return sorted(set(partners))



def _cp_product_state_event(facts: dict[str, Any], carbon_map: int, phosphorus_map: int) -> str | None:
    """Classify *product-state evidence* for a newly installed C--P bond.

    This does not claim Arbuzov, Michaelis--Becker, or Wittig mechanism.
    The electrophilic carbon must be sp3 on the reactant side. Unsupported
    P(+)-C(-) resonance representations intentionally fall back to topology.
    """
    if not _is_sp3_carbon(facts, "reactant", carbon_map):
        return None
    if _atom(facts, "product", phosphorus_map).get("element") != "P":
        return None
    pn = _attachments(facts, "product", phosphorus_map)
    cp = [
        a for a in pn
        if a.get("element") == "C" and _att_map(a) == carbon_map
    ]
    if len(cp) != 1:
        return None
    cp_order = float(cp[0].get("bond_order", 0.0))
    # Alkylidenephosphorane: a mapped P=C bond and exactly three other
    # singly bound carbon substituents on P. This is a product-state label,
    # not a claim that a particular ylide-generating base was used.
    if 1.8 <= cp_order <= 2.2:
        other = [a for a in pn if a is not cp[0]]
        if len(other) == 3 and all(
            a.get("element") == "C" and 0.9 <= float(a.get("bond_order", 0)) <= 1.2
            for a in other
        ):
            return "phosphorus_ylide_formation"
    # Alkylphosphonate: one C--P, one P=O and two P--O single bonds.
    # Require the complete tetracoordinate P motif, not just any P--O bond.
    if 0.9 <= cp_order <= 1.2 and len(pn) == 4:
        other = [a for a in pn if a is not cp[0]]
        if len(other) == 3 and all(a.get("element") == "O" for a in other):
            orders = sorted(float(a.get("bond_order", 0.0)) for a in other)
            if (0.9 <= orders[0] <= 1.2 and 0.9 <= orders[1] <= 1.2
                    and 1.8 <= orders[2] <= 2.2):
                return "alkylphosphonate_formation"
    return None

def _formal_carbon_alkylation_at_center(
    facts: dict[str, Any], center: int, lost: dict[str, Any], gained: dict[str, Any]
) -> bool:
    """Conservative structural carbon-alkylation criterion (not a mechanism call).

    Requires a mapped saturated C--X electrophile, a formerly C--H mapped
    carbon nucleophile, direct *intercomponent* C--C bond formation, and no
    ring construction or pi-bond change. Retains formal substitution evidence.
    """
    gm = _att_map(gained)
    if (gm is None or _att_element(lost) not in HALOGENS
            or _att_element(gained) != "C"
            or not _external_gain(facts, center, gained)
            or not _is_sp3_carbon(facts, "reactant", center)
            or not _is_sp3_carbon(facts, "product", center)
            or not _is_sp3_carbon(facts, "reactant", gm)
            or not _is_sp3_carbon(facts, "product", gm)
            or float(gained.get("bond_order", 0)) != 1.0):
        return False
    rdon = _atom(facts, "reactant", gm)
    pdon = _atom(facts, "product", gm)
    if (not rdon or not pdon or
            int(rdon.get("total_h", 0)) <= int(pdon.get("total_h", 0))):
        return False
    # A rearranged skeleton, ring contraction or formal cyclization is not
    # ordinary carbon alkylation merely because it has local C--X exchange.
    rings = facts.get("ring_lineage_v2") or {}
    if any(rings.get(k) for k in ("formations", "openings", "expansions",
                                  "contractions", "rewirings")):
        return False
    if facts.get("bond_order_changes"):
        return False
    # Mapped C--X cleavage is represented as lost attachment even when X is
    # retained elsewhere; disallow additional mapped skeletal bond cleavage.
    for edge in facts.get("broken_bonds", []):
        if (edge.get("a") in (center, gm) or edge.get("b") in (center, gm)):
            return False
    # No previously mapped heavy-atom neighbor of the donor may disappear.
    # This separates H substitution from C--C partner exchange.
    product_neighbors = {(a.get("map"), a.get("element"))
                         for a in _attachments(facts, "product", gm)}
    for att in _attachments(facts, "reactant", gm):
        if att.get("map") is not None and (att.get("map"), att.get("element")) not in product_neighbors:
            return False
    return True


def _sulfonyl_chloride_product_state_event(
    facts: dict[str, Any], center: int, lost: dict[str, Any], gained: dict[str, Any]
) -> str | None:
    """Conservative net product classifier; never infer a sulfonylation mechanism.

    R-S(=O)2-Cl + N -> R-S(=O)2-N: sulfonamide (not sulfamide).
    R-S(=O)2-Cl + O-C -> R-S(=O)2-O-C: sulfonate ester.
    The mapped S--C bond must survive and incoming mapped N/O must originate
    from a different reactant component. S-O-N and N-S(=O)2-N are excluded.
    """
    if _att_element(lost) != "Cl":
        return None
    incoming = _att_map(gained)
    if incoming is None or not _external_gain(facts, center, gained):
        return None
    r_c = {
        int(att["map"])
        for att in _attachments(facts, "reactant", center)
        if att.get("element") == "C" and att.get("map") is not None
    }
    p_c = {
        int(att["map"])
        for att in _attachments(facts, "product", center)
        if att.get("element") == "C" and att.get("map") is not None
    }
    if not (r_c & p_c):
        return None
    incoming_el = _att_element(gained)
    if incoming_el == "N":
        return "sulfonamide_formation"
    if incoming_el == "O" and any(
        att.get("element") == "C"
        and att.get("map") is not None
        and int(att["map"]) != center
        for att in _attachments(facts, "product", int(incoming))
    ):
        return "sulfonate_ester_formation"
    return None


def _acyl_chloride_ester_product_state(facts: dict[str, Any], center: int,
                                       lost: dict[str, Any], gained: dict[str, Any]) -> bool:
    """Net structural O-acylation by an acyl chloride; not a mechanism inference.

    Require mapped conserved carbonyl oxygen; newly attached external alcohol/enol
    O whose pre-existing O-C bond survives; an actual ester product state.
    Exclude anhydrides, intramolecular lactones and non-O substitution.
    """
    if _att_element(lost) != "Cl" or _att_element(gained) != "O":
        return False
    om = _att_map(gained)
    if om is None or not _external_gain(facts, center, gained):
        return False
    if _carbon_fg_state(facts, "reactant", center) != "acyl_halide":
        return False
    if _carbon_fg_state(facts, "product", center) != "ester":
        return False
    # Mapped carbonyl oxygen identity must persist; exclude migrating oxygen maps.
    def double_oxygen_maps(side: str) -> set[int]:
        return {
            int(a["map"]) for a in _attachments(facts, side, center)
            if a.get("element") == "O" and a.get("map") is not None
            and float(a.get("bond_order", 0.0)) >= 1.8
        }
    if not (double_oxygen_maps("reactant") & double_oxygen_maps("product")):
        return False
    if _bond_between(facts, "reactant", center, om) > 0.2:
        return False
    if _ring_formation_touches(facts, center, om):
        return False
    # A retained O-C bond from the oxygen donor is required (alcohol/enol O).
    for att in _attachments(facts, "product", om):
        cm = _att_map(att)
        if (cm is not None and cm != center and att.get("element") == "C"
            and 0.8 <= float(att.get("bond_order", 0.0)) <= 1.2
            and _bond_between(facts, "reactant", om, int(cm)) >= 0.8):
            return True
    return False


def detect_substitution(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Mapped-center exchange detector.

    Requirement: conserved mapped center. Unmapped neighbors may support the lost/gained
    attachment identity, but cannot create an event on their own.
    """
    out: list[dict[str, Any]] = []
    for delta in facts["attachment_deltas"]:
        if len(delta["lost"]) != 1 or len(delta["gained"]) != 1:
            continue
        center = int(delta["center_map"])
        lost, gained = delta["lost"][0], delta["gained"][0]
        center_el = delta["center_element"]
        lost_el, gained_el = _att_element(lost), _att_element(gained)
        lost_map, gained_map = _att_map(lost), _att_map(gained)

        if _internal_same_element_rewire(facts, center, lost, gained):
            out.append(
                {
                    "event": "mapping_rewire_candidate",
                    "role": "topology_evidence",
                    "subtype": f"{center_el}_{lost_el}_to_{gained_el}",
                    "center_map": center,
                    "lost_map": lost_map,
                    "gained_map": gained_map,
                    "support_maps": _support(center, lost_map, gained_map),
                    "mapping_scope": "same_origin_same_element_mapped_attachment_exchange",
                }
            )
            continue

        local = {
            "event": "same_center_replacement",
            "role": "local_event",
            "center_map": center,
            "center_element": center_el,
            "lost_element": lost_el,
            "gained_element": gained_el,
            "lost_map": lost_map,
            "gained_map": gained_map,
            "hydrogen_delta": delta.get("hydrogen_delta", 0),
            "support_maps": _support(center, lost_map, gained_map),
            "mapping_scope": "common_mapped_center_plus_one_sided_attachments",
        }
        out.append(local)

        # donor-transfer bookkeeping should not become substitution
        if (
            lost_el in DONOR_HANDLE_ELEMENTS
            or center_el in DONOR_HANDLE_ELEMENTS
            or center_el in HALOGENS
        ):
            out.append(
                {
                    "event": "transfer_center_replacement",
                    "role": "local_event",
                    "center_map": center,
                    "transfer_from": lost_el,
                    "transfer_to": gained_el,
                    "support_maps": _support(center, lost_map, gained_map),
                    "mapping_scope": "common_mapped_center_plus_one_sided_attachments",
                }
            )
            continue

        if _has_conserved_carbonyl(facts, center):
            acyl_chloride_ester = _acyl_chloride_ester_product_state(
                facts, center, lost, gained
            )
            out.append(
                {
                    "event": "esterification" if acyl_chloride_ester else "carbonyl_substitution",
                    "role": "specific_event",
                    "subtype": f"acyl_{lost_el}_to_{gained_el}",
                    "center_map": center,
                    "support_maps": _support(center, lost_map, gained_map),
                    "mapping_scope": (
                        "conserved_mapped_acyl_chloride_ester_external_oxygen_product_state"
                        if acyl_chloride_ester else "common_mapped_center"
                    ),
                }
            )
            continue

        if _has_conserved_sulfonyl(facts, center):
            sulfonyl_product = _sulfonyl_chloride_product_state_event(
                facts, center, lost, gained
            )
            out.append(
                {
                    "event": sulfonyl_product or "sulfonyl_substitution",
                    "role": "specific_event",
                    "subtype": f"sulfonyl_{lost_el}_to_{gained_el}",
                    "center_map": center,
                    "support_maps": _support(center, lost_map, gained_map),
                    "mapping_scope": (
                        "conserved_mapped_sulfonyl_product_state_external_n_or_o"
                        if sulfonyl_product else "common_mapped_center"
                    ),
                }
            )
            continue

        if center_el == "C":
            cp_product_state = (
                _cp_product_state_event(facts, center, int(gained_map))
                if (
                    gained_el == "P"
                    and gained_map is not None
                    and lost_el in (HALOGENS | {"O", "S"})
                    and _external_gain(facts, center, gained)
                )
                else None
            )
            # Composite structural events should win before generic replacement.
            # Mapping-limited imine formation: a one-sided carbonyl-bearing carbon
            # disappears while a mapped N forms a C=N bond at the conserved center.
            if cp_product_state is not None:
                event = cp_product_state
                subtype = f"C_{lost_el}_to_P"
            elif (
                lost_map is None
                and lost_el == "C"
                and bool(lost.get("carbonyl_like"))
                and gained_el == "N"
                and float(gained.get("bond_order", 0.0)) >= 1.8
            ):
                event = "condensation_like_topology"
                subtype = "mapping_limited_aza_condensation_like"
            # Mapping-limited oxidative alkene cleavage/reductive workup pattern:
            # loss of a one-sided alkene carbon and gain of O at the retained carbon.
            elif (
                lost_map is None
                and lost_el == "C"
                and float(lost.get("bond_order", 0.0)) >= 1.8
                and gained_el == "O"
            ):
                event = "alkene_cleavage_to_alcohol"
                subtype = "alkene_fragment_to_alcohol"
            elif _formal_carbon_alkylation_at_center(facts, center, lost, gained):
                event = "carbon_alkylation"
                subtype = f"C_{lost_el}_to_C"
            elif gained_el == "B":
                event = "borylation"
                subtype = f"C_{lost_el}_to_{gained_el}"
            elif gained_el == "Si":
                event = "silylation"
                subtype = f"C_{lost_el}_to_{gained_el}"
            elif gained_el in HALOGENS:
                event = "halogenation_or_halide_exchange"
                subtype = f"C_{lost_el}_to_{gained_el}"
            elif (
                lost_el == "O"
                and gained_el == "O"
                and _external_gain(facts, center, gained)
            ):
                event = "etherification_or_o_alkylation"
                subtype = "C_O_to_O"
            elif gained_el == "N" and lost_el in (HALOGENS | {"O", "N", "S"}):
                ca = _atom(facts, "reactant", center)
                if ca.get("aromatic"):
                    event = "aryl_or_heteroaryl_amination"
                    subtype = f"C_{lost_el}_to_N"
                elif _is_sp3_carbon(facts, "reactant", center):
                    event = "n_alkylation"
                    subtype = (
                        "n_methylation"
                        if int(ca.get("total_h", 0)) >= 3
                        else f"C_{lost_el}_to_N"
                    )
                else:
                    event = "c_n_substitution"
                    subtype = f"C_{lost_el}_to_N"
            elif gained_el == "O" and lost_el in (HALOGENS | {"N", "S"}):
                ca = _atom(facts, "reactant", center)
                if float(gained.get("bond_order", 0.0)) >= 1.8:
                    event = "heteroaryl_oxo_topology"
                    subtype = f"C_{lost_el}_to_O_double_bond"
                elif ca.get("aromatic"):
                    event = "o_arylation"
                    subtype = f"C_{lost_el}_to_O"
                elif _is_sp3_carbon(facts, "reactant", center):
                    event = "o_alkylation"
                    subtype = (
                        "o_benzylation"
                        if _is_benzylic_carbon(facts, "reactant", center)
                        else f"C_{lost_el}_to_O"
                    )
                else:
                    event = "c_o_substitution"
                    subtype = f"C_{lost_el}_to_O"
            elif gained_el == "S" and lost_el in (HALOGENS | {"O", "N"}):
                event = "c_s_substitution"
                subtype = f"C_{lost_el}_to_S"
            elif lost_el in {"O", "N", "S"} and gained_el in {"N", "O", "S", "C"}:
                event = "substitution"
                subtype = f"C_{lost_el}_to_{gained_el}"
            elif _external_gain(facts, center, gained):
                event = "substitution"
                subtype = f"C_{lost_el}_to_{gained_el}"
            else:
                event = "replacement_candidate"
                subtype = f"C_{lost_el}_to_{gained_el}"
            out.append(
                {
                    "event": event,
                    "role": (
                        "specific_event"
                        if event
                        not in {
                            "replacement_candidate",
                            "condensation_like_topology",
                            "heteroaryl_oxo_topology",
                        }
                        else "topology_evidence"
                    ),
                    "subtype": subtype,
                    "center_map": center,
                    "external_incoming": _external_gain(facts, center, gained),
                    "support_maps": _support(center, lost_map, gained_map),
                    "mapping_scope": "common_mapped_center_plus_one_sided_attachments",
                }
            )
        elif center_el in {"N", "O", "S", "P"}:
            out.append(
                {
                    "event": "heteroatom_attachment_exchange",
                    "role": "local_event",
                    "subtype": f"{center_el}_{lost_el}_to_{gained_el}",
                    "center_map": center,
                    "support_maps": _support(center, lost_map, gained_map),
                    "mapping_scope": "common_mapped_center_plus_one_sided_attachments",
                }
            )
    return out


def _ring_record_support(rec: dict[str, Any]) -> list[int]:
    return sorted({
        int(x)
        for key in (
            "support_maps", "ring_maps", "common_maps", "reactant_ring_maps",
            "product_ring_maps", "conserved_maps",
        )
        for x in (rec.get(key) or [])
        if x is not None
    })


def _ring_small_event(kind: str, direction: str) -> str | None:
    table = {
        ("epoxide", "formation"): "epoxide_formation",
        ("epoxide", "opening"): "epoxide_opening",
        ("aziridine", "formation"): "aziridine_formation",
        ("aziridine", "opening"): "aziridine_opening",
        ("thiirane", "formation"): "thiirane_formation",
        ("thiirane", "opening"): "thiirane_opening",
        ("cyclopropane", "formation"): "cyclopropane_formation",
        ("cyclopropane", "opening"): "cyclopropane_opening",
        ("oxetane", "formation"): "oxetane_formation",
        ("oxetane", "opening"): "oxetane_opening",
        ("azetidine", "formation"): "azetidine_formation",
        ("azetidine", "opening"): "azetidine_opening",
        ("beta_lactam", "formation"): "beta_lactam_formation",
        ("beta_lactam", "opening"): "beta_lactam_opening",
        ("cyclobutane", "formation"): "cyclobutane_formation",
        ("cyclobutane", "opening"): "cyclobutane_opening",
    }
    return table.get((str(kind), str(direction)))


def _lost_aziridine_ring_n_c_bond(facts: dict[str, Any], nitrogen: int, carbon: int | None) -> bool:
    """Veto dealkylation only for a graph-validated aziridine N-C ring scission.

    Both mapped atoms must remain in the product. The ring-lineage record is
    required: the presence of any unrelated aziridine elsewhere is insufficient.
    """
    if carbon is None:
        return False
    if (
        _atom(facts, "reactant", nitrogen).get("element") != "N"
        or _atom(facts, "product", nitrogen).get("element") != "N"
        or _atom(facts, "reactant", carbon).get("element") != "C"
        or _atom(facts, "product", carbon).get("element") != "C"
        or _bond_between(facts, "reactant", nitrogen, carbon) < 0.9
        or _bond_between(facts, "product", nitrogen, carbon) > 0.1
    ):
        return False
    if "ring_lineage_v2" in facts:
        rows = (facts.get("ring_lineage_v2") or {}).get("openings", []) or []
    else:
        rows = facts.get("ring_openings", []) or []
    for rec in rows:
        if str(rec.get("kind") or "") != "aziridine" and str(rec.get("subfamily") or "") != "aziridine_ring_opening":
            continue
        ring_maps = {int(m) for m in (
            rec.get("reactant_ring_maps") or rec.get("ring_maps") or rec.get("support_maps") or rec.get("conserved_maps") or []
        ) if m is not None}
        if {nitrogen, carbon} <= ring_maps:
            return True
    return False


def detect_ring_opening(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Full-lineage ring opening with small-ring specialization.

    The legacy conserved-common-ring facts remain available for QC/backward inspection,
    but semantic ownership uses the insertion/deletion-aware v2 lineage when present.
    """
    out: list[dict[str, Any]] = []
    lineage_present = "ring_lineage_v2" in facts
    lineage = facts.get("ring_lineage_v2") or {}
    rows = lineage.get("openings", []) or []
    if not lineage_present:
        # Compatibility fallback for old facts/fixtures.
        rows = [
            {
                "ring_maps": rec.get("ring_maps", []),
                "common_maps": rec.get("conserved_maps", []),
                "support_maps": rec.get("conserved_maps", []),
                "kind": {
                    "epoxide_ring_opening": "epoxide",
                    "aziridine_ring_opening": "aziridine",
                    "thiirane_ring_opening": "thiirane",
                }.get(rec.get("subfamily"), "ring"),
                "subfamily": rec.get("subfamily"),
                "mapping_scope": "conserved_mapped_ring_lineage_fallback",
            }
            for rec in facts.get("ring_openings", []) or []
        ]

    whole_delta = int(lineage.get("whole_cycle_rank_delta", facts.get("common_cycle_rank_delta", 0)))
    for rec in rows:
        support = _ring_record_support(rec)
        out.append(
            {
                "event": "ring_opening",
                "role": "specific_event" if whole_delta < 0 else "topology_evidence",
                "subtype": rec.get("subfamily") or "validated_ring_opening",
                "ring_maps": rec.get("ring_maps", support),
                "support_maps": support,
                "ring_size": rec.get("ring_size"),
                "ring_kind": rec.get("kind"),
                "mapping_scope": rec.get("mapping_scope", "full_ring_lineage_v2"),
            }
        )
        child = _ring_small_event(str(rec.get("kind") or ""), "opening")
        if child:
            out.append(
                {
                    "event": child,
                    "role": "specific_event",
                    "subtype": child,
                    "ring_maps": rec.get("ring_maps", support),
                    "support_maps": support,
                    "ring_size": rec.get("ring_size"),
                    "ring_kind": rec.get("kind"),
                    "mapping_scope": "full_ring_lineage_v2_small_ring_identity",
                }
            )
    return out


def _ring_v031_carbonyl_transfer_heterocyclization(
    facts: dict[str, Any], rec: dict[str, Any], ring_maps: list[int],
    size: int, whole_delta: int, is_hetero: bool,
) -> bool:
    """Conservative P6-preserving recovery of CDI-type product-ring closures.

    Does NOT infer formation from negative rank alone.  Requires a v2
    formation, two new mapped ring bonds, a 5/6-membered heterocycle,
    and a retained carbonyl-derived donor carbon from a different reactant
    component. Internal mapped bond cleavage excludes ring-expansion cases.
    Azole exchange, same-core rewiring and intramolecular
    single-bond closures are not recovered via this rule.
    """
    if (whole_delta >= 0 or not is_hetero or size not in (5, 6)
            or len(rec.get("formed_ring_bonds") or []) != 2):
        return False
    try:
        formed = [tuple(map(int, pair)) for pair in rec["formed_ring_bonds"]]
        if not all(len(x) == 2 and x[0] != x[1] for x in formed):
            return False
        ring = {int(m) for m in ring_maps}
        if len(ring) != size or not all(set(x) <= ring for x in formed):
            return False
        # A cyclopropane/other retained-core C--C cleavage followed by
        # acylative ring expansion is not the reagent-ring-loss closure
        # targeted here.  The CDI-type positives construct new ring bonds
        # without severing any bond between atoms of that product ring.
        # Keep this independent of overall cycle-rank and reagent names.
        if any(
            int(edge["a"]) in ring and int(edge["b"]) in ring
            for edge in (facts.get("broken_bonds") or [])
        ):
            return False
        origins = facts["reactant"]["map_to_component"]
        components = {origins.get(m) for m in ring}
        if None in components or len(components) != 2:
            return False
        for mn in ring:
            if (_atom(facts, "reactant", mn).get("element") != "C"
                    or _atom(facts, "product", mn).get("element") != "C"):
                continue
            if not any(mn in b for b in formed):
                continue
            if not any(
                a.get("element") == "O" and float(a.get("bond_order", 0)) >= 1.9
                for a in _attachments(facts, "reactant", mn)
            ):
                continue
            if not any(a.get("element") == "O" for a in _attachments(facts, "product", mn)):
                continue
            if any(origins.get(other) != origins.get(mn) for other in ring):
                return True
    except (TypeError, ValueError, KeyError, AttributeError):
        return False
    return False


def detect_ring_formation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Full-lineage ring construction and architecture specialization.

    Generic ring formation is retained as fallback evidence, while structurally
    identifiable small rings and annulation architectures emit more informative children.
    """
    out: list[dict[str, Any]] = []
    lineage_present = "ring_lineage_v2" in facts
    lineage = facts.get("ring_lineage_v2") or {}
    rows = lineage.get("formations", []) or []
    if not lineage_present:
        rows = [
            {
                "ring_maps": rec.get("ring_maps", []),
                "common_maps": rec.get("conserved_maps", []),
                "support_maps": rec.get("conserved_maps", []),
                "kind": "ring",
                "ring_size": len(rec.get("ring_maps", []) or []),
                "mapping_scope": "conserved_mapped_ring_lineage_fallback",
            }
            for rec in facts.get("ring_formations", []) or []
        ]

    whole_delta = int(lineage.get("whole_cycle_rank_delta", facts.get("common_cycle_rank_delta", 0)))
    has_component_split = bool(facts.get("component_splits"))

    # v0.3: mapped-core fusion annulations are deliberately separate from generic
    # ring formation.  They recover a narrow Pschorr/Scholl-like mapping regime in
    # which the conserved mapped core gains a cycle although the full molecular
    # cycle rank is unchanged because reactant ring atoms were left unmapped.
    for rec in lineage.get("mapped_core_annulations", []) or []:
        support = _ring_record_support(rec)
        kind = str(rec.get("kind") or "ring")
        is_hetero = kind in {
            "epoxide", "aziridine", "thiirane", "oxetane", "azetidine",
            "beta_lactam", "heterocycle",
        }
        out.append(
            {
                "event": "heterocycle_annulation" if is_hetero else "carbocycle_annulation",
                "role": "specific_event",
                "subtype": "mapped_core_fusion_annulation",
                "ring_maps": rec.get("ring_maps", support),
                "support_maps": support,
                "ring_size": rec.get("ring_size"),
                "mapping_scope": rec.get(
                    "mapping_scope",
                    "mapped_core_fusion_annulation_zero_full_rank_change",
                ),
            }
        )

    for rec in rows:
        support = _ring_record_support(rec)
        ring_maps = list(rec.get("ring_maps", support) or support)
        kind = str(rec.get("kind") or "ring")
        size = int(rec.get("ring_size") or len(ring_maps) or 0)
        out.append(
            {
                "event": "ring_formation",
                "role": "specific_event" if whole_delta > 0 else "topology_evidence",
                "subtype": "validated_ring_formation",
                "ring_maps": ring_maps,
                "support_maps": support,
                "ring_size": size,
                "ring_kind": kind,
                "mapping_scope": rec.get("mapping_scope", "full_ring_lineage_v2"),
            }
        )

        child = _ring_small_event(kind, "formation")
        if child:
            out.append(
                {
                    "event": child,
                    "role": "specific_event",
                    "subtype": child,
                    "ring_maps": ring_maps,
                    "support_maps": support,
                    "ring_size": size,
                    "ring_kind": kind,
                    "mapping_scope": "full_ring_lineage_v2_small_ring_identity",
                }
            )

        if size >= 10:
            out.append(
                {
                    "event": "macrocyclization",
                    "role": "specific_event",
                    "subtype": f"{size}_member_ring_formation",
                    "ring_maps": ring_maps,
                    "support_maps": support,
                    "ring_size": size,
                    "mapping_scope": "full_ring_lineage_v2_large_ring_formation",
                }
            )

        fused = bool(rec.get("fused_to_retained_ring"))
        spiro = bool(rec.get("spiro_to_retained_ring"))
        opening_supports = [
            set(x.get("ring_maps", []) or [])
            for x in (lineage.get("openings", []) or [])
        ]
        zero_rank_replacement_annulation = (
            whole_delta == 0
            and not has_component_split
            and len(rec.get("formed_ring_bonds", []) or []) >= 2
            and any(len(set(ring_maps) & x) >= 2 for x in opening_supports)
        )
        reactant_origins = {
            facts["reactant"]["map_to_component"].get(int(m))
            for m in (rec.get("common_maps") or [])
            if m is not None
        }
        reactant_origins.discard(None)
        is_hetero = kind in {
            "epoxide", "aziridine", "thiirane", "oxetane", "azetidine",
            "beta_lactam", "heterocycle",
        }
        # Ring v0.3.1: local heterocycle construction may coexist with net
        # ring loss when cyclic carbonyl-transfer reagents depart.
        negative_rank_carbonyl_closure = _ring_v031_carbonyl_transfer_heterocyclization(
            facts, rec, ring_maps, size, whole_delta, is_hetero,
        )
        multicomponent_multibond_construction = (
            not has_component_split
            and whole_delta >= 0
            and size >= 5
            and len(reactant_origins) >= 2
            and len(rec.get("formed_ring_bonds", []) or []) >= 2
        )
        # Carbocyclic multicomponent/multibond construction is structurally
        # annulative (Dötz/Pauson-Khand/Danheiser-like).  For heterocycles, keep
        # the v0.2 distinction between general heterocycle formation and annulation;
        # only recover the zero-rank multicomponent replacement regime where a new
        # heterocycle is proven but a leaving ring keeps whole cycle rank unchanged
        # (Kröhnke-like mapped architecture).
        intercomponent_annulation = bool(
            multicomponent_multibond_construction
            and (not is_hetero or whole_delta == 0)
        )
        architecture_primary = (
            whole_delta > 0 or (whole_delta == 0 and not has_component_split)
            or negative_rank_carbonyl_closure
        )
        retained_scaffold_polycycle = bool(rec.get("retained_scaffold_polycyclic_closure"))
        if child is None and spiro and architecture_primary:
            out.append(
                {
                    "event": "spirocyclization",
                    "role": "specific_event",
                    "subtype": "heterospirocyclization" if is_hetero else "carbospirocyclization",
                    "ring_maps": ring_maps,
                    "support_maps": support,
                    "ring_size": size,
                    "mapping_scope": "new_ring_shares_one_atom_with_retained_ring",
                }
            )
        elif child is None and (fused or zero_rank_replacement_annulation) and architecture_primary:
            out.append(
                {
                    "event": "heterocycle_annulation" if is_hetero else "carbocycle_annulation",
                    "role": "specific_event",
                    "subtype": (
                        "zero_rank_ring_replacement_annulation"
                        if zero_rank_replacement_annulation and not fused
                        else "fused_ring_annulation"
                    ),
                    "ring_maps": ring_maps,
                    "support_maps": support,
                    "ring_size": size,
                    "mapping_scope": (
                        "zero_rank_open_plus_multibond_ring_construction_without_component_split"
                        if zero_rank_replacement_annulation and not fused
                        else "new_ring_shares_mapped_edge_with_retained_ring"
                    ),
                }
            )
        elif child is None and retained_scaffold_polycycle and architecture_primary:
            out.append(
                {
                    "event": "carbocycle_annulation",
                    "role": "specific_event",
                    "subtype": "retained_scaffold_polycyclic_annulation",
                    "ring_maps": ring_maps,
                    "support_maps": support,
                    "ring_size": size,
                    "mapping_scope": rec.get(
                        "mapping_scope",
                        "full_ring_lineage_retained_scaffold_polycyclic_recovery",
                    ),
                }
            )
        elif child is None and intercomponent_annulation and architecture_primary:
            out.append(
                {
                    "event": "heterocycle_annulation" if is_hetero else "carbocycle_annulation",
                    "role": "specific_event",
                    "subtype": (
                        "intercomponent_heterocycle_annulation"
                        if is_hetero
                        else "intercomponent_carbocycle_annulation"
                    ),
                    "ring_maps": ring_maps,
                    "support_maps": support,
                    "ring_size": size,
                    "mapping_scope": "new_ring_joins_multiple_reactant_components_with_multibond_construction",
                }
            )
        elif is_hetero and child is None and (whole_delta > 0 or negative_rank_carbonyl_closure):
            formed_inside = [tuple(x) for x in rec.get("formed_ring_bonds", []) or []]
            intra = [
                b for b in formed_inside
                if _same_reactant_origin(facts, int(b[0]), int(b[1]))
            ]
            if len(formed_inside) == 1 and intra:
                out.append(
                    {
                        "event": "intramolecular_heterocyclization",
                        "role": "specific_event",
                        "subtype": "single_bond_intramolecular_heterocyclization",
                        "ring_maps": ring_maps,
                        "support_maps": support,
                        "ring_size": size,
                        "mapping_scope": "full_ring_lineage_v2_single_intracomponent_closure",
                    }
                )
            else:
                inserted_atom_closure = bool(rec.get("single_inserted_ring_atom_closure"))
                out.append(
                    {
                        "event": "heterocycle_formation",
                        "role": "specific_event",
                        "subtype": (
                            "single_inserted_ring_atom_heterocycle_formation"
                            if inserted_atom_closure
                            else "nonannulative_heterocycle_formation"
                        ),
                        "ring_maps": ring_maps,
                        "support_maps": support,
                        "ring_size": size,
                        "mapping_scope": (
                            rec.get("mapping_scope")
                            if inserted_atom_closure
                            else "full_ring_lineage_v2_heterocycle_formation"
                        ),
                    }
                )

    return out


def detect_ring_size_change(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    lineage = facts.get("ring_lineage_v2") or {}
    for event, key in (("ring_expansion", "expansions"), ("ring_contraction", "contractions")):
        for rec in lineage.get(key, []) or []:
            out.append(
                {
                    "event": event,
                    "role": "specific_event",
                    "subtype": f"{rec.get('reactant_ring_size')}_to_{rec.get('product_ring_size')}_member_ring",
                    "reactant_ring_size": rec.get("reactant_ring_size"),
                    "product_ring_size": rec.get("product_ring_size"),
                    "reactant_ring_maps": rec.get("reactant_ring_maps", []),
                    "product_ring_maps": rec.get("product_ring_maps", []),
                    "support_maps": _ring_record_support(rec),
                    "reactant_ring_kind": rec.get("reactant_kind"),
                    "product_ring_kind": rec.get("product_kind"),
                    "mapping_scope": "full_ring_lineage_v2_matched_ring_size_change",
                }
            )
    return out


def detect_fragmentation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Separate local fragment loss from true fragmentation.

    A disappearing unmapped branch is anchored to a conserved mapped center and may
    generate a specific local event, but generic fragmentation is primary only when
    mapped component splitting or stronger mapped multi-bond cleavage is present.
    """
    out: list[dict[str, Any]] = []
    if facts["component_splits"]:
        for split in facts["component_splits"]:
            out.append(
                {
                    "event": "fragmentation",
                    "role": "specific_event",
                    "subtype": "mapped_component_split",
                    "support_maps": split["maps"],
                    "mapping_scope": "common_mapped_component_lineage",
                }
            )

    if len(facts["broken_bonds"]) >= 2 and facts["component_splits"]:
        maps = []
        for b in facts["broken_bonds"]:
            maps.extend([b["a"], b["b"]])
        out.append(
            {
                "event": "multi_bond_fragmentation",
                "role": "specific_event",
                "support_maps": sorted(set(maps)),
                "mapping_scope": "common_mapped_bond_cleavage",
            }
        )

    for delta in facts["attachment_deltas"]:
        if not delta["lost"] or delta["gained"]:
            continue
        if int(delta.get("hydrogen_delta", 0)) <= 0:
            continue
        center = int(delta["center_map"])
        center_el = delta["center_element"]
        for lost in delta["lost"]:
            lost_el = _att_element(lost)
            subtype = "fragment_loss_to_h"
            specific = None
            if center_el == "O" and lost_el == "C":
                # S2_2A_ORIENTED_O_DEACYLATION
                #
                # Distinguish the two orientations of ester cleavage around
                # the conserved oxygen.
                #
                # O-deacylation / scaffold-side cleavage:
                #
                #     R-C(=O)-O*-scaffold -> HO*-scaffold
                #
                # The carbon LOST from O* is itself carbonyl-like, and the
                # product oxygen no longer has a carbonyl neighbor.
                #
                # Acid-side ester hydrolysis / de-esterification:
                #
                #     R-C(=O)-O*-R' -> R-C(=O)-O*H
                #
                # The carbonyl neighbor is retained on the conserved oxygen;
                # the lost carbon substituent is not itself the acyl carbon.
                #
                # Previous behavior tested only whether the reactant oxygen
                # had any carbonyl neighbor. That is also true for
                # O-deacylation, so the broad hydrolysis branch shadowed the
                # more specific lost-acyl branch.

                reactant_o_has_carbonyl = _oxygen_has_carbonyl_neighbor(
                    facts, "reactant", center
                )

                product_o_has_carbonyl = _oxygen_has_carbonyl_neighbor(
                    facts, "product", center
                )

                if (
                    bool(lost.get("carbonyl_like"))
                    and reactant_o_has_carbonyl
                    and not product_o_has_carbonyl
                ):
                    specific = "o_deacylation_or_ester_cleavage"

                elif reactant_o_has_carbonyl and product_o_has_carbonyl:
                    specific = "ester_hydrolysis_or_deesterification"

                else:
                    # A retained O-B bond denotes boronate-ester O-alkyl loss,
                    # not cleavage of an ordinary organic ether linkage.
                    if not _retained_mapped_oxygen_boron_bond(facts, center):
                        specific = "o_dealkylation"
            elif center_el == "C" and lost_el == "C" and lost.get("carboxyl_like"):
                specific = "decarboxylation"
            elif center_el == "C" and lost_el in HALOGENS:
                specific = "dehalogenation"
            elif center_el == "C" and lost_el == "S":
                specific = "desulfurization"
            elif (
                center_el == "N" and lost_el == "S" and bool(lost.get("sulfonyl_like"))
            ):
                specific = "n_desulfonylation"
            elif (
                center_el == "N"
                and lost_el == "C"
                and bool(lost.get("carbonyl_like"))
                and int(lost.get("degree", 0)) <= 2
            ):
                specific = "n_deformylation"
            elif (
                center_el == "C"
                and lost_el == "O"
                and float(lost.get("bond_order", 0)) >= 1.9
            ):
                specific = "carbonyl_deoxygenation_or_reduction"
            elif (
                center_el == "N"
                and lost_el == "O"
                and float(lost.get("bond_order", 0)) >= 1.9
            ):
                specific = "n_oxide_or_nitroso_reduction"

            out.append(
                {
                    "event": subtype,
                    "role": "local_event",
                    "center_map": center,
                    "lost_element": lost_el,
                    "hydrogen_delta": delta.get("hydrogen_delta", 0),
                    "support_maps": _support(center, _att_map(lost)),
                    "mapping_scope": "common_mapped_center_plus_one_sided_attachment",
                }
            )
            if specific:
                role = (
                    "topology_evidence"
                    if specific == "ester_hydrolysis_or_deesterification"
                    else "specific_event"
                )
                out.append(
                    {
                        "event": specific,
                        "role": role,
                        "center_map": center,
                        "support_maps": _support(center, _att_map(lost)),
                        "mapping_scope": (
                            "supporting_oxygen_side_evidence_for_acyl_center_conversion"
                            if specific == "ester_hydrolysis_or_deesterification"
                            else "common_mapped_center_plus_one_sided_attachment"
                        ),
                    }
                )
    return out


def _endpoint_has_handle(facts: dict[str, Any], mn: int, elements: set[str]) -> bool:
    return _has_neighbor_element(facts, "reactant", mn, elements)


def _acyl_cc_replacement_center(
    facts: dict[str, Any], center: int, incoming: int
) -> bool:
    """C-C formation at a conserved acyl carbon with simultaneous N/O/X displacement."""
    if not _has_conserved_carbonyl(facts, center):
        return False
    lost, gained = _lost_gained_for_center(facts, center)
    if not any(_att_map(g) == incoming and g.get("element") == "C" for g in gained):
        return False
    return any(x.get("element") in ({"N", "O"} | HALOGENS) for x in lost)


def _consumed_carbonyl_c(facts: dict[str, Any], center: int) -> bool:
    """Reacting carbonyl carbon loses its C=O; O is not a coupling leaving group.

    Evaluate the *same mapped bond-forming carbon*, not the entire reaction or
    another reactant. A persistent product carbonyl is deliberately not excluded
    (e.g. acyl substitutions have their own separate, pre-existing guard).
    """
    if _atom(facts, "reactant", center).get("element") != "C":
        return False
    return any(
        att.get("element") == "O" and float(att.get("bond_order", 0)) >= 1.9
        for att in _attachments(facts, "reactant", center)
    ) and not _has_conserved_carbonyl(facts, center)




def _formal_vinyl_alkynyl_coupling(
    facts: dict[str, Any], a: int, b: int
) -> tuple[int, int, int, int] | None:
    """Net terminal alkynyl-C(sp) + vinyl-C(sp2)-X bond formation.

    Describe preserved vinyl C=C and alkyne C#C, not a Sonogashira
    mechanism; the reaction record may lack catalytic evidence.
    """
    if (
        len(facts.get("formed_bonds") or []) != 1
        or len(facts.get("intercomponent_formed_bonds") or []) != 1
        or bool(facts.get("broken_bonds"))
        or bool(facts.get("bond_order_changes"))
    ):
        return None
    lineage = facts.get("ring_lineage_v2") or {}
    if any(lineage.get(k) for k in
           ("formations", "openings", "expansions", "contractions", "rewirings")):
        return None
    for vinyl, alkynyl in ((a, b), (b, a)):
        v_r = _atom(facts, "reactant", vinyl)
        v_p = _atom(facts, "product", vinyl)
        a_r = _atom(facts, "reactant", alkynyl)
        a_p = _atom(facts, "product", alkynyl)
        if not (
            v_r.get("element") == v_p.get("element") == "C"
            and a_r.get("element") == a_p.get("element") == "C"
            and not v_r.get("aromatic") and not v_p.get("aromatic")
            and v_r.get("hybridization") == v_p.get("hybridization") == "SP2"
            and a_r.get("hybridization") == a_p.get("hybridization") == "SP"
            and _endpoint_has_handle(facts, vinyl, HALOGENS)
            and _hydrogen_delta(facts, vinyl) == 0
            and int(a_r.get("total_h", 0)) == 1
            and _hydrogen_delta(facts, alkynyl) == -1
            and 0.9 <= _bond_between(facts, "product", vinyl, alkynyl) <= 1.2
        ):
            continue
        v_partners = [int(x["map"]) for x in _attachments(facts, "reactant", vinyl)
                      if x.get("map") is not None and x.get("element") == "C"
                      and 1.8 <= float(x.get("bond_order", 0)) <= 2.2
                      and 1.8 <= _bond_between(facts, "product", vinyl, int(x["map"])) <= 2.2]
        a_partners = [int(x["map"]) for x in _attachments(facts, "reactant", alkynyl)
                      if x.get("map") is not None and x.get("element") == "C"
                      and 2.8 <= float(x.get("bond_order", 0)) <= 3.2
                      and 2.8 <= _bond_between(facts, "product", alkynyl, int(x["map"])) <= 3.2]
        if len(v_partners) == 1 and len(a_partners) == 1:
            return vinyl, alkynyl, v_partners[0], a_partners[0]
    return None


def _formal_aryl_alkene_hydroarylation(
    facts: dict[str, Any], a: int, b: int
) -> tuple[int, int, int] | None:
    """Mapped net Ar--X + C=C -> Ar--C--C(H) without a mechanism claim.

    Only a single intermolecular C--C addition with one consumed C=C is
    admissible. Enol-ether-to-ketone couplings, arene C--H alkylations, Heck
    alkenylations retaining C=C, and cyclizations are excluded.
    """
    if (
        len(facts.get("formed_bonds") or []) != 1
        or len(facts.get("intercomponent_formed_bonds") or []) != 1
        or bool(facts.get("broken_bonds"))
        or len(facts.get("bond_order_changes") or []) != 1
    ):
        return None
    lineage = facts.get("ring_lineage_v2") or {}
    if any(lineage.get(k) for k in (
        "formations", "openings", "expansions", "contractions", "rewirings"
    )):
        return None
    change = facts["bond_order_changes"][0]
    if not (1.8 <= float(change.get("reactant_order", 0.0)) <= 2.2
            and 0.9 <= float(change.get("product_order", 0.0)) <= 1.2):
        return None
    for aryl, alkene in ((a, b), (b, a)):
        ar = _atom(facts, "reactant", aryl)
        ap = _atom(facts, "product", aryl)
        rn = _atom(facts, "reactant", alkene)
        if not (ar.get("element") == ap.get("element") == "C"
                and bool(ar.get("aromatic")) and bool(ap.get("aromatic"))
                and rn.get("element") == "C"
                and not rn.get("aromatic")
                and _is_sp3_carbon(facts, "product", alkene)
                and _endpoint_has_handle(facts, aryl, HALOGENS | {"O"})
                and int(rn.get("total_h", 0)) >= 2
                and _hydrogen_delta(facts, alkene) == 0):
            continue
        other = None
        if int(change.get("a", -1)) == alkene:
            other = int(change["b"])
        elif int(change.get("b", -1)) == alkene:
            other = int(change["a"])
        if other is None or not _is_sp3_carbon(facts, "product", other):
            continue
        if (
            _bond_between(facts, "reactant", alkene, other) >= 1.8
            and 0.9 <= _bond_between(facts, "product", alkene, other) <= 1.2
            and _hydrogen_delta(facts, other) == 1
            and 0.9 <= _bond_between(facts, "product", aryl, alkene) <= 1.2
        ):
            return aryl, alkene, other
    return None

def detect_coupling(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Intercomponent bond formation is a primitive, not automatically a family.

    Exact bond class is always recorded from mapped endpoints. Cross-coupling-compatible
    is promoted only when a C-C intercomponent bond has recognizable electrophile and
    donor handles on opposite mapped endpoints.
    """
    out: list[dict[str, Any]] = []
    for b in facts["intercomponent_formed_bonds"]:
        a, c = int(b["a"]), int(b["b"])
        ea = _atom(facts, "reactant", a).get("element", "?")
        ec = _atom(facts, "reactant", c).get("element", "?")
        pair = canonical_bond_class(ea, ec)
        out.append(
            {
                "event": "intercomponent_bond_formation",
                "role": "primitive",
                "formed_bond_class": pair,
                "a": a,
                "b": c,
                "support_maps": _support(a, c),
                "mapping_scope": "formed_bond_between_mapped_reactant_origin_components",
            }
        )

        if ea == "C" and ec == "C":
            vinyl_alkynyl = _formal_vinyl_alkynyl_coupling(facts, a, c)
            if vinyl_alkynyl is not None:
                vinyl, alkynyl, v_partner, a_partner = vinyl_alkynyl
                out.append({
                    "event": "vinyl_alkynyl_coupling",
                    "role": "specific_event",
                    "subtype": "vinyl_x_to_terminal_alkynyl_c_c_preserved_pi",
                    "center_map": vinyl,
                    "alkynyl_map": alkynyl,
                    "vinyl_pi_partner_map": v_partner,
                    "alkyne_pi_partner_map": a_partner,
                    "support_maps": _support(vinyl, alkynyl, v_partner, a_partner),
                    "mapping_scope": "intercomponent_vinyl_c_x_to_terminal_alkyne_c_h_preserved_c_c_pi_bonds",
                })
                continue
            formal_hydroarylation = _formal_aryl_alkene_hydroarylation(facts, a, c)
            if formal_hydroarylation is not None:
                aryl, alkene, partner = formal_hydroarylation
                out.append({
                    "event": "aryl_alkene_hydroarylation",
                    "role": "specific_event",
                    "subtype": "aryl_x_to_alkyl_c_c_with_pi_reduction_and_partner_h_gain",
                    "center_map": aryl,
                    "alkene_carbon_map": alkene,
                    "alkene_partner_map": partner,
                    "support_maps": _support(aryl, alkene, partner),
                    "mapping_scope": "mapped_aryl_c_x_replacement_and_alkene_pi_consumption_with_h_gain",
                })
                continue
            # Nitrile-carbon transfer is more informative than generic C-C formation.
            # It is intentionally structural (nitrile installation), not a named-mechanism call.
            nitrile_a = _is_nitrile_carbon(facts, a)
            nitrile_c = _is_nitrile_carbon(facts, c)
            if nitrile_a ^ nitrile_c:
                center = c if nitrile_a else a
                donor = a if nitrile_a else c
                out.append(
                    {
                        "event": "cyanation_or_nitrile_installation",
                        "role": "specific_event",
                        "subtype": "c_cn_installation",
                        "center_map": center,
                        "nitrile_carbon_map": donor,
                        "support_maps": _support(center, donor),
                        "mapping_scope": "mapped_intercomponent_c_c_bond_to_nitrile_carbon",
                    }
                )
            else:
                # A carbonyl C=O consumed during C-C bond construction is not
                # an electrophilic C-O leaving group for substitutive coupling.
                # Keep aryl/vinyl O-sulfonate handles and conserved-acyl guards.
                a_e = _endpoint_has_handle(
                    facts, a, LEAVING_GROUP_ELEMENTS | {"O"}
                ) and not _consumed_carbonyl_c(facts, a)
                c_e = _endpoint_has_handle(
                    facts, c, LEAVING_GROUP_ELEMENTS | {"O"}
                ) and not _consumed_carbonyl_c(facts, c)
                a_d = _endpoint_has_handle(facts, a, DONOR_HANDLE_ELEMENTS)
                c_d = _endpoint_has_handle(facts, c, DONOR_HANDLE_ELEMENTS)
                a_b = _endpoint_has_handle(facts, a, {"B"})
                c_b = _endpoint_has_handle(facts, c, {"B"})
                acyl_replacement = _acyl_cc_replacement_center(
                    facts, a, c
                ) or _acyl_cc_replacement_center(facts, c, a)

                # Alpha arylation: aryl/heteroaryl C-X couples to a carbon alpha to
                # carbonyl that loses H. This is not a generic C-Cl -> C substitution.
                a_alpha = (
                    _is_alpha_to_carbonyl(facts, "reactant", a)
                    and _hydrogen_delta(facts, a) < 0
                )
                c_alpha = (
                    _is_alpha_to_carbonyl(facts, "reactant", c)
                    and _hydrogen_delta(facts, c) < 0
                )
                a_aryl_x = bool(
                    _atom(facts, "reactant", a).get("aromatic")
                ) and _endpoint_has_handle(facts, a, LEAVING_GROUP_ELEMENTS)
                c_aryl_x = bool(
                    _atom(facts, "reactant", c).get("aromatic")
                ) and _endpoint_has_handle(facts, c, LEAVING_GROUP_ELEMENTS)
                if (a_alpha and c_aryl_x) or (c_alpha and a_aryl_x):
                    alpha, aryl = (a, c) if a_alpha else (c, a)
                    out.append(
                        {
                            "event": "carbonyl_alpha_arylation",
                            "role": "specific_event",
                            "subtype": "alpha_c_h_to_aryl_c_c",
                            "alpha_map": alpha,
                            "aryl_map": aryl,
                            "center_map": alpha,
                            "support_maps": _support(alpha, aryl),
                            "mapping_scope": "mapped_alpha_carbon_h_loss_plus_aryl_c_x_replacement_and_c_c_formation",
                        }
                    )
                elif ((a_e and c_b) or (c_e and a_b)) and not acyl_replacement:
                    electrophile, donor = (a, c) if a_e and c_b else (c, a)
                    out.append(
                        {
                            "event": "suzuki_type_cross_coupling",
                            "role": "specific_event",
                            "subtype": "organoboron_c_c_cross_coupling",
                            "center_map": electrophile,
                            "donor_map": donor,
                            "support_maps": _support(a, c),
                            "mapping_scope": "mapped_c_b_donor_plus_electrophilic_c_x_or_c_o_loss_and_intercomponent_c_c_formation",
                        }
                    )
                elif ((a_e and c_d) or (c_e and a_d)) and not acyl_replacement:
                    out.append(
                        {
                            "event": "cross_coupling_compatible",
                            "role": "specific_event",
                            "subtype": "c_c_cross_coupling",
                            "support_maps": _support(a, c),
                            "mapping_scope": "mapped_intercomponent_endpoints_plus_mapped_or_one_sided_handles",
                        }
                    )
                else:
                    out.append(
                        {
                            "event": "c_c_intercomponent_formation",
                            "role": "topology_evidence",
                            "support_maps": _support(a, c),
                            "mapping_scope": "mapped_intercomponent_endpoints",
                        }
                    )
        elif {ea, ec} & HALOGENS and "C" in {ea, ec}:
            carbon = a if ea == "C" else c
            halogen = c if ea == "C" else a
            x = ec if ea == "C" else ea
            pi_partners = _consumed_pi_bond_at(facts, carbon)
            if pi_partners:
                out.append(
                    {
                        "event": "addition_like_topology",
                        "role": "topology_evidence",
                        "subtype": "c_x_formation_at_consumed_pi_system",
                        "center_map": carbon,
                        "halogen": x,
                        "pi_partner_maps": pi_partners,
                        "support_maps": _support(carbon, halogen, *pi_partners),
                        "mapping_scope": "mapped_c_x_formation_plus_pi_bond_consumption_not_yet_validated_addition",
                    }
                )
            elif _hydrogen_delta(facts, carbon) < 0:
                out.append(
                    {
                        "event": "c_h_halogenation",
                        "role": "specific_event",
                        "subtype": f"C_H_to_{x}",
                        "center_map": carbon,
                        "halogen": x,
                        "support_maps": _support(carbon, halogen),
                        "mapping_scope": "mapped_c_x_formation_plus_hydrogen_loss",
                    }
                )
            else:
                out.append(
                    {
                        "event": f'{pair.lower().replace("-", "_")}_intercomponent_formation',
                        "role": "topology_evidence",
                        "support_maps": _support(a, c),
                        "mapping_scope": "mapped_intercomponent_endpoints",
                    }
                )
        else:
            out.append(
                {
                    "event": f'{pair.lower().replace("-", "_")}_intercomponent_formation',
                    "role": "topology_evidence",
                    "support_maps": _support(a, c),
                    "mapping_scope": "mapped_intercomponent_endpoints",
                }
            )
    return out


def detect_ring_rewiring(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Promote validated same-size ring lineage change; retain ambiguous legacy QC."""
    out: list[dict[str, Any]] = []
    lineage = facts.get("ring_lineage_v2") or {}

    for rec in lineage.get("rewirings", []) or []:
        out.append(
            {
                "event": "ring_rewiring",
                "role": "specific_event",
                "subtype": "validated_same_size_ring_lineage_change",
                "reactant_ring_maps": rec.get("reactant_ring_maps", []),
                "product_ring_maps": rec.get("product_ring_maps", []),
                "support_maps": _ring_record_support(rec),
                "reactant_ring_kind": rec.get("reactant_kind"),
                "product_ring_kind": rec.get("product_kind"),
                "mapping_scope": "full_ring_lineage_v2_same_size_membership_change",
            }
        )

    # Payne-type epoxide migration can share only one atom between old/new epoxide.
    # Treat one opening + one formation on the same mapped component and zero whole
    # cycle-rank change as a chemically specific epoxide rearrangement.
    ep_open = [x for x in lineage.get("openings", []) or [] if x.get("kind") == "epoxide"]
    ep_form = [x for x in lineage.get("formations", []) or [] if x.get("kind") == "epoxide"]
    if int(lineage.get("whole_cycle_rank_delta", facts.get("common_cycle_rank_delta", 0))) == 0:
        for old in ep_open:
            for new in ep_form:
                old_maps = set(old.get("ring_maps", []) or [])
                new_maps = set(new.get("ring_maps", []) or [])
                common = old_maps & new_maps
                if not common:
                    continue
                support = sorted(old_maps | new_maps)
                out.append(
                    {
                        "event": "epoxide_rearrangement",
                        "role": "specific_event",
                        "subtype": "epoxide_ring_migration",
                        "reactant_ring_maps": sorted(old_maps),
                        "product_ring_maps": sorted(new_maps),
                        "support_maps": support,
                        "mapping_scope": "epoxide_open_plus_epoxide_form_zero_net_cycle_rank_with_shared_mapped_atom",
                    }
                )

    if out:
        return out

    # Legacy ambiguous open+form state remains topology/QC evidence only.
    if int(facts.get("common_cycle_rank_delta", 0)) != 0:
        return []
    if not facts.get("ring_openings") or not facts.get("ring_formations"):
        return []
    maps = set()
    for rec in list(facts["ring_openings"]) + list(facts["ring_formations"]):
        maps.update(rec.get("conserved_maps", []))
    return [
        {
            "event": "ring_rewiring_candidate",
            "role": "topology_evidence",
            "support_maps": sorted(maps),
            "mapping_scope": "mapped_ring_lineage_open_plus_form_with_zero_net_cycle_rank",
        }
    ]


def _mapped_neighbor_signature(
    facts: dict[str, Any], side: str, center: int, exclude: set[int]
) -> set[tuple[int, str, float]]:
    out = set()
    for att in _attachments(facts, side, center):
        mn = _att_map(att)
        if mn is None or mn in exclude:
            continue
        out.add((mn, _att_element(att), round(float(att.get("bond_order", 0.0)), 2)))
    return out


def _has_one_sided_unmapped_change(facts: dict[str, Any], centers: set[int]) -> bool:
    for delta in facts["attachment_deltas"]:
        if int(delta["center_map"]) not in centers:
            continue
        for att in list(delta.get("lost", [])) + list(delta.get("gained", [])):
            if _att_map(att) is None:
                return True
    return False


def detect_rearrangement(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Positive conserved-skeleton migration detector.

    Generic rearrangement is promoted only when a mapped root changes attachment
    from mapped source S to mapped destination D inside the same reactant-origin
    skeleton, while the root's other mapped local environment is preserved.
    External incorporation, mapped atom loss/gain from that skeleton, one-sided
    unmapped leaving/entering attachments, component split/merge, or net cycle-rank
    change blocks promotion.
    """
    out: list[dict[str, Any]] = []
    if facts["intercomponent_formed_bonds"]:
        return out
    if facts["component_splits"] or facts["component_merges"]:
        return out
    if facts["common_cycle_rank_delta"] != 0:
        return out

    broken = facts["broken_bonds"]
    formed = facts["formed_bonds"]
    seen = set()
    rcomp = facts["reactant"]["map_to_component"]
    pcomp = facts["product"]["map_to_component"]
    r_only = set(facts["mapping"]["reactant_only_maps"])
    p_only = set(facts["mapping"]["product_only_maps"])

    for br in broken:
        for fo in formed:
            shared = {br["a"], br["b"]} & {fo["a"], fo["b"]}
            for root in shared:
                src = br["b"] if br["a"] == root else br["a"]
                dst = fo["b"] if fo["a"] == root else fo["a"]
                if src == dst:
                    continue
                origin_component = rcomp.get(root)
                if not (origin_component == rcomp.get(src) == rcomp.get(dst)):
                    continue

                # A pure rearrangement preserves membership of the mapped substrate
                # skeleton. Mapped atoms disappearing from, or appearing inside, the
                # root component indicate replacement/fragmentation/incorporation.
                if any(rcomp.get(m) == origin_component for m in r_only):
                    continue
                root_product_component = pcomp.get(root)
                if any(pcomp.get(m) == root_product_component for m in p_only):
                    continue

                centers = {root, src, dst}
                if _has_one_sided_unmapped_change(facts, centers):
                    continue

                # The migrating root must retain its other mapped local environment.
                r_sig = _mapped_neighbor_signature(facts, "reactant", root, {src, dst})
                p_sig = _mapped_neighbor_signature(facts, "product", root, {src, dst})
                if r_sig != p_sig:
                    continue

                key = (root, src, dst)
                if key in seen:
                    continue
                seen.add(key)

                aromatic_triplet = all(
                    _atom(facts, "reactant", m).get("aromatic") for m in key
                )
                out.append(
                    {
                        "event": "mapped_attachment_relocation",
                        "role": (
                            "topology_evidence"
                            if aromatic_triplet
                            else "specific_event"
                        ),
                        "root_map": root,
                        "source_map": src,
                        "destination_map": dst,
                        "root_environment_preserved": True,
                        "mapped_skeleton_membership_preserved": True,
                        "mapping_permutation_suspected": bool(aromatic_triplet),
                        "support_maps": _support(root, src, dst),
                        "mapping_scope": "conserved_mapped_intracomponent_triplet",
                    }
                )
                if not aromatic_triplet:
                    out.append(
                        {
                            "event": "rearrangement",
                            "role": "specific_event",
                            "subtype": "conserved_skeleton_intracomponent_migration",
                            "root_map": root,
                            "source_map": src,
                            "destination_map": dst,
                            "support_maps": _support(root, src, dst),
                            "mapping_scope": "conserved_mapped_intracomponent_triplet",
                        }
                    )
    return out


# =============================================================================
# Batch03: mapped directional-conversion / redox / elimination / condensation /
# pericyclic / aromatic-functionalization detectors.
# =============================================================================


def _bond_between(facts: dict[str, Any], side: str, a: int, b: int) -> float:
    return float(facts[side]["bonds"].get(tuple(sorted((int(a), int(b)))), 0.0) or 0.0)


def _element_attachments(
    facts: dict[str, Any],
    side: str,
    center: int,
    element: str,
    min_order: float = 0.0,
    max_order: float = 99.0,
) -> list[dict[str, Any]]:
    return [
        x
        for x in _attachments(facts, side, center)
        if x.get("element") == element
        and min_order <= float(x.get("bond_order", 0.0)) <= max_order
    ]


def _anhydride_partner_centers(
    facts: dict[str, Any], side: str, center: int
) -> list[int]:
    """Mapped acyl centers joined through a single-bond oxygen: C(=O)-O-C(=O)."""
    out: list[int] = []
    for oatt in _element_attachments(facts, side, center, "O", 0.8, 1.2):
        om = _att_map(oatt)
        if om is None:
            continue
        for catt in _attachments(facts, side, om):
            cm = _att_map(catt)
            if cm is None or cm == center or catt.get("element") != "C":
                continue
            if not (0.8 <= float(catt.get("bond_order", 0.0)) <= 1.2):
                continue
            if any(
                x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
                for x in _attachments(facts, side, cm)
            ):
                out.append(int(cm))
    return sorted(set(out))


def _carbon_fg_state(facts: dict[str, Any], side: str, mn: int) -> str:
    atom = _atom(facts, side, mn)
    if atom.get("element") != "C":
        return ""
    atts = _attachments(facts, side, mn)

    # Nitrile carbon: C#N with no requirement on the other substituent here.
    # Conversion-specific validators apply the stricter organic-nitrile guard.
    if any(
        x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 2.8 for x in atts
    ):
        return "nitrile"

    carbonyl_o = [
        x
        for x in atts
        if x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
    ]
    if carbonyl_o:
        single_o = [
            x
            for x in atts
            if x.get("element") == "O" and 0.8 <= float(x.get("bond_order", 0.0)) <= 1.2
        ]
        single_n = [
            x
            for x in atts
            if x.get("element") == "N" and float(x.get("bond_order", 0.0)) <= 1.2
        ]
        single_x = [
            x
            for x in atts
            if x.get("element") in HALOGENS and float(x.get("bond_order", 0.0)) <= 1.2
        ]
        if single_n:
            return "amide"
        # Chloro-/bromoformates contain both OR and X on the carbonyl carbon.
        # They must not be collapsed into generic ester state before amidation logic.
        if single_o and single_x:
            return "haloformate"
        if single_o:
            # Explicit anhydride state must precede generic ester recognition.
            if _anhydride_partner_centers(facts, side, mn):
                return "anhydride"
            # An O with another heavy neighbor is ester-like; an O with only the
            # carbonyl carbon as heavy neighbor is acid/hydroxyl-like.
            if any(int(x.get("degree", 0)) >= 2 for x in single_o):
                return "ester"
            return "carboxylic_acid"
        if single_x:
            return "acyl_halide"
        h = int(atom.get("total_h", 0))
        return "aldehyde" if h >= 1 else "ketone"

    # Halide-bearing carbon is a directional state even when aromatic.
    if any(
        x.get("element") in HALOGENS and float(x.get("bond_order", 0.0)) <= 1.2
        for x in atts
    ):
        return "aryl_halide" if atom.get("aromatic") else "alkyl_halide"

    # Alcohol carbon: exclude aromatic C-O systems; require an OH-like O.
    if atom.get("aromatic"):
        return ""
    oh = [
        x
        for x in atts
        if x.get("element") == "O"
        and 0.8 <= float(x.get("bond_order", 0.0)) <= 1.2
        and int(x.get("degree", 0)) <= 1
    ]
    if oh:
        h = int(atom.get("total_h", 0))
        if h >= 2:
            return "primary_alcohol"
        if h == 1:
            return "secondary_alcohol"
        return "tertiary_alcohol"

    single_n = [
        x
        for x in atts
        if x.get("element") == "N" and float(x.get("bond_order", 0.0)) <= 1.2
    ]
    if single_n and int(atom.get("total_h", 0)) >= 1:
        return "amine_carbon"

    return ""


def _is_azide_n(facts: dict[str, Any], side: str, mn: int) -> bool:
    """True for an N atom belonging to an N3 azide chain.

    Works for both free mapped azide and carbon-bound organic azide. It intentionally
    avoids assigning mechanism; this is only a structural functional-group predicate.
    """
    if _atom(facts, side, mn).get("element") != "N":
        return False
    for att in _attachments(facts, side, mn):
        if att.get("element") != "N":
            continue
        if any(str(x[0]) == "N" for x in (att.get("neighbor_shell") or [])):
            return True
    return False


def _is_intact_azide_n(facts: dict[str, Any], side: str, mn: int) -> bool:
    """Strict N3 structural recognition for retained azide installation.

    True for an N atom belonging to an N3 azide chain.

    Works for both free mapped azide and carbon-bound organic azide. It intentionally
    avoids assigning mechanism; this is only a structural functional-group predicate.
    """
    atom = _atom(facts, side, mn)
    if atom.get("element") != "N" or atom.get("aromatic"):
        return False
    for att in _attachments(facts, side, mn):
        if att.get("element") != "N":
            continue
        adjacent = _att_map(att)
        if adjacent is not None and _atom(facts, side, adjacent).get("aromatic"):
            continue
        # The N--N--N chain must be nonaromatic and contain azide-like pi bonding.
        # An aromatic benzotriazole N3 subgraph is not an azide reagent.
        shell_n = [
            x for x in (att.get("neighbor_shell") or [])
            if str(x[0]) == "N" and not bool(x[2])
        ]
        if any(
            max(float(att.get("bond_order", 0.0)), float(x[1])) >= 1.8
            for x in shell_n
        ):
            return True
    return False


def _nitrogen_fg_state(facts: dict[str, Any], side: str, mn: int) -> str:
    atom = _atom(facts, side, mn)
    if atom.get("element") != "N":
        return ""
    atts = _attachments(facts, side, mn)
    if _is_azide_n(facts, side, mn):
        return "azide_n"
    if any(
        x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 2.8 for x in atts
    ):
        return "nitrile_n"
    oxy = [x for x in atts if x.get("element") == "O"]
    if len(oxy) >= 2 and sum(float(x.get("bond_order", 0.0)) for x in oxy) >= 2.5:
        return "nitro"
    if len(oxy) == 1 and int(atom.get("charge", 0)) > 0:
        return "aza_n_oxide" if atom.get("aromatic") else "amine_n_oxide"
    if not oxy and int(atom.get("charge", 0)) <= 0:
        return "aza_nitrogen" if atom.get("aromatic") else "amine_nitrogen"
    return ""


def _sulfur_fg_state(facts: dict[str, Any], side: str, mn: int) -> str:
    atom = _atom(facts, side, mn)
    if atom.get("element") != "S":
        return ""
    n_oxo = sum(
        1
        for x in _attachments(facts, side, mn)
        if x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.5
    )
    if n_oxo >= 2:
        return "sulfone"
    if n_oxo == 1:
        return "sulfoxide"
    return "sulfide"


def _has_external_formed_at(facts: dict[str, Any], center: int) -> bool:
    return any(
        center in (int(x["a"]), int(x["b"]))
        for x in facts.get("intercomponent_formed_bonds", [])
    )


def _conversion_support(facts: dict[str, Any], center: int | None) -> list[int]:
    if center is None:
        return []
    maps = {int(center)}
    for side in ("reactant", "product"):
        for x in _attachments(facts, side, int(center)):
            if x.get("map") is not None:
                maps.add(int(x["map"]))
    return sorted(maps)


def _component_heavy_count(
    facts: dict[str, Any], side: str, component_id: int | None
) -> int:
    if component_id is None:
        return 0
    comps = facts[side].get("components", [])
    if component_id < 0 or component_id >= len(comps):
        return 0
    mol = comps[component_id].get("mol")
    if mol is None:
        return 0
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def _center_component_provenance(
    facts: dict[str, Any], center: int | None
) -> dict[str, Any]:
    if center is None:
        return {}
    rci = facts["reactant"]["map_to_component"].get(int(center))
    pci = facts["product"]["map_to_component"].get(int(center))
    common = set(facts["mapping"]["common_maps"])
    source_maps = sorted(
        m for m in common if facts["reactant"]["map_to_component"].get(m) == rci
    )
    product_maps = sorted(
        m for m in common if facts["product"]["map_to_component"].get(m) == pci
    )
    contributing_components: dict[int, int] = {}
    for m in product_maps:
        rc = facts["reactant"]["map_to_component"].get(m)
        if rc is not None:
            contributing_components[rc] = contributing_components.get(rc, 0) + 1
    source_common = len(source_maps)
    dominant_other = max(
        (n for rc, n in contributing_components.items() if rc != rci), default=0
    )
    r_heavy = _component_heavy_count(facts, "reactant", rci)
    merged = len(contributing_components) > 1
    role = "substrate"
    # Conservative reagent-side rule: small source component is merged into a product
    # whose mapped scaffold is dominated by another reactant component.
    if merged and r_heavy <= 8 and dominant_other >= max(4, source_common * 2):
        role = "incorporated_reagent"
    return {
        "source_component": rci,
        "product_component": pci,
        "source_component_heavy_atoms": r_heavy,
        "source_component_common_maps": source_common,
        "product_component_common_maps": len(product_maps),
        "contributing_reactant_components": contributing_components,
        "component_role": role,
    }


def _organic_nitrile_center(facts: dict[str, Any], side: str, center: int) -> bool:
    """Require R-C#N with a carbon substituent; reject O-C#N / N-C#N / S-C#N."""
    if _carbon_fg_state(facts, side, center) != "nitrile":
        return False
    atts = _attachments(facts, side, center)
    nitrile_n = [
        x
        for x in atts
        if x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 2.8
    ]
    other = [x for x in atts if x not in nitrile_n]
    return len(nitrile_n) == 1 and any(x.get("element") == "C" for x in other)


def _ester_alkyl_context(facts: dict[str, Any], side: str, center: int) -> str | None:
    """Return a small retrieval context for the ester O-alkyl group when mapped."""
    for oatt in _element_attachments(facts, side, center, "O", 0.8, 1.2):
        om = _att_map(oatt)
        if om is None:
            continue
        for x in _attachments(facts, side, om):
            cm = _att_map(x)
            if cm is None or cm == center or x.get("element") != "C":
                continue
            ca = _atom(facts, side, cm)
            if int(ca.get("total_h", 0)) >= 3 and int(ca.get("heavy_degree", 0)) <= 1:
                return "methyl"
            if int(ca.get("total_h", 0)) >= 2:
                return "alkyl"
    return None


def _conversion_records(
    name: str,
    center: int | None,
    support_maps: list[int],
    direction: str = "neutral",
    family: str = "functional_group_interconversion",
    source_state: str | None = None,
    target_state: str | None = None,
    primary_eligible: bool = True,
    specific_event_name: str | None = None,
    extra: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    base = {
        "event": name,
        "subtype": name,
        "center_map": center,
        "support_maps": sorted(set(support_maps)),
        "source_functional_group": source_state,
        "target_functional_group": target_state,
        "conversion_family": family,
        "redox_direction": direction,
        "primary_eligible": bool(primary_eligible),
        "mapping_scope": "same_mapped_functional_center_and_validated_local_environment",
    }
    if extra:
        base.update(extra)
    out = [dict(base, role="directional_conversion")]
    if direction in {"oxidation", "reduction"}:
        out.append(
            {
                "event": direction,
                "role": "axis_event",
                "subtype": name,
                "center_map": center,
                "support_maps": sorted(set(support_maps)),
                "primary_eligible": bool(primary_eligible),
                "mapping_scope": "derived_from_mapped_directional_conversion",
            }
        )
    if primary_eligible:
        spec = dict(base, role="specific_event", conversion_specific=True)
        if specific_event_name:
            spec["event"] = specific_event_name
            spec["subtype"] = name
            spec["directional_conversion"] = name
        out.append(spec)
    return out


def _has_other_heavy_edits(facts: dict[str, Any], excluded_maps: set[int]) -> bool:
    for key in ("formed_bonds", "broken_bonds"):
        for b in facts.get(key, []):
            if not ({int(b["a"]), int(b["b"])} <= excluded_maps):
                return True
    for b in facts.get("bond_order_changes", []):
        if not ({int(b["a"]), int(b["b"])} <= excluded_maps):
            return True
    return False


def _lost_gained_for_center(
    facts: dict[str, Any], center: int
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    for d in facts.get("attachment_deltas", []):
        if int(d["center_map"]) == int(center):
            return list(d.get("lost", [])), list(d.get("gained", []))
    return [], []


def _is_isocyanate_carbon(facts: dict[str, Any], side: str, center: int) -> bool:
    if _atom(facts, side, center).get("element") != "C":
        return False
    atts = _attachments(facts, side, center)
    return any(
        x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8 for x in atts
    ) and any(
        x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 1.8 for x in atts
    )


def _donor_class(facts: dict[str, Any], donor: int) -> str:
    a = _atom(facts, "reactant", donor)
    if a.get("element") != "C":
        return str(a.get("element", "unknown")).lower()
    if a.get("aromatic"):
        return "aryl_carbon"
    atts = _attachments(facts, "reactant", donor)
    if any(
        x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 2.8 for x in atts
    ):
        return "alkynyl_carbon"
    if any(
        x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 1.8 for x in atts
    ):
        return "alkenyl_carbon"
    return "alkyl_carbon"


def detect_functional_group_conversions(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Mapped directional conversions; broad redox/FGI parents are derived only.

    Batch04 adds strict source->target conversion contracts, component provenance,
    azide transformations, and conversion-driven oxidation/reduction. Formal local
    oxidation-state changes by themselves never create oxidation/reduction events.
    """
    out: list[dict[str, Any]] = []
    common = facts["mapping"]["common_maps"]
    seen: set[tuple[str, int | None]] = set()

    def add(
        name,
        center,
        direction="neutral",
        family="functional_group_interconversion",
        source=None,
        target=None,
        support=None,
        primary=True,
        specific_event_name=None,
        extra=None,
    ):
        key = (name, center)
        if key in seen:
            return
        seen.add(key)
        provenance = _center_component_provenance(facts, center)
        if (
            provenance.get("component_role") == "incorporated_reagent"
            and direction in {"oxidation", "reduction"}
            and center is not None
            and not _has_external_formed_at(facts, int(center))
        ):
            primary = False
        merged_extra = dict(provenance)
        if extra:
            merged_extra.update(extra)
        out.extend(
            _conversion_records(
                name,
                center,
                support or _conversion_support(facts, center),
                direction,
                family,
                source,
                target,
                primary,
                specific_event_name=specific_event_name,
                extra=merged_extra,
            )
        )

    for mn in common:
        el = _atom(facts, "reactant", mn).get("element")
        if el == "C":
            r = _carbon_fg_state(facts, "reactant", mn)
            p = _carbon_fg_state(facts, "product", mn)
            external = _has_external_formed_at(facts, mn)

            # External-heavy-atom incorporation veto applies only to *simple carbonyl
            # reductions*. Oxidations and directional conversions such as aldehyde ->
            # nitrile or alcohol -> acid may require mapped O/N incorporation.
            reduction_no_external = {
                ("aldehyde", "primary_alcohol"): "aldehyde_to_primary_alcohol",
                ("ketone", "secondary_alcohol"): "ketone_to_secondary_alcohol",
                ("ester", "aldehyde"): "ester_to_aldehyde",
                ("ester", "primary_alcohol"): "ester_to_primary_alcohol",
                (
                    "carboxylic_acid",
                    "primary_alcohol",
                ): "carboxylic_acid_to_primary_alcohol",
                ("nitrile", "amine_carbon"): "nitrile_to_amine",
                ("amide", "amine_carbon"): "amide_to_amine",
            }
            if (r, p) in reduction_no_external and not external:
                add(reduction_no_external[(r, p)], mn, "reduction", "reduction", r, p)

            unrestricted_table = {
                ("primary_alcohol", "aldehyde"): (
                    "primary_alcohol_to_aldehyde",
                    "oxidation",
                ),
                ("secondary_alcohol", "ketone"): (
                    "secondary_alcohol_to_ketone",
                    "oxidation",
                ),
                ("primary_alcohol", "carboxylic_acid"): (
                    "primary_alcohol_to_carboxylic_acid",
                    "oxidation",
                ),
                ("aldehyde", "carboxylic_acid"): (
                    "aldehyde_to_carboxylic_acid",
                    "oxidation",
                ),
                ("aldehyde", "nitrile"): ("aldehyde_to_nitrile", "oxidation"),
                ("amide", "nitrile"): ("amide_to_nitrile", "neutral"),
            }
            if (r, p) in unrestricted_table:
                name, direction = unrestricted_table[(r, p)]
                family = (
                    direction
                    if direction in {"oxidation", "reduction"}
                    else "functional_group_interconversion"
                )
                add(name, mn, direction, family, r, p)

            # Oxidative esterification: aldehyde carbon retains C=O and gains an
            # external alkoxy oxygen while losing the aldehydic hydrogen.
            if (
                r == "aldehyde"
                and p == "ester"
                and external
                and _hydrogen_delta(facts, mn) < 0
            ):
                add(
                    "aldehyde_to_ester",
                    mn,
                    "oxidation",
                    "oxidation",
                    r,
                    p,
                    specific_event_name="oxidative_esterification",
                )

            # Direct C-H oxidation to aldehyde from a mapped non-carbonyl carbon.
            # This covers benzylic/heteroaryl methyl oxidation without asserting reagent.
            if r == "" and p == "aldehyde":
                ra = _atom(facts, "reactant", mn)
                if (
                    not ra.get("aromatic")
                    and int(ra.get("total_h", 0)) >= 2
                    and _hydrogen_delta(facts, mn) <= -1
                    and any(
                        x.get("element") == "O"
                        and float(x.get("bond_order", 0.0)) >= 1.8
                        for x in _attachments(facts, "product", mn)
                    )
                ):
                    add(
                        "c_h_to_aldehyde",
                        mn,
                        "oxidation",
                        "oxidation",
                        "alkyl_c_h",
                        "aldehyde",
                        specific_event_name="c_h_to_aldehyde_oxidation",
                    )

            # Benzylic halide -> aldehyde.
            #
            # Covers the common structural transformation shared by Kornblum,
            # Sommelet, Krohnke-type, and related oxidations without encoding
            # the named method itself. Restrict to an sp3 halide-bearing carbon
            # directly attached to an aromatic carbon; generic alkyl-halide
            # chemistry must not be promoted to aldehyde oxidation.
            if r == "alkyl_halide" and p == "aldehyde":
                benzylic = any(
                    att.get("element") == "C"
                    and att.get("map") is not None
                    and bool(
                        _atom(
                            facts,
                            "reactant",
                            int(att["map"]),
                        ).get("aromatic")
                    )
                    for att in _attachments(facts, "reactant", mn)
                )
                product_carbonyl_o = any(
                    att.get("element") == "O"
                    and float(att.get("bond_order", 0.0)) >= 1.8
                    for att in _attachments(facts, "product", mn)
                )
                if benzylic and product_carbonyl_o:
                    add(
                        "benzylic_halide_to_aldehyde",
                        mn,
                        "oxidation",
                        "oxidation",
                        "benzylic_halide",
                        "aldehyde",
                    )

            # Heteroaryl halide hydrolysis/oxo tautomerization is not ordinary O-arylation.
            if r == "aryl_halide" and p in {"aldehyde", "ketone"}:
                lost, gained = _lost_gained_for_center(facts, mn)
                if any(x.get("element") in HALOGENS for x in lost) and any(
                    x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
                    for x in gained
                ):
                    add(
                        "heteroaryl_halide_to_oxo",
                        mn,
                        "neutral",
                        "functional_group_interconversion",
                        r,
                        p,
                        specific_event_name="heteroaryl_hydrolysis_to_oxo",
                    )

            # Nitrile -> aldehyde is a reduction but necessarily gains oxygen; validate
            # the organic nitrile center rather than vetoing external O incorporation.
            if (
                r == "nitrile"
                and p == "aldehyde"
                and _organic_nitrile_center(facts, "reactant", mn)
            ):
                add("nitrile_to_aldehyde", mn, "reduction", "reduction", r, p)

            # Neutral hydration/hydrolysis/esterification conversions.
            if r == "ester" and p == "carboxylic_acid":
                alkyl = _ester_alkyl_context(facts, "reactant", mn)
                extra = {"ester_alkyl_group": alkyl}
                if alkyl == "methyl":
                    extra["specific_conversion"] = "methyl_ester_to_carboxylic_acid"
                add(
                    "ester_to_carboxylic_acid",
                    mn,
                    "neutral",
                    "functional_group_interconversion",
                    r,
                    p,
                    extra=extra,
                )
            elif r == "carboxylic_acid" and p == "ester":
                ester_o, ester_donor = _ester_oxygen_and_donor(
                    facts, "product", mn, require_new_o_c=True
                )
                add(
                    "carboxylic_acid_to_ester",
                    mn,
                    "neutral",
                    "functional_group_interconversion",
                    r,
                    p,
                    specific_event_name="esterification",
                    extra={
                        "ester_alkyl_group": _ester_alkyl_context(facts, "product", mn),
                        "ester_oxygen_map": ester_o,
                        "ester_donor_map": ester_donor,
                        "mapping_scope_note": "same_acyl_center_not_mechanistic_oxygen_provenance",
                    },
                )
            elif (
                r == "nitrile"
                and p == "amide"
                and _organic_nitrile_center(facts, "reactant", mn)
            ):
                add(
                    "nitrile_to_amide",
                    mn,
                    "neutral",
                    "functional_group_interconversion",
                    r,
                    p,
                    specific_event_name="nitrile_hydration",
                )
            elif (
                r == "nitrile"
                and p == "carboxylic_acid"
                and _organic_nitrile_center(facts, "reactant", mn)
            ):
                add(
                    "nitrile_to_carboxylic_acid",
                    mn,
                    "neutral",
                    "functional_group_interconversion",
                    r,
                    p,
                )

            # Oxidative C-H -> carbonyl. A methylene carbon becoming a ketone is the
            # audited high-confidence subtype; keep c_h_to_carbonyl as its parent.
            if r == "" and p == "ketone":
                ra = _atom(facts, "reactant", mn)
                if (
                    not ra.get("aromatic")
                    and int(ra.get("total_h", 0)) >= 2
                    and any(
                        x.get("element") == "O"
                        and float(x.get("bond_order", 0.0)) >= 1.8
                        for x in _attachments(facts, "product", mn)
                    )
                ):
                    add(
                        "methylene_to_ketone",
                        mn,
                        "oxidation",
                        "oxidation",
                        "methylene",
                        "ketone",
                        extra={"conversion_parent": "c_h_to_carbonyl"},
                    )

            # Neutral halide/alcohol FGI at the same mapped carbon.
            if r in {
                "primary_alcohol",
                "secondary_alcohol",
                "tertiary_alcohol",
            } and p in {"alkyl_halide", "aryl_halide"}:
                add(
                    "alcohol_to_halide",
                    mn,
                    "neutral",
                    "functional_group_interconversion",
                    r,
                    p,
                )
            elif r in {"alkyl_halide", "aryl_halide"} and p in {
                "primary_alcohol",
                "secondary_alcohol",
                "tertiary_alcohol",
            }:
                add(
                    "halide_to_alcohol",
                    mn,
                    "neutral",
                    "functional_group_interconversion",
                    r,
                    p,
                )

            # Hydrodehalogenation / hydrogenolysis at the same mapped carbon.
            if (
                r in {"alkyl_halide", "aryl_halide"}
                and p == ""
                and _hydrogen_delta(facts, mn) > 0
                and not external
            ):
                name = (
                    "aryl_halide_to_hydrogen"
                    if _atom(facts, "reactant", mn).get("aromatic")
                    else "alkyl_halide_to_hydrogen"
                )
                add(name, mn, "reduction", "reduction", r, "hydrocarbon")

            # Reductive amination: carbonyl -> C-N/H with external N incorporation.
            if r in {"aldehyde", "ketone"} and external:
                formed_n = []
                for b in facts.get("intercomponent_formed_bonds", []):
                    if mn not in (b["a"], b["b"]):
                        continue
                    other = int(b["b"] if int(b["a"]) == mn else b["a"])
                    if _atom(facts, "reactant", other).get("element") == "N":
                        formed_n.append(other)
                single_n = [
                    n
                    for n in formed_n
                    if 0.8 <= _bond_between(facts, "product", mn, n) <= 1.2
                ]
                if (
                    single_n
                    and not _element_attachments(facts, "product", mn, "O", 1.8)
                    and _hydrogen_delta(facts, mn) > 0
                ):
                    add(
                        "carbonyl_to_amine_reductive_amination",
                        mn,
                        "reduction",
                        "reductive_amination",
                        r,
                        "amine_carbon",
                        _support(mn, *single_n),
                        specific_event_name="reductive_amination",
                    )

            # Azide installation at carbon. The proximal azide N must be structurally
            # part of an N3 chain in the reactant, so generic C-N coupling is excluded.
            lost, gained = _lost_gained_for_center(facts, mn)
            azide_gains = [
                g
                for g in gained
                if _att_map(g) is not None
                and _is_intact_azide_n(facts, "reactant", int(_att_map(g)))
                and _is_intact_azide_n(facts, "product", int(_att_map(g)))
                and 0.8 <= _bond_between(
                    facts, "product", mn, int(_att_map(g))
                ) <= 1.2
            ]
            if azide_gains:
                azn = int(_att_map(azide_gains[0]))
                # Epoxide azidolysis is a ring-opening event first.
                epoxide_open = any(
                    rec.get("subfamily") == "epoxide_ring_opening"
                    and mn in set(rec.get("conserved_maps", []))
                    for rec in facts.get("ring_openings", [])
                )
                if epoxide_open:
                    add(
                        "epoxide_to_azido_alcohol",
                        mn,
                        "neutral",
                        "functional_group_interconversion",
                        "epoxide",
                        "azido_alcohol",
                        _support(mn, azn),
                        primary=False,
                        extra={"specific_event": "epoxide_azidolysis"},
                    )
                else:
                    lost_hal = [x for x in lost if x.get("element") in HALOGENS]
                    lost_oh = [
                        x
                        for x in lost
                        if x.get("element") == "O" and int(x.get("degree", 0)) <= 1
                    ]
                    lost_sulf = [
                        x
                        for x in lost
                        if x.get("element") == "O"
                        and any(
                            str(y[0]) == "S" for y in (x.get("neighbor_shell") or [])
                        )
                    ]
                    if lost_hal:
                        aromatic = bool(_atom(facts, "reactant", mn).get("aromatic"))
                        name = (
                            "aryl_halide_to_azide"
                            if aromatic
                            else "alkyl_halide_to_azide"
                        )
                        add(
                            name,
                            mn,
                            "neutral",
                            "functional_group_interconversion",
                            "aryl_halide" if aromatic else "alkyl_halide",
                            "organic_azide",
                            _support(mn, azn),
                            specific_event_name="azide_substitution",
                            extra={"leaving_group": lost_hal[0].get("element")},
                        )
                    elif lost_sulf:
                        add(
                            "sulfonate_to_azide",
                            mn,
                            "neutral",
                            "functional_group_interconversion",
                            "sulfonate",
                            "organic_azide",
                            _support(mn, azn),
                            specific_event_name="azide_substitution",
                        )
                    elif (
                        lost_oh
                        and r == "carboxylic_acid"
                        and _has_conserved_carbonyl(facts, mn)
                    ):
                        add(
                            "carboxylic_acid_to_acyl_azide",
                            mn,
                            "neutral",
                            "functional_group_interconversion",
                            "carboxylic_acid",
                            "acyl_azide",
                            _support(mn, azn),
                            specific_event_name="azide_substitution",
                        )
                    elif (
                        lost_oh
                        and r in {
                            "primary_alcohol", "secondary_alcohol", "tertiary_alcohol"
                        }
                        and _is_sp3_carbon(facts, "reactant", mn)
                    ):
                        add(
                            "alcohol_to_azide",
                            mn,
                            "neutral",
                            "functional_group_interconversion",
                            r,
                            "organic_azide",
                            _support(mn, azn),
                            specific_event_name="azide_substitution",
                        )

        elif el == "N":
            r = _nitrogen_fg_state(facts, "reactant", mn)
            p = _nitrogen_fg_state(facts, "product", mn)
            if r == "nitro" and p in {"amine_nitrogen", "aza_nitrogen"}:
                carbon_neighbors = [
                    x
                    for x in _attachments(facts, "reactant", mn)
                    if x.get("element") == "C"
                ]
                retained_carbon_maps = {
                    int(_att_map(x))
                    for x in carbon_neighbors
                    if _att_map(x) is not None
                } & {
                    int(_att_map(x))
                    for x in _attachments(facts, "product", mn)
                    if x.get("element") == "C" and _att_map(x) is not None
                }
                product_noncarbon_heavy = [
                    x
                    for x in _attachments(facts, "product", mn)
                    if x.get("element") not in {"C", "H"}
                ]
                # True nitro reduction keeps the carbon skeleton attachment, removes
                # both nitro oxygens and does not create N-N/new heteroatom connectivity.
                if (
                    carbon_neighbors
                    and retained_carbon_maps
                    and not product_noncarbon_heavy
                ):
                    context = (
                        "aryl"
                        if any(
                            _att_map(x) is not None
                            and _atom(facts, "reactant", int(_att_map(x))).get(
                                "aromatic"
                            )
                            for x in carbon_neighbors
                        )
                        else "alkyl"
                    )
                    add(
                        "nitro_to_amine",
                        mn,
                        "reduction",
                        "reduction",
                        r,
                        p,
                        extra={"nitro_context": context},
                    )
            elif r == "azide_n" and p in {"amine_nitrogen", "aza_nitrogen"}:
                add("azide_to_amine", mn, "reduction", "reduction", r, p)
            elif r == "aza_n_oxide" and p == "aza_nitrogen":
                support = set(_conversion_support(facts, mn))
                primary = not _has_other_heavy_edits(facts, support)
                add(
                    "aza_n_oxide_to_aza_nitrogen",
                    mn,
                    "reduction",
                    "reduction",
                    r,
                    p,
                    primary=primary,
                    extra={"specific_event": "aza_n_oxide_deoxygenation"},
                )
            elif r == "amine_n_oxide" and p == "amine_nitrogen":
                support = set(_conversion_support(facts, mn))
                primary = not _has_other_heavy_edits(facts, support)
                add(
                    "tertiary_amine_n_oxide_to_amine",
                    mn,
                    "reduction",
                    "reduction",
                    r,
                    p,
                    primary=primary,
                )
            elif r == "aza_nitrogen" and p == "aza_n_oxide":
                add("aza_nitrogen_to_n_oxide", mn, "oxidation", "oxidation", r, p)
            elif r == "amine_nitrogen" and p == "amine_n_oxide":
                add("tertiary_amine_to_n_oxide", mn, "oxidation", "oxidation", r, p)

        elif el == "S":
            r = _sulfur_fg_state(facts, "reactant", mn)
            p = _sulfur_fg_state(facts, "product", mn)
            table = {
                ("sulfide", "sulfoxide"): ("sulfide_to_sulfoxide", "oxidation"),
                ("sulfide", "sulfone"): ("sulfide_to_sulfone", "oxidation"),
                ("sulfoxide", "sulfone"): ("sulfoxide_to_sulfone", "oxidation"),
                ("sulfone", "sulfoxide"): ("sulfone_to_sulfoxide", "reduction"),
                ("sulfoxide", "sulfide"): ("sulfoxide_to_sulfide", "reduction"),
                ("sulfone", "sulfide"): ("sulfone_to_sulfide", "reduction"),
            }
            if (r, p) in table:
                name, direction = table[(r, p)]
                add(name, mn, direction, direction, r, p)

    # Pair-centered C-C bond-order conversions.
    for bc in facts.get("bond_order_changes", []):
        a, b = int(bc["a"]), int(bc["b"])
        if (
            _atom(facts, "reactant", a).get("element") != "C"
            or _atom(facts, "reactant", b).get("element") != "C"
        ):
            continue
        ro, po = float(bc["reactant_order"]), float(bc["product_order"])

        # Dehydrogenation: C-C -> C=C with loss of H at both endpoints and no heavy
        # attachment exchange. This is distinct from beta-elimination.
        if (
            ro <= 1.2
            and po >= 1.8
            and _hydrogen_delta(facts, a) < 0
            and _hydrogen_delta(facts, b) < 0
        ):
            la, ga = _lost_gained_for_center(facts, a)
            lb, gb = _lost_gained_for_center(facts, b)
            heavy_loss = [x for x in la + lb if x.get("element") != "H"]
            heavy_gain = [x for x in ga + gb if x.get("element") != "H"]
            if not heavy_loss and not heavy_gain:
                name = "alkane_to_alkene"
                extra = {}
                # More specific audited subtype when the new alkene is alpha,beta to a carbonyl.
                carbonyl_adj = any(
                    any(
                        att.get("element") == "C" and bool(att.get("carbonyl_like"))
                        for att in _attachments(facts, "product", x)
                    )
                    for x in (a, b)
                )
                spec_name = (
                    "carbonyl_alpha_beta_dehydrogenation"
                    if carbonyl_adj
                    else "dehydrogenation"
                )
                add(
                    name,
                    a,
                    "oxidation",
                    "oxidation",
                    "alkane",
                    "alkene",
                    _support(a, b),
                    specific_event_name=spec_name,
                    extra=extra,
                )

        if ro >= 1.8 and po <= 1.2:
            # Epoxide: same O becomes a single-bond neighbor of both alkene carbons.
            oa = {
                x.get("map")
                for x in _element_attachments(facts, "product", a, "O", 0.8, 1.2)
                if x.get("map") is not None
            }
            ob = {
                x.get("map")
                for x in _element_attachments(facts, "product", b, "O", 0.8, 1.2)
                if x.get("map") is not None
            }
            shared_o = sorted(oa & ob)
            if shared_o and ro < 2.8:
                add(
                    "alkene_to_epoxide",
                    a,
                    "oxidation",
                    "oxidation",
                    "alkene",
                    "epoxide",
                    _support(a, b, *shared_o),
                )
                continue

            # Dihydroxylation: both endpoints gain OH-like O attachments.
            def gained_oh(center: int) -> list[int]:
                outm = []
                for d in facts["attachment_deltas"]:
                    if int(d["center_map"]) != center:
                        continue
                    for x in d.get("gained", []):
                        if (
                            x.get("element") == "O"
                            and int(x.get("degree", 0)) <= 1
                            and x.get("map") is not None
                        ):
                            outm.append(int(x["map"]))
                return outm

            ga, gb = gained_oh(a), gained_oh(b)
            if ga and gb and ro < 2.8:
                add(
                    "alkene_to_diol",
                    a,
                    "oxidation",
                    "oxidation",
                    "alkene",
                    "vicinal_diol",
                    _support(a, b, *ga, *gb),
                )
                continue

            # Hydrogenation: no new heavy attachment at either endpoint and H increases.
            if (
                _hydrogen_delta(facts, a) > 0
                and _hydrogen_delta(facts, b) > 0
                and not _has_external_formed_at(facts, a)
                and not _has_external_formed_at(facts, b)
            ):
                if ro >= 2.8 and po <= 1.2:
                    name = "alkyne_to_alkane"
                elif ro >= 2.8 and po >= 1.8:
                    name = "alkyne_to_alkene"
                else:
                    name = "alkene_to_alkane"
                spec_name = (
                    "complete_alkyne_hydrogenation"
                    if name == "alkyne_to_alkane"
                    else (
                        "alkene_hydrogenation"
                        if name == "alkene_to_alkane"
                        else "partial_alkyne_hydrogenation"
                    )
                )
                add(
                    name,
                    a,
                    "reduction",
                    "reduction",
                    "pi_bond",
                    "reduced_pi_bond",
                    _support(a, b),
                    specific_event_name=spec_name,
                )

    # Oxidative cleavage where one mapped alkene carbon survives as an aldehyde and
    # the partner is unmapped/removed. This is a mapping-limited but chemically exact
    # directional pattern seen repeatedly in the oxidation audit.
    for d in facts.get("attachment_deltas", []):
        c = int(d["center_map"])
        if _atom(facts, "reactant", c).get("element") != "C":
            continue
        lost_pi_c = [
            x
            for x in d.get("lost", [])
            if x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 1.8
        ]
        gained_co = [
            x
            for x in d.get("gained", [])
            if x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
        ]
        if (
            lost_pi_c
            and gained_co
            and _carbon_fg_state(facts, "product", c) == "aldehyde"
        ):
            add(
                "alkene_oxidative_cleavage_to_aldehyde",
                c,
                "oxidation",
                "oxidation",
                "alkene_carbon",
                "aldehyde",
            )

    return out


def detect_addition_events(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Mapped addition validation separated from generic pi-bond consumption.

    `addition_like_topology` is evidence only. Specific events require new sigma bonds
    to the same mapped pi-system endpoints that undergo bond-order loss.
    """
    out: list[dict[str, Any]] = []
    consumed: list[tuple[int, int, float, float]] = []
    for bc in facts.get("bond_order_changes", []):
        ro, po = float(bc["reactant_order"]), float(bc["product_order"])
        if ro >= 1.8 and po < ro - 0.7:
            a, b = int(bc["a"]), int(bc["b"])
            consumed.append((a, b, ro, po))
            out.append(
                {
                    "event": "addition_like_topology",
                    "role": "topology_evidence",
                    "subtype": "mapped_pi_bond_consumption",
                    "support_maps": _support(a, b),
                    "mapping_scope": "mapped_pi_bond_order_decrease_only_not_synthetic_event",
                }
            )

    # Carbonyl addition: C=O -> C-O and exactly one new heavy sigma bond at the same
    # carbonyl C. No leaving-group exchange at that center.
    for c in facts["mapping"]["common_maps"]:
        if _carbon_fg_state(facts, "reactant", c) not in {"aldehyde", "ketone"}:
            continue
        r_o = [x for x in _element_attachments(facts, "reactant", c, "O", 1.8)]
        if not r_o:
            continue
        om = _att_map(r_o[0])
        if om is None or not (0.8 <= _bond_between(facts, "product", c, om) <= 1.2):
            continue
        formed = [
            b for b in facts.get("formed_bonds", []) if c in (int(b["a"]), int(b["b"]))
        ]
        new_heavy = []
        for b in formed:
            other = int(b["b"] if int(b["a"]) == c else b["a"])
            if other != om:
                new_heavy.append(other)
        broken_at_c = [
            b for b in facts.get("broken_bonds", []) if c in (int(b["a"]), int(b["b"]))
        ]
        if len(new_heavy) != 1 or broken_at_c:
            continue
        donor = new_heavy[0]
        donor_el = _atom(facts, "reactant", donor).get("element", "?")
        subtype = f"carbonyl_{str(donor_el).lower()}_addition"
        if donor_el == "C":
            subtype = "carbonyl_c_c_addition"
        elif donor_el == "N":
            subtype = "carbonyl_c_n_addition"
        elif donor_el == "O":
            subtype = "carbonyl_c_o_addition"
        elif donor_el == "S":
            subtype = "carbonyl_c_s_addition"
        out.append(
            {
                "event": "carbonyl_addition",
                "role": "specific_event",
                "subtype": subtype,
                "center_map": c,
                "incoming_map": donor,
                "donor_class": _donor_class(facts, donor),
                "support_maps": _support(c, om, donor),
                "mapping_scope": "same_mapped_carbonyl_c_o_reduced_plus_one_new_sigma_bond",
            }
        )

    # Heterocumulene addition: amine + isocyanate -> urea linkage.
    for c in facts["mapping"]["common_maps"]:
        if not _is_isocyanate_carbon(facts, "reactant", c):
            continue
        if not _element_attachments(facts, "product", c, "O", 1.8):
            continue
        formed_n = []
        for b in facts.get("intercomponent_formed_bonds", []):
            if c not in (int(b["a"]), int(b["b"])):
                continue
            other = int(b["b"] if int(b["a"]) == c else b["a"])
            if _atom(facts, "reactant", other).get("element") == "N":
                formed_n.append(other)
        if formed_n:
            out.append(
                {
                    "event": "heterocumulene_addition",
                    "role": "specific_event",
                    "subtype": "amine_to_isocyanate_urea_formation",
                    "center_map": c,
                    "support_maps": _support(c, *formed_n),
                    "mapping_scope": "mapped_isocyanate_c_n_o_plus_external_amine_c_n_formation",
                }
            )

    # Alkene/alkyne additions with mapped incoming atoms tied to consumed endpoints.
    for a, b, ro, po in consumed:
        if (
            _atom(facts, "reactant", a).get("element") != "C"
            or _atom(facts, "reactant", b).get("element") != "C"
        ):
            continue
        formed_a = []
        formed_b = []
        for x in facts.get("formed_bonds", []):
            aa, bb = int(x["a"]), int(x["b"])
            if a in (aa, bb):
                formed_a.append(bb if aa == a else aa)
            if b in (aa, bb):
                formed_b.append(bb if aa == b else aa)
        # Specific redox conversions already own hydrogenation/dihydroxylation/epoxidation.
        # Here we add reagent-specific sigma additions not represented by those conversions.
        for center, other_center, formed_here in ((a, b, formed_a), (b, a, formed_b)):
            for inc in formed_here:
                el = _atom(facts, "reactant", inc).get("element")
                if el == "Si" and _hydrogen_delta(facts, other_center) > 0:
                    out.append(
                        {
                            "event": "alkene_hydrosilylation",
                            "role": "specific_event",
                            "subtype": "hydrosilylation",
                            "center_map": center,
                            "support_maps": _support(a, b, inc),
                            "mapping_scope": "mapped_c_c_pi_consumption_plus_c_si_formation_and_h_gain",
                        }
                    )
                if el in HALOGENS and _hydrogen_delta(facts, other_center) > 0:
                    out.append(
                        {
                            "event": "alkene_hydrohalogenation",
                            "role": "specific_event",
                            "subtype": {
                                "F": "hydrofluorination",
                                "Cl": "hydrochlorination",
                                "Br": "hydrobromination",
                                "I": "hydroiodination",
                            }[str(el)],
                            "center_map": center,
                            "support_maps": _support(a, b, inc),
                            "mapping_scope": "mapped_c_c_pi_consumption_plus_c_x_formation_and_h_gain",
                        }
                    )

        # Halocyclization / haloetherification: one endpoint gains halogen and the
        # other gains an intramolecular heteroatom bond while a ring is formed.
        hal = []
        het = []
        for center, formed_here in ((a, formed_a), (b, formed_b)):
            for inc in formed_here:
                el = _atom(facts, "reactant", inc).get("element")
                if el in HALOGENS:
                    hal.append((center, inc, el))
                elif el in {"O", "N", "S"} and _same_reactant_origin(
                    facts, center, inc
                ):
                    het.append((center, inc, el))
        if hal and het and facts.get("ring_formations"):
            hc, hm, hel = hal[0]
            xc, xm, xel = het[0]
            subtype = {
                "O": "haloetherification",
                "N": "haloamination",
                "S": "halothioetherification",
            }.get(xel, "halocyclization")
            out.append(
                {
                    "event": "halocyclization",
                    "role": "specific_event",
                    "subtype": subtype,
                    "halogen": hel,
                    "support_maps": _support(a, b, hm, xm),
                    "mapping_scope": "mapped_pi_consumption_plus_c_x_and_intramolecular_c_hetero_formation_with_ring_gain",
                }
            )

    return out



def detect_conjugate_acyl_addition(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Net Stetter-compatible aldehyde-derived 1,4-acyl incorporation.

    Requires an enone C=C -> C-C on the same mapped alpha/beta pair,
    a new beta-C--aldehyde-C bond from a separate reactant component,
    unchanged electrophile and donor carbonyl oxygen identities, and
    complementary mapped H loss/gain.  No catalytic mechanism is inferred.
    """
    out: list[dict[str, Any]] = []
    for bc in facts.get("bond_order_changes", []):
        if not (1.8 <= float(bc["reactant_order"]) <= 2.2
                and 0.8 <= float(bc["product_order"]) <= 1.2):
            continue
        a, b = int(bc["a"]), int(bc["b"])
        if any(_atom(facts, "reactant", m).get("element") != "C" for m in (a, b)):
            continue
        for alpha, beta in ((a, b), (b, a)):
            # Alpha carbon is bonded to a conserved ketone carbonyl.
            carbonyls = [
                int(att["map"])
                for att in _attachments(facts, "reactant", alpha)
                if att.get("element") == "C"
                and att.get("map") is not None
                and 0.8 <= float(att.get("bond_order", 0)) <= 1.2
                and _carbon_fg_state(facts, "reactant", int(att["map"])) == "ketone"
                and _has_conserved_carbonyl(facts, int(att["map"]))
            ]
            if not carbonyls or _hydrogen_delta(facts, alpha) <= 0:
                continue
            # Require exactly one heavy-atom addition to either original alkene atom.
            endpoint_formed = [
                bond for bond in facts.get("formed_bonds", [])
                if a in (int(bond["a"]), int(bond["b"]))
                or b in (int(bond["a"]), int(bond["b"]))
            ]
            if len(endpoint_formed) != 1:
                continue
            new_bond = endpoint_formed[0]
            if beta not in (int(new_bond["a"]), int(new_bond["b"])):
                continue
            donor = int(new_bond["b"] if int(new_bond["a"]) == beta else new_bond["a"])
            if _carbon_fg_state(facts, "reactant", donor) != "aldehyde":
                continue
            if _carbon_fg_state(facts, "product", donor) != "ketone":
                continue
            if _hydrogen_delta(facts, donor) >= 0:
                continue
            if not any(
                {int(bond["a"]), int(bond["b"])} == {beta, donor}
                and 0.8 <= float(bond.get("order", 0)) <= 1.2
                for bond in facts.get("intercomponent_formed_bonds", [])
            ):
                continue
            donor_oxygens = [
                int(oa["map"]) for oa in _element_attachments(facts, "reactant", donor, "O", 1.8)
                if oa.get("map") is not None
                and _bond_between(facts, "product", donor, int(oa["map"])) >= 1.8
            ]
            if not donor_oxygens:
                continue
            carbonyl = carbonyls[0]
            out.append({
                "event": "conjugate_acyl_addition",
                "role": "specific_event",
                "subtype": "aldehyde_derived_acyl_1_4_addition",
                "center_map": beta,
                "alpha_map": alpha,
                "donor_carbonyl_map": donor,
                "acceptor_carbonyl_map": carbonyl,
                "support_maps": _support(alpha, beta, carbonyl, donor, donor_oxygens[0]),
                "mapping_scope": "mapped_enone_beta_c_acyl_cc_formation_with_conserved_aldehyde_carbonyl",
            })
            break
    return out

def _ring_formation_touches(facts: dict[str, Any], *maps: int | None) -> bool:
    targets = {int(m) for m in maps if m is not None}
    if not targets:
        return False
    for rec in facts.get("ring_formations", []):
        if targets & set(rec.get("conserved_maps") or rec.get("ring_maps") or []):
            return True
    return False


def _is_weinreb_amide_center(facts: dict[str, Any], center: int) -> bool:
    """Structural Weinreb-amide signature: acyl N bears both O and C substituents."""
    if _carbon_fg_state(facts, "reactant", center) != "amide":
        return False
    for natt in _element_attachments(facts, "reactant", center, "N", 0.8, 1.6):
        shell = list(natt.get("neighbor_shell") or [])
        elems = {str(x[0]) for x in shell}
        if "O" in elems and "C" in elems:
            return True
    return False


def _new_acyl_attachment(
    facts: dict[str, Any], center: int, element: str, max_order: float = 1.6
) -> list[int]:
    out: list[int] = []
    for x in _attachments(facts, "product", center):
        if x.get("element") != element or x.get("map") is None:
            continue
        bo = float(x.get("bond_order", 0.0))
        if not (0.8 <= bo <= max_order):
            continue
        m = int(x["map"])
        if _bond_between(facts, "reactant", center, m) < 0.2:
            out.append(m)
    return sorted(set(out))


def _ester_oxygen_and_donor(
    facts: dict[str, Any], side: str, center: int, require_new_o_c: bool = False
) -> tuple[int | None, int | None]:
    """Mapped ester O and its non-acyl carbon substituent; optionally require new O-C bond."""
    for oatt in _element_attachments(facts, side, center, "O", 0.8, 1.2):
        om = _att_map(oatt)
        if om is None:
            continue
        for x in _attachments(facts, side, om):
            cm = _att_map(x)
            if cm is None or cm == center or x.get("element") != "C":
                continue
            if require_new_o_c and _bond_between(facts, "reactant", om, int(cm)) > 0.2:
                continue
            return int(om), int(cm)
    return None, None


def _new_ester_oxygen_and_donor(
    facts: dict[str, Any], center: int
) -> tuple[int | None, int | None]:
    """Return newly attached ester O and its non-acyl carbon substituent when mapped."""
    for om in _new_acyl_attachment(facts, center, "O", 1.2):
        for x in _attachments(facts, "product", om):
            cm = _att_map(x)
            if cm is None or cm == center or x.get("element") != "C":
                continue
            return om, int(cm)
    return None, None


def detect_acyl_conversions(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Batch05 acyl semantics above broad `carbonyl_substitution` bookkeeping.

    The detector is net-structural: it distinguishes acyl exchange, amidation,
    lactam/lactone formation, anhydride formation, Weinreb ketone formation and
    aromatic acylation without asserting a detailed mechanism.
    """
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, int | None]] = set()

    def add_conversion(
        name: str,
        center: int,
        source: str,
        target: str,
        event: str | None = None,
        support: list[int] | None = None,
        extra: dict[str, Any] | None = None,
    ):
        key = (name, center)
        if key in seen:
            return
        seen.add(key)
        merged = _center_component_provenance(facts, center)
        if extra:
            merged.update(extra)
        out.extend(
            _conversion_records(
                name,
                center,
                support or _conversion_support(facts, center),
                "neutral",
                "acyl_conversion",
                source,
                target,
                True,
                specific_event_name=event,
                extra=merged,
            )
        )

    # Pair cyclic/bridged anhydrides once, not as two esterifications.
    anhydride_pairs: set[tuple[int, int]] = set()
    for c in facts["mapping"]["common_maps"]:
        if _carbon_fg_state(facts, "reactant", c) != "carboxylic_acid":
            continue
        if _carbon_fg_state(facts, "product", c) != "anhydride":
            continue
        for other in _anhydride_partner_centers(facts, "product", c):
            if other not in facts["mapping"]["common_maps"]:
                continue
            if _carbon_fg_state(facts, "reactant", other) != "carboxylic_acid":
                continue
            pair = tuple(sorted((int(c), int(other))))
            if pair in anhydride_pairs:
                continue
            anhydride_pairs.add(pair)
            intra = _same_reactant_origin(facts, *pair)
            out.append(
                {
                    "event": (
                        "intramolecular_anhydride_formation"
                        if intra
                        else "anhydride_formation"
                    ),
                    "role": "specific_event",
                    "subtype": "carboxylic_acids_to_anhydride",
                    "center_map": pair[0],
                    "acyl_center_maps": list(pair),
                    "support_maps": list(pair),
                    "mapping_scope": "paired_mapped_acyl_centers_joined_by_product_anhydride_oxygen",
                }
            )

    for c in facts["mapping"]["common_maps"]:
        if _atom(facts, "reactant", c).get(
            "element"
        ) != "C" or not _has_conserved_carbonyl(facts, c):
            continue
        r = _carbon_fg_state(facts, "reactant", c)
        p = _carbon_fg_state(facts, "product", c)
        lost, gained = _lost_gained_for_center(facts, c)
        lost_elements = [x.get("element") for x in lost]
        gained_elements = [x.get("element") for x in gained]

        # O-for-O and N-for-N exchange on a conserved acyl center.
        if (
            r == "ester"
            and p == "ester"
            and "O" in lost_elements
            and "O" in gained_elements
        ):
            gained_o = next((x for x in gained if x.get("element") == "O"), None)
            if gained_o is not None and _external_gain(facts, c, gained_o):
                add_conversion(
                    "transesterification",
                    c,
                    "ester",
                    "ester",
                    event="transesterification",
                )

        if (
            r == "amide"
            and p == "amide"
            and "N" in lost_elements
            and "N" in gained_elements
        ):
            gained_n = next((x for x in gained if x.get("element") == "N"), None)
            if gained_n is not None and _external_gain(facts, c, gained_n):
                add_conversion(
                    "transamidation", c, "amide", "amide", event="transamidation"
                )

        # Acid/ester/acyl-halide -> amide. Aromaticized lactams can carry a 1.5 C-N
        # bond, so use direct product attachment instead of relying only on FG state.
        new_ns = [
            n for n in _new_acyl_attachment(facts, c, "N", 1.6)
            if not _is_intact_azide_n(facts, "product", n)
        ]
        if r in {"carboxylic_acid", "ester", "acyl_halide", "haloformate"} and new_ns:
            nmap = new_ns[0]
            cyclic = int(
                facts.get("common_cycle_rank_delta", 0)
            ) > 0 and _ring_formation_touches(facts, c, nmap)
            # A cyclic carbamate has O=C(N)O, with C, N and O on the
            # newly constructed ring; a conventional lactam lacks that O.
            cyclic_carbamate = cyclic and any(
                _att_map(oatt) is not None
                and any(
                    {c, nmap, int(_att_map(oatt))}.issubset(
                        {int(m) for m in (
                            ring.get("ring_maps") or ring.get("conserved_maps") or []
                        )}
                    )
                    for ring in (facts.get("ring_formations") or [])
                )
                for oatt in _element_attachments(facts, "product", c, "O", 0.8, 1.2)
            )
            if cyclic_carbamate:
                direction, event = f"{r}_to_cyclic_carbamate", "carbamate_formation"
            elif r == "carboxylic_acid":
                direction, event = "carboxylic_acid_to_amide", (
                    "lactamization" if cyclic else "amidation"
                )
            elif r == "ester":
                direction, event = "ester_to_amide", (
                    "lactamization" if cyclic else "amidation"
                )
            elif r == "acyl_halide":
                direction, event = "acyl_halide_to_amide", (
                    "lactamization" if cyclic else "amidation"
                )
            else:
                direction, event = "haloformate_to_carbamate", "carbamate_formation"
            add_conversion(
                direction,
                c,
                r,
                "carbamate" if (r == "haloformate" or cyclic_carbamate) else "amide",
                event=event,
                support=_support(c, nmap),
                extra={"nucleophile_map": nmap, "cyclic": cyclic},
            )

        # Carboxylic acid -> acyl halide is a canonical activation event, not generic
        # carbonyl substitution. Preserve the directional state explicitly.
        if r == "carboxylic_acid" and p == "acyl_halide":
            gained_x = next((x for x in gained if x.get("element") in HALOGENS), None)
            if gained_x is not None:
                add_conversion(
                    "carboxylic_acid_to_acyl_halide",
                    c,
                    r,
                    p,
                    event="acyl_halide_formation",
                    support=_support(c, _att_map(gained_x)),
                    extra={"halogen": gained_x.get("element")},
                )

        # Acid -> ester in a newly formed ring is lactonization. The orthogonal
        # carboxylic_acid_to_ester directional conversion is emitted elsewhere.
        if r == "carboxylic_acid" and p == "ester":
            om, donor = _new_ester_oxygen_and_donor(facts, c)
            if (
                om is not None
                and int(facts.get("common_cycle_rank_delta", 0)) > 0
                and _ring_formation_touches(facts, c, om)
            ):
                out.append(
                    {
                        "event": "lactonization",
                        "role": "specific_event",
                        "subtype": "intramolecular_carboxylic_acid_to_ester",
                        "center_map": c,
                        "oxygen_map": om,
                        "donor_map": donor,
                        "support_maps": _support(c, om, donor),
                        "mapping_scope": "mapped_acid_to_ester_plus_validated_ring_formation",
                    }
                )

        # Weinreb amide -> ketone: C-N(OR)R is displaced by a newly installed carbon.
        if _is_weinreb_amide_center(facts, c) and p == "ketone":
            new_cs = _new_acyl_attachment(facts, c, "C", 1.2)
            if new_cs:
                add_conversion(
                    "weinreb_amide_to_ketone",
                    c,
                    "weinreb_amide",
                    "ketone",
                    event="weinreb_amide_to_ketone",
                    support=_support(c, new_cs[0]),
                    extra={
                        "incoming_carbon_map": new_cs[0],
                        "incoming_carbon_class": _donor_class(facts, new_cs[0]),
                    },
                )

        # Acyl halide + aromatic carbon -> aryl ketone. This is acyl C-C formation,
        # not generic cross-coupling.
        if r == "acyl_halide" and p == "ketone":
            new_cs = _new_acyl_attachment(facts, c, "C", 1.2)
            for cm in new_cs:
                if _atom(facts, "reactant", cm).get("aromatic"):
                    add_conversion(
                        "acyl_halide_to_aryl_ketone",
                        c,
                        "acyl_halide",
                        "ketone",
                        event="aromatic_acylation",
                        support=_support(c, cm),
                        extra={"incoming_aromatic_carbon_map": cm},
                    )
                    break

    return out


def detect_carbonyl_olefination(facts: dict[str, Any]) -> list[dict[str, Any]]:
    # v6.7: specific olefination requires LOSS of a direct P/S/Si donor handle.
    out: list[dict[str, Any]] = []

    for c in facts.get("mapping", {}).get("common_maps", []):
        c = int(c)
        if _carbon_fg_state(facts, "reactant", c) not in {"aldehyde", "ketone"}:
            continue
        if _element_attachments(facts, "product", c, "O", 1.8):
            continue

        formed_double_c = [
            x for x in _attachments(facts, "product", c)
            if x.get("element") == "C"
            and float(x.get("bond_order", 0.0)) >= 1.8
            and x.get("map") is not None
            and _bond_between(facts, "reactant", c, int(x["map"])) < 0.2
        ]
        if not formed_double_c:
            continue

        donor = int(formed_double_c[0]["map"])
        lost_handles = _cnd_lost_direct_handle_elements(facts, donor)

        if not lost_handles:
            out.append(
                {
                    "event": "carbonyl_to_alkene_topology",
                    "role": "topology_evidence",
                    "subtype": "unvalidated_or_retained_handle_carbonyl_to_alkene",
                    "center_map": c,
                    "donor_map": donor,
                    "retained_handles": sorted(_cnd_retained_direct_handle_elements(facts, donor)),
                    "support_maps": _support(c, donor),
                    "mapping_scope": "v67_carbonyl_to_alkene_without_lost_p_s_si_handle",
                    "promotion_blockers": ["missing_lost_P_S_Si_donor_handle"],
                }
            )
            continue

        subtype = (
            "wittig_or_hwe_like_olefination" if "P" in lost_handles
            else "sulfur_reagent_olefination" if "S" in lost_handles
            else "peterson_like_olefination"
        )

        out.append(
            {
                "event": "carbonyl_olefination",
                "role": "specific_event",
                "subtype": subtype,
                "center_map": c,
                "donor_map": donor,
                "donor_handle": sorted(lost_handles),
                "support_maps": _support(c, donor),
                "mapping_scope": "v67_carbonyl_to_alkene_plus_direct_p_s_si_handle_loss",
            }
        )

    return out


def detect_oxidative_structural_conversions(
    facts: dict[str, Any]
) -> list[dict[str, Any]]:
    """Specific oxidation motifs not captured by simple source/target FG states."""
    out: list[dict[str, Any]] = []

    # Oxidative deformylation: aromatic carbon loses a carbonyl-like C substituent and gains O.
    for d in facts.get("attachment_deltas", []):
        c = int(d["center_map"])
        if not _atom(facts, "reactant", c).get("aromatic"):
            continue
        lost_formyl = [
            x
            for x in d.get("lost", [])
            if x.get("element") == "C" and bool(x.get("carbonyl_like"))
        ]
        gained_o = [x for x in d.get("gained", []) if x.get("element") == "O"]
        if lost_formyl and gained_o:
            out.append(
                {
                    "event": "oxidative_deformylation",
                    "role": "specific_event",
                    "subtype": "aryl_formyl_to_phenol",
                    "center_map": c,
                    "support_maps": _support(
                        c, *[_att_map(x) for x in lost_formyl + gained_o]
                    ),
                    "mapping_scope": "mapped_aromatic_center_formyl_loss_plus_oxygen_gain",
                    "redox_direction": "oxidation",
                }
            )

    # Quinone formation: two or more aromatic C-O single bonds become C=O while the
    # ring loses aromaticity. Preserve one composite event rather than local oxidations.
    qcenters = []
    for c in facts["mapping"]["common_maps"]:
        ra, pa = _atom(facts, "reactant", c), _atom(facts, "product", c)
        if ra.get("element") != "C" or not ra.get("aromatic") or pa.get("aromatic"):
            continue
        ro = [x for x in _element_attachments(facts, "reactant", c, "O", 0.8, 1.2)]
        po = [x for x in _element_attachments(facts, "product", c, "O", 1.8)]
        if ro and po:
            qcenters.append(c)
    if len(qcenters) >= 2:
        out.append(
            {
                "event": "quinone_formation",
                "role": "specific_event",
                "subtype": "aromatic_di_oxygenation_to_quinone",
                "support_maps": sorted(qcenters),
                "mapping_scope": "multiple_mapped_aromatic_c_o_single_to_c_o_double_centers",
                "redox_direction": "oxidation",
            }
        )

    # Oxidative aromatization: many atoms in one conserved ring become aromatic with no
    # external heavy-atom incorporation controlling the changed ring.
    newly_aromatic = [
        m
        for m in facts["mapping"]["common_maps"]
        if not _atom(facts, "reactant", m).get("aromatic")
        and _atom(facts, "product", m).get("aromatic")
    ]
    if len(newly_aromatic) >= 4 and not facts.get("intercomponent_formed_bonds"):
        out.append(
            {
                "event": "oxidative_aromatization",
                "role": "specific_event",
                "subtype": "dehydrogenative_aromatization",
                "support_maps": sorted(newly_aromatic),
                "mapping_scope": "conserved_mapped_ring_gains_aromaticity_without_external_incorporation",
                "redox_direction": "oxidation",
            }
        )

    return out


def _leaving_losses_at(facts: dict[str, Any], center: int) -> list[dict[str, Any]]:
    for d in facts["attachment_deltas"]:
        if int(d["center_map"]) == center:
            return list(d.get("lost", []))
    return []


def detect_elimination(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Strict mapped alpha/beta elimination detector.

    A validated elimination requires a bond-order increase on a specific mapped
    alpha-beta pair, leaving-group loss at one endpoint, and H/attachment loss at the
    other endpoint. Unrelated bond-order changes do not promote elimination.
    """
    out: list[dict[str, Any]] = []
    for bc in facts.get("bond_order_changes", []):
        a, b = int(bc["a"]), int(bc["b"])
        ro, po = float(bc["reactant_order"]), float(bc["product_order"])
        if not ((ro <= 1.2 and po >= 1.8) or (1.8 <= ro < 2.8 and po >= 2.8)):
            continue
        if (
            _atom(facts, "reactant", a).get("element") != "C"
            or _atom(facts, "reactant", b).get("element") != "C"
        ):
            continue
        # Carbonyl olefination is not generic elimination.
        if _element_attachments(facts, "reactant", a, "O", 1.8) or _element_attachments(
            facts, "reactant", b, "O", 1.8
        ):
            continue

        candidates = []
        for alpha, beta in ((a, b), (b, a)):
            losses = _leaving_losses_at(facts, alpha)
            lg = [
                x
                for x in losses
                if x.get("element") in (HALOGENS | {"O", "N", "S", "P"})
            ]
            beta_h_loss = _hydrogen_delta(facts, beta) < 0
            beta_other_loss = bool(_leaving_losses_at(facts, beta))
            if lg:
                candidates.append((alpha, beta, lg, beta_h_loss or beta_other_loss))

        if not candidates:
            continue

        # Like-topology is still mapped and pair-local, but lacks complete beta evidence.
        for alpha, beta, lg, beta_support in candidates:
            subtype_map = {
                "F": "dehydrohalogenation",
                "Cl": "dehydrohalogenation",
                "Br": "dehydrohalogenation",
                "I": "dehydrohalogenation",
                "O": "dehydration_or_dealkoxy_elimination",
                "N": "deamination_elimination",
                "S": "desulfurative_elimination",
                "P": "phosphorus_leaving_group_elimination",
            }
            subtype = subtype_map.get(str(lg[0].get("element")), "beta_elimination")
            support = _support(alpha, beta, *[_att_map(x) for x in lg])
            if not beta_support:
                out.append(
                    {
                        "event": "elimination_like_topology",
                        "role": "topology_evidence",
                        "subtype": subtype,
                        "alpha_map": alpha,
                        "beta_map": beta,
                        "support_maps": support,
                        "promotion_blockers": [
                            "beta_h_or_attachment_loss_not_confirmed"
                        ],
                        "mapping_scope": "same_mapped_alpha_beta_pair",
                    }
                )
                continue
            # External incorporation at the elimination pair argues for substitution/addition.
            if _has_external_formed_at(facts, alpha) or _has_external_formed_at(
                facts, beta
            ):
                out.append(
                    {
                        "event": "elimination_like_topology",
                        "role": "topology_evidence",
                        "subtype": subtype,
                        "alpha_map": alpha,
                        "beta_map": beta,
                        "support_maps": support,
                        "promotion_blockers": [
                            "external_heavy_atom_incorporation_at_elimination_pair"
                        ],
                        "mapping_scope": "same_mapped_alpha_beta_pair",
                    }
                )
                continue
            out.append(
                {
                    "event": "elimination",
                    "role": "specific_event",
                    "subtype": subtype,
                    "alpha_map": alpha,
                    "beta_map": beta,
                    "support_maps": support,
                    "mapping_scope": "validated_mapped_alpha_beta_elimination_unit",
                }
            )
    return out


def _nucleophile_condensation_subtype(facts: dict[str, Any], nmap: int) -> str:
    r_atts = _attachments(facts, "reactant", nmap)
    if any(x.get("element") == "N" for x in r_atts):
        return "hydrazone_condensation"
    if any(x.get("element") == "O" for x in r_atts):
        return "oxime_condensation"
    return "imine_condensation"


def _cnd_map(att: dict[str, Any]) -> int | None:
    m = att.get("map")
    return int(m) if m is not None else None


def _cnd_atom_h(facts: dict[str, Any], side: str, mn: int) -> int:
    return int(_atom(facts, side, mn).get("total_h", 0) or 0)


def _cnd_carbonyl_states() -> set[str]:
    return {
        "aldehyde", "ketone", "ester", "carboxylic_acid", "amide",
        "anhydride", "acyl_halide", "haloformate",
    }


def _cnd_is_carbonyl_carbon(facts: dict[str, Any], side: str, mn: int) -> bool:
    return _carbon_fg_state(facts, side, mn) in _cnd_carbonyl_states()


def _cnd_nitrile_pair(facts: dict[str, Any], side: str, carbon: int) -> tuple[int, float] | None:
    if _atom(facts, side, carbon).get("element") != "C":
        return None
    for x in _attachments(facts, side, carbon):
        if (
            x.get("element") == "N"
            and x.get("map") is not None
            and float(x.get("bond_order", 0.0)) >= 2.8
        ):
            return int(x["map"]), float(x.get("bond_order", 0.0))
    return None


def _cnd_is_isocyanide_carbon(facts: dict[str, Any], carbon: int) -> bool:
    pair = _cnd_nitrile_pair(facts, "reactant", carbon)
    if pair is None:
        return False
    nmap, _ = pair
    heavy_other = [
        x for x in _attachments(facts, "reactant", nmap)
        if x.get("map") is not None
        and int(x["map"]) != carbon
        and x.get("element") != "H"
    ]
    return bool(heavy_other)


def _cnd_is_cyanide_carbon(facts: dict[str, Any], carbon: int) -> bool:
    return _cnd_nitrile_pair(facts, "reactant", carbon) is not None and not _cnd_is_isocyanide_carbon(facts, carbon)


def _cnd_is_amine_like_n(facts: dict[str, Any], nmap: int, *, allow_acyl_nh: bool = True) -> bool:
    if _atom(facts, "reactant", nmap).get("element") != "N":
        return False

    for x in _attachments(facts, "reactant", nmap):
        if x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 2.8:
            return False

    attached_to_carbonyl = False
    for x in _attachments(facts, "reactant", nmap):
        if x.get("element") != "C" or x.get("map") is None:
            continue
        if _cnd_is_carbonyl_carbon(facts, "reactant", int(x["map"])):
            attached_to_carbonyl = True
            break

    if not attached_to_carbonyl:
        return True
    return bool(allow_acyl_nh and _cnd_atom_h(facts, "reactant", nmap) > 0)


def _cnd_is_alcohol_or_thiol_hetero(facts: dict[str, Any], hmap: int) -> bool:
    return (
        _atom(facts, "reactant", hmap).get("element") in {"O", "S"}
        and _cnd_atom_h(facts, "reactant", hmap) > 0
    )


def _cnd_original_carbonyl_oxygen(facts: dict[str, Any], center: int) -> int | None:
    for x in _attachments(facts, "reactant", center):
        if (
            x.get("element") == "O"
            and x.get("map") is not None
            and float(x.get("bond_order", 0.0)) >= 1.8
        ):
            return int(x["map"])
    return None


def _cnd_oxygen_fate(facts: dict[str, Any], center: int, oxygen: int | None) -> str:
    if oxygen is None:
        return "OXYGEN_UNMAPPED"
    if _bond_between(facts, "product", center, oxygen) >= 1.8:
        return "RETAINED_AS_ORIGINAL_CARBONYL"

    p_atom = _atom(facts, "product", oxygen)
    if not p_atom:
        return "OXYGEN_MISSING"

    heavy = [x for x in _attachments(facts, "product", oxygen) if x.get("element") != "H"]
    if not heavy:
        return "FREE_O_WATER_LIKE"

    for x in heavy:
        if x.get("element") != "C" or x.get("map") is None:
            continue
        cm = int(x["map"])
        if _bond_between(facts, "product", oxygen, cm) >= 1.8 and cm != center:
            return "TRANSFERRED_TO_OTHER_CARBONYL"

    return "OTHER_RETAINED_OXYGEN"


def _cnd_gained_at_center(facts: dict[str, Any], center: int) -> list[dict[str, Any]]:
    out = []
    for x in _attachments(facts, "product", center):
        m = _cnd_map(x)
        if m is not None and _bond_between(facts, "reactant", center, m) < 0.2:
            out.append(x)
    return out


def _cnd_lost_direct_handle_elements(facts: dict[str, Any], center: int) -> set[str]:
    lost: set[str] = set()
    for d in facts.get("attachment_deltas", []) or []:
        if int(d.get("center_map", -1)) != int(center):
            continue
        for x in d.get("lost", []) or []:
            el = str(x.get("element") or "")
            if el in {"P", "S", "Si"}:
                lost.add(el)
    return lost


def _cnd_retained_direct_handle_elements(facts: dict[str, Any], center: int) -> set[str]:
    out: set[str] = set()
    for x in _attachments(facts, "reactant", center):
        if x.get("element") not in {"P", "S", "Si"} or x.get("map") is None:
            continue
        hm = int(x["map"])
        if _bond_between(facts, "product", center, hm) > 0.2:
            out.add(str(x["element"]))
    return out


def _cnd_neighbor_ewg_type(facts: dict[str, Any], carbon: int) -> str | None:
    state = _carbon_fg_state(facts, "reactant", carbon)
    if state in {"aldehyde", "ketone"}:
        return "ALDEHYDE_KETONE"
    if state in {"ester", "carboxylic_acid", "amide", "anhydride", "acyl_halide", "haloformate"}:
        return "ACYL_CARBONYL"
    if state == "nitrile":
        return "NITRILE"
    if any(
        x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 1.7
        for x in _attachments(facts, "reactant", carbon)
    ):
        return "IMINE"
    if any(
        x.get("element") == "S" and float(x.get("bond_order", 0.0)) >= 1.7
        for x in _attachments(facts, "reactant", carbon)
    ):
        return "THIOCARBONYL"
    return None


def _cnd_is_nitro_n(facts: dict[str, Any], nmap: int) -> bool:
    if _atom(facts, "reactant", nmap).get("element") != "N":
        return False
    oxy = [x for x in _attachments(facts, "reactant", nmap) if x.get("element") == "O"]
    return len(oxy) >= 2 and sum(float(x.get("bond_order", 0.0)) for x in oxy) >= 2.5


def _cnd_vinylogous_ewg(facts: dict[str, Any], donor: int, max_depth: int = 3) -> bool:
    frontier = [(donor, 0)]
    seen = {donor}
    while frontier:
        node, depth = frontier.pop(0)
        if depth >= max_depth:
            continue
        for x in _attachments(facts, "reactant", node):
            if x.get("element") != "C" or x.get("map") is None:
                continue
            nxt = int(x["map"])
            if nxt in seen:
                continue
            seen.add(nxt)
            nd = depth + 1
            if nd >= 2 and _cnd_neighbor_ewg_type(facts, nxt) is not None:
                return True
            frontier.append((nxt, nd))
    return False


def _cnd_conjugated_active_carbon(facts: dict[str, Any], donor: int) -> bool:
    branches = []
    for x in _attachments(facts, "reactant", donor):
        if x.get("element") != "C" or x.get("map") is None:
            continue
        cm = int(x["map"])
        branch_pi = False
        for q in _attachments(facts, "reactant", cm):
            qm = _cnd_map(q)
            if qm is None or qm == donor or q.get("element") != "C":
                continue
            if float(q.get("bond_order", 0.0)) >= 1.7 or bool(_atom(facts, "reactant", qm).get("aromatic")):
                branch_pi = True
                break
        if branch_pi:
            branches.append(cm)
    return len(set(branches)) >= 2


def _cnd_f1_activation_reasons(facts: dict[str, Any], donor: int) -> list[str]:
    reasons: list[str] = []
    for x in _attachments(facts, "reactant", donor):
        m = _cnd_map(x)
        el = x.get("element")
        if m is None:
            continue
        if el == "C":
            t = _cnd_neighbor_ewg_type(facts, m)
            if t is not None:
                reasons.append(f"ALPHA_TO_{t}")
            if bool(_atom(facts, "reactant", m).get("aromatic")):
                reasons.append("BENZYLIC_OR_HETEROARYL_ACTIVE_CARBON")
        elif el == "N" and _cnd_is_nitro_n(facts, m):
            reasons.append("NITRO_ACTIVE_CARBON")
        elif el in {"O", "N", "S"}:
            reasons.append("HETEROATOM_ADJACENT_ACTIVE_CARBON")

    if _cnd_retained_direct_handle_elements(facts, donor):
        reasons.append("RETAINED_EXPLICIT_EWG_ACTIVE_CARBON")
    if _cnd_vinylogous_ewg(facts, donor):
        reasons.append("VINYLOGOUS_EWG_ACTIVE_CARBON")
    if _cnd_conjugated_active_carbon(facts, donor):
        reasons.append("CONJUGATED_CYCLOPENTADIENE_FULVENE_LIKE")
    return list(dict.fromkeys(reasons))


def _cnd_allenyl_product(facts: dict[str, Any], donor: int) -> bool:
    dbl_c = [
        x for x in _attachments(facts, "product", donor)
        if x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 1.7
    ]
    return len(dbl_c) >= 2


def _cnd_ring_opening_without_pair_reclosure(facts: dict[str, Any], a: int, b: int) -> bool:
    touched_open = False
    for rec in facts.get("ring_openings", []) or []:
        maps = {int(x) for x in (rec.get("conserved_maps") or rec.get("ring_maps") or [])}
        if a in maps or b in maps:
            touched_open = True
            break
    if not touched_open:
        return False
    for rec in facts.get("ring_formations", []) or []:
        maps = {int(x) for x in (rec.get("conserved_maps") or rec.get("ring_maps") or [])}
        if a in maps and b in maps:
            return False
    return True


def _cnd_same_product_ring(facts: dict[str, Any], maps: set[int], size: int | None = None) -> bool:
    for rec in facts.get("ring_formations", []) or []:
        ring = {int(x) for x in (rec.get("ring_maps") or rec.get("conserved_maps") or [])}
        if size is not None and len(ring) != size:
            continue
        if maps.issubset(ring):
            return True
    return False


def _cnd_carbonic_acid_like(facts: dict[str, Any], c: int) -> bool:
    if _atom(facts, "reactant", c).get("element") != "C":
        return False
    atts = _attachments(facts, "reactant", c)
    oxy = [x for x in atts if x.get("element") == "O"]
    carbons = [x for x in atts if x.get("element") == "C"]
    return len(oxy) >= 2 and not carbons and any(float(x.get("bond_order", 0.0)) >= 1.8 for x in oxy)


def _cnd_product_carbonyl_like(facts: dict[str, Any], c: int) -> bool:
    return bool(_element_attachments(facts, "product", c, "O", 1.8))


def _cnd_emit(
    out: list[dict[str, Any]],
    *,
    family: str,
    subtype: str,
    support: list[int],
    mapping_scope: str,
    **extra: Any,
) -> None:
    rec = {
        "event": "condensation",
        "role": "specific_event",
        "family": family,
        "subtype": subtype,
        "support_maps": support,
        "mapping_scope": mapping_scope,
        "derived_parent": "condensation",
    }
    rec.update(extra)
    out.append(rec)


def detect_condensation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    # Runtime v6.7 audited condensation detector.
    out: list[dict[str, Any]] = []
    common = [int(x) for x in facts.get("mapping", {}).get("common_maps", [])]

    # F3: aldehyde/ketone + heteronucleophile.
    for c in common:
        if _carbon_fg_state(facts, "reactant", c) not in {"aldehyde", "ketone"}:
            continue
        if _element_attachments(facts, "product", c, "O", 1.8):
            continue

        for n in _attachments(facts, "product", c):
            if (
                n.get("element") == "N"
                and n.get("map") is not None
                and 1.7 <= float(n.get("bond_order", 0.0)) <= 2.2
            ):
                nmap = int(n["map"])
                if _bond_between(facts, "reactant", c, nmap) <= 0.2:
                    _cnd_emit(
                        out, family="F3",
                        subtype=_nucleophile_condensation_subtype(facts, nmap),
                        center_map=c, nucleophile_map=nmap,
                        support=_support(c, nmap),
                        mapping_scope="v67_f3_aldehyde_ketone_to_new_c_n_double_bond",
                    )

        single_os = [
            x for x in _attachments(facts, "product", c)
            if x.get("element") in {"O", "S"}
            and x.get("map") is not None
            and 0.8 <= float(x.get("bond_order", 0.0)) <= 1.2
            and _cnd_atom_h(facts, "product", int(x["map"])) == 0
        ]
        if len(single_os) >= 2:
            elems = [str(x.get("element")) for x in single_os]
            hmaps = [int(x["map"]) for x in single_os]
            subtype = (
                "acetal_or_ketal_condensation" if elems.count("O") >= 2
                else "thioacetal_or_thioketal_condensation" if elems.count("S") >= 2
                else "mixed_o_s_acetal_condensation"
            )
            _cnd_emit(
                out, family="F3", subtype=subtype,
                center_map=c, support=_support(c, *hmaps),
                mapping_scope="v67_f3_completed_o_s_acetal_only",
            )

    # F1: activated-carbon -> aldehyde/ketone alkene condensation.
    for c in common:
        if _carbon_fg_state(facts, "reactant", c) not in {"aldehyde", "ketone"}:
            continue
        if _element_attachments(facts, "product", c, "O", 1.8):
            continue
        for x in _attachments(facts, "product", c):
            if (
                x.get("element") != "C"
                or x.get("map") is None
                or float(x.get("bond_order", 0.0)) < 1.7
            ):
                continue
            donor = int(x["map"])
            if _bond_between(facts, "reactant", c, donor) > 0.2:
                continue
            if _carbon_fg_state(facts, "reactant", donor) in _cnd_carbonyl_states() | {"nitrile"}:
                continue
            if _cnd_lost_direct_handle_elements(facts, donor):
                continue
            if _cnd_allenyl_product(facts, donor):
                continue
            if _cnd_ring_opening_without_pair_reclosure(facts, c, donor):
                continue

            activation = _cnd_f1_activation_reasons(facts, donor)
            if not activation:
                continue
            _cnd_emit(
                out, family="F1",
                subtype="activated_carbon_aldehyde_ketone_condensation",
                center_map=c, donor_map=donor, donor_activation=activation,
                support=_support(c, donor),
                mapping_scope="v67_f1_consumed_carbonyl_new_c_c_pi_activated_retained_donor",
            )

    # F2A: enolate-driven acyl condensation.
    for acyl in common:
        rstate = _carbon_fg_state(facts, "reactant", acyl)
        if rstate not in {"ester", "carboxylic_acid", "acyl_halide", "anhydride"}:
            continue
        if not _element_attachments(facts, "product", acyl, "O", 1.8):
            continue

        lost_single = []
        for d in facts.get("attachment_deltas", []) or []:
            if int(d.get("center_map", -1)) != acyl:
                continue
            for q in d.get("lost", []) or []:
                if q.get("element") in {"O", "Cl", "Br", "I", "F"} and float(q.get("bond_order", 0.0)) <= 1.2:
                    lost_single.append(q)
        if not lost_single:
            continue

        for q in _cnd_gained_at_center(facts, acyl):
            if q.get("element") != "C" or q.get("map") is None:
                continue
            donor = int(q["map"])
            alpha_carbonyl = any(
                a.get("element") == "C"
                and a.get("map") is not None
                and int(a["map"]) != acyl
                and _cnd_is_carbonyl_carbon(facts, "reactant", int(a["map"]))
                for a in _attachments(facts, "reactant", donor)
            )
            if not alpha_carbonyl:
                continue

            inter = any(
                {int(b["a"]), int(b["b"])} == {acyl, donor}
                for b in facts.get("intercomponent_formed_bonds", []) or []
            )
            ring = _cnd_same_product_ring(facts, {acyl, donor})
            if not inter and not ring:
                continue

            _cnd_emit(
                out, family="F2A",
                subtype="dieckmann_enolate_acyl_cyclization" if ring and not inter else "enolate_acyl_condensation",
                center_map=acyl, donor_map=donor,
                support=_support(acyl, donor),
                mapping_scope="v67_f2a_retained_acyl_carbonyl_lg_loss_alpha_carbon_c_c_formation",
            )

    # F2B: Thorpe / Thorpe-Ziegler.
    for acceptor in common:
        pair = _cnd_nitrile_pair(facts, "reactant", acceptor)
        if pair is None or _cnd_is_isocyanide_carbon(facts, acceptor):
            continue
        nmap, _ = pair
        p_cn = _bond_between(facts, "product", acceptor, nmap)
        if p_cn >= 2.5:
            continue

        for q in _cnd_gained_at_center(facts, acceptor):
            if q.get("element") != "C" or q.get("map") is None:
                continue
            donor = int(q["map"])
            alpha_to_nitrile = any(
                a.get("element") == "C"
                and a.get("map") is not None
                and _carbon_fg_state(facts, "reactant", int(a["map"])) == "nitrile"
                for a in _attachments(facts, "reactant", donor)
            )
            if not alpha_to_nitrile:
                continue

            p_cc = _bond_between(facts, "product", acceptor, donor)
            if 1.7 <= p_cn <= 2.2 and p_cc <= 1.2:
                tautomer = "imine_tautomer"
            elif p_cn <= 1.2 and p_cc >= 1.7:
                tautomer = "enamine_tautomer"
            else:
                continue

            inter = any(
                {int(b["a"]), int(b["b"])} == {acceptor, donor}
                for b in facts.get("intercomponent_formed_bonds", []) or []
            )
            ring = _cnd_same_product_ring(facts, {acceptor, donor})
            if not inter and not ring:
                continue

            _cnd_emit(
                out, family="F2B",
                subtype=f"thorpe_ziegler_{tautomer}" if ring and not inter else f"thorpe_intercomponent_{tautomer}",
                center_map=acceptor, donor_map=donor, nitrile_n_map=nmap,
                support=_support(acceptor, donor, nmap),
                mapping_scope="v67_f2b_consumed_nitrile_alpha_nitrile_carbon_c_c_formation",
            )

    # F4 local: N/O/S-assisted carbonyl capture.
    consumed_centers: list[dict[str, Any]] = []
    for c in common:
        if _carbon_fg_state(facts, "reactant", c) not in {"aldehyde", "ketone"}:
            continue
        original_o = _cnd_original_carbonyl_oxygen(facts, c)
        if _element_attachments(facts, "product", c, "O", 1.8):
            continue

        gained = _cnd_gained_at_center(facts, c)
        fate = _cnd_oxygen_fate(facts, c, original_o)
        consumed_centers.append({"center": c, "oxygen": original_o, "fate": fate, "gained": gained})

        if fate not in {"FREE_O_WATER_LIKE", "OXYGEN_MISSING", "OXYGEN_UNMAPPED"}:
            continue

        gained_n = [
            q for q in gained
            if q.get("element") == "N"
            and q.get("map") is not None
            and _cnd_is_amine_like_n(facts, int(q["map"]), allow_acyl_nh=True)
        ]
        gained_c = [q for q in gained if q.get("element") == "C" and q.get("map") is not None]
        gained_p = [q for q in gained if q.get("element") == "P" and q.get("map") is not None]

        if gained_n and (gained_c or gained_p):
            nmap = int(gained_n[0]["map"])
            partner = int((gained_c or gained_p)[0]["map"])
            _cnd_emit(
                out, family="F4",
                subtype="dehydrative_n_nucleophile_assisted_carbonyl_capture",
                center_map=c, nucleophile_map=nmap, capture_partner_map=partner,
                oxygen_fate=fate, support=_support(c, nmap, partner),
                mapping_scope="v67_f4_dehydrative_n_plus_c_or_p_capture",
            )

        gained_os = [
            q for q in gained
            if q.get("element") in {"O", "S"}
            and q.get("map") is not None
            and _cnd_is_alcohol_or_thiol_hetero(facts, int(q["map"]))
        ]
        if gained_os and gained_c:
            hmap = int(gained_os[0]["map"])
            cmap = int(gained_c[0]["map"])
            _cnd_emit(
                out, family="F4",
                subtype="dehydrative_o_s_assisted_carbonyl_capture",
                center_map=c, nucleophile_map=hmap, capture_partner_map=cmap,
                oxygen_fate=fate, support=_support(c, hmap, cmap),
                mapping_scope="v67_f4_dehydrative_o_s_plus_c_capture",
            )

    # F4C2 Ugi-like.
    for ctr in consumed_centers:
        c = int(ctr["center"])
        gained = ctr["gained"]
        amines = [
            int(q["map"]) for q in gained
            if q.get("element") == "N" and q.get("map") is not None
            and _cnd_is_amine_like_n(facts, int(q["map"]), allow_acyl_nh=False)
        ]
        isoc = [
            int(q["map"]) for q in gained
            if q.get("element") == "C" and q.get("map") is not None
            and _cnd_is_isocyanide_carbon(facts, int(q["map"]))
        ]

        for nmap in amines:
            for iso_c in isoc:
                if not _cnd_product_carbonyl_like(facts, iso_c):
                    continue
                acid_maps = []
                for q in _attachments(facts, "product", nmap):
                    if q.get("element") != "C" or q.get("map") is None:
                        continue
                    am = int(q["map"])
                    if _bond_between(facts, "reactant", nmap, am) <= 0.2 and _carbon_fg_state(facts, "reactant", am) == "carboxylic_acid":
                        acid_maps.append(am)
                if not acid_maps:
                    continue

                _cnd_emit(
                    out, family="F4",
                    subtype="mcr_supported_n_nucleophile_assisted_carbonyl_capture",
                    center_map=c, nucleophile_map=nmap, capture_partner_map=iso_c,
                    oxygen_fate=ctr["fate"], support=_support(c, nmap, iso_c, *acid_maps),
                    mapping_scope="v67_f4_local_capture_supported_by_ugi",
                )
                _cnd_emit(
                    out, family="F4C2", subtype="ugi_like_isocyanide_mcr",
                    center_map=c, amine_n_map=nmap, isocyanide_c_map=iso_c,
                    acid_maps=sorted(acid_maps), architecture="multicomponent_reaction",
                    support=_support(c, nmap, iso_c, *acid_maps),
                    mapping_scope="v67_f4c2_carbonyl_amine_isocyanide_acid_architecture",
                )

    # F4C3 Bucherer-Bergs-like; same five-member ring invariant.
    cyanides: list[tuple[int, int]] = []
    amines_all: list[int] = []
    carbonic: list[int] = []
    for m in common:
        pair = _cnd_nitrile_pair(facts, "reactant", m)
        if pair is not None and _cnd_is_cyanide_carbon(facts, m):
            cyanides.append((m, int(pair[0])))
        if _atom(facts, "reactant", m).get("element") == "N" and _cnd_is_amine_like_n(facts, m, allow_acyl_nh=False):
            amines_all.append(m)
        if _cnd_carbonic_acid_like(facts, m):
            carbonic.append(m)

    for ctr in consumed_centers:
        c = int(ctr["center"])
        for cy_c, cy_n in cyanides:
            if not _cnd_product_carbonyl_like(facts, cy_c):
                continue
            for am_n in amines_all:
                for c1 in carbonic:
                    participants = {c, cy_c, cy_n, am_n, c1}
                    if not _cnd_product_carbonyl_like(facts, c1):
                        continue
                    if not _cnd_same_product_ring(facts, participants, size=5):
                        continue

                    _cnd_emit(
                        out, family="F4",
                        subtype="mcr_supported_n_nucleophile_assisted_carbonyl_capture",
                        center_map=c, nucleophile_map=am_n, capture_partner_map=cy_c,
                        oxygen_fate=ctr["fate"], support=_support(*sorted(participants)),
                        mapping_scope="v67_f4_local_capture_supported_by_bucherer_bergs",
                    )
                    _cnd_emit(
                        out, family="F4C3", subtype="bucherer_bergs_like_hydantoin_mcr",
                        center_map=c, amine_n_map=am_n, cyanide_c_map=cy_c,
                        cyanide_n_map=cy_n, carbonic_c_map=c1,
                        architecture="multicomponent_ring_condensation",
                        support=_support(*sorted(participants)),
                        mapping_scope="v67_f4c3_five_participant_same_ring_invariant",
                    )

    # F4C1 Biginelli-like; deliberately strict.
    dehydrative_centers = [
        x for x in consumed_centers
        if x["fate"] in {"FREE_O_WATER_LIKE", "OXYGEN_MISSING", "OXYGEN_UNMAPPED"}
    ]
    f4a_centers = []
    for x in dehydrative_centers:
        gained = x["gained"]
        has_n = any(
            q.get("element") == "N" and q.get("map") is not None
            and _cnd_is_amine_like_n(facts, int(q["map"]), allow_acyl_nh=True)
            for q in gained
        )
        has_cp = any(q.get("element") in {"C", "P"} for q in gained)
        if has_n and has_cp:
            f4a_centers.append(int(x["center"]))

    six_hetero_rings = []
    for rec in facts.get("ring_formations", []) or []:
        ring = [int(x) for x in (rec.get("ring_maps") or rec.get("conserved_maps") or [])]
        if len(ring) == 6 and any(_atom(facts, "product", m).get("element") == "N" for m in ring):
            six_hetero_rings.append(ring)

    if len(dehydrative_centers) >= 2 and f4a_centers and six_hetero_rings:
        ring = six_hetero_rings[0]
        _cnd_emit(
            out, family="F4C1", subtype="biginelli_like_multicenter_cyclocondensation",
            architecture="multicenter_ring_condensation", ring_maps=ring,
            center_maps=sorted(int(x["center"]) for x in dehydrative_centers),
            support=_support(*ring),
            mapping_scope="v67_f4c1_multi_dehydrative_centers_n_capture_six_member_n_ring",
        )

    return out


def _changed_maps(facts: dict[str, Any]) -> set[int]:
    maps: set[int] = set()
    for k in ("formed_bonds", "broken_bonds", "bond_order_changes"):
        for x in facts.get(k, []):
            maps.update([int(x["a"]), int(x["b"])])
    return maps


def _mapped_adj(
    facts: dict[str, Any], side: str, nodes: set[int] | None = None
) -> dict[int, set[int]]:
    nodes = nodes or set(facts["mapping"]["common_maps"])
    adj = {n: set() for n in nodes}
    for (a, b), _bo in facts[side]["bonds"].items():
        if a in nodes and b in nodes:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def _shortest_path(
    adj: dict[int, set[int]], start: int, end: int, max_edges: int = 7
) -> list[int] | None:
    q = [(start, [start])]
    seen = {start}
    while q:
        cur, path = q.pop(0)
        if len(path) - 1 > max_edges:
            continue
        if cur == end:
            return path
        for nb in adj.get(cur, set()):
            if nb not in seen:
                seen.add(nb)
                q.append((nb, path + [nb]))
    return None


def _alternating_pi_path(facts: dict[str, Any], path: list[int]) -> bool:
    bos = [
        _bond_between(facts, "reactant", path[i], path[i + 1])
        for i in range(len(path) - 1)
    ]
    if len(path) not in {4, 6, 8}:
        return False
    # A 4/6/8-atom 4pi/6pi/8pi electron chain must carry the terminal
    # pi bonds explicitly (2-1-2, 2-1-2-1-2, ...).  A 1-2-1 path is
    # geometrically alternating but does NOT identify a 4pi diene.
    return all(
        (1.8 <= bo <= 2.2 if i % 2 == 0 else 0.8 <= bo <= 1.2)
        for i, bo in enumerate(bos)
    )


def _pericyclic_pi_shift_matches(facts: dict[str, Any], path: list[int]) -> bool:
    """Demand classical alternating pi migration on this *mapped* chain.

    This deliberately excludes aromatized/composite endpoint records whose
    intermediate pericyclic step cannot be recovered uniquely.
    """
    if not _alternating_pi_path(facts, path):
        return False
    for i in range(len(path) - 1):
        order = _bond_between(facts, "product", path[i], path[i + 1])
        if i % 2 == 0:
            if not 0.8 <= order <= 1.2:
                return False
        elif not 1.8 <= order <= 2.2:
            return False
    return True


def _pericyclic_pi_consumed(facts: dict[str, Any], pair: list[int], *, aromatic_ok: bool = False) -> bool:
    if len(pair) != 2:
        return False
    react = _bond_between(facts, "reactant", *pair)
    prod = _bond_between(facts, "product", *pair)
    return (react >= 1.8 and prod >= 0.8 and react - prod >= (0.4 if aromatic_ok else 0.7))


def _pericyclic_preformed_dipole(facts: dict[str, Any], path: list[int]) -> bool:
    """Whitelist explicit preformed 1,3-dipolar precursors, not generated-in-situ synthon chemistry.

    This is deliberately incomplete; other dipoles remain deferred until
    curated positive and negative controls establish their graph contract.
    """
    if len(path) != 3:
        return False
    states = [_atom(facts, "reactant", n) for n in path]
    elements = [a.get("element") for a in states]
    charges = [int(a.get("charge", 0)) for a in states]
    orders = [_bond_between(facts, "reactant", path[i], path[i + 1]) for i in range(2)]
    # Conventional covalent azide R-N=N+=N-, either map direction.
    if elements == ["N", "N", "N"]:
        return charges[1] > 0 and (charges[0] < 0 or charges[2] < 0) and all(o >= 1.8 for o in orders)
    # An *already present* nitrile oxide C#N+-O- (not an oxime precursor).
    if elements in (["C", "N", "O"], ["O", "N", "C"]):
        neg = 0 if elements[0] == "O" else 2
        cn_bond_index = 0 if elements[0] == "C" else 1
        return charges[1] > 0 and charges[neg] < 0 and orders[cn_bond_index] >= 2.8
    return False


def detect_pericyclic(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Recognize only explicit mapped *formal pericyclic-compatible* topologies.

    This does not assert a concerted mechanism. Generic ring formation + pi addition is
    insufficient. The mapped changes must form a recognized cycloaddition,
    electrocyclization, or [3,3]-sigmatropic edit pattern.
    """
    out: list[dict[str, Any]] = []
    if not facts["mapping"]["valid"]:
        return out

    # Formal cycloaddition: exactly two new intercomponent sigma bonds between
    # the same two reactant components, mapped cycle-rank gain, and no mapped bond
    # cleavage.  Class is inferred from the mapped reactant paths between the two
    # connection endpoints on each component (2+2, 3+2, 4+2), not from generic
    # "ring formation + pi addition".
    groups: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for fb in facts.get("intercomponent_formed_bonds", []):
        key = tuple(
            sorted((int(fb["reactant_component_a"]), int(fb["reactant_component_b"])))
        )
        groups.setdefault(key, []).append(fb)
    if int(facts.get("common_cycle_rank_delta", 0)) > 0 and not facts.get(
        "broken_bonds"
    ):
        radj = _mapped_adj(facts, "reactant")
        for comps, fbs in groups.items():
            if len(fbs) != 2:
                continue
            endpoints = {comps[0]: [], comps[1]: []}
            for fb in fbs:
                for m, ci in (
                    (int(fb["a"]), int(fb["reactant_component_a"])),
                    (int(fb["b"]), int(fb["reactant_component_b"])),
                ):
                    if ci in endpoints:
                        endpoints[ci].append(m)
            if any(len(set(v)) != 2 for v in endpoints.values()):
                continue
            paths = {}
            for ci, vals in endpoints.items():
                a, b = sorted(set(vals))
                path = _shortest_path(radj, a, b, max_edges=5)
                if not path or any(
                    facts["reactant"]["map_to_component"].get(x) != ci for x in path
                ):
                    paths = {}
                    break
                paths[ci] = path
            if len(paths) != 2:
                continue
            lens = sorted(len(p) for p in paths.values())
            subtype = None
            # [2+2]: each reacting component contributes a mapped alkene pair.
            if lens == [2, 2] and all(
                _pericyclic_pi_consumed(facts, p)
                for p in paths.values()
            ) and all(0.8 <= float(fb["order"]) <= 1.2 for fb in fbs):
                subtype = "formal_2plus2_cycloaddition_compatible"
            # [3+2]: one alkene pair + a 3-atom heteroatom-containing conjugated/dipolar path.
            elif lens == [2, 3]:
                p2 = next(p for p in paths.values() if len(p) == 2)
                p3 = next(p for p in paths.values() if len(p) == 3)
                if (
                    _pericyclic_pi_consumed(facts, p2, aromatic_ok=True)
                    and _pericyclic_preformed_dipole(facts, p3)
                    and any(
                        _bond_between(facts, "reactant", p3[i], p3[i + 1])
                        > _bond_between(facts, "product", p3[i], p3[i + 1]) + 0.4
                        for i in range(2)
                    )
                ):
                    subtype = "formal_3plus2_cycloaddition_compatible"
            # [4+2]: alkene pair + alternating four-atom diene path.
            elif lens == [2, 4]:
                p2 = next(p for p in paths.values() if len(p) == 2)
                p4 = next(p for p in paths.values() if len(p) == 4)
                if _bond_between(
                    facts, "reactant", p2[0], p2[1]
                ) >= 1.8 and _pericyclic_pi_consumed(facts, p2) and _pericyclic_pi_shift_matches(facts, p4) and all(
                    0.8 <= float(fb["order"]) <= 1.2 for fb in fbs
                ):
                    subtype = "formal_4plus2_cycloaddition_compatible"
            if subtype:
                support = sorted({m for p in paths.values() for m in p})
                out.append(
                    {
                        "event": "pericyclic_compatible",
                        "role": "specific_event",
                        "subtype": subtype,
                        "support_maps": support,
                        "mechanism_claim": "formal_topology_only_not_concertedness",
                        "mapping_scope": "mapped_two_component_two_bond_cycloaddition_circuit",
                    }
                )

    # Electrocyclization: intracomponent terminal bond closes a mapped alternating
    # 4pi/6pi/8pi chain with coherent bond-order reorganization.
    if (
        int(facts.get("common_cycle_rank_delta", 0)) > 0
        and not facts.get("intercomponent_formed_bonds")
        and not facts.get("broken_bonds")
    ):
        adj = _mapped_adj(facts, "reactant")
        for fb in facts.get("intracomponent_formed_bonds", []):
            a, b = int(fb["a"]), int(fb["b"])
            path = _shortest_path(adj, a, b, max_edges=7)
            if (
                path and _pericyclic_pi_shift_matches(facts, path)
                and 0.8 <= float(fb["order"]) <= 1.2
            ):
                changed = _changed_maps(facts)
                if len(changed & set(path)) >= 4:
                    out.append(
                        {
                            "event": "pericyclic_compatible",
                            "role": "specific_event",
                            "subtype": f"electrocyclization_{len(path)}_atom_conjugated_chain",
                            "support_maps": path,
                            "mechanism_claim": "formal_topology_only_not_concertedness",
                            "mapping_scope": "mapped_alternating_pi_chain_terminal_ring_closure",
                        }
                    )

    # [3,3]-sigmatropic-compatible: exactly one mapped sigma bond broken and one
    # formed, no component exchange/cycle change, and a connected six-atom changed
    # framework with at least two pi-bond shifts.
    if (
        int(facts.get("common_cycle_rank_delta", 0)) == 0
        and not facts.get("intercomponent_formed_bonds")
        and not facts.get("component_splits")
        and not facts.get("component_merges")
        and len(facts.get("broken_bonds", [])) == 1
        and len(facts.get("formed_bonds", [])) == 1
    ):
        br = facts["broken_bonds"][0]
        fo = facts["formed_bonds"][0]
        if float(br.get("order", 0.0)) <= 1.2 and float(fo.get("order", 0.0)) <= 1.2:
            bo_changes = [
                x
                for x in facts.get("bond_order_changes", [])
                if abs(float(x["product_order"]) - float(x["reactant_order"])) >= 0.7
            ]
            changed = _changed_maps(facts)
            if len(bo_changes) >= 2 and 5 <= len(changed) <= 7:
                adj = _mapped_adj(facts, "reactant", changed)
                # connected changed framework
                start = next(iter(changed))
                stack = [start]
                seen = {start}
                while stack:
                    cur = stack.pop()
                    for nb in adj.get(cur, set()):
                        if nb not in seen:
                            seen.add(nb)
                            stack.append(nb)
                if seen == changed:
                    out.append(
                        {
                            "event": "pericyclic_compatible",
                            "role": "specific_event",
                            "subtype": "sigmatropic_3_3_compatible",
                            "support_maps": sorted(changed),
                            "mechanism_claim": "formal_topology_only_not_concertedness",
                            "mapping_scope": "mapped_sigma_break_form_plus_connected_pi_shift_framework",
                        }
                    )
    return out


def detect_aromatic_functionalization(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Mapped aromatic C-H functionalization on a conserved aromatic carbon."""
    out: list[dict[str, Any]] = []
    for center in facts["mapping"]["common_maps"]:
        ra = _atom(facts, "reactant", center)
        pa = _atom(facts, "product", center)
        if ra.get("element") != "C" or not ra.get("aromatic") or not pa.get("aromatic"):
            continue
        if int(ra.get("total_h", 0)) <= int(pa.get("total_h", 0)):
            continue
        gains = []
        for d in facts["attachment_deltas"]:
            if int(d["center_map"]) == center:
                gains = list(d.get("gained", []))
                break
        for g in gains:
            gm = _att_map(g)
            if gm is None or not _external_gain(facts, center, g):
                continue
            el = _att_element(g)
            subtype = "aromatic_c_h_functionalization"
            if el == "B":
                subtype = "aromatic_borylation"
            elif el in HALOGENS:
                subtype = "aromatic_halogenation"
            elif el == "N":
                subtype = "aromatic_amination"
            elif el == "O":
                subtype = "aromatic_oxygenation"
            elif el == "S":
                subtype = "aromatic_sulfur_functionalization"
            elif el == "C":
                gst = _carbon_fg_state(facts, "product", gm)
                if gst == "aldehyde":
                    subtype = "aromatic_formylation"
                elif gst == "carboxylic_acid":
                    subtype = "aromatic_carboxylation"
                else:
                    subtype = "aromatic_c_c_functionalization"
            out.append(
                {
                    "event": "aromatic_functionalization",
                    "role": "specific_event",
                    "subtype": subtype,
                    "center_map": center,
                    "incoming_map": gm,
                    "support_maps": _support(center, gm),
                    "mapping_scope": "conserved_mapped_aromatic_c_h_to_heavy_atom_attachment",
                }
            )
    return out


def detect_pi_reorganizations(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """High-confidence bond-order reorganizations missing from earlier rule packs."""
    out: list[dict[str, Any]] = []

    # Local C=C migration: one C-C double bond is lost while an adjacent C-C double
    # bond is formed, without atom incorporation, cleavage, or attachment exchange.
    changes = list(facts.get("bond_order_changes", []))
    if (
        len(changes) == 2
        and not facts.get("formed_bonds")
        and not facts.get("broken_bonds")
        and not facts.get("attachment_deltas")
    ):
        inc = [
            b
            for b in changes
            if float(b["reactant_order"]) <= 1.2
            and 1.8 <= float(b["product_order"]) < 2.8
        ]
        dec = [
            b
            for b in changes
            if 1.8 <= float(b["reactant_order"]) < 2.8
            and float(b["product_order"]) <= 1.2
        ]
        if len(inc) == 1 and len(dec) == 1:
            im = {int(inc[0]["a"]), int(inc[0]["b"])}
            dm = {int(dec[0]["a"]), int(dec[0]["b"])}
            support = sorted(im | dm)
            if len(im & dm) == 1 and all(
                _atom(facts, "reactant", m).get("element") == "C" for m in support
            ):
                out.append(
                    {
                        "event": "double_bond_migration",
                        "role": "specific_event",
                        "subtype": "alkene_isomerization",
                        "support_maps": support,
                        "mapping_scope": "adjacent_mapped_c_c_single_double_bond_exchange_without_heavy_atom_edits",
                    }
                )

    # Partial hydrogenation/dearomatization of a heteroarene: multiple aromatic bonds
    # become single, mapped atoms lose aromaticity, and the same scaffold gains H.
    aromatic_reductions = [
        b
        for b in changes
        if 1.4 <= float(b["reactant_order"]) <= 1.6 and float(b["product_order"]) <= 1.2
    ]
    if (
        len(aromatic_reductions) >= 2
        and not facts.get("formed_bonds")
        and not facts.get("broken_bonds")
        and not facts.get("attachment_deltas")
    ):
        support = sorted(
            {int(m) for b in aromatic_reductions for m in (b["a"], b["b"])}
        )
        lost_aromatic = [
            m
            for m in support
            if _atom(facts, "reactant", m).get("aromatic")
            and not _atom(facts, "product", m).get("aromatic")
        ]
        h_gain = sum(max(0, _hydrogen_delta(facts, m)) for m in support)
        hetero = any(
            _atom(facts, "reactant", m).get("element") in {"N", "O", "S"}
            for m in support
        )
        if len(lost_aromatic) >= 2 and h_gain >= 2 and hetero:
            out.append(
                {
                    "event": "heteroarene_partial_hydrogenation",
                    "role": "specific_event",
                    "subtype": "partial_dearomatizing_hydrogenation",
                    "support_maps": support,
                    "mapping_scope": "multiple_mapped_aromatic_bonds_to_single_plus_h_gain_on_conserved_heteroarene",
                    "redox_direction": "reduction",
                }
            )

    return out


def detect_isothiocyanate_formation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Amine + thiocarbonyl dihalide-like carbon -> mapped N=C=S product."""
    out: list[dict[str, Any]] = []
    for b in facts.get("formed_bonds", []):
        a, c = int(b["a"]), int(b["b"])
        if float(b.get("order", 0.0)) < 1.8:
            continue
        ea = _atom(facts, "reactant", a).get("element")
        ec = _atom(facts, "reactant", c).get("element")
        if {ea, ec} != {"C", "N"}:
            continue
        carbon, nitrogen = (a, c) if ea == "C" else (c, a)
        r_atts = _attachments(facts, "reactant", carbon)
        p_atts = _attachments(facts, "product", carbon)
        has_cs_r = any(
            x.get("element") == "S" and float(x.get("bond_order", 0.0)) >= 1.8
            for x in r_atts
        )
        has_cs_p = any(
            x.get("element") == "S" and float(x.get("bond_order", 0.0)) >= 1.8
            for x in p_atts
        )
        halides = [x for x in r_atts if x.get("element") in HALOGENS]
        has_cn_p = any(
            _att_map(x) == nitrogen and float(x.get("bond_order", 0.0)) >= 1.8
            for x in p_atts
        )
        if (
            has_cs_r
            and has_cs_p
            and len(halides) >= 1
            and has_cn_p
            and _hydrogen_delta(facts, nitrogen) <= -1
        ):
            out.append(
                {
                    "event": "isothiocyanate_formation",
                    "role": "specific_event",
                    "subtype": "amine_to_isothiocyanate",
                    "center_map": carbon,
                    "nitrogen_map": nitrogen,
                    "support_maps": _support(carbon, nitrogen),
                    "mapping_scope": "mapped_thiocarbonyl_halide_carbon_plus_amine_n_to_n_c_s",
                }
            )
    return out


# =============================================================================
# S2_3A_HETEROCUMULENE_CLOSURE
#
# Product-state / local-geometry closure for:
#
#   amine -> isocyanate
#   amine -> isothiocyanate
#   amine + isocyanate -> urea
#   amine + cyanate-equivalent N#C-O -> urea
#   amine + thiocyanate-equivalent N#C-S -> thiourea
#
# This detector is additive.  Existing heterocumulene detectors remain intact.
# =============================================================================


def _s23a_is_cyanate_equivalent_carbon(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:

    if (
        _atom(
            facts,
            side,
            center,
        ).get("element")
        != "C"
    ):
        return False

    atts = (
        _attachments(
            facts,
            side,
            center,
        )
        or []
    )

    n_triple = [
        x
        for x in atts
        if (
            x.get("element") == "N"
            and float(
                x.get(
                    "bond_order",
                    0.0,
                )
            )
            >= 2.8
        )
    ]

    o_single = [
        x
        for x in atts
        if (
            x.get("element") == "O"
            and 0.8
            <= float(
                x.get(
                    "bond_order",
                    0.0,
                )
            )
            <= 1.2
        )
    ]

    # Conservative N#C-O local state.
    return len(n_triple) == 1 and len(o_single) == 1 and len(atts) == 2


def _s23a_is_thiocyanate_equivalent_carbon(
    facts: dict[str, Any],
    side: str,
    center: int,
) -> bool:

    if (
        _atom(
            facts,
            side,
            center,
        ).get("element")
        != "C"
    ):
        return False

    atts = (
        _attachments(
            facts,
            side,
            center,
        )
        or []
    )

    n_triple = [
        x
        for x in atts
        if (
            x.get("element") == "N"
            and float(
                x.get(
                    "bond_order",
                    0.0,
                )
            )
            >= 2.8
        )
    ]

    s_single = [
        x
        for x in atts
        if (
            x.get("element") == "S"
            and 0.8
            <= float(
                x.get(
                    "bond_order",
                    0.0,
                )
            )
            <= 1.2
        )
    ]

    # Conservative N#C-S local state.
    return len(n_triple) == 1 and len(s_single) == 1 and len(atts) == 2


def _s23a_product_urea_like_center(
    facts: dict[str, Any],
    carbon: int,
    incoming_n: int,
    chalcogen: str,
) -> bool:

    p_atts = (
        _attachments(
            facts,
            "product",
            carbon,
        )
        or []
    )

    has_chalcogen_double = any(
        x.get("element") == chalcogen
        and float(
            x.get(
                "bond_order",
                0.0,
            )
        )
        >= 1.8
        for x in p_atts
    )

    incoming_attached = any(
        x.get("element") == "N"
        and _att_map(x) is not None
        and int(_att_map(x)) == int(incoming_n)
        and float(
            x.get(
                "bond_order",
                0.0,
            )
        )
        <= 1.2
        for x in p_atts
    )

    n_single = [
        x
        for x in p_atts
        if (
            x.get("element") == "N"
            and float(
                x.get(
                    "bond_order",
                    0.0,
                )
            )
            <= 1.2
        )
    ]

    return has_chalcogen_double and incoming_attached and len(n_single) >= 2


def _s23a_local_support(
    facts: dict[str, Any],
    carbon: int,
    incoming: int,
) -> list[int]:

    maps = {
        int(carbon),
        int(incoming),
    }

    for side in (
        "reactant",
        "product",
    ):
        for x in (
            _attachments(
                facts,
                side,
                carbon,
            )
            or []
        ):
            m = _att_map(x)

            if m is not None:
                maps.add(int(m))

    return sorted(maps)


def detect_s23a_heterocumulene_closure(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []

    # -------------------------------------------------------------
    # Existing narrow isothiocyanate detector owns anything it
    # already recognizes.  The S2.3A detector only broadens coverage.
    # -------------------------------------------------------------

    existing_itc_centers = {
        int(r["center_map"])
        for r in (detect_isothiocyanate_formation(facts) or [])
        if r.get("center_map") is not None
    }

    # =============================================================
    # A. Heterocumulene FORMATION
    #
    # Require an intercomponent C=N double bond and an exact
    # N=C=O / N=C=S product state.
    # =============================================================

    for b in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        if (
            float(
                b.get(
                    "order",
                    0.0,
                )
            )
            < 1.8
        ):
            continue

        a = int(b["a"])

        c = int(b["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        ec = _atom(
            facts,
            "reactant",
            c,
        ).get("element")

        if {
            ea,
            ec,
        } != {
            "C",
            "N",
        }:
            continue

        carbon, nitrogen = (a, c) if ea == "C" else (c, a)

        if (
            _hydrogen_delta(
                facts,
                nitrogen,
            )
            > -1
        ):
            continue

        support = _s23a_local_support(
            facts,
            carbon,
            nitrogen,
        )

        # ---------------------------------------------------------
        # Isocyanate formation
        # ---------------------------------------------------------

        if _is_isocyanate_carbon(
            facts,
            "product",
            carbon,
        ) and not _is_isocyanate_carbon(
            facts,
            "reactant",
            carbon,
        ):

            out.append(
                {
                    "event": "isocyanate_formation",
                    "role": "specific_event",
                    "subtype": "amine_to_isocyanate",
                    "center_map": carbon,
                    "nitrogen_map": nitrogen,
                    "support_maps": support,
                    "mapping_scope": (
                        "mapped_intercomponent_c_n_double_formation_"
                        "to_exact_n_c_o_heterocumulene_product"
                    ),
                }
            )

            continue

        # ---------------------------------------------------------
        # Broaden isothiocyanate formation beyond the historical
        # thiocarbonyl-halide precursor restriction.
        # ---------------------------------------------------------

        if (
            _is_isothiocyanate_carbon(
                facts,
                "product",
                carbon,
            )
            and not _is_isothiocyanate_carbon(
                facts,
                "reactant",
                carbon,
            )
            and carbon not in existing_itc_centers
        ):

            out.append(
                {
                    "event": "isothiocyanate_formation",
                    "role": "specific_event",
                    "subtype": "amine_to_isothiocyanate",
                    "center_map": carbon,
                    "nitrogen_map": nitrogen,
                    "support_maps": support,
                    "mapping_scope": (
                        "mapped_intercomponent_c_n_double_formation_"
                        "to_exact_n_c_s_heterocumulene_product"
                    ),
                }
            )

    # =============================================================
    # B. Heterocumulene / heterocumulene-equivalent ADDITION
    #
    # Incoming N forms a single bond to the mapped heterocumulene
    # carbon and the product is an exact urea/thiourea-like center.
    # =============================================================

    for b in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        bo = float(
            b.get(
                "order",
                0.0,
            )
        )

        if not (0.8 <= bo <= 1.2):
            continue

        a = int(b["a"])

        c = int(b["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        ec = _atom(
            facts,
            "reactant",
            c,
        ).get("element")

        if {
            ea,
            ec,
        } != {
            "C",
            "N",
        }:
            continue

        carbon, nitrogen = (a, c) if ea == "C" else (c, a)

        if (
            _hydrogen_delta(
                facts,
                nitrogen,
            )
            > -1
        ):
            continue

        support = _s23a_local_support(
            facts,
            carbon,
            nitrogen,
        )

        # ---------------------------------------------------------
        # Actual isocyanate + amine -> urea
        # ---------------------------------------------------------

        if _is_isocyanate_carbon(
            facts,
            "reactant",
            carbon,
        ) and _s23a_product_urea_like_center(
            facts,
            carbon,
            nitrogen,
            "O",
        ):

            out.append(
                {
                    "event": "heterocumulene_addition",
                    "role": "specific_event",
                    "subtype": "amine_to_isocyanate_urea_formation",
                    "center_map": carbon,
                    "nucleophile_map": nitrogen,
                    "support_maps": support,
                    "mapping_scope": (
                        "mapped_n_addition_to_isocyanate_carbon_"
                        "with_exact_urea_product_state"
                    ),
                }
            )

            continue

        # ---------------------------------------------------------
        # N#C-O mapped equivalent -> urea
        # ---------------------------------------------------------

        if _s23a_is_cyanate_equivalent_carbon(
            facts,
            "reactant",
            carbon,
        ) and _s23a_product_urea_like_center(
            facts,
            carbon,
            nitrogen,
            "O",
        ):

            out.append(
                {
                    "event": "heterocumulene_addition",
                    "role": "specific_event",
                    "subtype": "amine_to_cyanate_equivalent_urea_formation",
                    "center_map": carbon,
                    "nucleophile_map": nitrogen,
                    "support_maps": support,
                    "mapping_scope": (
                        "mapped_n_addition_to_n_c_o_cyanate_"
                        "equivalent_with_exact_urea_product_state"
                    ),
                }
            )

            continue

        # ---------------------------------------------------------
        # N#C-S mapped equivalent -> thiourea
        # ---------------------------------------------------------

        if _s23a_is_thiocyanate_equivalent_carbon(
            facts,
            "reactant",
            carbon,
        ) and _s23a_product_urea_like_center(
            facts,
            carbon,
            nitrogen,
            "S",
        ):

            out.append(
                {
                    "event": "heterocumulene_addition",
                    "role": "specific_event",
                    "subtype": "amine_to_thiocyanate_equivalent_thiourea_formation",
                    "center_map": carbon,
                    "nucleophile_map": nitrogen,
                    "support_maps": support,
                    "mapping_scope": (
                        "mapped_n_addition_to_n_c_s_thiocyanate_"
                        "equivalent_with_exact_thiourea_product_state"
                    ),
                }
            )

    return out


# =============================================================================
# Batch07: holdout-driven semantic closure and event-specialization detectors.
# =============================================================================


def _bond_change(facts: dict[str, Any], a: int, b: int) -> tuple[float, float] | None:
    pair = {int(a), int(b)}
    for rec in facts.get("bond_order_changes", []) or []:
        if {int(rec["a"]), int(rec["b"])} == pair:
            return float(rec.get("reactant_order", 0.0)), float(
                rec.get("product_order", 0.0)
            )
    return None


def _attachment_delta(facts: dict[str, Any], center: int) -> dict[str, Any] | None:
    for rec in facts.get("attachment_deltas", []) or []:
        if int(rec["center_map"]) == int(center):
            return rec
    return None


def _neighbor_is_sulfonyl_s(facts: dict[str, Any], side: str, mn: int) -> bool:
    if _atom(facts, side, mn).get("element") != "S":
        return False
    return (
        sum(
            1
            for x in _attachments(facts, side, mn)
            if x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
        )
        >= 2
    )


def _alkene_activated(facts: dict[str, Any], a: int, b: int) -> bool:
    """Conservative Michael-acceptor context on either alkene carbon."""
    for c, other in ((a, b), (b, a)):
        for x in _attachments(facts, "reactant", c):
            m = _att_map(x)
            if m is None or int(m) == int(other):
                continue
            if x.get("element") == "C" and bool(x.get("carbonyl_like")):
                return True
            if x.get("element") == "S" and _neighbor_is_sulfonyl_s(
                facts, "reactant", int(m)
            ):
                return True
            if (
                x.get("element") in {"N", "P"}
                and float(x.get("bond_order", 0.0)) >= 1.8
            ):
                return True
    return False


# =============================================================================
# S2_3D_ALPHA_HYDROXYIMINO_AMIDE
#
# Composite mapped transformation:
#
#   trihalomethyl-C -- C(OH)2
#          + amine-N
#          + N-O nucleophile
#
#                 ->
#
#       N-C(=O)-C=N-O
#
# The mapped O_C_to_C relocation commonly observed in this family is mapping /
# topology evidence only.  Oxygen provenance is deliberately NOT part of the
# semantic signature.
# =============================================================================


def detect_s23d_alpha_hydroxyimino_amide(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []

    def bond_order(
        side: str,
        a: int,
        b: int,
    ) -> float:

        for att in (
            _attachments(
                facts,
                side,
                a,
            )
            or []
        ):
            mn = _att_map(att)

            if mn is not None and int(mn) == int(b):
                return float(
                    att.get(
                        "bond_order",
                        0.0,
                    )
                )

        return 0.0

    def element_atts(
        side: str,
        center: int,
        element: str,
        low: float,
        high: float,
    ) -> list[dict[str, Any]]:

        return [
            x
            for x in (
                _attachments(
                    facts,
                    side,
                    center,
                )
                or []
            )
            if (
                x.get("element") == element
                and low
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= high
            )
        ]

    formed_cn_single = []
    formed_cn_double = []

    # ---------------------------------------------------------
    # Identify the two intercomponent C-N formations.
    # ---------------------------------------------------------

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        a = int(fb["a"])

        b = int(fb["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element")

        if {ea, eb} != {"C", "N"}:
            continue

        carbon, nitrogen = (a, b) if ea == "C" else (b, a)

        bo = float(
            fb.get(
                "order",
                0.0,
            )
        )

        if 0.8 <= bo <= 1.2:
            formed_cn_single.append(
                (
                    carbon,
                    nitrogen,
                )
            )

        elif 1.8 <= bo <= 2.2:
            formed_cn_double.append(
                (
                    carbon,
                    nitrogen,
                )
            )

    # ---------------------------------------------------------
    # Composite ownership requires both centers.
    # ---------------------------------------------------------

    for (
        amide_c,
        amide_n,
    ) in formed_cn_single:

        for (
            oxime_c,
            oxime_n,
        ) in formed_cn_double:

            if amide_c == oxime_c:
                continue

            # Retained adjacent C-C framework.
            if not (
                0.8
                <= bond_order(
                    "reactant",
                    amide_c,
                    oxime_c,
                )
                <= 1.2
            ):
                continue

            if not (
                0.8
                <= bond_order(
                    "product",
                    amide_c,
                    oxime_c,
                )
                <= 1.2
            ):
                continue

            # -------------------------------------------------
            # Reactant amide precursor:
            # exactly three singly attached halogens.
            # -------------------------------------------------

            r_hal = [
                x
                for x in (
                    _attachments(
                        facts,
                        "reactant",
                        amide_c,
                    )
                    or []
                )
                if (
                    x.get("element")
                    in {
                        "F",
                        "Cl",
                        "Br",
                        "I",
                    }
                    and 0.8
                    <= float(
                        x.get(
                            "bond_order",
                            0.0,
                        )
                    )
                    <= 1.2
                )
            ]

            if len(r_hal) != 3:
                continue

            # -------------------------------------------------
            # Adjacent carbon is hydrate/geminal O,O-like.
            # -------------------------------------------------

            r_hydrate_o = element_atts(
                "reactant",
                oxime_c,
                "O",
                0.8,
                1.2,
            )

            if len(r_hydrate_o) != 2:
                continue

            # -------------------------------------------------
            # Product amide center: C(=O)-N
            # -------------------------------------------------

            p_amide_o = element_atts(
                "product",
                amide_c,
                "O",
                1.8,
                2.2,
            )

            if len(p_amide_o) != 1:
                continue

            p_amide_n = [
                x
                for x in element_atts(
                    "product",
                    amide_c,
                    "N",
                    0.8,
                    1.2,
                )
                if (_att_map(x) is not None and int(_att_map(x)) == amide_n)
            ]

            if len(p_amide_n) != 1:
                continue

            # -------------------------------------------------
            # Product hydroxyimino center: C=N-O
            # -------------------------------------------------

            p_oxime_n = [
                x
                for x in element_atts(
                    "product",
                    oxime_c,
                    "N",
                    1.8,
                    2.2,
                )
                if (_att_map(x) is not None and int(_att_map(x)) == oxime_n)
            ]

            if len(p_oxime_n) != 1:
                continue

            p_n_o = element_atts(
                "product",
                oxime_n,
                "O",
                0.8,
                1.2,
            )

            if not p_n_o:
                continue

            # Both incoming nitrogens lose hydrogen.
            if (
                _hydrogen_delta(
                    facts,
                    amide_n,
                )
                > -1
            ):
                continue

            if (
                _hydrogen_delta(
                    facts,
                    oxime_n,
                )
                > -1
            ):
                continue

            carbonyl_o_map = _att_map(p_amide_o[0])

            oxime_o_map = _att_map(p_n_o[0])

            support = _support(
                amide_c,
                amide_n,
                oxime_c,
                oxime_n,
                carbonyl_o_map,
                oxime_o_map,
            )

            out.append(
                {
                    "event": "alpha_hydroxyimino_amide_formation",
                    "role": "specific_event",
                    "subtype": ("trihalomethyl_hydrate_to_" "alpha_hydroxyimino_amide"),
                    "center_map": amide_c,
                    "amide_carbon_map": amide_c,
                    "amide_nitrogen_map": amide_n,
                    "oxime_carbon_map": oxime_c,
                    "oxime_nitrogen_map": oxime_n,
                    "carbonyl_oxygen_map": carbonyl_o_map,
                    "oxime_oxygen_map": oxime_o_map,
                    "support_maps": support,
                    "mapping_scope": (
                        "retained_adjacent_carbon_pair_with_"
                        "trihalomethyl_to_amide_carbonyl_and_"
                        "hydrate_carbon_to_c_n_o_hydroxyimino_"
                        "construction"
                    ),
                }
            )

    return out


def _p2_alpha_hydroxyimino_nitrile_support(
    facts: dict[str, Any], alpha_carbon: int, nitroso_nitrogen: int
) -> tuple[int, int, int] | None:
    """Qualify C(H)-C#N -> C(=N-OH)-C#N mapped product-state chemistry.

    Requires the same mapped nitrile C#N and adjacent C-C bond on both sides.
    The incoming N=O must become C=N-OH. No named nitrosation mechanism
    or untracked atom provenance is asserted.
    """
    if any(
        _atom(facts, side, alpha_carbon).get("element") != "C"
        or _atom(facts, side, nitroso_nitrogen).get("element") != "N"
        for side in ("reactant", "product")
    ):
        return None
    if _hydrogen_delta(facts, alpha_carbon) >= 0:
        return None
    product_n = _attachments(facts, "product", nitroso_nitrogen)
    if any(a.get("element") == "N" for a in product_n):
        return None
    oxygens = [
        int(_att_map(a)) for a in product_n
        if a.get("element") == "O"
        and _att_map(a) is not None
        and 0.8 <= float(a.get("bond_order", 0)) <= 1.2
        and _atom(facts, "product", int(_att_map(a))).get("total_h", 0) >= 1
        and 1.8 <= _bond_between(facts, "reactant", nitroso_nitrogen, int(_att_map(a)))
        and 0.8 <= _bond_between(facts, "product", nitroso_nitrogen, int(_att_map(a))) <= 1.2
    ]
    if len(oxygens) != 1 or sum(a.get("element") == "O" for a in product_n) != 1:
        return None
    for a in _attachments(facts, "reactant", alpha_carbon):
        nitrile_c = _att_map(a)
        if (
            a.get("element") != "C" or nitrile_c is None
            or not 0.8 <= float(a.get("bond_order", 0)) <= 1.2
            or not 0.8 <= _bond_between(facts, "product", alpha_carbon, nitrile_c) <= 1.2
        ):
            continue
        for n in _attachments(facts, "reactant", nitrile_c):
            nitrile_n = _att_map(n)
            if (
                n.get("element") == "N" and nitrile_n is not None
                and float(n.get("bond_order", 0)) >= 2.8
                and _atom(facts, "product", nitrile_c).get("element") == "C"
                and _atom(facts, "product", nitrile_n).get("element") == "N"
                and _bond_between(facts, "product", nitrile_c, nitrile_n) >= 2.8
            ):
                return nitrile_c, nitrile_n, oxygens[0]
    return None


def detect_alpha_oximation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Alpha C-H to C=N-OH beside a retained mapped carbonyl or nitrile."""
    out = []
    for b in facts.get("formed_bonds", []) or []:
        if float(b.get("order", 0.0)) < 1.8:
            continue
        a, c = int(b["a"]), int(b["b"])
        ea, ec = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", c
        ).get("element")
        if {ea, ec} != {"C", "N"}:
            continue
        carbon, nitrogen = (a, c) if ea == "C" else (c, a)
        if not _is_alpha_to_carbonyl(facts, "reactant", carbon):
            # A separate, narrower substrate class: alpha C-H next to a
            # *retained mapped nitrile*, rather than a retained carbonyl.
            nitrile = _p2_alpha_hydroxyimino_nitrile_support(facts, carbon, nitrogen)
            if nitrile is None:
                continue
            nitrile_c, nitrile_n, oxime_o = nitrile
            out.append({
                "event": "alpha_oximation",
                "role": "specific_event",
                "subtype": "alpha_c_h_to_hydroxyimino_alpha_nitrile",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "oxygen_map": oxime_o,
                "adjacent_nitrile_map": nitrile_c,
                "nitrile_nitrogen_map": nitrile_n,
                "support_maps": _support(carbon, nitrogen, oxime_o, nitrile_c, nitrile_n),
                "mapping_scope": "mapped_alpha_c_h_to_c_n_oh_with_retained_adjacent_c_n_triple_bond",
            })
            continue
        # Nitrogen must remain bonded to oxygen in the product; N-N composites are excluded.
        p_n_atts = _attachments(facts, "product", nitrogen)
        if not any(x.get("element") == "O" for x in p_n_atts):
            continue
        if any(x.get("element") == "N" for x in p_n_atts):
            continue
        if _hydrogen_delta(facts, carbon) > -1:
            continue
        # Adjacent carbonyl must remain a carbonyl.
        retained = []
        for x in _attachments(facts, "reactant", carbon):
            cm = _att_map(x)
            if (
                cm is not None
                and x.get("element") == "C"
                and bool(x.get("carbonyl_like"))
                and _has_conserved_carbonyl(facts, int(cm))
            ):
                retained.append(int(cm))
        if not retained:
            continue
        out.append(
            {
                "event": "alpha_oximation",
                "role": "specific_event",
                "subtype": "alpha_c_h_to_oxime",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "adjacent_carbonyl_map": retained[0],
                "support_maps": _support(carbon, nitrogen, retained[0]),
                "mapping_scope": "mapped_alpha_carbon_gains_c_n_double_bond_to_n_o_h_with_retained_adjacent_carbonyl",
            }
        )
    return out


def detect_alkene_sulfenyl_halogenation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Vicinal sulfenyl-halogen addition across a mapped alkene."""
    out = []
    for oc in facts.get("bond_order_changes", []) or []:
        if not (
            float(oc.get("reactant_order", 0.0)) >= 1.8
            and float(oc.get("product_order", 0.0)) <= 1.2
        ):
            continue
        a, b = int(oc["a"]), int(oc["b"])
        if (
            _atom(facts, "reactant", a).get("element") != "C"
            or _atom(facts, "reactant", b).get("element") != "C"
        ):
            continue
        s_links, x_links = [], []
        for fb in facts.get("formed_bonds", []) or []:
            x, y = int(fb["a"]), int(fb["b"])
            if float(fb.get("order", 0.0)) > 1.2:
                continue
            if a in (x, y) or b in (x, y):
                alk = a if a in (x, y) else b
                other = y if alk == x else x
                el = _atom(facts, "product", other).get("element")
                if el == "S":
                    s_links.append((alk, other))
                if el in HALOGENS:
                    x_links.append((alk, other, el))
        if not s_links or not x_links or s_links[0][0] == x_links[0][0]:
            continue
        sm = s_links[0][1]
        xm, xel = x_links[0][1], x_links[0][2]
        sx_broken = any(
            {int(bb["a"]), int(bb["b"])} == {sm, xm}
            for bb in facts.get("broken_bonds", []) or []
        )
        if not sx_broken:
            continue
        out.append(
            {
                "event": "alkene_sulfenyl_halogenation",
                "role": "specific_event",
                "subtype": (
                    "chlorosulfenylation"
                    if xel == "Cl"
                    else f"{xel.lower()}_sulfenylation"
                ),
                "center_map": a,
                "alkene_maps": sorted((a, b)),
                "sulfur_map": sm,
                "halogen_map": xm,
                "support_maps": _support(a, b, sm, xm),
                "mapping_scope": "mapped_alkene_consumption_plus_vicinal_c_s_and_c_x_formation_and_s_x_cleavage",
            }
        )
    return out


# =============================================================================
# S2_3C_UNMAPPED_LG_N_CYCLOPROPYLATION
#
# Representation-boundary extension of n_alkylation.
#
# Structural contract:
#
#   cyclopropyl C bearing two unmapped O substituents
#                     +
#                  mapped N-H
#
#                     ->
#
#             mapped cyclopropyl C-N
#
# Both O substituents disappear, C gains H, and N loses H.
#
# This is graph-defined. Conditions/reagents are not part of the semantic
# detector contract.
# =============================================================================


def detect_unmapped_lg_n_cyclopropylation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []

    def bond_order(
        side: str,
        a: int,
        b: int,
    ) -> float:

        for att in (
            _attachments(
                facts,
                side,
                a,
            )
            or []
        ):

            mn = _att_map(att)

            if mn is not None and int(mn) == int(b):
                return float(
                    att.get(
                        "bond_order",
                        0.0,
                    )
                )

        return 0.0

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        # New C-N single bond.
        bo = float(
            fb.get(
                "order",
                0.0,
            )
        )

        if not (0.8 <= bo <= 1.2):
            continue

        a = int(fb["a"])

        b = int(fb["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element")

        if {
            ea,
            eb,
        } != {
            "C",
            "N",
        }:
            continue

        carbon, nitrogen = (a, b) if ea == "C" else (b, a)

        # ---------------------------------------------------------
        # Donor center must be aliphatic sp3 carbon.
        # ---------------------------------------------------------

        ca = _atom(
            facts,
            "reactant",
            carbon,
        )

        if ca.get("aromatic"):
            continue

        if not _is_sp3_carbon(
            facts,
            "reactant",
            carbon,
        ):
            continue

        if _has_conserved_carbonyl(
            facts,
            carbon,
        ):
            continue

        # ---------------------------------------------------------
        # This is not ring opening or pi-addition chemistry.
        # ---------------------------------------------------------

        if (
            int(
                facts.get(
                    "common_cycle_rank_delta",
                    0,
                )
            )
            != 0
        ):
            continue

        if _consumed_pi_bond_at(
            facts,
            carbon,
        ):
            continue

        # ---------------------------------------------------------
        # Incoming N-H loses hydrogen.
        # ---------------------------------------------------------

        if (
            _hydrogen_delta(
                facts,
                nitrogen,
            )
            > -1
        ):
            continue

        # ---------------------------------------------------------
        # Exact donor-center attachment change:
        #
        #     lost:   O(unmapped), O(unmapped)
        #     gained: N(mapped)
        #     C Hdelta >= +1
        #
        # Ordinary detect_substitution() cannot represent this
        # because it requires one lost + one gained attachment.
        # ---------------------------------------------------------

        d = _attachment_delta(
            facts,
            carbon,
        )

        if not d:
            continue

        lost = (
            d.get(
                "lost",
                [],
            )
            or []
        )

        gained = (
            d.get(
                "gained",
                [],
            )
            or []
        )

        if len(lost) != 2:
            continue

        if len(gained) != 1:
            continue

        lost_o = [
            x
            for x in lost
            if (
                x.get("element") == "O"
                and x.get("map") is None
                and 0.8
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= 1.2
            )
        ]

        if len(lost_o) != 2:
            continue

        gained_n = gained[0]

        if gained_n.get("element") != "N":
            continue

        gained_n_map = _att_map(gained_n)

        if gained_n_map is None or int(gained_n_map) != nitrogen:
            continue

        if not (
            0.8
            <= float(
                gained_n.get(
                    "bond_order",
                    0.0,
                )
            )
            <= 1.2
        ):
            continue

        if (
            int(
                d.get(
                    "hydrogen_delta",
                    0,
                )
            )
            < 1
        ):
            continue

        # ---------------------------------------------------------
        # Cyclopropyl geometry:
        #
        #       C1
        #      /  \
        # donor--C2
        #
        # Donor carbon must have exactly two mapped carbon
        # neighbours and those neighbours must themselves be bonded.
        # ---------------------------------------------------------

        r_atts = (
            _attachments(
                facts,
                "reactant",
                carbon,
            )
            or []
        )

        ring_c = [
            int(_att_map(x))
            for x in r_atts
            if (
                x.get("element") == "C"
                and _att_map(x) is not None
                and 0.8
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= 1.2
            )
        ]

        if len(ring_c) != 2:
            continue

        c1, c2 = ring_c

        if not (
            0.8
            <= bond_order(
                "reactant",
                c1,
                c2,
            )
            <= 1.2
        ):
            continue

        # ---------------------------------------------------------
        # Product must preserve the cyclopropane and contain the
        # new C-N bond, with no oxygen remaining on donor carbon.
        # ---------------------------------------------------------

        p_atts = (
            _attachments(
                facts,
                "product",
                carbon,
            )
            or []
        )

        if any(x.get("element") == "O" for x in p_atts):
            continue

        product_n = [
            x
            for x in p_atts
            if (
                x.get("element") == "N"
                and _att_map(x) is not None
                and int(_att_map(x)) == nitrogen
                and 0.8
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= 1.2
            )
        ]

        if len(product_n) != 1:
            continue

        retained_ring_c = {
            int(_att_map(x))
            for x in p_atts
            if (
                x.get("element") == "C"
                and _att_map(x) is not None
                and 0.8
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= 1.2
            )
        }

        if (
            not {
                c1,
                c2,
            }
            <= retained_ring_c
        ):
            continue

        if not (
            0.8
            <= bond_order(
                "product",
                c1,
                c2,
            )
            <= 1.2
        ):
            continue

        out.append(
            {
                "event": "n_alkylation",
                "role": "specific_event",
                "subtype": ("n_cyclopropylation_from_" "geminal_o_o_donor"),
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "cyclopropyl_maps": sorted(
                    [
                        carbon,
                        c1,
                        c2,
                    ]
                ),
                "support_maps": _support(
                    carbon,
                    nitrogen,
                    c1,
                    c2,
                ),
                "mapping_scope": (
                    "mapped_intercomponent_c_n_single_"
                    "formation_at_cyclopropyl_carbon_"
                    "with_loss_of_two_unmapped_o_"
                    "substituents_and_c_h_gain"
                ),
            }
        )

    return out


def detect_unmapped_lg_o_alkylation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """O-alkylation when the leaving group is absent/unmapped but net substitution is clear."""
    out = []
    for b in facts.get("intercomponent_formed_bonds", []) or []:
        if float(b.get("order", 0.0)) > 1.2:
            continue
        a, c = int(b["a"]), int(b["b"])
        ea, ec = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", c
        ).get("element")
        if {ea, ec} != {"C", "O"}:
            continue
        carbon, oxygen = (a, c) if ea == "C" else (c, a)
        if not _is_sp3_carbon(facts, "reactant", carbon):
            continue
        if _hydrogen_delta(facts, oxygen) > -1:
            continue
        if _has_conserved_carbonyl(facts, carbon) or _atom(
            facts, "reactant", carbon
        ).get("aromatic"):
            continue
        # A mapped leaving-group loss is handled by the ordinary substitution detector.
        d = _attachment_delta(facts, carbon)
        if d and any(
            x.get("element") in (HALOGENS | {"O", "N", "S"})
            for x in d.get("lost", []) or []
        ):
            continue
        # Avoid epoxide/ring-opening and other pi-addition chemistry.
        if int(facts.get("common_cycle_rank_delta", 0)) < 0 or _consumed_pi_bond_at(
            facts, carbon
        ):
            continue
        out.append(
            {
                "event": "o_alkylation",
                "role": "specific_event",
                "subtype": "unmapped_leaving_group_sp3",
                "center_map": carbon,
                "oxygen_map": oxygen,
                "support_maps": _support(carbon, oxygen),
                "mapping_scope": "intercomponent_sp3_c_o_formation_plus_o_h_loss_without_mapped_leaving_group",
            }
        )
    return out


def detect_n_glycosylation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Anomeric/sugar-like C-O -> C-N substitution with ring-oxygen context."""
    out = []
    for d in facts.get("attachment_deltas", []) or []:
        center = int(d["center_map"])
        if _atom(facts, "reactant", center).get("element") != "C" or not _atom(
            facts, "reactant", center
        ).get("in_ring"):
            continue
        # Cyclic anhydride/lactone/imide/carbonate carbonyls are acyl
        # centers, never glycosylation anomeric centers.
        if any(
            a.get("element") == "O" and float(a.get("bond_order", 0.0)) >= 1.8
            for a in _attachments(facts, "reactant", center)
        ):
            continue
        lost_o = [x for x in d.get("lost", []) or [] if x.get("element") == "O"]
        gained_n = [
            x
            for x in d.get("gained", []) or []
            if x.get("element") == "N" and _att_map(x) is not None
        ]
        if not lost_o or not gained_n:
            continue
        # Anomeric carbon should retain another O in the ring/sugar framework.
        retained_o = [
            x
            for x in _attachments(facts, "reactant", center)
            if x.get("element") == "O" and x not in lost_o
        ]
        if not retained_o:
            continue
        nmap = int(_att_map(gained_n[0]))
        if not _external_gain(facts, center, gained_n[0]):
            continue
        out.append(
            {
                "event": "n_glycosylation",
                "role": "specific_event",
                "subtype": "anomeric_c_o_to_n",
                "center_map": center,
                "nitrogen_map": nmap,
                "support_maps": _support(center, nmap, _att_map(lost_o[0])),
                "mapping_scope": "mapped_sugar_ring_carbon_loses_o_attachment_and_gains_external_n_with_retained_ring_oxygen",
            }
        )
    return out


def detect_vinylic_amination(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """C-O -> C-N substitution at a carbon of a retained C=C framework."""
    out = []
    for d in facts.get("attachment_deltas", []) or []:
        center = int(d["center_map"])
        if _atom(facts, "reactant", center).get("element") != "C" or _atom(
            facts, "reactant", center
        ).get("aromatic"):
            continue
        lost_o = [x for x in d.get("lost", []) or [] if x.get("element") == "O"]
        gained_n = [x for x in d.get("gained", []) or [] if x.get("element") == "N"]
        if not lost_o or not gained_n:
            continue
        retained_alkene = False
        partner = None
        for x in _attachments(facts, "reactant", center):
            m = _att_map(x)
            if (
                m is None
                or x.get("element") != "C"
                or float(x.get("bond_order", 0.0)) < 1.8
            ):
                continue
            if _bond_between(facts, "product", center, int(m)) >= 1.8:
                retained_alkene = True
                partner = int(m)
                break
        if not retained_alkene:
            continue
        nmap = _att_map(gained_n[0])
        out.append(
            {
                "event": "vinylic_amination",
                "role": "specific_event",
                "subtype": "vinylic_c_o_to_n",
                "center_map": center,
                "nitrogen_map": nmap,
                "alkene_partner_map": partner,
                "support_maps": _support(center, nmap, partner, _att_map(lost_o[0])),
                "mapping_scope": "mapped_vinylic_carbon_retains_c_c_pi_bond_while_o_attachment_is_replaced_by_n",
            }
        )
    return out


def detect_terminal_alkyne_deprotection(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    """Carbon-cap removal that unmasks a terminal alkyne.

    S2_2D1_TERMINAL_ALKYNE_DEPROTECTION

    Validated structural contract:

        cap-C-C#C-R  ->  H-C#C-R

    Require:
      * conserved mapped carbon center;
      * H gain at that center;
      * loss of exactly one one-sided/unmapped carbon attachment;
      * no heavy-atom gain;
      * the same mapped carbon remains triple-bonded to the center
        on both sides.

    The retained C#C requirement prevents generic C-C fragmentation
    from being promoted to terminal-alkyne deprotection.
    """
    out = []

    for d in facts.get("attachment_deltas", []) or []:

        c = int(d["center_map"])

        if (
            _atom(
                facts,
                "reactant",
                c,
            ).get("element")
            != "C"
        ):
            continue

        if (
            _atom(
                facts,
                "product",
                c,
            ).get("element")
            != "C"
        ):
            continue

        if (
            int(
                d.get(
                    "hydrogen_delta",
                    0,
                )
            )
            <= 0
        ):
            continue

        lost = list(d.get("lost", []) or [])

        gained = list(d.get("gained", []) or [])

        gained_heavy = [x for x in gained if x.get("element") != "H"]

        if gained_heavy:
            continue

        lost_c = [x for x in lost if x.get("element") == "C"]

        if len(lost) != 1 or len(lost_c) != 1:
            continue

        cap = lost_c[0]

        # Protecting cap must be one-sided/unmapped.
        if _att_map(cap) is not None:
            continue

        # Cap leaves through an ordinary single C-C bond.
        if (
            float(
                cap.get(
                    "bond_order",
                    0.0,
                )
                or 0.0
            )
            > 1.2
        ):
            continue

        # Do not absorb acyl cleavage into this family.
        if bool(cap.get("carbonyl_like")):
            continue

        react_triple = {
            int(_att_map(x))
            for x in _attachments(
                facts,
                "reactant",
                c,
            )
            if (
                x.get("element") == "C"
                and _att_map(x) is not None
                and float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                >= 2.8
            )
        }

        prod_triple = {
            int(_att_map(x))
            for x in _attachments(
                facts,
                "product",
                c,
            )
            if (
                x.get("element") == "C"
                and _att_map(x) is not None
                and float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                >= 2.8
            )
        }

        retained_partners = react_triple & prod_triple

        if not retained_partners:
            continue

        out.append(
            {
                "event": "terminal_alkyne_deprotection",
                "role": "specific_event",
                "subtype": "carbon_cap_to_terminal_alkyne",
                "center_map": c,
                "alkyne_partner_maps": sorted(retained_partners),
                "support_maps": _support(
                    c,
                    *sorted(retained_partners),
                ),
                "mapping_scope": "mapped_alkynyl_carbon_retains_c_c_triple_bond_loses_one_sided_carbon_cap_and_gains_h",
            }
        )

    return out


def detect_s_deprotection(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Detect thiol unmasking by loss of a carbon protecting attachment.

    S2_2B2_S_DEPROTECTION

    Required net geometry:

        R-S-C(protecting group) -> R-SH

    The surviving sulfur is mapped.  The departing carbon may be
    one-sided/unmapped, but attachment metadata preserves whether it
    is carbonyl-like.

    This intentionally does not classify arbitrary S-C cleavage:
      * sulfur must survive on both sides;
      * sulfur must gain H;
      * no new heavy attachment may be gained;
      * exactly one carbon attachment must be lost.

    Subtypes:
      s_deacylation     lost carbon is carbonyl-like
      s_c_deprotection  other carbon-based sulfur protecting group
    """
    out = []

    for d in facts.get("attachment_deltas", []) or []:
        s = int(d["center_map"])

        if _atom(facts, "reactant", s).get("element") != "S":
            continue

        if _atom(facts, "product", s).get("element") != "S":
            continue

        lost = list(d.get("lost", []) or [])
        gained = list(d.get("gained", []) or [])

        lost_c = [x for x in lost if x.get("element") == "C"]

        gained_heavy = [x for x in gained if x.get("element") != "H"]

        if len(lost_c) != 1:
            continue

        if len(lost) != 1:
            continue

        if gained_heavy:
            continue

        if int(d.get("hydrogen_delta", 0)) <= 0:
            continue

        x = lost_c[0]

        subtype = (
            "s_deacylation" if bool(x.get("carbonyl_like")) else "s_c_deprotection"
        )

        out.append(
            {
                "event": "s_deprotection",
                "role": "specific_event",
                "subtype": subtype,
                "center_map": s,
                "support_maps": _support(
                    s,
                    _att_map(x),
                ),
                "mapping_scope": "mapped_s_loses_one_carbon_protecting_attachment_and_gains_h",
            }
        )

    return out


def detect_o_heteroatom_deprotection(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    """O-deprotection by one-sided phosphorus or sulfonyl loss.

    S2_2C_O_HETEROATOM_DEPROTECTION

    Validated geometries:

        R-O-P(...)   -> R-OH
        R-O-S(=O)2R -> R-OH

    The surviving oxygen is mapped; the departing heteroatom is often
    unmapped but its attachment descriptor retains local environment.
    """
    out = []

    for d in facts.get("attachment_deltas", []) or []:

        o = int(d["center_map"])

        if (
            _atom(
                facts,
                "reactant",
                o,
            ).get("element")
            != "O"
        ):
            continue

        if (
            _atom(
                facts,
                "product",
                o,
            ).get("element")
            != "O"
        ):
            continue

        lost = list(d.get("lost", []) or [])

        gained = list(d.get("gained", []) or [])

        gained_heavy = [x for x in gained if x.get("element") != "H"]

        if gained_heavy:
            continue

        if (
            int(
                d.get(
                    "hydrogen_delta",
                    0,
                )
            )
            <= 0
        ):
            continue

        # ----------------------------------------------------
        # O-dephosphorylation
        #
        # Require exactly one lost P attachment and a
        # phosphate-like local environment with at least
        # three oxygen neighbors around phosphorus.
        # ----------------------------------------------------
        lost_p = [x for x in lost if x.get("element") == "P"]

        if len(lost) == 1 and len(lost_p) == 1:
            x = lost_p[0]

            shell = list(x.get("neighbor_shell") or [])

            n_o = sum(1 for y in shell if str(y[0]) == "O")

            n_poxo = sum(1 for y in shell if (str(y[0]) == "O" and float(y[1]) >= 1.5))

            if n_o >= 3 and n_poxo >= 1:
                out.append(
                    {
                        "event": "o_dephosphorylation",
                        "role": "specific_event",
                        "subtype": "phosphate_ester_to_alcohol",
                        "center_map": o,
                        "support_maps": _support(
                            o,
                            _att_map(x),
                        ),
                        "mapping_scope": "mapped_o_loses_one_sided_phosphate_attachment_and_gains_h",
                    }
                )

                continue

        # ----------------------------------------------------
        # O-desulfonylation
        #
        # Require the lost sulfur itself to be sulfonyl-like.
        # This prevents arbitrary O-S cleavage from being
        # promoted to protecting-group chemistry.
        # ----------------------------------------------------
        lost_s = [x for x in lost if x.get("element") == "S"]

        if len(lost) == 1 and len(lost_s) == 1:
            x = lost_s[0]

            if bool(x.get("sulfonyl_like")):
                out.append(
                    {
                        "event": "o_desulfonylation",
                        "role": "specific_event",
                        "subtype": "o_sulfonyl_deprotection",
                        "center_map": o,
                        "support_maps": _support(
                            o,
                            _att_map(x),
                        ),
                        "mapping_scope": "mapped_o_loses_one_sided_sulfonyl_attachment_and_gains_h",
                    }
                )

    return out


def detect_n_deprotection(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Net N-protecting-group removal classified by local lost-carbon context.

    Mixed protecting groups on the same N are emitted as separate semantic events
    (e.g. carbamate deprotection + N-debenzylation), which is a deliberate Batch07
    test case for multiple SemanticEvents within one connected edit center.
    """
    out = []
    for d in facts.get("attachment_deltas", []) or []:
        n = int(d["center_map"])
        if _atom(facts, "reactant", n).get("element") != "N":
            continue
        lost_c = [x for x in d.get("lost", []) or [] if x.get("element") == "C"]
        lost_s = [x for x in d.get("lost", []) or [] if x.get("element") == "S"]
        gained_heavy = [
            x for x in d.get("gained", []) or [] if x.get("element") not in {"H"}
        ]

        # S2_2B1_N_S_DEPROTECTION
        #
        # N-protecting-group removal may cross the mapped/unmapped
        # boundary through sulfur rather than carbon:
        #
        #   R-N-S(=O)-R'     -> R-NH       sulfinyl deprotection
        #   R-N-S(=O)2-R'    -> R-NH       sulfonyl deprotection
        #
        # The sulfur atom is commonly one-sided/unmapped but its local
        # neighbor shell remains available in the attachment descriptor.
        if not (lost_c or lost_s) or gained_heavy:
            continue

        hdelta = int(d.get("hydrogen_delta", 0))
        aromatic_n = bool(_atom(facts, "reactant", n).get("aromatic"))
        degree_drop = int(_atom(facts, "product", n).get("heavy_degree", 0)) < int(
            _atom(facts, "reactant", n).get("heavy_degree", 0)
        )
        tautomeric_h_gain = any(
            _atom(facts, "reactant", m).get("element") == "N"
            and _atom(facts, "product", m).get("element") == "N"
            and _hydrogen_delta(facts, int(m)) > 0
            for m in facts["mapping"]["common_maps"]
            if int(m) != n
            and facts["reactant"]["map_to_component"].get(int(m))
            == facts["reactant"]["map_to_component"].get(n)
        )
        if hdelta <= 0 and not (aromatic_n and degree_drop and tautomeric_h_gain):
            continue

        # S2_2B1_N_S_DEPROTECTION
        #
        # Preserve the existing ontology:
        #
        #   sulfonyl loss  -> specific child n_desulfonylation
        #   sulfinyl loss -> generic parent n_deprotection,
        #                    subtype n_sulfinyl_deprotection
        #
        # Do not generalize arbitrary N-S cleavage to deprotection.
        # The lost sulfur must carry at least one oxo group.
        #
        # neighbor_shell entries follow the existing attachment
        # descriptor convention:
        #     (element, bond_order, aromatic, ...)
        #
        # sulfinyl:  one S=O
        # sulfonyl: two S=O, also normally sulfonyl_like=True

        sulfur_owned = False

        for x in lost_s:
            shell = list(x.get("neighbor_shell") or [])

            oxo_count = sum(1 for y in shell if str(y[0]) == "O" and float(y[1]) >= 1.5)

            if bool(x.get("sulfonyl_like")) or oxo_count >= 2:
                out.append(
                    {
                        "event": "n_desulfonylation",
                        "role": "specific_event",
                        "subtype": "n_sulfonyl_deprotection",
                        "center_map": n,
                        "support_maps": _support(
                            n,
                            _att_map(x),
                        ),
                        "mapping_scope": "mapped_n_loses_one_sided_sulfonyl_protecting_attachment",
                    }
                )
                sulfur_owned = True

            elif oxo_count == 1:
                out.append(
                    {
                        "event": "n_deprotection",
                        "role": "specific_event",
                        "subtype": "n_sulfinyl_deprotection",
                        "center_map": n,
                        "support_maps": _support(
                            n,
                            _att_map(x),
                        ),
                        "mapping_scope": "mapped_n_loses_one_sided_sulfinyl_protecting_attachment",
                    }
                )
                sulfur_owned = True

        # A pure N-S deprotection is fully owned above.  If the same N
        # independently loses carbon protecting groups as well, retain
        # processing below so multiple deprotections may coexist.
        if sulfur_owned and not lost_c:
            continue

        carbonyl = [x for x in lost_c if bool(x.get("carbonyl_like"))]
        noncarbonyl = [x for x in lost_c if not bool(x.get("carbonyl_like"))]

        # Two carbonyl attachments are characteristic of imide/phthalimide release.
        if len(carbonyl) >= 2:
            out.append(
                {
                    "event": "n_deprotection",
                    "role": "specific_event",
                    "subtype": "phthalimide_or_imide_deprotection",
                    "center_map": n,
                    "support_maps": _support(n, *[_att_map(x) for x in carbonyl]),
                    "mapping_scope": "mapped_n_loses_two_carbonyl_protecting_attachments",
                }
            )
            continue

        # Carbonyl-derived protecting group.
        for x in carbonyl:
            shell = list(x.get("neighbor_shell") or [])
            has_single_o = any(str(y[0]) == "O" and float(y[1]) <= 1.2 for y in shell)
            subtype = "carbamate_deprotection" if has_single_o else "n_deacylation"
            out.append(
                {
                    "event": "n_deprotection",
                    "role": "specific_event",
                    "subtype": subtype,
                    "center_map": n,
                    "support_maps": _support(n, _att_map(x)),
                    "mapping_scope": "mapped_n_loses_carbonyl_derived_protecting_attachment",
                }
            )

        # Non-carbonyl carbon protecting groups, especially benzyl.
        # Aziridine N-C ring scission is ring opening, not N-deprotection.
        for x in noncarbonyl:
            if _lost_aziridine_ring_n_c_bond(facts, n, _att_map(x)):
                continue
            benzylic = any(
                str(y[0]) == "C" and bool(y[2]) for y in (x.get("neighbor_shell") or [])
            )
            subtype = "n_debenzylation" if benzylic else "n_dealkylation"
            out.append(
                {
                    "event": "n_deprotection",
                    "role": "specific_event",
                    "subtype": subtype,
                    "center_map": n,
                    "support_maps": _support(n, _att_map(x)),
                    "mapping_scope": "mapped_n_loses_noncarbonyl_carbon_protecting_attachment",
                }
            )
    return out


def detect_n_n_bond_formation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for b in facts.get("intercomponent_formed_bonds", []) or []:
        if float(b.get("order", 0.0)) > 1.2:
            continue
        a, c = int(b["a"]), int(b["b"])
        if (
            _atom(facts, "reactant", a).get("element") != "N"
            or _atom(facts, "reactant", c).get("element") != "N"
        ):
            continue
        # Incoming hydrazine-like N retains another N neighbor; chloramine/amino-transfer is n-amination.
        hydrazine_like = any(
            any(
                x.get("element") == "N" and _att_map(x) not in {a, c}
                for x in _attachments(facts, "reactant", n)
            )
            for n in (a, c)
        )
        subtype = "hydrazination" if hydrazine_like else "n_amination"
        out.append(
            {
                "event": "n_n_bond_formation",
                "role": "specific_event",
                "subtype": subtype,
                "center_map": a,
                "nitrogen_maps": sorted((a, c)),
                "support_maps": _support(a, c),
                "mapping_scope": "mapped_intercomponent_n_n_single_bond_formation",
            }
        )
    return out


def _is_isothiocyanate_carbon(facts: dict[str, Any], side: str, center: int) -> bool:
    if _atom(facts, side, center).get("element") != "C":
        return False
    atts = _attachments(facts, side, center)
    return any(
        x.get("element") == "S" and float(x.get("bond_order", 0.0)) >= 1.8 for x in atts
    ) and any(
        x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 1.8 for x in atts
    )


def detect_expanded_heterocumulene_addition(
    facts: dict[str, Any]
) -> list[dict[str, Any]]:
    out = []
    for b in facts.get("intercomponent_formed_bonds", []) or []:
        a, c = int(b["a"]), int(b["b"])
        ea, ec = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", c
        ).get("element")
        carbon = a if ea == "C" else (c if ec == "C" else None)
        nuc = c if carbon == a else a if carbon == c else None
        if carbon is None or nuc is None or float(b.get("order", 0.0)) > 1.2:
            continue
        nel = _atom(facts, "reactant", nuc).get("element")
        if _is_isocyanate_carbon(facts, "reactant", carbon) and nel == "O":
            # Addition should reduce C=N while keeping C=O.
            if any(
                _bond_change(facts, carbon, int(_att_map(x)))
                and _bond_change(facts, carbon, int(_att_map(x)))[0] >= 1.8
                for x in _attachments(facts, "reactant", carbon)
                if x.get("element") == "N" and _att_map(x) is not None
            ):
                out.append(
                    {
                        "event": "heterocumulene_addition",
                        "role": "specific_event",
                        "subtype": "alcohol_to_isocyanate_carbamate_formation",
                        "center_map": carbon,
                        "nucleophile_map": nuc,
                        "support_maps": _support(carbon, nuc),
                        "mapping_scope": "mapped_o_addition_to_isocyanate_carbon_with_heterocumulene_pi_reduction",
                    }
                )
        elif _is_isothiocyanate_carbon(facts, "reactant", carbon) and nel == "N":
            out.append(
                {
                    "event": "heterocumulene_addition",
                    "role": "specific_event",
                    "subtype": "amine_to_isothiocyanate_thiourea_formation",
                    "center_map": carbon,
                    "nucleophile_map": nuc,
                    "support_maps": _support(carbon, nuc),
                    "mapping_scope": "mapped_n_addition_to_isothiocyanate_carbon",
                }
            )

    # Thiocarbonyl dihalide -> thiourea-like substitution by N nucleophiles.
    for c in facts["mapping"]["common_maps"]:
        if _atom(facts, "reactant", c).get("element") != "C":
            continue
        r_atts = _attachments(facts, "reactant", c)
        p_atts = _attachments(facts, "product", c)
        if not any(
            x.get("element") == "S" and float(x.get("bond_order", 0.0)) >= 1.8
            for x in r_atts
        ):
            continue
        if sum(1 for x in r_atts if x.get("element") in HALOGENS) < 1:
            continue
        p_ns = [
            x
            for x in p_atts
            if x.get("element") == "N" and float(x.get("bond_order", 0.0)) <= 1.2
        ]
        if len(p_ns) >= 2:
            out.append(
                {
                    "event": "thiourea_formation",
                    "role": "specific_event",
                    "subtype": "thiocarbonyl_dihalide_to_thiourea",
                    "center_map": int(c),
                    "support_maps": _support(int(c), *[_att_map(x) for x in p_ns]),
                    "mapping_scope": "mapped_thiocarbonyl_carbon_replaces_halides_with_two_n_nucleophiles",
                }
            )
    return out


def detect_carbonyl_deprotection(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for c in facts["mapping"]["common_maps"]:
        if _atom(facts, "reactant", c).get("element") != "C":
            continue
        r_o_single = _element_attachments(facts, "reactant", c, "O", 0.8, 1.2)
        p_o_double = _element_attachments(facts, "product", c, "O", 1.8, 2.2)
        if len(r_o_single) < 2 or not p_o_double:
            continue
        if _carbon_fg_state(facts, "reactant", c):
            continue
        lost, _gained = _lost_gained_for_center(facts, c)
        lost_o = [x for x in lost if x.get("element") == "O"]
        if not lost_o and not any(
            int(c) in (int(x["a"]), int(x["b"]))
            and _bond_class_from_pair(facts, int(x["a"]), int(x["b"])) == "C-O"
            for x in facts.get("bond_order_changes", []) or []
        ):
            continue
        target = (
            "aldehyde"
            if int(_atom(facts, "product", c).get("total_h", 0)) >= 1
            else "ketone"
        )
        out.append(
            {
                "event": "carbonyl_deprotection",
                "role": "specific_event",
                "subtype": f"acetal_or_ketal_to_{target}",
                "center_map": int(c),
                "support_maps": _conversion_support(facts, int(c)),
                "mapping_scope": "mapped_carbon_with_two_single_o_attachments_restores_c_o_double_bond",
            }
        )
    return out


def _bond_class_from_pair(facts: dict[str, Any], a: int, b: int) -> str:
    return canonical_bond_class(
        _atom(facts, "reactant", a).get("element")
        or _atom(facts, "product", a).get("element"),
        _atom(facts, "reactant", b).get("element")
        or _atom(facts, "product", b).get("element"),
    )


def detect_heteroaryl_alkoxy_to_oxo(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for c in facts["mapping"]["common_maps"]:
        if _atom(facts, "reactant", c).get("element") != "C" or not _atom(
            facts, "reactant", c
        ).get("aromatic"):
            continue
        if not any(
            x.get("element") == "O" and float(x.get("bond_order", 0.0)) <= 1.2
            for x in _attachments(facts, "reactant", c)
        ):
            continue
        if not any(
            x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
            for x in _attachments(facts, "product", c)
        ):
            continue
        # Require a ring nitrogen neighbor on either side to avoid generic phenol/quinone chemistry.
        n_neighbor = any(
            x.get("element") == "N" for x in _attachments(facts, "reactant", c)
        ) or any(x.get("element") == "N" for x in _attachments(facts, "product", c))
        if not n_neighbor:
            continue
        out.append(
            {
                "event": "heteroaryl_hydrolysis_to_oxo",
                "role": "specific_event",
                "subtype": "heteroaryl_alkoxy_to_oxo",
                "center_map": int(c),
                "support_maps": _conversion_support(facts, int(c)),
                "mapping_scope": "mapped_heteroaryl_c_o_single_to_c_o_double_with_ring_n_context",
            }
        )
    return out


def detect_aza_michael_addition(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for oc in facts.get("bond_order_changes", []) or []:
        if not (
            float(oc.get("reactant_order", 0.0)) >= 1.8
            and float(oc.get("product_order", 0.0)) <= 1.2
        ):
            continue
        a, b = int(oc["a"]), int(oc["b"])
        if (
            _atom(facts, "reactant", a).get("element") != "C"
            or _atom(facts, "reactant", b).get("element") != "C"
        ):
            continue
        if not _alkene_activated(facts, a, b):
            continue
        formed = []
        for fb in facts.get("intercomponent_formed_bonds", []) or []:
            x, y = int(fb["a"]), int(fb["b"])
            if float(fb.get("order", 0.0)) > 1.2:
                continue
            if a in (x, y) or b in (x, y):
                carbon = a if a in (x, y) else b
                other = y if carbon == x else x
                if _atom(facts, "reactant", other).get("element") == "N":
                    formed.append((carbon, other))
        if not formed:
            continue
        carbon, nitrogen = formed[0]
        out.append(
            {
                "event": "aza_michael_addition",
                "role": "specific_event",
                "subtype": "amine_conjugate_addition",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "alkene_maps": sorted((a, b)),
                "support_maps": _support(a, b, nitrogen),
                "mapping_scope": "mapped_activated_c_c_pi_bond_consumption_plus_external_c_n_formation",
            }
        )
    return out


def detect_imine_reduction(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for oc in facts.get("bond_order_changes", []) or []:
        if not (
            float(oc.get("reactant_order", 0.0)) >= 1.8
            and float(oc.get("product_order", 0.0)) <= 1.2
        ):
            continue
        a, b = int(oc["a"]), int(oc["b"])
        if {
            _atom(facts, "reactant", a).get("element"),
            _atom(facts, "reactant", b).get("element"),
        } != {"C", "N"}:
            continue
        carbon, nitrogen = (
            (a, b) if _atom(facts, "reactant", a).get("element") == "C" else (b, a)
        )
        # Avoid heterocumulene/amide resonance changes and reactions with new heavy-atom bonds at the pair.
        if (
            _has_conserved_carbonyl(facts, carbon)
            or _is_isocyanate_carbon(facts, "reactant", carbon)
            or _is_isothiocyanate_carbon(facts, "reactant", carbon)
        ):
            continue
        if any(
            carbon in (int(x["a"]), int(x["b"]))
            or nitrogen in (int(x["a"]), int(x["b"]))
            for x in facts.get("formed_bonds", []) or []
        ):
            continue
        if _hydrogen_delta(facts, carbon) + _hydrogen_delta(facts, nitrogen) <= 0:
            continue
        out.append(
            {
                "event": "imine_reduction",
                "role": "specific_event",
                "subtype": "c_n_double_to_single",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "support_maps": _support(carbon, nitrogen),
                "mapping_scope": "mapped_c_n_double_to_single_with_net_h_gain_and_no_new_heavy_bond",
                "redox_direction": "reduction",
            }
        )
    return out


def detect_v67_carbon_to_nitrogen_migration(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """C(source)-C(migrator) -> N(target)-C(migrator).

    Cross-side bonds require retained atom maps. Functional context inspects
    all local attachments, including neighbors without atom maps. This is a
    formal mapped transformation classifier, not proof of reaction mechanism.
    """
    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int, int]] = set()

    def atoms(side: str, mn: int) -> dict[str, Any]:
        return _atom(facts, side, mn)

    def attached(side: str, mn: int) -> list[dict[str, Any]]:
        return _attachments(facts, side, mn)

    def order(side: str, a: int, b: int) -> float:
        return float(_bond_between(facts, side, a, b) or 0.0)

    def double_oxygen(side: str, carbon: int) -> bool:
        return any(a.get('element') == 'O'
                   and float(a.get('bond_order', 0.0) or 0.0) >= 1.8
                   for a in attached(side, carbon))

    def is_nitro_or_nitrile_n(n: int) -> bool:
        rn = atoms('reactant', n)
        if rn.get('element') != 'N':
            return True
        na = attached('reactant', n)
        oxygens = sum(q.get('element') == 'O' for q in na)
        if int(rn.get('charge', 0) or 0) > 0 and oxygens >= 2:
            return True
        if any(q.get('element') == 'C'
               and float(q.get('bond_order', 0.0) or 0.0) >= 2.8
               for q in na):
            return True
        return False

    # Iterating common maps prevents unmapped context atoms from becoming
    # cross-side correspondence anchors.
    common = {int(m) for m in facts.get('mapping', {}).get('common_maps', [])}
    for source in sorted(common):
        if atoms('reactant', source).get('element') != 'C':
            continue
        if atoms('product', source).get('element') != 'C':
            continue
        if not (
            double_oxygen('product', source)
            or any(q.get('element') == 'N'
                   and float(q.get('bond_order', 0.0) or 0.0) >= 1.8
                   for q in attached('reactant', source))
        ):
            continue

        # Source--migrator must be genuinely broken.
        for rroot in attached('reactant', source):
            if rroot.get('element') != 'C' or rroot.get('map') is None:
                continue
            root = int(rroot['map'])
            if root not in common or order('product', source, root) > 0.2:
                continue
            if atoms('product', root).get('element') != 'C':
                continue

            # Product: migrator--nitrogen new; source--same-N present.
            for ptarget in attached('product', root):
                if ptarget.get('element') != 'N' or ptarget.get('map') is None:
                    continue
                n = int(ptarget['map'])
                if n not in common or is_nitro_or_nitrile_n(n):
                    continue
                if order('reactant', root, n) > 0.2:
                    continue
                if order('product', source, n) <= 0.2:
                    continue

                rn = attached('reactant', n)
                n_at_source = order('reactant', source, n) > 0.2
                if not n_at_source and not any(
                    a.get('element') in {'N','O'} for a in rn
                ):
                    continue

                if n_at_source:
                    if order('reactant', source, n) >= 1.8:
                        if not any(a.get('element') == 'O' for a in rn):
                            continue
                        family, subtype = (
                          'RB1_OXIME_TO_AMIDE_MIGRATION', 'oxime_to_amide_migration'
                        )
                    elif double_oxygen('reactant', source):
                        if any(a.get('element') == 'N' for a in rn):
                            family, subtype = (
                              'RC1_ACYL_AZIDE_DERIVED_C_TO_N_MIGRATION',
                              'acyl_azide_derived_c_to_n_migration'
                            )
                        elif any(a.get('element') == 'O' for a in rn):
                            family, subtype = (
                              'RC2_HYDROXAMATE_DERIVED_C_TO_N_MIGRATION',
                              'hydroxamate_derived_c_to_n_migration'
                            )
                        else:
                            family, subtype = (
                              'RC3_AMIDE_DERIVED_C_TO_N_MIGRATION',
                              'amide_derived_c_to_n_migration'
                            )
                    else:
                        continue
                else:
                    if not double_oxygen('reactant', source):
                        continue
                    acyl_oxygen = any(
                        a.get('element') == 'O'
                        and 0.8 <= float(a.get('bond_order', 0.0) or 0.0) <= 1.2
                        for a in attached('reactant', source)
                    )
                    acyl_halide = any(
                        a.get('element') in {'F','Cl','Br','I'}
                        for a in attached('reactant', source)
                    )
                    if acyl_oxygen or acyl_halide:
                        family, subtype = (
                          'RC4_EXTERNAL_AZIDE_ACYL_C_TO_N_MIGRATION',
                          'external_azide_acyl_c_to_n_migration'
                        )
                    else:
                        family, subtype = (
                          'RB2_EXTERNAL_N_CARBONYL_TO_AMIDE_INSERTION',
                          'external_n_carbonyl_to_amide_insertion'
                        )

                key = (source, root, n)
                if key in seen:
                    continue
                seen.add(key)
                out.append({
                    'event': 'carbon_to_nitrogen_migration',
                    'role': 'specific_event',
                    'family': family,
                    'subtype': subtype,
                    'derived_parent': 'rearrangement',
                    'center_map': source,
                    'source_center_map': source,
                    'migration_carbon_map': root,
                    'migrating_root_map': root,
                    'target_n_map': n,
                    'support_maps': _support(source, root, n),
                    'mapping_scope': 'v67_broken_c_c_new_c_n_with_context_validated_target_n',
                })
    return out


def detect_hofmann_rearrangement(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Amide N migrates from a lost carbonyl carbon to its carbon substituent."""
    out = []
    for d in facts.get("attachment_deltas", []) or []:
        n = int(d["center_map"])
        if _atom(facts, "reactant", n).get("element") != "N":
            continue
        lost_carbonyl = [
            x
            for x in d.get("lost", []) or []
            if x.get("element") == "C" and bool(x.get("carbonyl_like"))
        ]
        gained_c = [
            x
            for x in d.get("gained", []) or []
            if x.get("element") == "C" and _att_map(x) is not None
        ]
        if not lost_carbonyl or not gained_c:
            continue
        cm = int(_att_map(gained_c[0]))
        cd = _attachment_delta(facts, cm)
        if cd is None or not any(
            x.get("element") == "C" and bool(x.get("carbonyl_like"))
            for x in cd.get("lost", []) or []
        ):
            continue
        out.append(
            {
                "event": "hofmann_rearrangement",
                "role": "specific_event",
                "subtype": "amide_to_amine_one_carbon_loss",
                "center_map": n,
                "migration_carbon_map": cm,
                "support_maps": _support(n, cm),
                "mapping_scope": "mapped_amide_n_and_carbon_substituent_both_lose_same_carbonyl_attachment_then_form_c_n",
            }
        )
    return out


def detect_nitrile_thioamide_and_ester(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for c in facts["mapping"]["common_maps"]:
        if not _organic_nitrile_center(facts, "reactant", int(c)):
            continue
        p_atts = _attachments(facts, "product", int(c))
        has_s_double = any(
            x.get("element") == "S" and float(x.get("bond_order", 0.0)) >= 1.8
            for x in p_atts
        )
        has_n_single = any(
            x.get("element") == "N" and float(x.get("bond_order", 0.0)) <= 1.2
            for x in p_atts
        )
        if has_s_double and has_n_single:
            out.extend(
                _conversion_records(
                    "nitrile_to_thioamide",
                    int(c),
                    _conversion_support(facts, int(c)),
                    "neutral",
                    "functional_group_interconversion",
                    "nitrile",
                    "thioamide",
                    True,
                    specific_event_name="nitrile_to_thioamide",
                )
            )
            continue
        if _carbon_fg_state(facts, "product", int(c)) == "ester":
            out.extend(
                _conversion_records(
                    "nitrile_to_ester",
                    int(c),
                    _conversion_support(facts, int(c)),
                    "neutral",
                    "functional_group_interconversion",
                    "nitrile",
                    "ester",
                    True,
                    specific_event_name="nitrile_to_ester",
                )
            )
    return out


def detect_guanylation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for fb in facts.get("intercomponent_formed_bonds", []) or []:
        a, b = int(fb["a"]), int(fb["b"])
        ea, eb = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", b
        ).get("element")
        if {ea, eb} != {"C", "N"}:
            continue
        carbon, nitrogen = (a, b) if ea == "C" else (b, a)
        r_atts = _attachments(facts, "reactant", carbon)
        if not any(
            x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 2.8
            for x in r_atts
        ):
            continue
        if not any(
            x.get("element") == "N" and float(x.get("bond_order", 0.0)) <= 1.2
            for x in r_atts
        ):
            continue
        out.append(
            {
                "event": "guanylation",
                "role": "specific_event",
                "subtype": "amine_to_cyanamide_guanylation",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "support_maps": _conversion_support(facts, carbon)
                + (
                    [nitrogen]
                    if nitrogen not in _conversion_support(facts, carbon)
                    else []
                ),
                "mapping_scope": "mapped_amine_addition_to_cyanamide_carbon_with_c_n_triple_reorganization",
            }
        )
    return out


def detect_chan_lam_n_arylation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for fb in facts.get("intercomponent_formed_bonds", []) or []:
        if float(fb.get("order", 0.0)) > 1.2:
            continue
        a, b = int(fb["a"]), int(fb["b"])
        ea, eb = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", b
        ).get("element")
        if {ea, eb} != {"C", "N"}:
            continue
        carbon, nitrogen = (a, b) if ea == "C" else (b, a)
        if not any(
            x.get("element") == "B" for x in _attachments(facts, "reactant", carbon)
        ):
            continue
        if _hydrogen_delta(facts, nitrogen) > -1:
            continue
        out.append(
            {
                "event": "chan_lam_n_arylation",
                "role": "specific_event",
                "subtype": "organoboron_n_arylation",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "support_maps": _support(carbon, nitrogen),
                "mapping_scope": "mapped_organoboron_carbon_loses_b_and_forms_intercomponent_c_n_bond",
            }
        )
    return out



def detect_n_phosphorus_functionalization(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """External N-H + P-Cl -> P-N product-state formation.

    Restrict to a single mapped intercomponent bond, a conserved N-H loss,
    and depletion of a P-Cl leaving group. Composite P-N/C-P ring formation
    is intentionally not classified by this rule. No mechanism is inferred.
    """
    if (len(facts.get("formed_bonds", [])) != 1
            or facts.get("broken_bonds")
            or facts.get("bond_order_changes")):
        return []
    out = []
    for fb in facts.get("intercomponent_formed_bonds", []):
        if not (0.8 <= float(fb.get("order", 0.0)) <= 1.2):
            continue
        a, b = int(fb["a"]), int(fb["b"])
        ea, eb = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", b
        ).get("element")
        if {ea, eb} != {"P", "N"}:
            continue
        phosphorus, nitrogen = (a, b) if ea == "P" else (b, a)
        if _hydrogen_delta(facts, nitrogen) != -1:
            continue
        if not (0.8 <= _bond_between(facts, "product", phosphorus, nitrogen) <= 1.2):
            continue
        def count_pcl(side: str) -> int:
            return sum(
                x.get("element") == "Cl"
                and 0.8 <= float(x.get("bond_order", 0.0)) <= 1.2
                for x in _attachments(facts, side, phosphorus)
            )
        if count_pcl("reactant") < 1 or count_pcl("product") >= count_pcl("reactant"):
            continue
        p_shell = _attachments(facts, "product", phosphorus)
        if any(
            x.get("element") == "O"
            and float(x.get("bond_order", 0.0)) >= 1.8
            for x in p_shell
        ):
            # Distinguish phosphoryl amides with ester O from chlorophosphoramides.
            subtype = (
                "phosphoramidate_formation"
                if any(x.get("element") == "O"
                       and float(x.get("bond_order", 0.0)) < 1.8
                       for x in p_shell)
                else "phosphoramidic_chloride_formation"
            )
        elif sum(
            x.get("element") == "O"
            and float(x.get("bond_order", 0.0)) <= 1.2
            for x in p_shell
        ) >= 2:
            subtype = "phosphoramidite_formation"
        elif sum(
            x.get("element") == "C"
            and float(x.get("bond_order", 0.0)) <= 1.2
            for x in p_shell
        ) >= 2:
            subtype = "aminophosphine_formation"
        else:
            continue  # unknown P state; no catch-all assignment
        out.append({
            "event": "n_phosphorus_functionalization",
            "role": "specific_event",
            "subtype": subtype,
            "center_map": phosphorus,
            "support_maps": _support(phosphorus, nitrogen),
            "mapping_scope": "mapped_external_p_cl_n_h_to_p_n_single_product_state",
        })
    return out

def _p2_oxo_heteroarene_phosphorothioyl_support(
    facts: dict[str, Any], oxygen: int, phosphorus: int
) -> list[int]:
    """Structural O-P(=S) endpoint from a conserved oxo heteroarene.

    Only the *mapped product state* is asserted. The oxo O becomes an O-P
    single-bond substituent on the same aromatic carbon as adjacent N-H
    loses H; the reagent is a P(=S)-Cl electrophile. No mechanism assigned.
    """
    if _hydrogen_delta(facts, oxygen) != 0:
        return []
    if not any(
        a.get("element") == "S" and float(a.get("bond_order", 0)) >= 1.8
        for a in _attachments(facts, "reactant", phosphorus)
    ):
        return []
    if not any(
        a.get("element") == "S" and float(a.get("bond_order", 0)) >= 1.8
        for a in _attachments(facts, "product", phosphorus)
    ):
        return []
    def n_pcl(side: str) -> int:
        return sum(
            a.get("element") == "Cl" and 0.8 <= float(a.get("bond_order", 0)) <= 1.2
            for a in _attachments(facts, side, phosphorus)
        )
    if n_pcl("reactant") < 1 or n_pcl("product") >= n_pcl("reactant"):
        return []
    for change in facts.get("bond_order_changes", []) or []:
        a, b = int(change["a"]), int(change["b"])
        if oxygen not in (a, b):
            continue
        carbon = b if a == oxygen else a
        if not (
            float(change.get("reactant_order", 0)) >= 1.8
            and 0.8 <= float(change.get("product_order", 0)) <= 1.2
        ):
            continue
        if not all(
            _atom(facts, side, carbon).get("element") == "C"
            and bool(_atom(facts, side, carbon).get("aromatic"))
            for side in ("reactant", "product")
        ):
            continue
        old_n = {
            _att_map(a) for a in _attachments(facts, "reactant", carbon)
            if a.get("element") == "N" and _att_map(a) is not None
        }
        new_n = {
            _att_map(a) for a in _attachments(facts, "product", carbon)
            if a.get("element") == "N" and _att_map(a) is not None
        }
        for nitrogen in sorted(old_n & new_n):
            if _hydrogen_delta(facts, nitrogen) == -1:
                return [carbon, nitrogen]
    return []


def detect_o_phosphorylation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for fb in facts.get("intercomponent_formed_bonds", []) or []:
        if float(fb.get("order", 0.0)) > 1.2:
            continue
        a, b = int(fb["a"]), int(fb["b"])
        ea, eb = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", b
        ).get("element")
        if {ea, eb} != {"O", "P"}:
            continue
        oxygen, phosphorus = (a, b) if ea == "O" else (b, a)
        oxo_support = []
        if _hydrogen_delta(facts, oxygen) > -1:
            oxo_support = _p2_oxo_heteroarene_phosphorothioyl_support(
                facts, oxygen, phosphorus
            )
            if not oxo_support:
                continue
        subtype = (
            "phosphorothioate_formation"
            if any(
                x.get("element") == "S" and float(x.get("bond_order", 0.0)) >= 1.8
                for x in _attachments(facts, "reactant", phosphorus)
            )
            else "phosphate_formation"
        )
        out.append(
            {
                "event": "o_phosphorylation",
                "role": "specific_event",
                "subtype": subtype,
                "center_map": phosphorus,
                "oxygen_map": oxygen,
                "support_maps": _support(phosphorus, oxygen, *oxo_support),
                "mapping_scope": (
                    "mapped_oxo_heteroarene_o_p_phosphorothioate_product_state"
                    if oxo_support else "mapped_intercomponent_o_p_formation_with_o_h_loss"
                ),
            }
        )
    return out


def detect_organozinc_c_s_coupling(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for fb in facts.get("intercomponent_formed_bonds", []) or []:
        if float(fb.get("order", 0.0)) > 1.2:
            continue
        a, b = int(fb["a"]), int(fb["b"])
        ea, eb = _atom(facts, "reactant", a).get("element"), _atom(
            facts, "reactant", b
        ).get("element")
        if {ea, eb} != {"C", "S"}:
            continue
        carbon, sulfur = (a, b) if ea == "C" else (b, a)
        if not any(
            x.get("element") == "Zn" for x in _attachments(facts, "reactant", carbon)
        ):
            continue
        out.append(
            {
                "event": "organozinc_c_s_coupling",
                "role": "specific_event",
                "subtype": "organozinc_thiol_c_s_coupling",
                "center_map": carbon,
                "sulfur_map": sulfur,
                "support_maps": _support(carbon, sulfur),
                "mapping_scope": "mapped_c_zn_donor_carbon_forms_intercomponent_c_s_bond",
            }
        )
    return out


def detect_hydrodealkoxylation(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for d in facts.get("attachment_deltas", []) or []:
        c = int(d["center_map"])
        if _atom(facts, "reactant", c).get("element") != "C" or not _atom(
            facts, "reactant", c
        ).get("aromatic"):
            continue
        lost_o = [x for x in d.get("lost", []) or [] if x.get("element") == "O"]
        if not lost_o or d.get("gained") or int(d.get("hydrogen_delta", 0)) <= 0:
            continue
        out.append(
            {
                "event": "hydrodealkoxylation",
                "role": "specific_event",
                "subtype": "aryl_or_heteroaryl_c_o_to_c_h",
                "center_map": c,
                "support_maps": _support(c, *[_att_map(x) for x in lost_o]),
                "mapping_scope": "mapped_aromatic_carbon_loses_o_attachment_and_gains_h",
                "redox_direction": "reduction",
            }
        )
    return out


def detect_s1_semantic_closure(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Conservative semantic-closure fallbacks discovered from the S1 named-reaction audit.

    These rules encode mapped structural transformations only. They do not use
    reaction names, reagents, solvents, catalysts, or method knowledge.

    The detector intentionally abstains on complex rearrangements and
    non-diagnostic/no-edit benchmark representations.
    """
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, int | None]] = set()

    def atom(side: str, mn: int) -> dict[str, Any]:
        try:
            return _atom(facts, side, int(mn)) or {}
        except Exception:
            return {}

    def atts(side: str, mn: int) -> list[dict[str, Any]]:
        try:
            return list(_attachments(facts, side, int(mn)) or [])
        except Exception:
            return []

    def bond(side: str, a: int, b: int) -> float:
        try:
            return float(_bond_between(facts, side, int(a), int(b)) or 0.0)
        except Exception:
            return 0.0

    def fg(side: str, mn: int) -> str:
        try:
            return str(_carbon_fg_state(facts, side, int(mn)) or "")
        except Exception:
            return ""

    def emit_specific(
        event: str,
        center: int | None,
        subtype: str | None = None,
        support=None,
        redox: str | None = None,
        extra: dict[str, Any] | None = None,
    ) -> None:
        key = (event, center)
        if key in seen:
            return
        seen.add(key)

        rec = {
            "event": event,
            "role": "specific_event",
            "subtype": subtype or event,
            "center_map": center,
            "support_maps": list(
                support or ([] if center is None else _support(center))
            ),
            "mapping_scope": "s1_mapped_structural_semantic_closure",
        }
        if redox:
            rec["redox_direction"] = redox
        if extra:
            rec.update(extra)
        out.append(rec)

    def emit_conversion(
        name: str,
        center: int,
        direction: str,
        source: str,
        target: str,
        specific: str | None = None,
        support=None,
    ) -> None:
        key = (specific or name, center)
        if key in seen:
            return
        seen.add(key)

        family = (
            direction
            if direction in {"oxidation", "reduction"}
            else "functional_group_interconversion"
        )

        out.extend(
            _conversion_records(
                name,
                center,
                support or _conversion_support(facts, center),
                direction,
                family,
                source,
                target,
                True,
                specific_event_name=specific,
                extra={"mapping_scope": "s1_mapped_structural_semantic_closure"},
            )
        )

    # ------------------------------------------------------------------
    # Collect mapped centers touched by the reaction.
    # ------------------------------------------------------------------
    deltas: dict[int, dict[str, Any]] = {}
    maps: set[int] = set()

    for d in facts.get("attachment_deltas", []) or []:
        cm = d.get("center_map")
        if cm is None:
            continue
        c = int(cm)
        maps.add(c)
        deltas[c] = d

        for x in list(d.get("lost", []) or []) + list(d.get("gained", []) or []):
            m = _att_map(x)
            if m is not None:
                maps.add(int(m))

    for fb in facts.get("intercomponent_formed_bonds", []) or []:
        for k in ("a", "b"):
            if fb.get(k) is not None:
                maps.add(int(fb[k]))

    # Different graph-build versions have used slightly different names.
    # Reading absent keys is harmless.
    for key in (
        "bond_order_changes",
        "bond_order_deltas",
        "changed_bonds",
        "bond_changes",
    ):
        for rec in facts.get(key, []) or []:
            if not isinstance(rec, dict):
                continue
            for k in ("a", "b", "atom1", "atom2", "map1", "map2"):
                v = rec.get(k)
                if isinstance(v, int):
                    maps.add(int(v))

    def elems(xs):
        return [str(x.get("element") or "") for x in xs]

    def lost(c):
        return list(deltas.get(c, {}).get("lost", []) or [])

    def gained(c):
        return list(deltas.get(c, {}).get("gained", []) or [])

    def hdelta(c):
        try:
            return int(
                deltas.get(c, {}).get("hydrogen_delta", _hydrogen_delta(facts, c))
            )
        except Exception:
            return 0

    # S1_CONFLICT_FIX_V2
    #
    # These helpers make the closure predicates topology-conserving.
    # "Exact" semantic conversions must prove that the relevant mapped
    # skeleton is preserved rather than infer reaction identity from a
    # single local bond-order/attachment change.

    def mapped_heavy_neighbors(
        side: str,
        center: int,
        *,
        exclude_maps=(),
        exclude_elements=(),
    ) -> set[int]:
        excluded_maps = {int(x) for x in exclude_maps if x is not None}
        excluded_elements = set(exclude_elements)

        out = set()
        for x in atts(side, center):
            el = str(x.get("element") or "")
            if not el or el == "H" or el in excluded_elements:
                continue
            m = _att_map(x)
            if m is None:
                continue
            m = int(m)
            if m not in excluded_maps:
                out.add(m)
        return out

    def mapped_element_neighbors(
        side: str,
        center: int,
        element: str,
    ) -> set[int]:
        out = set()
        for x in atts(side, center):
            if x.get("element") != element:
                continue
            m = _att_map(x)
            if m is not None:
                out.add(int(m))
        return out

    def intercomponent_bond_was_formed(a: int, b: int) -> bool:
        target = {int(a), int(b)}
        for fb in facts.get("intercomponent_formed_bonds", []) or []:
            if {int(fb["a"]), int(fb["b"])} == target:
                return True
        return False

    # ------------------------------------------------------------------
    # 1. Generic C-H halogenation fallback.
    #
    # C-H -> C-X requires genuine halogen installation. A mapped halogen
    # migrating from elsewhere in the same substrate is not halogenation.
    # The non-halogen heavy skeleton at the carbon must remain unchanged.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if atom("reactant", c).get("element") != "C":
            continue

        gx = [x for x in gained(c) if x.get("element") in HALOGENS]

        if not gx or hdelta(c) >= 0:
            continue

        external_x = []

        for x in gx:
            xm = _att_map(x)
            if xm is None:
                continue

            xm = int(xm)

            existed_on_reactant = bool(atom("reactant", xm).get("element"))

            if not existed_on_reactant or intercomponent_bond_was_formed(c, xm):
                external_x.append(x)

        if not external_x:
            continue

        # A true C-H halogenation does not simultaneously replace or
        # cleave another non-halogen heavy-atom attachment at this carbon.
        lost_nonhalogen_heavy = [
            x
            for x in lost(c)
            if x.get("element") not in {"", "H"} and x.get("element") not in HALOGENS
        ]

        if lost_nonhalogen_heavy:
            continue

        r_skeleton = mapped_heavy_neighbors(
            "reactant",
            c,
            exclude_elements=HALOGENS,
        )

        p_skeleton = mapped_heavy_neighbors(
            "product",
            c,
            exclude_elements=HALOGENS,
        )

        if r_skeleton != p_skeleton:
            continue

        subtype = "c_h_halogenation"

        alpha_acid = False
        for a in atts("reactant", c):
            m = _att_map(a)
            if m is not None and fg("reactant", int(m)) == "carboxylic_acid":
                alpha_acid = True
                break

        if alpha_acid:
            subtype = "alpha_carboxylic_acid_halogenation"

        emit_specific(
            "c_h_halogenation",
            c,
            subtype=subtype,
            support=_support(
                c,
                *[_att_map(x) for x in external_x],
            ),
        )

    # ------------------------------------------------------------------
    # 2. Aliphatic single-bond C-O -> C-H deoxygenation.
    #
    # S1_CONFLICT_FIX_V3
    #
    # This event must represent an isolated replacement of a
    # single-bond oxygen substituent by hydrogen.
    #
    # Reject:
    #   * carbonyl reductions
    #   * oxygen migration
    #   * coupled oxygen reorganization elsewhere in the mapped substrate
    #
    # Barton-McCombie / Marko-Lam preserve the carbon skeleton and
    # simply remove the oxygen functionality at the reacting carbon.
    # ------------------------------------------------------------------

    def oxygen_attachment_signature(
        side: str,
        center: int,
    ) -> tuple[tuple[int | None, float], ...]:

        sig = []

        for x in atts(side, center):
            if x.get("element") != "O":
                continue

            m = _att_map(x)

            sig.append(
                (
                    None if m is None else int(m),
                    round(
                        float(x.get("bond_order", 0.0)),
                        2,
                    ),
                )
            )

        return tuple(
            sorted(
                sig,
                key=lambda z: (
                    -1 if z[0] is None else z[0],
                    z[1],
                ),
            )
        )

    for c in sorted(maps):

        ra = atom("reactant", c)

        if ra.get("element") != "C" or ra.get("aromatic"):
            continue

        lo = [x for x in lost(c) if x.get("element") == "O"]

        reactant_o_single = [
            x
            for x in atts("reactant", c)
            if (
                x.get("element") == "O"
                and 0.8 <= float(x.get("bond_order", 0.0)) <= 1.2
            )
        ]

        reactant_o_double = [
            x
            for x in atts("reactant", c)
            if (x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8)
        ]

        product_o = [x for x in atts("product", c) if x.get("element") == "O"]

        heavy_gained = [x for x in gained(c) if x.get("element") not in {"", "H"}]

        lost_non_o_heavy = [
            x for x in lost(c) if x.get("element") not in {"", "H", "O"}
        ]

        r_skeleton = mapped_heavy_neighbors(
            "reactant",
            c,
            exclude_elements={"O"},
        )

        p_skeleton = mapped_heavy_neighbors(
            "product",
            c,
            exclude_elements={"O"},
        )

        # ----------------------------------------------------------
        # Oxygen provenance.
        #
        # If an oxygen lost from the center is mapped and survives
        # in the product attached elsewhere, this is migration /
        # rearrangement rather than deoxygenation.
        # ----------------------------------------------------------

        lost_o_maps = {int(m) for x in lo if (m := _att_map(x)) is not None}

        oxygen_migrated = False

        for om in lost_o_maps:

            po = atom("product", om)

            if po.get("element") != "O":
                continue

            # The mapped oxygen survives in the product but is no
            # longer attached to the original carbon.
            if (
                bond(
                    "product",
                    c,
                    om,
                )
                <= 0.1
            ):
                oxygen_migrated = True
                break

        if oxygen_migrated:
            continue

        # ----------------------------------------------------------
        # Isolation test.
        #
        # Simple deoxygenation should not coexist with a separate
        # mapped C-O installation, removal, or bond-order change at
        # another carbon.
        #
        # This rejects Pinacol-like and Ferrier-like global oxygen
        # reorganization while preserving simple C-O -> C-H events.
        # ----------------------------------------------------------

        other_oxygen_reorganization = False

        for d in sorted(maps):

            if d == c:
                continue

            if (
                atom("reactant", d).get("element") != "C"
                and atom("product", d).get("element") != "C"
            ):
                continue

            rsig = oxygen_attachment_signature(
                "reactant",
                d,
            )

            psig = oxygen_attachment_signature(
                "product",
                d,
            )

            if rsig != psig:
                other_oxygen_reorganization = True
                break

        if other_oxygen_reorganization:
            continue

        if (
            lo
            and reactant_o_single
            and not reactant_o_double
            and not product_o
            and not heavy_gained
            and not lost_non_o_heavy
            and hdelta(c) > 0
            and r_skeleton == p_skeleton
        ):
            emit_conversion(
                "carbon_oxygen_to_hydrogen_deoxygenation",
                c,
                "reduction",
                "c_o_substituted_carbon",
                "hydrocarbon",
            )

    # ------------------------------------------------------------------
    # 3. Fleming-Tamao-type WHAT: organosilane carbon -> alcohol carbon.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if atom("reactant", c).get("element") != "C":
            continue

        lsi = [x for x in lost(c) if x.get("element") == "Si"]
        go = [x for x in gained(c) if x.get("element") == "O"]

        if lsi and go:
            emit_conversion(
                "organosilane_to_alcohol",
                c,
                "oxidation",
                "organosilane_carbon",
                "alcohol",
            )

    # ------------------------------------------------------------------
    # 4. Same-center high-value directional conversions.
    #
    # Acid/amide conversions require preservation of the mapped carbon
    # substituent skeleton. This prevents donor carbonyls participating
    # in aromatic formylation or malformed reaction pairs from being
    # interpreted as simple same-center FG conversions.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if atom("reactant", c).get("element") != "C":
            continue

        r = fg("reactant", c)
        pfg = fg("product", c)

        r_carbons = mapped_element_neighbors("reactant", c, "C")
        p_carbons = mapped_element_neighbors("product", c, "C")

        carbon_skeleton_conserved = r_carbons == p_carbons

        if r == "carboxylic_acid" and pfg == "nitrile" and carbon_skeleton_conserved:
            emit_conversion(
                "carboxylic_acid_to_nitrile",
                c,
                "neutral",
                r,
                pfg,
            )

        if (
            r == "amide"
            and pfg == "aldehyde"
            and carbon_skeleton_conserved
            and hdelta(c) > 0
        ):
            emit_conversion(
                "amide_to_aldehyde",
                c,
                "reduction",
                r,
                pfg,
            )

        if (
            r == "carboxylic_acid"
            and pfg == "aldehyde"
            and carbon_skeleton_conserved
            and hdelta(c) > 0
        ):
            emit_conversion(
                "carboxylic_acid_to_aldehyde",
                c,
                "reduction",
                r,
                pfg,
            )

        # Carbonyl -> alkyne.
        if r in {"aldehyde", "ketone"}:
            pc = [
                x
                for x in atts("product", c)
                if x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 2.8
            ]
            ro = [
                x
                for x in atts("reactant", c)
                if x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
            ]

            if pc and ro:
                emit_specific(
                    "carbonyl_to_alkyne",
                    c,
                    subtype="deoxygenative_carbonyl_to_alkyne",
                    support=_support(
                        c,
                        *[_att_map(x) for x in pc],
                    ),
                )

        # Acyl halide -> ketene.
        if r == "acyl_halide":
            pc_double = [
                x
                for x in atts("product", c)
                if x.get("element") == "C" and float(x.get("bond_order", 0.0)) >= 1.8
            ]
            po_double = [
                x
                for x in atts("product", c)
                if x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
            ]

            if pc_double and po_double:
                emit_conversion(
                    "acyl_halide_to_ketene",
                    c,
                    "neutral",
                    "acyl_halide",
                    "ketene",
                )

    # ------------------------------------------------------------------
    # 5. Pinner-family structural products from an organic nitrile.
    #
    # A nitrile carbon incorporated into a newly formed/aromatic ring
    # must not be classified as a simple amidine/imidate/orthoester FGI.
    # Preserve the mapped carbon substituent skeleton.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if not _organic_nitrile_center(
            facts,
            "reactant",
            c,
        ):
            continue

        if _ring_formation_touches(facts, c):
            continue

        if atom("product", c).get("aromatic"):
            continue

        r_carbons = mapped_element_neighbors("reactant", c, "C")
        p_carbons = mapped_element_neighbors("product", c, "C")

        if r_carbons != p_carbons:
            continue

        patts = atts("product", c)

        n_atts = [x for x in patts if x.get("element") == "N"]

        o_atts = [x for x in patts if x.get("element") == "O"]

        # Amidine = C attached to exactly two nitrogens,
        # with one imine-like and one single-bond N.
        if len(n_atts) == 2:
            n_orders = sorted(float(x.get("bond_order", 0.0)) for x in n_atts)

            if n_orders[0] <= 1.2 and n_orders[1] >= 1.5:
                emit_conversion(
                    "nitrile_to_amidine",
                    c,
                    "neutral",
                    "nitrile",
                    "amidine",
                )

        elif len(o_atts) >= 3 and not n_atts:
            emit_conversion(
                "nitrile_to_orthoester",
                c,
                "neutral",
                "nitrile",
                "orthoester",
            )

        elif n_atts and o_atts:
            if any(float(x.get("bond_order", 0.0)) >= 1.5 for x in n_atts):
                emit_conversion(
                    "nitrile_to_imidate",
                    c,
                    "neutral",
                    "nitrile",
                    "imidate",
                )

    # ------------------------------------------------------------------
    # 6. Amine -> diazonium structural conversion.
    # ------------------------------------------------------------------
    for n in sorted(maps):
        if atom("reactant", n).get("element") != "N":
            continue

        r_c = [x for x in atts("reactant", n) if x.get("element") == "C"]
        p_n = [
            x
            for x in atts("product", n)
            if x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 2.5
        ]

        if r_c and p_n:
            aromatic_attachment = False
            for x in r_c:
                m = _att_map(x)
                if m is not None and atom("reactant", int(m)).get("aromatic"):
                    aromatic_attachment = True
                    break

            if aromatic_attachment:
                emit_specific(
                    "amine_to_diazonium",
                    n,
                    subtype="aryl_amine_to_diazonium",
                    support=_support(n, *[_att_map(x) for x in p_n]),
                )

    # ------------------------------------------------------------------
    # 7. Primary amine + polyhalomethane -> isocyanide topology.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if atom("reactant", c).get("element") != "C":
            continue

        r_hal = [x for x in atts("reactant", c) if x.get("element") in HALOGENS]
        p_n = [
            x
            for x in atts("product", c)
            if x.get("element") == "N" and float(x.get("bond_order", 0.0)) >= 2.8
        ]

        if len(r_hal) >= 2 and p_n:
            emit_specific(
                "primary_amine_to_isocyanide",
                c,
                subtype="isocyanide_formation",
                support=_support(c, *[_att_map(x) for x in p_n]),
            )

    # ------------------------------------------------------------------
    # 8. Quaternary/tertiary N-dealkylation fallback.
    # ------------------------------------------------------------------
    for n in sorted(maps):
        if atom("reactant", n).get("element") != "N":
            continue

        rc = [x for x in atts("reactant", n) if x.get("element") == "C"]
        pc = [x for x in atts("product", n) if x.get("element") == "C"]
        # P1_STAGE2_N_ACYL_PRECISION:
        # Removal of an N-acyl / N-carbamate carbonyl carbon is not
        # ordinary alkyl cleavage. Keep the N-C bond loss as graph evidence
        # and leave N-Boc/Cbz/deacylation to the protecting-group detectors.
        lc = [
            x for x in lost(n)
            if x.get("element") == "C"
            and not bool(x.get("carbonyl_like"))
            and not _lost_aziridine_ring_n_c_bond(facts, n, _att_map(x))
        ]

        if (
            len(rc) >= 3
            and len(pc) < len(rc)
            and lc
            and not any(float(x.get("bond_order", 0.0)) >= 2.8 for x in pc)
        ):
            emit_specific(
                "n_dealkylation",
                n,
                subtype="mapped_n_c_cleavage_dealkylation",
                support=_support(n, *[_att_map(x) for x in lc]),
                redox="reduction",
            )

    # ------------------------------------------------------------------
    # 9. Amino-acid alpha carbon -> aldehyde.
    #
    # The reacting carbon is an sp3 alpha carbon bearing N, not an
    # amide/acyl carbon. The benchmark may omit mapping on the departing
    # carboxyl carbon, so use explicit N attachment + lost N/C evidence
    # while vetoing any reactant C=O at the center.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if atom("reactant", c).get("element") != "C":
            continue

        if fg("product", c) != "aldehyde":
            continue

        reactant_double_o = any(
            x.get("element") == "O" and float(x.get("bond_order", 0.0)) >= 1.8
            for x in atts("reactant", c)
        )

        if reactant_double_o:
            continue

        reactant_n = any(
            x.get("element") == "N" and float(x.get("bond_order", 0.0)) <= 1.2
            for x in atts("reactant", c)
        )

        le = elems(lost(c))

        if reactant_n and "N" in le and "C" in le:
            emit_conversion(
                "alpha_amino_acid_to_aldehyde",
                c,
                "neutral",
                "alpha_amino_acid_carbon",
                "aldehyde",
            )

    # ------------------------------------------------------------------
    # 10. Pi-bond reductions / hydration / oxidation.
    # ------------------------------------------------------------------
    carbon_maps = sorted(
        m
        for m in maps
        if atom("reactant", m).get("element") == "C"
        and atom("product", m).get("element") == "C"
    )

    for i, a in enumerate(carbon_maps):
        for b in carbon_maps[i + 1 :]:
            rb = bond("reactant", a, b)
            pb = bond("product", a, b)

            ga = gained(a)
            gb = gained(b)

            ga_o = [x for x in ga if x.get("element") == "O"]
            gb_o = [x for x in gb if x.get("element") == "O"]
            ga_n = [x for x in ga if x.get("element") == "N"]
            gb_n = [x for x in gb if x.get("element") == "N"]

            # Alkyne -> alkene hydrogenation.
            #
            # Bond-order reduction alone is insufficient: rearrangements,
            # cyclizations and additions may also contain C#C -> C=C.
            # Require hydrogen gain at both alkyne carbons and preservation
            # of all external mapped heavy-atom neighbors.
            if rb >= 2.8 and 1.7 <= pb <= 2.2:

                def external_heavy_delta(
                    center: int,
                    partner: int,
                ) -> bool:
                    for x in lost(center) + gained(center):
                        el = str(x.get("element") or "")
                        if not el or el == "H":
                            continue

                        m = _att_map(x)

                        # Ignore the mutual alkyne/alkene bond-order
                        # representation itself.
                        if m is not None and int(m) == int(partner):
                            continue

                        return True

                    return False

                r_a = mapped_heavy_neighbors(
                    "reactant",
                    a,
                    exclude_maps={b},
                )
                p_a = mapped_heavy_neighbors(
                    "product",
                    a,
                    exclude_maps={b},
                )
                r_b = mapped_heavy_neighbors(
                    "reactant",
                    b,
                    exclude_maps={a},
                )
                p_b = mapped_heavy_neighbors(
                    "product",
                    b,
                    exclude_maps={a},
                )

                if (
                    hdelta(a) > 0
                    and hdelta(b) > 0
                    and not external_heavy_delta(a, b)
                    and not external_heavy_delta(b, a)
                    and r_a == p_a
                    and r_b == p_b
                ):
                    emit_specific(
                        "partial_alkyne_hydrogenation",
                        a,
                        subtype="alkyne_to_alkene",
                        support=_support(a, b),
                        redox="reduction",
                    )

            # Alkyne hydration -> ketone.
            if rb >= 2.8 and 0.8 <= pb <= 1.2:
                if fg("product", a) == "ketone" or fg("product", b) == "ketone":
                    kc = a if fg("product", a) == "ketone" else b
                    emit_specific(
                        "alkyne_hydration_to_ketone",
                        kc,
                        subtype="alkyne_to_ketone",
                        support=_support(a, b),
                    )

            # Alkene consumption.
            if 1.7 <= rb <= 2.2 and 0.8 <= pb <= 1.2:

                # Wacker-like structural outcome: alkene -> ketone.
                if fg("product", a) == "ketone" or fg("product", b) == "ketone":
                    kc = a if fg("product", a) == "ketone" else b
                    emit_specific(
                        "alkene_to_ketone_oxidation",
                        kc,
                        subtype="alkene_to_ketone",
                        support=_support(a, b),
                        redox="oxidation",
                    )

                # Aminohydroxylation: N and O added across opposite endpoints.
                if (ga_n and gb_o) or (gb_n and ga_o):
                    emit_specific(
                        "alkene_aminohydroxylation",
                        a,
                        subtype="vicinal_aminohydroxylation",
                        support=_support(
                            a,
                            b,
                            *[_att_map(x) for x in ga_n + gb_n + ga_o + gb_o],
                        ),
                        redox="oxidation",
                    )

                # Simple alkene hydration.
                #
                # S2_4_PI_SYSTEM_CLOSURE
                #
                # A single O installation is not automatically hydration:
                # ROH/ArOH addition gives the same local C=C -> C-C + C-O
                # geometry. Require a single product C-O bond and an oxygen
                # that is absent from the reactant or water-like there.
                total_o = ga_o + gb_o
                total_n = ga_n + gb_n

                single_o = [
                    x
                    for x in total_o
                    if 0.8
                    <= float(
                        x.get(
                            "bond_order",
                            x.get(
                                "order",
                                0.0,
                            ),
                        )
                        or 0.0
                    )
                    <= 1.2
                ]

                if len(single_o) == 1 and not total_n:

                    ox = single_o[0]
                    om = _att_map(ox)

                    hydration_oxygen = False

                    if om is None:

                        hydration_oxygen = True

                    else:

                        om = int(om)

                        r_ox = atom(
                            "reactant",
                            om,
                        )

                        if not r_ox.get("element"):

                            # Product-only mapped oxygen.
                            hydration_oxygen = True

                        else:

                            r_ox_heavy = [
                                y
                                for y in atts(
                                    "reactant",
                                    om,
                                )
                                if y.get("element")
                                not in {
                                    "",
                                    "H",
                                }
                            ]

                            # Mapped water / hydroxide-like oxygen,
                            # not an alcohol/phenol/ether oxygen.
                            hydration_oxygen = len(r_ox_heavy) == 0

                    if hydration_oxygen:

                        emit_specific(
                            "alkene_hydration",
                            a if ga_o else b,
                            subtype="alkene_to_alcohol",
                            support=_support(
                                a,
                                b,
                                *[_att_map(x) for x in single_o],
                            ),
                        )

    # ------------------------------------------------------------------
    # 11. Carbocyclic arene partial hydrogenation.
    # ------------------------------------------------------------------
    dearomatized = [
        m
        for m in maps
        if atom("reactant", m).get("element") == "C"
        and atom("reactant", m).get("aromatic")
        and atom("product", m)
        and not atom("product", m).get("aromatic")
    ]

    if len(dearomatized) >= 4:
        # Do not claim this for heteroarenes; existing detector owns those.
        emit_specific(
            "arene_partial_hydrogenation",
            dearomatized[0],
            subtype="partial_dearomatizing_hydrogenation",
            support=_support(*dearomatized),
            redox="reduction",
        )

    # ------------------------------------------------------------------
    # 12. Kabachnik-Fields structural pattern:
    #     carbonyl carbon gains both N and P.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if fg("reactant", c) not in {"aldehyde", "ketone"}:
            continue

        ge = set(elems(gained(c)))
        if {"N", "P"}.issubset(ge):
            emit_specific(
                "aminophosphonate_multicomponent_addition",
                c,
                subtype="carbonyl_aminophosphonate_addition",
                support=_support(c, *[_att_map(x) for x in gained(c)]),
            )

    # ------------------------------------------------------------------
    # 13. Nitroso condensation fallback.
    #
    # Require a genuine nitroso N=O center: exactly one oxygen
    # attachment in the reactant and no oxygen retained on product N.
    # This excludes nitro -> oxime mappings such as Angeli-Rimini.
    # ------------------------------------------------------------------
    for c in sorted(maps):
        if atom("reactant", c).get("element") != "C":
            continue

        for x in atts("product", c):
            if x.get("element") != "N" or float(x.get("bond_order", 0.0)) < 1.8:
                continue

            nm = _att_map(x)

            if nm is None or bond("reactant", c, int(nm)) > 0.1:
                continue

            nm = int(nm)

            r_o = [y for y in atts("reactant", nm) if y.get("element") == "O"]

            p_o = [y for y in atts("product", nm) if y.get("element") == "O"]

            genuine_nitroso = (
                len(r_o) == 1
                and float(r_o[0].get("bond_order", 0.0)) >= 1.8
                and not p_o
            )

            if genuine_nitroso:
                emit_specific(
                    "condensation",
                    c,
                    subtype="nitroso_condensation",
                    support=_support(c, nm),
                )

    # ------------------------------------------------------------------
    # 14. Intercomponent N-N formation fallback.
    #
    # S1_CONFLICT_FIX_V3
    #
    # The S1 recovery is specifically a newly assembled azide N3 chain.
    # A new N-N bond incorporated into an aromatic N-heterocycle must
    # not be interpreted as azide formation.
    # ------------------------------------------------------------------

    def product_n_neighbors(
        nmap: int,
    ) -> set[int]:

        out_n = set()

        for x in atts(
            "product",
            nmap,
        ):
            if x.get("element") != "N":
                continue

            m = _att_map(x)

            if m is not None:
                out_n.add(int(m))

        return out_n

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        a = int(fb["a"])
        b = int(fb["b"])

        if (
            atom("reactant", a).get("element") != "N"
            or atom("reactant", b).get("element") != "N"
        ):
            continue

        # Aromatic N-N construction belongs to heterocycle
        # annulation / aromatic ring formation, not azide formation.
        if atom("product", a).get("aromatic") or atom("product", b).get("aromatic"):
            continue

        nset = {a, b}
        nset.update(product_n_neighbors(a))
        nset.update(product_n_neighbors(b))

        # Exactly three mapped non-aromatic nitrogens must form
        # the local product N3 chain.
        if len(nset) != 3:
            continue

        if any(
            atom("product", n).get("element") != "N"
            or atom("product", n).get("aromatic")
            for n in nset
        ):
            continue

        # Linear N3 connectivity gives N-neighbor degrees 1,2,1.
        n_degrees = sorted(len(product_n_neighbors(n) & nset) for n in nset)

        if n_degrees != [1, 1, 2]:
            continue

        # Retain the existing azide recognition as an independent
        # chemical sanity check.
        if not any(
            _is_azide_n(
                facts,
                "product",
                n,
            )
            for n in nset
        ):
            continue

        emit_specific(
            "n_n_bond_formation",
            a,
            subtype="azide_n_n_formation",
            support=_support(
                *sorted(nset),
            ),
        )

    # ------------------------------------------------------------------
    # 15. Allylic C-H -> C-O installation fallback.
    #
    # S2_5_OWNERSHIP_REFINEMENT:
    # the reacting carbon itself must be a non-aromatic sp3 C-H center.
    # This excludes carbonyl carbons in alpha,beta-unsaturated
    # aldehydes/ketones from masquerading as allylic C-H oxidation.
    # ------------------------------------------------------------------
    for c in sorted(maps):

        ra = atom(
            "reactant",
            c,
        )

        if ra.get("element") != "C":
            continue

        if ra.get("aromatic"):
            continue

        if (
            int(
                ra.get(
                    "total_h",
                    0,
                )
                or 0
            )
            <= 0
        ):
            continue

        # Allylic C-H center must itself be sp3:
        # no multiple heavy-atom bond at the reacting carbon.
        if any(
            (
                x.get("element")
                not in {
                    None,
                    "",
                    "H",
                }
            )
            and float(
                x.get(
                    "bond_order",
                    0.0,
                )
                or 0.0
            )
            >= 1.5
            for x in atts(
                "reactant",
                c,
            )
        ):
            continue

        go = [x for x in gained(c) if x.get("element") == "O"]
        if not go or hdelta(c) >= 0:
            continue

        allylic = False
        for n in atts("reactant", c):
            nm = _att_map(n)
            if nm is None or n.get("element") != "C":
                continue

            for nn in atts("reactant", int(nm)):
                mm = _att_map(nn)
                if mm is None or int(mm) == c or nn.get("element") != "C":
                    continue

                if bond("reactant", int(nm), int(mm)) >= 1.8:
                    allylic = True
                    break

            if allylic:
                break

        if allylic:
            emit_specific(
                "allylic_oxidative_acyloxylation",
                c,
                subtype="allylic_c_o_formation",
                support=_support(c, *[_att_map(x) for x in go]),
                redox="oxidation",
            )

    return out


# ======================================================================
# S2_1_SILICON_TRANSFORMATIONS
#
# Reaction-derived silicon protection / deprotection semantics.
#
# Scope:
#   * O-H  -> O-Si          silylation
#   * O-Si -> O-H           desilylation
#   * N-Si -> N-H           desilylation
#   * C-Si -> C-H           desilylation
#       - terminal alkynyl
#       - aryl
#       - generic carbon
#   * C-Si -> C-X           desilylative halogenation
#   * C-O-Si -> C=O         silyl-ether oxidation to carbonyl
#
# Deliberately NOT handled here:
#   * alkene hydrosilylation
#   * existing C-X -> C-Si substitution/silylation
#   * Fleming-Tamao organosilane -> alcohol
#
# Those already have validated owners elsewhere.
# ======================================================================



def detect_sub05_o_silylation_product_state(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Two structurally verified O-Si product states without O-H reactant loss.

    Only isolated external mapped O-Si construction is eligible.  These
    subtypes do not assert mechanistic Si-transfer pathways, and preserve
    silicon migrations, composite condensations, and conventional silylation.
    """
    formed = list(facts.get("formed_bonds", []) or [])
    if len(formed) != 1 or facts.get("broken_bonds"):
        return []
    fb = formed[0]
    if not 0.8 <= float(fb.get("order", 0)) <= 1.2:
        return []
    a, b = int(fb["a"]), int(fb["b"])
    ea = _atom(facts, "reactant", a).get("element")
    eb = _atom(facts, "reactant", b).get("element")
    if {ea, eb} != {"O", "Si"}:
        return []
    if not any({int(x["a"]), int(x["b"])} == {a, b}
               for x in (facts.get("intercomponent_formed_bonds", []) or [])):
        return []
    o, si = (a, b) if ea == "O" else (b, a)
    if _atom(facts, "reactant", o).get("total_h") != 0 or _atom(facts, "product", o).get("total_h") != 0:
        return []
    if not 0.8 <= _bond_between(facts, "product", o, si) <= 1.2:
        return []
    before_o = _attachments(facts, "reactant", o)
    after_o = _attachments(facts, "product", o)
    changes = facts.get("bond_order_changes", []) or []
    subtype = None
    support = [o, si]
    scope = None

    # Type A: carbamate O-alkyl -> O-silyl group exchange, with retained
    # mapped O-C(=O)-N core and an unmapped departing alkyl carbon.
    if not changes:
        departed_one_sided_c = any(
            x.get("element") == "C" and x.get("map") is None
            and 0.8 <= float(x.get("bond_order", 0)) <= 1.2
            for x in before_o
        )
        no_new_one_sided_c = not any(x.get("element") == "C" and x.get("map") is None for x in after_o)
        if departed_one_sided_c and no_new_one_sided_c:
            for x in before_o:
                c = x.get("map")
                if x.get("element") != "C" or c is None:
                    continue
                c = int(c)
                if not (0.8 <= _bond_between(facts, "reactant", o, c) <= 1.2
                        and 0.8 <= _bond_between(facts, "product", o, c) <= 1.2):
                    continue
                def carbamate_core(side):
                    shell = _attachments(facts, side, c)
                    return (any(z.get("element") == "O" and float(z.get("bond_order", 0)) >= 1.8
                                and z.get("map") != o for z in shell)
                            and any(z.get("element") == "N" and 0.8 <= float(z.get("bond_order", 0)) <= 1.2
                                for z in shell))
                if carbamate_core("reactant") and carbamate_core("product"):
                    subtype = "silyloxy_carbamate_o_alkyl_exchange"
                    support = [o, si, c]
                    scope = "mapped_carbamate_o_c_to_o_si_exchange_with_one_sided_alkyl_loss"
                    break
    # Type B: heteroaryl oxo tautomer -> silylated lactim, with adjacent
    # mapped aromatic N-H loss and exactly one C=O -> C-O order change.
    elif len(changes) == 1:
        ch = changes[0]
        x, y = int(ch["a"]), int(ch["b"])
        c = y if x == o else x if y == o else None
        if (c is not None and _atom(facts, "reactant", c).get("element") == "C"
                and bool(_atom(facts, "reactant", c).get("aromatic"))
                and bool(_atom(facts, "product", c).get("aromatic"))
                and float(ch.get("reactant_order", 0)) >= 1.8
                and 0.8 <= float(ch.get("product_order", 0)) <= 1.2):
            for n_att in _attachments(facts, "reactant", c):
                n = n_att.get("map")
                if n is None or n_att.get("element") != "N":
                    continue
                n = int(n)
                nr, np = _atom(facts, "reactant", n), _atom(facts, "product", n)
                if (nr.get("aromatic") and np.get("aromatic")
                        and int(nr.get("total_h", 0)) == 1 and int(np.get("total_h", 0)) == 0
                        and _bond_between(facts, "product", c, n) > 0):
                    subtype = "oxo_heteroaryl_tautomer_o_silylation"
                    support = [o, si, c, n]
                    scope = "mapped_oxo_heteroaryl_c_o_to_siloxy_with_adjacent_n_h_loss"
                    break
    if subtype is None:
        return []
    return [{
        "event": "silylation", "role": "specific_event", "subtype": subtype,
        "center_map": o, "support_maps": sorted(set(support)), "mapping_scope": scope,
    }]

def detect_silicon_transformations(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[str, str, int | None]] = set()

    HALOGENS = {"F", "Cl", "Br", "I"}

    def total_h(side: str, mn: int) -> int:
        try:
            return int(
                _atom(
                    facts,
                    side,
                    mn,
                ).get("total_h", 0)
                or 0
            )
        except Exception:
            return 0

    def hydrogen_delta(mn: int) -> int:
        return total_h("product", mn) - total_h("reactant", mn)

    def add(
        event: str,
        subtype: str,
        center: int | None,
        support: list[int],
        mapping_scope: str,
        *,
        redox_direction: str | None = None,
        extra: dict[str, Any] | None = None,
    ) -> None:

        key = (
            event,
            subtype,
            center,
        )

        if key in seen:
            return

        seen.add(key)

        rec: dict[str, Any] = {
            "event": event,
            "role": "specific_event",
            "subtype": subtype,
            "center_map": center,
            "support_maps": sorted({int(x) for x in support if x is not None}),
            "mapping_scope": mapping_scope,
        }

        if redox_direction:
            rec["redox_direction"] = redox_direction

        if extra:
            rec.update(extra)

        out.append(rec)

    def has_triple_c_neighbor(
        side: str,
        center: int,
    ) -> bool:

        for x in _attachments(
            facts,
            side,
            center,
        ):
            if (
                x.get("element") == "C"
                and float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                >= 2.8
            ):
                return True

        return False

    # ------------------------------------------------------------------
    # S2_1B_INTERNAL_SILICON_MIGRATION_GUARD
    #
    # A mapped Si atom that is lost from one conserved center and gained
    # by another conserved center is an internal silicon migration, not
    # an independent desilylation + silylation pair.
    #
    # Brook-type rearrangements are the canonical example.
    # ------------------------------------------------------------------

    si_lost_maps: set[int] = set()
    si_gained_maps: set[int] = set()

    for d0 in facts.get("attachment_deltas", []) or []:

        for x in d0.get("lost", []) or []:
            if x.get("element") != "Si":
                continue
            mn = _att_map(x)
            if mn is not None:
                si_lost_maps.add(int(mn))

        for x in d0.get("gained", []) or []:
            if x.get("element") != "Si":
                continue
            mn = _att_map(x)
            if mn is not None:
                si_gained_maps.add(int(mn))

    migrating_si_maps = si_lost_maps & si_gained_maps

    # ------------------------------------------------------------------
    # A. Attachment-delta silicon chemistry
    # ------------------------------------------------------------------

    for d in (
        facts.get(
            "attachment_deltas",
            [],
        )
        or []
    ):

        if d.get("center_map") is None:
            continue

        center = int(d["center_map"])

        react_atom = _atom(
            facts,
            "reactant",
            center,
        )

        prod_atom = _atom(
            facts,
            "product",
            center,
        )

        center_el = react_atom.get("element") or prod_atom.get("element")

        if center_el not in {
            "C",
            "N",
            "O",
        }:
            continue

        lost = d.get("lost", []) or []

        gained = d.get("gained", []) or []

        lost_si = [
            x
            for x in lost
            if (
                x.get("element") == "Si"
                and (_att_map(x) is None or int(_att_map(x)) not in migrating_si_maps)
            )
        ]

        gained_si = [
            x
            for x in gained
            if (
                x.get("element") == "Si"
                and (_att_map(x) is None or int(_att_map(x)) not in migrating_si_maps)
            )
        ]

        lost_si_maps = [_att_map(x) for x in lost_si]

        gained_si_maps = [_att_map(x) for x in gained_si]

        hdelta = hydrogen_delta(center)

        gained_heavy_non_si = [
            x
            for x in gained
            if (
                x.get("element")
                and x.get("element")
                not in {
                    "H",
                    "Si",
                }
            )
        ]

        # --------------------------------------------------------------
        # O-H -> O-Si
        #
        # Require:
        #   conserved oxygen
        #   new Si attachment
        #   hydrogen loss / alcohol-like reactant O
        #
        # This is deliberately narrower than generic O-Si bond formation.
        # --------------------------------------------------------------

        if center_el == "O" and gained_si and not lost_si:

            react_h = total_h(
                "reactant",
                center,
            )

            prod_h = total_h(
                "product",
                center,
            )

            if hdelta < 0 or (react_h >= 1 and prod_h == 0):
                add(
                    "silylation",
                    "O_H_to_Si",
                    center,
                    _support(
                        center,
                        *gained_si_maps,
                    ),
                    "conserved_oxygen_loses_hydrogen_and_gains_silicon_attachment",
                )

        # --------------------------------------------------------------
        # S2_2D2_N_SILYLATION
        # N-H -> N-Si
        #
        # Directional sibling of the validated O-H -> O-Si branch.
        #
        # Require:
        #   conserved nitrogen
        #   gain of Si
        #   no loss of Si
        #   no simultaneous non-Si heavy-atom gain
        #   explicit hydrogen loss / N-H reactant state
        #
        # Internal mapped Si migration has already been excluded above
        # when gained_si was constructed.
        # --------------------------------------------------------------

        if center_el == "N" and gained_si and not lost_si and not gained_heavy_non_si:

            react_h = total_h(
                "reactant",
                center,
            )

            prod_h = total_h(
                "product",
                center,
            )

            if hdelta < 0 or (react_h >= 1 and prod_h == 0):
                add(
                    "silylation",
                    "N_H_to_Si",
                    center,
                    _support(
                        center,
                        *gained_si_maps,
                    ),
                    "conserved_nitrogen_loses_hydrogen_and_gains_silicon_attachment",
                )

        # --------------------------------------------------------------
        # O-Si / N-Si -> O-H / N-H
        #
        # Require:
        #   loss of Si
        #   no new heavy-atom replacement
        #   hydrogen gain
        #
        # This prevents ordinary substitution from being mislabeled as
        # deprotection.
        # --------------------------------------------------------------

        if (
            center_el in {"O", "N"}
            and lost_si
            and not gained_heavy_non_si
            and hdelta > 0
        ):

            subtype = "o_desilylation" if center_el == "O" else "n_desilylation"

            add(
                "desilylation",
                subtype,
                center,
                _support(
                    center,
                    *lost_si_maps,
                ),
                ("conserved_heteroatom_loses_silicon_" "attachment_and_gains_hydrogen"),
            )

        # --------------------------------------------------------------
        # C-Si -> C-X
        #
        # Desilylative halogenation is chemically distinct from simple
        # desilylation and must not be collapsed into C-Si -> C-H.
        # --------------------------------------------------------------

        if center_el == "C" and lost_si:

            gained_halogen = [x for x in gained if x.get("element") in HALOGENS]

            if gained_halogen:

                halogens = sorted({str(x.get("element")) for x in gained_halogen})

                halogen = halogens[0] if len(halogens) == 1 else "halogen"

                add(
                    "desilylative_halogenation",
                    f"c_si_to_c_{halogen.lower()}",
                    center,
                    _support(
                        center,
                        *lost_si_maps,
                        *[_att_map(x) for x in gained_halogen],
                    ),
                    "conserved_carbon_loses_silicon_and_gains_halogen",
                )

                # Do not also call this simple desilylation.
                continue

            # ----------------------------------------------------------
            # C-Si -> C-H
            #
            # Require hydrogen gain and no heavy-atom replacement.
            # ----------------------------------------------------------

            if hdelta > 0 and not gained_heavy_non_si:

                aromatic = bool(react_atom.get("aromatic") or prod_atom.get("aromatic"))

                alkynyl = has_triple_c_neighbor(
                    "reactant",
                    center,
                ) or has_triple_c_neighbor(
                    "product",
                    center,
                )

                if alkynyl:
                    subtype = "terminal_alkynyl_desilylation"
                elif aromatic:
                    subtype = "aryl_desilylation"
                else:
                    subtype = "c_desilylation"

                add(
                    "desilylation",
                    subtype,
                    center,
                    _support(
                        center,
                        *lost_si_maps,
                    ),
                    "conserved_carbon_loses_silicon_attachment_and_gains_hydrogen",
                )

    # ------------------------------------------------------------------
    # B. Silyl ether -> carbonyl
    #
    # Pattern:
    #
    #       C-O-Si        C=O
    #          \    ->     |
    #
    # The conserved oxygen loses its Si attachment while the conserved
    # C-O bond changes single -> double.
    #
    # This is kept separate from generic desilylation because the
    # defining transformation is oxidation to carbonyl.
    # ------------------------------------------------------------------

    oxygen_losing_si: dict[int, list[int | None]] = {}

    for d in (
        facts.get(
            "attachment_deltas",
            [],
        )
        or []
    ):

        if d.get("center_map") is None:
            continue

        o = int(d["center_map"])

        if (
            _atom(
                facts,
                "reactant",
                o,
            ).get("element")
            != "O"
        ):
            continue

        lost_si = [x for x in (d.get("lost", []) or []) if x.get("element") == "Si"]

        if lost_si:
            oxygen_losing_si[o] = [_att_map(x) for x in lost_si]

    for bc in (
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    ):

        if bc.get("a") is None or bc.get("b") is None:
            continue

        a = int(bc["a"])
        b = int(bc["b"])

        rorder = float(
            bc.get(
                "reactant_order",
                0.0,
            )
            or 0.0
        )

        porder = float(
            bc.get(
                "product_order",
                0.0,
            )
            or 0.0
        )

        if not (rorder <= 1.2 and porder >= 1.8):
            continue

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element") or _atom(
            facts,
            "product",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element") or _atom(
            facts,
            "product",
            b,
        ).get("element")

        if {ea, eb} != {"C", "O"}:
            continue

        carbon = a if ea == "C" else b

        oxygen = a if ea == "O" else b

        if oxygen not in oxygen_losing_si:
            continue

        # Require an enol-silane reactant state rather than a generic
        # silyl ether: the carbon becoming carbonyl must have had a
        # C=C bond in the reactant.
        # S2_1C_MUKAIYAMA_ENOL_PARTNER
        #
        # Preserve the mapped carbon partner of the original enol C=C.
        # The resolver can then determine whether carbonyl restoration is
        # subordinate to a Mukaiyama-type carbonyl addition rather than an
        # independent defining event.
        reactant_alkene_partners = [
            x
            for x in _attachments(
                facts,
                "reactant",
                carbon,
            )
            if (
                x.get("element") == "C"
                and float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                >= 1.8
            )
        ]

        # Presence of the enol C=C is the detector requirement.
        # Mapping of the partner is optional.
        #
        # A mapped partner is used only by the resolver when deciding
        # whether the silicon-carbonyl event is subordinate to a
        # Mukaiyama-type carbonyl addition.
        # S2_1D_SILYLATED_ALCOHOL_TO_CARBONYL
        #
        # Two distinct structural families share O-Si loss + C=O
        # formation:
        #
        #   1. C=C-O-Si -> C=O
        #      silyl enol ether carbonyl restoration
        #
        #   2. saturated C(H)-O-Si -> C=O
        #      oxidation of a silylated alcohol
        #
        # Keep them as sibling semantic events rather than merging
        # both into a generic silyl-ether conversion.

        reactant_alkene_partner_maps = sorted(
            {
                int(_att_map(x))
                for x in reactant_alkene_partners
                if _att_map(x) is not None
            }
        )

        if reactant_alkene_partners:

            add(
                "silyl_enol_ether_to_carbonyl",
                "silyl_enol_ether_to_carbonyl",
                carbon,
                _support(
                    carbon,
                    oxygen,
                    *oxygen_losing_si[oxygen],
                ),
                (
                    "mapped_silyl_enol_ether_restores_carbonyl_"
                    "while_oxygen_loses_silicon"
                ),
                extra={
                    "enol_partner_maps": reactant_alkene_partner_maps,
                },
            )

        elif hydrogen_delta(carbon) < 0:

            add(
                "silyl_ether_to_carbonyl",
                "silylated_alcohol_to_carbonyl",
                carbon,
                _support(
                    carbon,
                    oxygen,
                    *oxygen_losing_si[oxygen],
                ),
                (
                    "mapped_saturated_c_o_si_center_becomes_"
                    "carbonyl_with_carbon_hydrogen_loss"
                ),
                redox_direction="oxidation",
            )

    return out


# =============================================================================
# S2_3B_FORMAMIDINE_FORMATION
#
# Structural contract:
#
#   formamide-acetal-like C(H)(OR)2-N
#                +
#             incoming N-H
#
#                  ->
#
#             N-C(H)=N
#
# No reagent-name dependence and no template dependence.
# =============================================================================


def detect_s23b_formamidine_formation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        # Exact new C=N-like double bond.
        bo = float(
            fb.get(
                "order",
                0.0,
            )
        )

        if not (1.8 <= bo <= 2.2):
            continue

        a = int(fb["a"])

        b = int(fb["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element")

        if {
            ea,
            eb,
        } != {
            "C",
            "N",
        }:
            continue

        carbon, incoming_n = (a, b) if ea == "C" else (b, a)

        # ---------------------------------------------------------
        # Formamidine rather than generic amidine:
        # donor center is C(H).
        # ---------------------------------------------------------

        carbon_atom = _atom(
            facts,
            "reactant",
            carbon,
        )

        if (
            int(
                carbon_atom.get(
                    "total_h",
                    0,
                )
            )
            < 1
        ):
            continue

        # ---------------------------------------------------------
        # Reactant donor state:
        #
        #            OR
        #            |
        #      N --- C(H) --- OR
        #
        # one retained mapped N + at least two single O attachments
        # ---------------------------------------------------------

        r_atts = (
            _attachments(
                facts,
                "reactant",
                carbon,
            )
            or []
        )

        donor_n = [
            x
            for x in r_atts
            if (
                x.get("element") == "N"
                and _att_map(x) is not None
                and 0.8
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= 1.2
            )
        ]

        r_o_single = [
            x
            for x in r_atts
            if (
                x.get("element") == "O"
                and 0.8
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= 1.2
            )
        ]

        if len(donor_n) != 1:
            continue

        if len(r_o_single) < 2:
            continue

        donor_n_map = int(_att_map(donor_n[0]))

        # ---------------------------------------------------------
        # Carbon attachment delta:
        # both O substituents disappear and are unmapped.
        # ---------------------------------------------------------

        deltas = [
            d
            for d in (
                facts.get(
                    "attachment_deltas",
                    [],
                )
                or []
            )
            if int(d["center_map"]) == carbon
        ]

        if len(deltas) != 1:
            continue

        delta = deltas[0]

        lost_o_unmapped = [
            x
            for x in (
                delta.get(
                    "lost",
                    [],
                )
                or []
            )
            if (
                x.get("element") == "O"
                and x.get("map") is None
                and 0.8
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                )
                <= 1.2
            )
        ]

        if len(lost_o_unmapped) < 2:
            continue

        # ---------------------------------------------------------
        # Product:
        #
        # donor N remains single-bonded;
        # incoming N is double-bonded;
        # no oxygen remains attached to the center.
        # ---------------------------------------------------------

        p_atts = (
            _attachments(
                facts,
                "product",
                carbon,
            )
            or []
        )

        incoming_double = any(
            x.get("element") == "N"
            and _att_map(x) is not None
            and int(_att_map(x)) == incoming_n
            and 1.8
            <= float(
                x.get(
                    "bond_order",
                    0.0,
                )
            )
            <= 2.2
            for x in p_atts
        )

        donor_retained = any(
            x.get("element") == "N"
            and _att_map(x) is not None
            and int(_att_map(x)) == donor_n_map
            and 0.8
            <= float(
                x.get(
                    "bond_order",
                    0.0,
                )
            )
            <= 1.2
            for x in p_atts
        )

        product_o = [x for x in p_atts if x.get("element") == "O"]

        if not incoming_double:
            continue

        if not donor_retained:
            continue

        if product_o:
            continue

        # Incoming nitrogen loses hydrogen during imine formation.
        if (
            _hydrogen_delta(
                facts,
                incoming_n,
            )
            > -1
        ):
            continue

        out.append(
            {
                "event": "formamidine_formation",
                "role": "specific_event",
                "subtype": "formamide_acetal_to_formamidine",
                "center_map": carbon,
                "incoming_n_map": incoming_n,
                "donor_n_map": donor_n_map,
                "support_maps": _support(
                    carbon,
                    incoming_n,
                    donor_n_map,
                ),
                "mapping_scope": (
                    "mapped_formamide_acetal_carbon_"
                    "loses_two_unmapped_o_substituents_"
                    "and_forms_intercomponent_c_n_double_bond"
                ),
            }
        )

    return out


# =============================================================================
# S2_3E1_PRODUCT_CENTERED_N_METHYLATION
#
# Representation-boundary extension of the existing:
#
#     n_alkylation | n_methylation
#
# Problem addressed:
# Atom mapping may assign a carboxylic-acid carbon to the product methyl carbon
# of an N-methylated amine.  The existing generic substitution detector requires
# an sp3 reactant carbon and therefore correctly does not fire on that mapped
# reactant representation.
#
# Semantic ownership is product-centered:
#
#     mapped N-H + mapped carbon  ->  N-CH3
#
# The reactant carboxylic-acid state is used only to isolate the known
# representation boundary.  It is NOT interpreted mechanistically as an acid
# being the physical methylating reagent.
#
# Conditions/reagents are deliberately excluded from this structural detector.
# =============================================================================


def detect_s23e1_product_centered_n_methylation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        if not (
            0.8
            <= float(
                fb.get(
                    "order",
                    0.0,
                )
            )
            <= 1.2
        ):
            continue

        a = int(fb["a"])

        b = int(fb["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element")

        if {
            ea,
            eb,
        } != {
            "C",
            "N",
        }:
            continue

        if ea == "C":
            carbon = a
            nitrogen = b
        else:
            carbon = b
            nitrogen = a

        # Representation-boundary trigger:
        # mapped source carbon appears as a carboxylic-acid carbon.
        if (
            _carbon_fg_state(
                facts,
                "reactant",
                carbon,
            )
            != "carboxylic_acid"
        ):
            continue

        # Product carbon must be an exact terminal methyl attached to N.
        pcarbon = _atom(
            facts,
            "product",
            carbon,
        )

        if pcarbon.get("element") != "C":
            continue

        if (
            int(
                pcarbon.get(
                    "total_h",
                    0,
                )
                or 0
            )
            != 3
        ):
            continue

        patts = (
            _attachments(
                facts,
                "product",
                carbon,
            )
            or []
        )

        # Carbon has exactly one heavy-atom attachment.
        if len(patts) != 1:
            continue

        patt = patts[0]

        if patt.get("element") != "N":
            continue

        if not (
            0.8
            <= float(
                patt.get(
                    "bond_order",
                    0.0,
                )
            )
            <= 1.2
        ):
            continue

        nmap = _att_map(patt)

        if nmap is None or int(nmap) != nitrogen:
            continue

        # N-H → N-C methylation evidence.
        if (
            _hydrogen_delta(
                facts,
                nitrogen,
            )
            > -1
        ):
            continue

        key = (
            carbon,
            nitrogen,
        )

        if key in seen:
            continue

        seen.add(key)

        out.append(
            {
                "event": "n_alkylation",
                "role": "specific_event",
                "subtype": "n_methylation",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "support_maps": _support(
                    carbon,
                    nitrogen,
                ),
                "mapping_scope": (
                    "product_centered_n_methylation_with_"
                    "carboxylic_acid_mapped_carbon_"
                    "representation_boundary"
                ),
                "representation_boundary": (
                    "reactant_carboxylic_acid_carbon_"
                    "provenance_not_interpreted_mechanistically"
                ),
            }
        )

    return out


# ======================================================================
# S2_3E2_CARBOXYLIC_ACID_TO_AMINOMETHYLENE
#
# Formal mapped transformation:
#
#     R-C(=O)OH  ->  R-CH2-N
#
# Contract:
#   * formed intercomponent C-N single bond
#   * mapped carbon is carboxylic-acid carbon in reactant
#   * exactly one mapped carbon substituent retained
#   * product center is exact R-CH2-N
#   * carbon H delta = +2
#   * nitrogen H delta = -1
#
# This is a formal graph-semantic transformation. It does not claim
# that the carboxylic acid itself is mechanistically the alkylating
# reagent.
# ======================================================================


def detect_s23e2_carboxylic_acid_to_aminomethylene(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    def bond_order(x: dict[str, Any]) -> float:

        return float(
            x.get(
                "order",
                x.get(
                    "bond_order",
                    0.0,
                ),
            )
            or 0.0
        )

    def mapped_element_neighbors(
        side: str,
        center: int,
        element: str,
    ) -> set[int]:

        vals: set[int] = set()

        for att in (
            _attachments(
                facts,
                side,
                center,
            )
            or []
        ):
            if att.get("element") != element:
                continue

            m = _att_map(att)

            if m is not None:
                vals.add(int(m))

        return vals

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        order = float(
            fb.get(
                "order",
                fb.get(
                    "bond_order",
                    0.0,
                ),
            )
            or 0.0
        )

        if not (0.8 <= order <= 1.2):
            continue

        a = int(fb["a"])
        b = int(fb["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element")

        if {
            ea,
            eb,
        } != {
            "C",
            "N",
        }:
            continue

        carbon = a if ea == "C" else b

        nitrogen = b if ea == "C" else a

        # --------------------------------------------------------------
        # Reactant center = substituted carboxylic acid
        # --------------------------------------------------------------

        if (
            _carbon_fg_state(
                facts,
                "reactant",
                carbon,
            )
            != "carboxylic_acid"
        ):
            continue

        ra = (
            _attachments(
                facts,
                "reactant",
                carbon,
            )
            or []
        )

        r_o_single = [
            x for x in ra if (x.get("element") == "O" and 0.8 <= bond_order(x) <= 1.2)
        ]

        r_o_double = [
            x for x in ra if (x.get("element") == "O" and 1.8 <= bond_order(x) <= 2.2)
        ]

        r_c_single = [
            x
            for x in ra
            if (
                x.get("element") == "C"
                and 0.8 <= bond_order(x) <= 1.2
                and _att_map(x) is not None
            )
        ]

        if not (len(r_o_single) == 1 and len(r_o_double) == 1 and len(r_c_single) == 1):
            continue

        # --------------------------------------------------------------
        # Product center = exact R-CH2-N
        # --------------------------------------------------------------

        pa = (
            _attachments(
                facts,
                "product",
                carbon,
            )
            or []
        )

        p_c_single = [
            x
            for x in pa
            if (
                x.get("element") == "C"
                and 0.8 <= bond_order(x) <= 1.2
                and _att_map(x) is not None
            )
        ]

        p_n_single = [
            x
            for x in pa
            if (
                x.get("element") == "N"
                and 0.8 <= bond_order(x) <= 1.2
                and _att_map(x) is not None
            )
        ]

        p_other = [
            x
            for x in pa
            if not (
                (
                    x.get("element") == "C"
                    and 0.8 <= bond_order(x) <= 1.2
                    and _att_map(x) is not None
                )
                or (
                    x.get("element") == "N"
                    and 0.8 <= bond_order(x) <= 1.2
                    and _att_map(x) is not None
                )
            )
        ]

        if not (len(p_c_single) == 1 and len(p_n_single) == 1 and len(p_other) == 0):
            continue

        if int(_att_map(p_n_single[0])) != nitrogen:
            continue

        # --------------------------------------------------------------
        # Retained carbon skeleton
        # --------------------------------------------------------------

        r_c_maps = mapped_element_neighbors(
            "reactant",
            carbon,
            "C",
        )

        p_c_maps = mapped_element_neighbors(
            "product",
            carbon,
            "C",
        )

        if not (len(r_c_maps) == 1 and r_c_maps == p_c_maps):
            continue

        retained_carbon = next(iter(r_c_maps))

        # --------------------------------------------------------------
        # Hydrogen geometry
        # --------------------------------------------------------------

        if (
            _hydrogen_delta(
                facts,
                carbon,
            )
            != 2
        ):
            continue

        if (
            _hydrogen_delta(
                facts,
                nitrogen,
            )
            != -1
        ):
            continue

        product_carbon = _atom(
            facts,
            "product",
            carbon,
        )

        if (
            int(
                product_carbon.get(
                    "total_h",
                    0,
                )
                or 0
            )
            != 2
        ):
            continue

        key = (
            carbon,
            nitrogen,
        )

        if key in seen:
            continue

        seen.add(key)

        out.append(
            {
                "event": "n_alkylation",
                "role": "specific_event",
                "subtype": "carboxylic_acid_to_aminomethylene",
                "center_map": carbon,
                "nitrogen_map": nitrogen,
                "retained_carbon_map": retained_carbon,
                "support_maps": _support(
                    carbon,
                    nitrogen,
                    retained_carbon,
                ),
                "redox_direction": "reduction",
                "mapping_scope": (
                    "mapped_carboxylic_acid_carbon_to_"
                    "retained_skeleton_aminomethylene"
                ),
            }
        )

    return out


# ======================================================================
# S2_3E3_PRODUCT_CENTERED_UREA_FORMATION
#
# Product-centered formal urea formation.
#
# Structural contract:
#
#   reactant mapped center:
#       C-O(mapped) single
#       optionally one unmapped C single attachment
#       no other heavy attachment
#
#   product:
#       exact O=C(N)N center
#
#   both C-N bonds:
#       newly formed intercomponent single bonds
#
#   both incoming N atoms:
#       hydrogen_delta == -1
#
#   retained mapped oxygen:
#       C-O single -> C=O
#
#   product carbon:
#       total_h == 0
#
#   carbon hydrogen delta:
#       -3 or -2
#
# This is deliberately a formal product-centered semantic event.
# It makes no mechanistic claim about the precursor representation.
# ======================================================================


def detect_s23e3_product_centered_urea_formation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []

    def bond_order(x: dict[str, Any]) -> float:

        return float(
            x.get(
                "order",
                x.get(
                    "bond_order",
                    0.0,
                ),
            )
            or 0.0
        )

    # --------------------------------------------------------------
    # Group newly formed intercomponent C-N single bonds
    # by mapped carbon center.
    # --------------------------------------------------------------

    formed_by_carbon: dict[int, set[int]] = {}

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        order = float(
            fb.get(
                "order",
                fb.get(
                    "bond_order",
                    0.0,
                ),
            )
            or 0.0
        )

        if not (0.8 <= order <= 1.2):
            continue

        a = int(fb["a"])

        b = int(fb["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element")

        if {
            ea,
            eb,
        } != {
            "C",
            "N",
        }:
            continue

        carbon = a if ea == "C" else b

        nitrogen = b if ea == "C" else a

        formed_by_carbon.setdefault(
            carbon,
            set(),
        ).add(nitrogen)

    # --------------------------------------------------------------
    # Exact E3 geometry
    # --------------------------------------------------------------

    for carbon, formed_ns in formed_by_carbon.items():

        if len(formed_ns) != 2:
            continue

        formed_ns = set(int(x) for x in formed_ns)

        if any(
            _hydrogen_delta(
                facts,
                n,
            )
            != -1
            for n in formed_ns
        ):
            continue

        # ----------------------------------------------------------
        # Product must be exact O=C(N)N
        # ----------------------------------------------------------

        pa = (
            _attachments(
                facts,
                "product",
                carbon,
            )
            or []
        )

        p_o_double = [
            x
            for x in pa
            if (
                x.get("element") == "O"
                and 1.8 <= bond_order(x) <= 2.2
                and _att_map(x) is not None
            )
        ]

        p_n_single = [
            x
            for x in pa
            if (
                x.get("element") == "N"
                and 0.8 <= bond_order(x) <= 1.2
                and _att_map(x) is not None
            )
        ]

        allowed_product_ids = {id(x) for x in (p_o_double + p_n_single)}

        p_other = [x for x in pa if id(x) not in allowed_product_ids]

        if not (len(p_o_double) == 1 and len(p_n_single) == 2 and len(p_other) == 0):
            continue

        product_o = int(_att_map(p_o_double[0]))

        product_ns = {int(_att_map(x)) for x in p_n_single}

        if product_ns != formed_ns:
            continue

        if (
            int(
                _atom(
                    facts,
                    "product",
                    carbon,
                ).get(
                    "total_h",
                    0,
                )
                or 0
            )
            != 0
        ):
            continue

        # ----------------------------------------------------------
        # Reactant mapped center
        # ----------------------------------------------------------

        ra = (
            _attachments(
                facts,
                "reactant",
                carbon,
            )
            or []
        )

        r_o_single = [
            x
            for x in ra
            if (
                x.get("element") == "O"
                and 0.8 <= bond_order(x) <= 1.2
                and _att_map(x) is not None
            )
        ]

        r_unmapped_c_single = [
            x
            for x in ra
            if (
                x.get("element") == "C"
                and 0.8 <= bond_order(x) <= 1.2
                and _att_map(x) is None
            )
        ]

        allowed_reactant_ids = {id(x) for x in (r_o_single + r_unmapped_c_single)}

        r_other = [x for x in ra if id(x) not in allowed_reactant_ids]

        if not (
            len(r_o_single) == 1 and len(r_unmapped_c_single) <= 1 and len(r_other) == 0
        ):
            continue

        if int(_att_map(r_o_single[0])) != product_o:
            continue

        carbon_hd = _hydrogen_delta(
            facts,
            carbon,
        )

        if carbon_hd not in {
            -3,
            -2,
        }:
            continue

        ns = sorted(formed_ns)

        out.append(
            {
                "event": "urea_formation",
                "role": "specific_event",
                "subtype": "product_centered_urea_formation",
                "center_map": carbon,
                "oxygen_map": product_o,
                "nitrogen_maps": ns,
                "support_maps": _support(
                    carbon,
                    product_o,
                    *ns,
                ),
                "mapping_scope": (
                    "mapped_product_centered_urea_" "formation_with_retained_oxygen"
                ),
            }
        )

    return out


# ======================================================================
# S2_4_PI_SYSTEM_CLOSURE
# ======================================================================


def detect_s24_pi_system_closure(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []

    def bo(rec: dict[str, Any]) -> float:

        return float(
            rec.get(
                "order",
                rec.get(
                    "bond_order",
                    0.0,
                ),
            )
            or 0.0
        )

    def external_mapped_heavy_neighbors(
        side: str,
        center: int,
        support: set[int],
    ) -> set[tuple[int, str, float]]:

        ans = set()

        for x in (
            _attachments(
                facts,
                side,
                center,
            )
            or []
        ):

            el = str(x.get("element") or "")

            if not el or el == "H":
                continue

            m = _att_map(x)

            if m is None:
                continue

            m = int(m)

            if m in support:
                continue

            ans.add(
                (
                    m,
                    el,
                    round(
                        float(
                            x.get(
                                "bond_order",
                                x.get(
                                    "order",
                                    0.0,
                                ),
                            )
                            or 0.0
                        ),
                        2,
                    ),
                )
            )

        return ans

    # ------------------------------------------------------------------
    # A. Formal H-X addition to C=C / C#C.
    # ------------------------------------------------------------------

    for bc in (
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    ):

        a = int(bc["a"])

        b = int(bc["b"])

        if (
            _atom(
                facts,
                "reactant",
                a,
            ).get("element")
            != "C"
            or _atom(
                facts,
                "reactant",
                b,
            ).get("element")
            != "C"
        ):
            continue

        ro = float(
            bc.get(
                "reactant_order",
                0.0,
            )
            or 0.0
        )

        po = float(
            bc.get(
                "product_order",
                0.0,
            )
            or 0.0
        )

        if 1.8 <= ro <= 2.2 and po <= 1.2:

            pi_kind = "alkene"

        elif ro >= 2.8 and 1.8 <= po <= 2.2:

            pi_kind = "alkyne"

        else:

            continue

        # All newly formed heavy bonds touching either endpoint.
        touching = []

        for fb in (
            facts.get(
                "formed_bonds",
                [],
            )
            or []
        ):

            x = int(fb["a"])

            y = int(fb["b"])

            if a not in {
                x,
                y,
            } and b not in {
                x,
                y,
            }:
                continue

            touching.append(fb)

        # Formal H-X addition must install exactly one heavy atom
        # at the pi pair. Multi-additions are deliberately deferred.
        if len(touching) != 1:

            continue

        fb = touching[0]

        if not (0.8 <= bo(fb) <= 1.2):
            continue

        x = int(fb["a"])

        y = int(fb["b"])

        if x in {
            a,
            b,
        }:

            carbon = x
            hetero = y

        elif y in {
            a,
            b,
        }:

            carbon = y
            hetero = x

        else:

            continue

        if carbon not in {
            a,
            b,
        }:

            continue

        element = str(
            _atom(
                facts,
                "reactant",
                hetero,
            ).get("element")
            or ""
        )

        if element not in {
            "O",
            "N",
            "S",
        }:

            continue

        # Must be genuine intercomponent incorporation.
        is_intercomponent = any(
            {
                int(z["a"]),
                int(z["b"]),
            }
            == {
                carbon,
                hetero,
            }
            for z in (
                facts.get(
                    "intercomponent_formed_bonds",
                    [],
                )
                or []
            )
        )

        if not is_intercomponent:
            continue

        # Donor X-H -> X-C.
        r_hetero = _atom(
            facts,
            "reactant",
            hetero,
        )

        if (
            int(
                r_hetero.get(
                    "total_h",
                    0,
                )
                or 0
            )
            < 1
        ):

            continue

        if (
            _hydrogen_delta(
                facts,
                hetero,
            )
            >= 0
        ):

            continue

        partner = b if carbon == a else a

        # Complementary formal hydrogen addition.
        if (
            _hydrogen_delta(
                facts,
                partner,
            )
            <= 0
        ):

            continue

        # Water-like oxygen remains the existing alkene-hydration
        # ontology rather than becoming generic O-H hydroaddition.
        if element == "O":

            r_heavy = [
                z
                for z in (
                    _attachments(
                        facts,
                        "reactant",
                        hetero,
                    )
                    or []
                )
                if z.get("element")
                not in {
                    "",
                    "H",
                }
            ]

            if not r_heavy:
                continue

        # Existing activated-alkene N-addition has a more specific
        # aza-Michael owner.
        if (
            element == "N"
            and pi_kind == "alkene"
            and _alkene_activated(
                facts,
                a,
                b,
            )
        ):
            continue

        subtype = f"{pi_kind}_" f"{element.lower()}_h_addition"

        out.append(
            {
                "event": "pi_heteroatom_hydroaddition",
                "role": "specific_event",
                "subtype": subtype,
                "center_map": carbon,
                "pi_maps": sorted(
                    [
                        a,
                        b,
                    ]
                ),
                "heteroatom_map": hetero,
                "partner_map": partner,
                "heteroatom_element": element,
                "support_maps": _support(
                    a,
                    b,
                    hetero,
                ),
                "mapping_scope": (
                    "mapped_c_c_pi_bond_reduction_"
                    "plus_intercomponent_c_x_single_"
                    "formation_plus_x_h_loss_and_"
                    "partner_h_gain"
                ),
            }
        )

    # ------------------------------------------------------------------
    # B. Clean adjacent mapped double-bond migration.
    # ------------------------------------------------------------------

    changes = list(
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    )

    if (
        len(changes) == 2
        and not (
            facts.get(
                "formed_bonds",
                [],
            )
            or []
        )
        and not (
            facts.get(
                "broken_bonds",
                [],
            )
            or []
        )
    ):

        inc = [
            z
            for z in changes
            if (
                float(
                    z.get(
                        "reactant_order",
                        0.0,
                    )
                    or 0.0
                )
                <= 1.2
                and 1.8
                <= float(
                    z.get(
                        "product_order",
                        0.0,
                    )
                    or 0.0
                )
                < 2.8
            )
        ]

        dec = [
            z
            for z in changes
            if (
                1.8
                <= float(
                    z.get(
                        "reactant_order",
                        0.0,
                    )
                    or 0.0
                )
                < 2.8
                and float(
                    z.get(
                        "product_order",
                        0.0,
                    )
                    or 0.0
                )
                <= 1.2
            )
        ]

        if len(inc) == 1 and len(dec) == 1:

            im = {
                int(inc[0]["a"]),
                int(inc[0]["b"]),
            }

            dm = {
                int(dec[0]["a"]),
                int(dec[0]["b"]),
            }

            support = im | dm

            if (
                len(im & dm) == 1
                and len(support) == 3
                and all(
                    _atom(
                        facts,
                        "reactant",
                        m,
                    ).get("element")
                    == "C"
                    for m in support
                )
            ):

                skeleton_ok = all(
                    external_mapped_heavy_neighbors(
                        "reactant",
                        m,
                        support,
                    )
                    == external_mapped_heavy_neighbors(
                        "product",
                        m,
                        support,
                    )
                    for m in support
                )

                hd = [
                    _hydrogen_delta(
                        facts,
                        m,
                    )
                    for m in support
                ]

                hydrogen_shift = sum(hd) == 0 and min(hd) < 0 and max(hd) > 0

                if skeleton_ok and hydrogen_shift:

                    out.append(
                        {
                            "event": "double_bond_migration",
                            "role": "specific_event",
                            "subtype": "alkene_isomerization",
                            "support_maps": sorted(support),
                            "mapping_scope": (
                                "adjacent_mapped_c_c_"
                                "single_double_exchange_"
                                "with_conserved_external_"
                                "heavy_skeleton_and_net_"
                                "hydrogen_shift"
                            ),
                        }
                    )

    return out


# ======================================================================
# S2_5_DIRECTIONAL_FGI
#
# Final directional-FGI closure.
#
# Structural contracts only:
#
#   1. C=C -> C-C-OH outcome with product hydroxyl O installation
#   2. aryl/heteroaryl C-NH2 -> C-H
#   3. aldoxime C=N-OH -> nitrile C#N
#   4. O-alkyl enol ether -> aldehyde/ketone
#   5. mapping-limited ethyl ester -> methyl ester
#
# No template identities, reagent names, or method names are used.
# ======================================================================


def detect_s25_directional_fgi(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []

    common = {
        int(x)
        for x in (
            facts.get(
                "mapping",
                {},
            ).get(
                "common_maps",
                [],
            )
            or []
        )
    }

    seen = set()

    def atom(
        side: str,
        mn: int,
    ) -> dict[str, Any]:

        return _atom(
            facts,
            side,
            int(mn),
        )

    def atts(
        side: str,
        mn: int,
    ) -> list[dict[str, Any]]:

        return (
            _attachments(
                facts,
                side,
                int(mn),
            )
            or []
        )

    def bo(
        side: str,
        a: int,
        b: int,
    ) -> float:

        return float(
            _bond_between(
                facts,
                side,
                int(a),
                int(b),
            )
            or 0.0
        )

    def heavy(
        rows,
    ):

        return [
            x
            for x in rows
            if str(x.get("element") or "")
            not in {
                "",
                "H",
            }
        ]

    def mapped_heavy_neighbors(
        side: str,
        center: int,
        exclude=None,
    ) -> set[int]:

        exclude = {int(x) for x in (exclude or [])}

        ans = set()

        for x in atts(
            side,
            center,
        ):

            if x.get("element") == "H":
                continue

            m = _att_map(x)

            if m is None:
                continue

            m = int(m)

            if m in exclude:
                continue

            ans.add(m)

        return ans

    def emit(
        directional_name: str,
        specific_event: str,
        center: int,
        support,
        source: str,
        target: str,
        *,
        direction: str = "neutral",
        family: str = "functional_group_interconversion",
        note: str,
    ):

        support = _support(*support)

        key = (
            directional_name,
            specific_event,
            int(center),
            tuple(sorted(support)),
        )

        if key in seen:
            return

        seen.add(key)

        out.extend(
            _conversion_records(
                directional_name,
                int(center),
                support,
                direction,
                family,
                source,
                target,
                True,
                specific_event_name=specific_event,
                extra={
                    "s2_5_directional_fgi": True,
                    "mapping_scope_note": note,
                },
            )
        )

    # ------------------------------------------------------------------
    # 1. Alkene -> alcohol outcome.
    #
    # Existing alkene_hydration owns water-like/product-only oxygen.
    # This rule covers the remaining product-centered C=C -> C-C-OH
    # pattern where the mapped product oxygen already had a heavy-atom
    # environment in the reactant.
    #
    # This deliberately does NOT infer hydroboration as a mechanism.
    # Mapped O provenance may reflect peroxide, solvent, or mapper choice.
    # ------------------------------------------------------------------

    for bc in (
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    ):

        a = int(bc["a"])

        b = int(bc["b"])

        if (
            atom(
                "reactant",
                a,
            ).get("element")
            != "C"
            or atom(
                "reactant",
                b,
            ).get("element")
            != "C"
        ):
            continue

        ro = float(
            bc.get(
                "reactant_order",
                0.0,
            )
            or 0.0
        )

        po = float(
            bc.get(
                "product_order",
                0.0,
            )
            or 0.0
        )

        if not (1.7 <= ro <= 2.2 and 0.8 <= po <= 1.2):
            continue

        for center, partner in (
            (
                a,
                b,
            ),
            (
                b,
                a,
            ),
        ):

            lost_c, gained_c = _lost_gained_for_center(
                facts,
                center,
            )

            lost_p, gained_p = _lost_gained_for_center(
                facts,
                partner,
            )

            gained_o = [
                x
                for x in gained_c
                if (
                    x.get("element") == "O"
                    and 0.8
                    <= float(
                        x.get(
                            "bond_order",
                            0.0,
                        )
                        or 0.0
                    )
                    <= 1.2
                    and _att_map(x) is not None
                )
            ]

            if len(gained_o) != 1:

                continue

            om = int(_att_map(gained_o[0]))

            if om not in common:
                continue

            if not (
                0.8
                <= bo(
                    "product",
                    center,
                    om,
                )
                <= 1.2
            ):

                continue

            # Product O must be hydroxyl-like:
            # one heavy neighbor, namely the alkene carbon.
            p_o_heavy = heavy(
                atts(
                    "product",
                    om,
                )
            )

            if len(p_o_heavy) != 1:

                continue

            p_o_map = _att_map(p_o_heavy[0])

            if p_o_map is None or int(p_o_map) != center:

                continue

            # Skip ordinary hydration oxygen, already owned by
            # alkene_hydration.  Here the mapped O must have existed
            # in a non-water-like heavy-atom environment.
            if (
                atom(
                    "reactant",
                    om,
                ).get("element")
                != "O"
            ):

                continue

            r_o_heavy = heavy(
                atts(
                    "reactant",
                    om,
                )
            )

            if not r_o_heavy:
                continue

            # No other heavy attachment exchange at the former
            # alkene endpoints.
            other_center_gain = [x for x in heavy(gained_c) if _att_map(x) != om]

            if heavy(lost_c) or heavy(lost_p) or other_center_gain or heavy(gained_p):

                continue

            if (
                _hydrogen_delta(
                    facts,
                    partner,
                )
                <= 0
            ):

                continue

            pfg = _carbon_fg_state(
                facts,
                "product",
                center,
            )

            if pfg not in {
                "primary_alcohol",
                "secondary_alcohol",
                "tertiary_alcohol",
            }:

                continue

            emit(
                "alkene_to_alcohol",
                "alkene_to_alcohol_oxygenation",
                center,
                [
                    a,
                    b,
                    om,
                ],
                "alkene",
                pfg,
                note=(
                    "mapped_c_c_double_to_single_plus_"
                    "product_hydroxyl_oxygen_installation;"
                    "oxygen_provenance_not_interpreted_"
                    "mechanistically"
                ),
            )

    # ------------------------------------------------------------------
    # 2. Aryl / heteroaryl amine -> C-H deamination.
    #
    # A conserved aromatic carbon loses exactly one N attachment,
    # gains no heavy atom, remains aromatic, and gains H.
    # ------------------------------------------------------------------

    for d in (
        facts.get(
            "attachment_deltas",
            [],
        )
        or []
    ):

        c = int(d["center_map"])

        if d.get("center_element") != "C":

            continue

        ra = atom(
            "reactant",
            c,
        )

        pa = atom(
            "product",
            c,
        )

        if not (ra.get("aromatic") and pa.get("aromatic")):

            continue

        lost_h = heavy(
            d.get(
                "lost",
                [],
            )
            or []
        )

        gained_h = heavy(
            d.get(
                "gained",
                [],
            )
            or []
        )

        if len(lost_h) != 1 or lost_h[0].get("element") != "N" or gained_h:

            continue

        if (
            int(
                d.get(
                    "hydrogen_delta",
                    0,
                )
                or 0
            )
            <= 0
        ):

            continue

        nm = _att_map(lost_h[0])

        emit(
            "aryl_amine_to_hydrogen",
            "aryl_amine_deamination",
            c,
            [
                c,
                nm,
            ],
            "aryl_or_heteroaryl_amine",
            "aryl_or_heteroaryl_hydrogen",
            note=("conserved_aromatic_carbon_loses_" "single_n_attachment_and_gains_h"),
        )

    # ------------------------------------------------------------------
    # 3. Aldoxime -> nitrile.
    #
    #     R-CH=N-OH  ->  R-C#N
    #
    # Same mapped C and N; C=N becomes C#N; N loses O; the carbon
    # skeleton around the oxime carbon is retained.
    # ------------------------------------------------------------------

    for c in sorted(common):

        ca = atom(
            "reactant",
            c,
        )

        if ca.get("element") != "C":

            continue

        if (
            int(
                ca.get(
                    "total_h",
                    0,
                )
                or 0
            )
            < 1
        ):

            # Restrict to aldoxime-like centers.
            continue

        rn = [
            x
            for x in atts(
                "reactant",
                c,
            )
            if (
                x.get("element") == "N"
                and 1.7
                <= float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                <= 2.2
                and _att_map(x) is not None
            )
        ]

        if len(rn) != 1:

            continue

        n = int(_att_map(rn[0]))

        if not (
            bo(
                "product",
                c,
                n,
            )
            >= 2.8
        ):

            continue

        if (
            _carbon_fg_state(
                facts,
                "product",
                c,
            )
            != "nitrile"
        ):

            continue

        r_n_o = [
            x
            for x in atts(
                "reactant",
                n,
            )
            if x.get("element") == "O"
        ]

        p_n_o = [
            x
            for x in atts(
                "product",
                n,
            )
            if x.get("element") == "O"
        ]

        if not r_n_o or p_n_o:

            continue

        r_skeleton = mapped_heavy_neighbors(
            "reactant",
            c,
            exclude={
                n,
            },
        )

        p_skeleton = mapped_heavy_neighbors(
            "product",
            c,
            exclude={
                n,
            },
        )

        if r_skeleton != p_skeleton:

            continue

        emit(
            "oxime_to_nitrile",
            "oxime_to_nitrile",
            c,
            [
                c,
                n,
                *[_att_map(x) for x in r_n_o],
            ],
            "aldoxime",
            "nitrile",
            note=(
                "same_mapped_c_n_pair_changes_" "c_n_double_to_triple_while_n_loses_o"
            ),
        )

    # ------------------------------------------------------------------
    # 4. O-alkyl enol ether -> carbonyl.
    #
    #     R-C(OR)=C  ->  R-C(=O)-C
    #
    # Excludes silyl enol ethers: reactant O must have another carbon
    # neighbor, not Si.
    # ------------------------------------------------------------------

    for bc in (
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    ):

        a = int(bc["a"])

        b = int(bc["b"])

        if (
            atom(
                "reactant",
                a,
            ).get("element")
            != "C"
            or atom(
                "reactant",
                b,
            ).get("element")
            != "C"
        ):

            continue

        ro = float(
            bc.get(
                "reactant_order",
                0.0,
            )
            or 0.0
        )

        po = float(
            bc.get(
                "product_order",
                0.0,
            )
            or 0.0
        )

        if not (1.7 <= ro <= 2.2 and 0.8 <= po <= 1.2):

            continue

        for c, partner in (
            (
                a,
                b,
            ),
            (
                b,
                a,
            ),
        ):

            pfg = _carbon_fg_state(
                facts,
                "product",
                c,
            )

            if pfg not in {
                "aldehyde",
                "ketone",
            }:

                continue

            r_oxy = [
                x
                for x in atts(
                    "reactant",
                    c,
                )
                if (
                    x.get("element") == "O"
                    and 0.8
                    <= float(
                        x.get(
                            "bond_order",
                            0.0,
                        )
                        or 0.0
                    )
                    <= 1.2
                    and _att_map(x) is not None
                )
            ]

            for ox in r_oxy:

                om = int(_att_map(ox))

                if not (
                    bo(
                        "product",
                        c,
                        om,
                    )
                    >= 1.8
                ):

                    continue

                # Reactant oxygen must be alkoxy/enol-ether O:
                # another carbon neighbor besides the enol carbon.
                other_r_c = [
                    x
                    for x in atts(
                        "reactant",
                        om,
                    )
                    if (
                        x.get("element") == "C"
                        and (_att_map(x) is None or int(_att_map(x)) != c)
                    )
                ]

                if not other_r_c:
                    continue

                # Product O is now the carbonyl O only.
                p_o_heavy = heavy(
                    atts(
                        "product",
                        om,
                    )
                )

                if len(p_o_heavy) != 1:

                    continue

                pm = _att_map(p_o_heavy[0])

                if pm is None or int(pm) != c:

                    continue

                emit(
                    "enol_ether_to_carbonyl",
                    "enol_ether_to_carbonyl",
                    c,
                    [
                        c,
                        partner,
                        om,
                    ],
                    "o_alkyl_enol_ether",
                    pfg,
                    note=(
                        "mapped_enol_c_c_double_reduces_"
                        "while_same_enol_oxygen_becomes_"
                        "carbonyl_oxygen_and_loses_alkyl_context"
                    ),
                )

    # ------------------------------------------------------------------
    # 5. Mapping-limited ethyl ester -> methyl ester.
    #
    # The same ester oxygen and alpha alkoxy carbon are retained:
    #
    #     -C(=O)-O-CH2-CH3  ->  -C(=O)-O-CH3
    #
    # The mapped CH2 becomes CH3 and loses its distal carbon.
    # This is the common atom-mapped representation of an ester alkyl
    # exchange for which the acyl-center O-for-O detector cannot fire.
    # ------------------------------------------------------------------

    for c in sorted(common):

        if (
            _carbon_fg_state(
                facts,
                "reactant",
                c,
            )
            != "ester"
            or _carbon_fg_state(
                facts,
                "product",
                c,
            )
            != "ester"
        ):

            continue

        r_os = [
            x
            for x in _element_attachments(
                facts,
                "reactant",
                c,
                "O",
                0.8,
                1.2,
            )
            if _att_map(x) is not None
        ]

        for ro in r_os:

            om = int(_att_map(ro))

            if not (
                0.8
                <= bo(
                    "product",
                    c,
                    om,
                )
                <= 1.2
            ):

                continue

            r_carbons = {
                int(_att_map(x))
                for x in atts(
                    "reactant",
                    om,
                )
                if (
                    x.get("element") == "C"
                    and _att_map(x) is not None
                    and int(_att_map(x)) != c
                )
            }

            p_carbons = {
                int(_att_map(x))
                for x in atts(
                    "product",
                    om,
                )
                if (
                    x.get("element") == "C"
                    and _att_map(x) is not None
                    and int(_att_map(x)) != c
                )
            }

            for alkyl_c in sorted(r_carbons & p_carbons):

                ra = atom(
                    "reactant",
                    alkyl_c,
                )

                pa = atom(
                    "product",
                    alkyl_c,
                )

                if ra.get("element") != "C" or pa.get("element") != "C":

                    continue

                if not (
                    int(
                        ra.get(
                            "total_h",
                            0,
                        )
                        or 0
                    )
                    == 2
                    and int(
                        pa.get(
                            "total_h",
                            0,
                        )
                        or 0
                    )
                    >= 3
                ):

                    continue

                lost_k, gained_k = _lost_gained_for_center(
                    facts,
                    alkyl_c,
                )

                lost_carbons = [x for x in heavy(lost_k) if x.get("element") == "C"]

                if not lost_carbons or heavy(gained_k):

                    continue

                emit(
                    "ethyl_ester_to_methyl_ester",
                    "transesterification",
                    c,
                    [
                        c,
                        om,
                        alkyl_c,
                        *[_att_map(x) for x in lost_carbons],
                    ],
                    "ethyl_ester",
                    "methyl_ester",
                    note=(
                        "mapping_limited_same_ester_o_"
                        "retained_while_alkoxy_ch2_"
                        "becomes_ch3"
                    ),
                )

    return out


# =============================================================================
# POST_V655_N_THP_PROTECTION
#
# Residual semantic closure discovered after the frozen v6.5.5 release.
#
# Structural ownership is ring-local and tautomer-aware:
#
#   N-H + six-membered cyclic enol ether
#       -> N-tetrahydropyranyl product
#
# Exact mapped C=C endpoint identity is NOT required.
# Exact identity of the proton-bearing tautomeric N is NOT required.
# =============================================================================


def detect_n_tetrahydropyranyl_protection(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    """Detect N-THP protection from a six-membered cyclic enol ether."""

    common = {int(x) for x in facts["mapping"]["common_maps"]}

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    def nbrs(
        side: str,
        m: int,
    ) -> set[int]:

        return {
            int(_att_map(a))
            for a in _attachments(
                facts,
                side,
                m,
            )
            if _att_map(a) is not None
        }

    def bond_order(
        side: str,
        a: int,
        b: int,
    ) -> float:

        return next(
            (
                float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                for x in _attachments(
                    facts,
                    side,
                    a,
                )
                if (_att_map(x) is not None and int(_att_map(x)) == b)
            ),
            0.0,
        )

    def cycles(
        side: str,
        start: int,
        size: int,
        *,
        aromatic: bool = False,
    ) -> list[tuple[int, ...]]:
        """Enumerate small mapped simple cycles containing start."""

        found: dict[
            frozenset[int],
            tuple[int, ...],
        ] = {}

        def walk(
            path: list[int],
        ) -> None:

            cur = path[-1]

            if len(path) == size:

                if start in nbrs(
                    side,
                    cur,
                ):
                    found.setdefault(
                        frozenset(path),
                        tuple(path),
                    )

                return

            for nxt in sorted(
                nbrs(
                    side,
                    cur,
                )
            ):

                if nxt == start or nxt in path or nxt not in common:
                    continue

                if aromatic and not _atom(
                    facts,
                    side,
                    nxt,
                ).get("aromatic"):
                    continue

                walk(path + [nxt])

        if start in common and (
            not aromatic
            or _atom(
                facts,
                side,
                start,
            ).get("aromatic")
        ):
            walk([start])

        return list(found.values())

    def same_product_cycle(
        ring: set[int],
        start: int,
        *,
        aromatic: bool = False,
    ) -> bool:

        return any(
            set(cycle) == ring
            for cycle in cycles(
                "product",
                start,
                len(ring),
                aromatic=aromatic,
            )
        )

    def ring_cc_double_count(
        side: str,
        ring: set[int],
    ) -> int:

        pairs: set[tuple[int, int]] = set()

        for a in ring:

            for b in (
                nbrs(
                    side,
                    a,
                )
                & ring
            ):

                pairs.add(
                    tuple(
                        sorted(
                            (
                                a,
                                b,
                            )
                        )
                    )
                )

        return sum(
            1
            for a, b in pairs
            if (
                _atom(
                    facts,
                    side,
                    a,
                ).get("element")
                == "C"
                and _atom(
                    facts,
                    side,
                    b,
                ).get("element")
                == "C"
                and bond_order(
                    side,
                    a,
                    b,
                )
                >= 1.8
            )
        )

    def nh_support(
        nitrogen: int,
    ) -> tuple[
        bool,
        str,
        list[int],
    ]:

        r_n = _atom(
            facts,
            "reactant",
            nitrogen,
        )

        p_n = _atom(
            facts,
            "product",
            nitrogen,
        )

        r_h = int(
            r_n.get(
                "total_h",
                0,
            )
            or 0
        )

        p_h = int(
            p_n.get(
                "total_h",
                0,
            )
            or 0
        )

        # --------------------------------------------------------
        # Mode A:
        # exact substituted N bears H in reactant.
        # --------------------------------------------------------

        if r_h > 0 and p_h < r_h:

            return (
                True,
                "exact_n_h",
                [nitrogen],
            )

        # --------------------------------------------------------
        # Mode B:
        # [nH] belongs to another mapped N in the same
        # aromatic tautomeric ring.
        # --------------------------------------------------------

        if not (r_n.get("aromatic") and p_n.get("aromatic")):

            return (
                False,
                "",
                [],
            )

        for size in (
            5,
            6,
        ):

            for cycle in cycles(
                "reactant",
                nitrogen,
                size,
                aromatic=True,
            ):

                ring = set(cycle)

                if not same_product_cycle(
                    ring,
                    nitrogen,
                    aromatic=True,
                ):
                    continue

                n_maps = [
                    m
                    for m in ring
                    if _atom(
                        facts,
                        "reactant",
                        m,
                    ).get("element")
                    == "N"
                ]

                if len(n_maps) < 2:
                    continue

                reactant_nh = [
                    m
                    for m in n_maps
                    if int(
                        _atom(
                            facts,
                            "reactant",
                            m,
                        ).get(
                            "total_h",
                            0,
                        )
                        or 0
                    )
                    > 0
                ]

                product_nh = [
                    m
                    for m in n_maps
                    if int(
                        _atom(
                            facts,
                            "product",
                            m,
                        ).get(
                            "total_h",
                            0,
                        )
                        or 0
                    )
                    > 0
                ]

                if reactant_nh and len(product_nh) < len(reactant_nh):

                    return (
                        True,
                        "tautomeric_ring_n_h",
                        sorted(ring),
                    )

        return (
            False,
            "",
            [],
        )

    # ------------------------------------------------------------
    # New N-C bond is the entry point.
    # ------------------------------------------------------------

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        if (
            float(
                fb.get(
                    "order",
                    1.0,
                )
                or 1.0
            )
            > 1.2
        ):
            continue

        a = int(fb["a"])

        b = int(fb["b"])

        ea = _atom(
            facts,
            "reactant",
            a,
        ).get("element")

        eb = _atom(
            facts,
            "reactant",
            b,
        ).get("element")

        if {
            ea,
            eb,
        } != {
            "C",
            "N",
        }:
            continue

        carbon, nitrogen = (a, b) if ea == "C" else (b, a)

        nh_ok, nh_mode, nh_ring = nh_support(nitrogen)

        if not nh_ok:
            continue

        # --------------------------------------------------------
        # Find conserved 6-membered 5C/1O cyclic enol ether.
        # --------------------------------------------------------

        for cycle in cycles(
            "reactant",
            carbon,
            6,
        ):

            ring = set(cycle)

            atoms = [
                _atom(
                    facts,
                    "reactant",
                    m,
                )
                for m in ring
            ]

            elements = [atom.get("element") for atom in atoms]

            if not (elements.count("C") == 5 and elements.count("O") == 1):
                continue

            if any(atom.get("aromatic") for atom in atoms):
                continue

            if not same_product_cycle(
                ring,
                carbon,
            ):
                continue

            oxygen = next(
                m
                for m in ring
                if _atom(
                    facts,
                    "reactant",
                    m,
                ).get("element")
                == "O"
            )

            # N attaches to an O-alpha carbon.
            if oxygen not in nbrs(
                "reactant",
                carbon,
            ) or oxygen not in nbrs(
                "product",
                carbon,
            ):
                continue

            # Ring-local invariant:
            #
            # DHP:  one C=C
            # THP:  zero C=C
            #
            # Do NOT demand that carbon itself is the mapped
            # C=C endpoint; this is the mapping-robust feature.
            if (
                ring_cc_double_count(
                    "reactant",
                    ring,
                )
                != 1
                or ring_cc_double_count(
                    "product",
                    ring,
                )
                != 0
            ):
                continue

            key = (
                carbon,
                nitrogen,
            )

            if key in seen:
                continue

            seen.add(key)

            out.append(
                {
                    "event": "n_tetrahydropyranyl_protection",
                    "role": "specific_event",
                    "subtype": "dihydropyran_n_protection",
                    "center_map": carbon,
                    "nitrogen_map": nitrogen,
                    "thp_oxygen_map": oxygen,
                    "thp_ring_maps": sorted(ring),
                    "nh_support_mode": nh_mode,
                    "tautomer_ring_maps": nh_ring,
                    "support_maps": _support(
                        carbon,
                        nitrogen,
                        oxygen,
                        *sorted(ring),
                        *nh_ring,
                    ),
                    "mapping_scope": (
                        "intercomponent_n_c_formation_at_"
                        "o_alpha_carbon_of_conserved_"
                        "six_membered_cyclic_enol_ether_"
                        "with_ring_pi_consumption_and_"
                        "exact_or_tautomeric_n_h_support"
                    ),
                }
            )

            break

    return out



# =============================================================================
# POST_V655_FINAL_RESIDUAL_CLOSURE
#
# Final bounded semantic closure after frozen Runtime v6.5.5.
#
# Families:
#
#   1. disulfide_formation
#   2. n_o_acetal_formation
#   3. nitro_to_hydroxylamine
#
# These detectors are reaction-structural only:
#
#   * no reagent-name dependence
#   * no condition dependence
#   * no template-ID dependence
#   * no Mapping Correction changes
#
# Allyl/decarboxylative-transfer residuals are deliberately excluded because
# their observed atom-map identity maps an ester carbonyl carbon onto the
# product allylic methylene.  That family remains representation-limited.
# =============================================================================


def detect_post_v655_disulfide_formation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    """
    Detect exact intercomponent construction of one S-S single bond.

    Two supported structural modes are retained under one product-level
    semantic owner:

        thiol + thiol/thiolate
            -> disulfide

        thiol/thiolate + sulfur-transfer precursor
            -> disulfide

    The latter includes mapped representations in which the electrophilic
    sulfur initially carries an unmapped S/Cl/etc. leaving attachment.
    Therefore the family event is deliberately `disulfide_formation`
    rather than universally claiming oxidative thiol coupling.
    """

    formed = (
        facts.get(
            "formed_bonds",
            [],
        )
        or []
    )

    broken = (
        facts.get(
            "broken_bonds",
            [],
        )
        or []
    )

    order_changed = (
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    )

    # Residual contract established during diagnosis:
    #
    #     F[S-S] | B[] | O[]
    #
    # Keep this detector intentionally narrow.
    if (
        len(formed) != 1
        or broken
        or order_changed
    ):
        return []

    out: list[dict[str, Any]] = []

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):

        try:
            a = int(fb["a"])
            b = int(fb["b"])
        except Exception:
            continue

        if float(
            fb.get(
                "order",
                1.0,
            )
            or 1.0
        ) > 1.2:
            continue

        ra = _atom(
            facts,
            "reactant",
            a,
        )

        rb = _atom(
            facts,
            "reactant",
            b,
        )

        pa = _atom(
            facts,
            "product",
            a,
        )

        pb = _atom(
            facts,
            "product",
            b,
        )

        if (
            ra.get("element") != "S"
            or rb.get("element") != "S"
            or pa.get("element") != "S"
            or pb.get("element") != "S"
        ):
            continue

        a_thiol_like = (
            int(
                ra.get(
                    "total_h",
                    0,
                )
                or 0
            )
            > 0
            or
            int(
                ra.get(
                    "formal_charge",
                    0,
                )
                or 0
            )
            < 0
        )

        b_thiol_like = (
            int(
                rb.get(
                    "total_h",
                    0,
                )
                or 0
            )
            > 0
            or
            int(
                rb.get(
                    "formal_charge",
                    0,
                )
                or 0
            )
            < 0
        )

        # At least one endpoint must be a thiol/thiolate sulfur.
        # This preserves the validated residual family while avoiding
        # arbitrary S-S construction between two unrelated sulfur centers.
        if not (
            a_thiol_like
            or b_thiol_like
        ):
            continue

        if (
            a_thiol_like
            and b_thiol_like
        ):
            subtype = (
                "thiol_pair_disulfide_formation"
            )

            source_mode = (
                "thiol_pair"
            )

        else:
            subtype = (
                "thiol_sulfur_transfer_"
                "disulfide_formation"
            )

            source_mode = (
                "thiol_plus_sulfur_transfer_center"
            )

        sulfur_maps = sorted(
            (
                a,
                b,
            )
        )

        out.append(
            {
                "event":
                    "disulfide_formation",

                "role":
                    "specific_event",

                "subtype":
                    subtype,

                "center_map":
                    sulfur_maps[0],

                "sulfur_maps":
                    sulfur_maps,

                "source_mode":
                    source_mode,

                "support_maps":
                    sulfur_maps,

                "mapping_scope":
                    (
                        "exact_intercomponent_"
                        "mapped_s_s_single_bond_formation"
                    ),
            }
        )

    return out


def detect_post_v655_n_o_acetal_formation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    """
    Detect the diagnosed formaldehyde-derived N,O-acetal center:

        CH2=O + O-nucleophile + N-nucleophile
            ->
        N-CH2-O

    The O nucleophile may be an alcohol oxygen or another O-H center
    (for example N-hydroxy oxygen); therefore the structural owner is
    N,O-acetal formation rather than the narrower aminomethyl ether name.
    """

    formed = (
        facts.get(
            "formed_bonds",
            [],
        )
        or []
    )

    broken = (
        facts.get(
            "broken_bonds",
            [],
        )
        or []
    )

    order_changed = (
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    )

    # Diagnosed residual contract:
    #
    #     F[C-N,C-O] | B[] | O[]
    #
    if (
        len(formed) != 2
        or broken
        or order_changed
    ):
        return []

    inter_pairs = set()

    for fb in (
        facts.get(
            "intercomponent_formed_bonds",
            [],
        )
        or []
    ):
        try:
            a = int(fb["a"])
            b = int(fb["b"])
        except Exception:
            continue

        if float(
            fb.get(
                "order",
                1.0,
            )
            or 1.0
        ) <= 1.2:

            inter_pairs.add(
                frozenset(
                    (
                        a,
                        b,
                    )
                )
            )

    out: list[dict[str, Any]] = []

    for c in (
        facts.get(
            "mapping",
            {},
        ).get(
            "common_maps",
            [],
        )
        or []
    ):

        c = int(c)

        ra = _atom(
            facts,
            "reactant",
            c,
        )

        pa = _atom(
            facts,
            "product",
            c,
        )

        if (
            ra.get("element") != "C"
            or pa.get("element") != "C"
        ):
            continue

        if (
            bool(
                ra.get(
                    "aromatic"
                )
            )
            or
            bool(
                pa.get(
                    "aromatic"
                )
            )
        ):
            continue

        # Formaldehyde-like source carbon.
        if int(
            ra.get(
                "total_h",
                0,
            )
            or 0
        ) != 2:
            continue

        r_atts = list(
            _attachments(
                facts,
                "reactant",
                c,
            )
            or []
        )

        r_heavy = [
            x
            for x in r_atts
            if x.get(
                "element"
            )
            not in {
                "",
                "H",
                None,
            }
        ]

        if len(
            r_heavy
        ) != 1:
            continue

        if not (
            r_heavy[0].get(
                "element"
            )
            == "O"
            and
            float(
                r_heavy[0].get(
                    "bond_order",
                    0.0,
                )
                or 0.0
            )
            >= 1.8
        ):
            continue

        # Product carbon remains CH2 and is attached to exactly one O
        # and one N through single bonds.
        if int(
            pa.get(
                "total_h",
                0,
            )
            or 0
        ) != 2:
            continue

        p_atts = list(
            _attachments(
                facts,
                "product",
                c,
            )
            or []
        )

        p_heavy = [
            x
            for x in p_atts
            if x.get(
                "element"
            )
            not in {
                "",
                "H",
                None,
            }
        ]

        if len(
            p_heavy
        ) != 2:
            continue

        o_atts = [
            x
            for x in p_heavy
            if (
                x.get(
                    "element"
                )
                == "O"
                and
                float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                <= 1.2
                and
                _att_map(
                    x
                )
                is not None
            )
        ]

        n_atts = [
            x
            for x in p_heavy
            if (
                x.get(
                    "element"
                )
                == "N"
                and
                float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                <= 1.2
                and
                _att_map(
                    x
                )
                is not None
            )
        ]

        if (
            len(o_atts) != 1
            or
            len(n_atts) != 1
        ):
            continue

        o = int(
            _att_map(
                o_atts[0]
            )
        )

        n = int(
            _att_map(
                n_atts[0]
            )
        )

        if (
            frozenset(
                (
                    c,
                    o,
                )
            )
            not in inter_pairs
            or
            frozenset(
                (
                    c,
                    n,
                )
            )
            not in inter_pairs
        ):
            continue

        out.append(
            {
                "event":
                    "n_o_acetal_formation",

                "role":
                    "specific_event",

                "subtype":
                    (
                        "formaldehyde_to_n_o_acetal"
                    ),

                "center_map":
                    c,

                "oxygen_map":
                    o,

                "nitrogen_map":
                    n,

                "support_maps":
                    _support(
                        c,
                        o,
                        n,
                    ),

                "mapping_scope":
                    (
                        "mapped_formaldehyde_carbon_"
                        "forms_intercomponent_c_o_"
                        "and_c_n_single_bonds"
                    ),
            }
        )

    return out


def detect_post_v655_nitro_to_hydroxylamine(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    """
    Detect exact mapped nitro -> N-hydroxylamine reduction.

    Validated residual geometry:

        Ar-N(+)(=O)-O(-)
              ->
        Ar-NH-OH

    The mapped N=O oxygen is retained as the hydroxylamine oxygen.
    The second nitro oxygen may be unmapped and therefore is not required
    to participate in the mapped edit graph.
    """

    formed = (
        facts.get(
            "formed_bonds",
            [],
        )
        or []
    )

    broken = (
        facts.get(
            "broken_bonds",
            [],
        )
        or []
    )

    order_changed = (
        facts.get(
            "bond_order_changes",
            [],
        )
        or []
    )

    if (
        formed
        or broken
        or len(
            order_changed
        ) != 1
    ):
        return []

    out: list[dict[str, Any]] = []

    for n in (
        facts.get(
            "mapping",
            {},
        ).get(
            "common_maps",
            [],
        )
        or []
    ):

        n = int(n)

        rn = _atom(
            facts,
            "reactant",
            n,
        )

        pn = _atom(
            facts,
            "product",
            n,
        )

        if (
            rn.get("element") != "N"
            or pn.get("element") != "N"
        ):
            continue

        if int(
            rn.get(
                "charge",
                0,
            )
            or 0
        ) != 1:
            continue

        if int(
            pn.get(
                "charge",
                0,
            )
            or 0
        ) != 0:
            continue

        if int(
            rn.get(
                "total_h",
                0,
            )
            or 0
        ) != 0:
            continue

        if int(
            pn.get(
                "total_h",
                0,
            )
            or 0
        ) < 1:
            continue

        r_atts = list(
            _attachments(
                facts,
                "reactant",
                n,
            )
            or []
        )

        p_atts = list(
            _attachments(
                facts,
                "product",
                n,
            )
            or []
        )

        r_oxygens = [
            x
            for x in r_atts
            if x.get(
                "element"
            )
            == "O"
        ]

        if len(
            r_oxygens
        ) < 2:
            continue

        r_double_o = [
            x
            for x in r_oxygens
            if (
                float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                >= 1.8
                and
                _att_map(
                    x
                )
                is not None
            )
        ]

        if len(
            r_double_o
        ) != 1:
            continue

        o = int(
            _att_map(
                r_double_o[0]
            )
        )

        p_oxygens = [
            x
            for x in p_atts
            if (
                x.get(
                    "element"
                )
                == "O"
                and
                _att_map(
                    x
                )
                is not None
                and
                int(
                    _att_map(
                        x
                    )
                )
                == o
                and
                float(
                    x.get(
                        "bond_order",
                        0.0,
                    )
                    or 0.0
                )
                <= 1.2
            )
        ]

        if len(
            p_oxygens
        ) != 1:
            continue

        po = _atom(
            facts,
            "product",
            o,
        )

        if int(
            po.get(
                "total_h",
                0,
            )
            or 0
        ) < 1:
            continue

        # Same mapped C-N attachment must survive.
        r_c = {
            int(
                _att_map(
                    x
                )
            )
            for x in r_atts
            if (
                x.get(
                    "element"
                )
                == "C"
                and
                _att_map(
                    x
                )
                is not None
            )
        }

        p_c = {
            int(
                _att_map(
                    x
                )
            )
            for x in p_atts
            if (
                x.get(
                    "element"
                )
                == "C"
                and
                _att_map(
                    x
                )
                is not None
            )
        }

        retained_c = sorted(
            r_c
            &
            p_c
        )

        if not retained_c:
            continue

        # The one mapped order change must be N=O -> N-O
        # at the same N/O pair.
        ch = order_changed[0]

        try:
            pair = {
                int(
                    ch["a"]
                ),
                int(
                    ch["b"]
                ),
            }
        except Exception:
            continue

        if pair != {
            n,
            o,
        }:
            continue

        if not (
            float(
                ch.get(
                    "reactant_order",
                    0.0,
                )
                or 0.0
            )
            >= 1.8
            and
            float(
                ch.get(
                    "product_order",
                    0.0,
                )
                or 0.0
            )
            <= 1.2
        ):
            continue

        c = retained_c[0]

        out.extend(
            _conversion_records(
                "nitro_to_hydroxylamine",
                n,
                _support(
                    n,
                    o,
                    c,
                ),
                "reduction",
                "reduction",
                "nitro",
                "hydroxylamine",
                True,
                specific_event_name=(
                    "nitro_to_hydroxylamine"
                ),
                extra={
                    "mapping_scope":
                        (
                            "mapped_nitro_n_o_"
                            "to_n_hydroxylamine"
                        ),
                },
            )
        )

    return out


# =============================================================================
# POST_V655_FINAL_FIVE_SEMANTIC_FAMILIES
#
# Final bounded post-v6.5.5 semantic families.
#
# Structural contracts only:
#
#   1. arylboron_o_coupling
#   2. arylboronic_acid_to_phenol
#   3. o_difluoromethylation
#   4. nitrate_ester_formation
#   5. boekelheide_rearrangement
#
# These rules do not use conditions/reagents for semantic ownership.
# =============================================================================


_POST_V655_ORGANIC_O_NEIGHBORS = {
    "C",
    "N",
    "S",
    "P",
    "Si",
}


def _post_v655_formed_c_o_pairs(
    facts: dict[str, Any],
) -> list[tuple[int, int]]:
    """
    Return newly formed mapped C-O bonds as (carbon_map, oxygen_map).
    """

    out: list[tuple[int, int]] = []

    for rec in facts.get("formed_bonds", []) or []:

        a = int(rec["a"])
        b = int(rec["b"])

        ae = _atom(
            facts,
            "product",
            a,
        ).get("element")

        be = _atom(
            facts,
            "product",
            b,
        ).get("element")

        if ae == "C" and be == "O":
            out.append((a, b))

        elif ae == "O" and be == "C":
            out.append((b, a))

    return out


def _post_v655_reactant_oh_or_oxide(
    facts: dict[str, Any],
    oxygen_map: int,
) -> bool:

    a = _atom(
        facts,
        "reactant",
        oxygen_map,
    )

    if a.get("element") != "O":
        return False

    return (
        int(
            a.get(
                "total_h",
                0,
            )
            or 0
        )
        > 0
        or
        int(
            a.get(
                "charge",
                0,
            )
            or 0
        )
        < 0
    )


def _post_v655_aryl_c_has_boron(
    facts: dict[str, Any],
    carbon_map: int,
) -> bool:

    c = _atom(
        facts,
        "reactant",
        carbon_map,
    )

    if (
        c.get("element") != "C"
        or not bool(
            c.get("aromatic")
        )
    ):
        return False

    return any(
        x.get("element") == "B"
        for x in _attachments(
            facts,
            "reactant",
            carbon_map,
        )
    )


def _post_v655_retained_organic_o_neighbor(
    facts: dict[str, Any],
    oxygen_map: int,
    incoming_carbon_map: int,
) -> list[int]:

    before = {
        int(x["map"])
        for x in _attachments(
            facts,
            "reactant",
            oxygen_map,
        )
        if (
            x.get("map") is not None
            and
            x.get("element")
            in _POST_V655_ORGANIC_O_NEIGHBORS
        )
    }

    after = {
        int(x["map"])
        for x in _attachments(
            facts,
            "product",
            oxygen_map,
        )
        if (
            x.get("map") is not None
            and
            int(x["map"])
            != int(
                incoming_carbon_map
            )
            and
            x.get("element")
            in _POST_V655_ORGANIC_O_NEIGHBORS
        )
    }

    return sorted(
        before & after
    )


# -----------------------------------------------------------------------------
# 1. Arylboron O-coupling
#
# Ar-B + R-OH/O(-) -> Ar-O-R
#
# The recipient oxygen must retain an independently mapped organic substituent.
# This separates ether formation from Ar-B -> Ar-OH hydroxylation.
# -----------------------------------------------------------------------------

def detect_post_v655_arylboron_o_coupling(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    for c, o in _post_v655_formed_c_o_pairs(
        facts
    ):

        if not _post_v655_aryl_c_has_boron(
            facts,
            c,
        ):
            continue

        if not _post_v655_reactant_oh_or_oxide(
            facts,
            o,
        ):
            continue

        retained = (
            _post_v655_retained_organic_o_neighbor(
                facts,
                o,
                c,
            )
        )

        if not retained:
            continue

        key = (
            int(c),
            int(o),
        )

        if key in seen:
            continue

        seen.add(key)

        out.append(
            {
                "event":
                    "arylboron_o_coupling",

                "role":
                    "specific_event",

                "subtype":
                    "arylboron_to_aryl_ether",

                "center_map":
                    int(c),

                "oxygen_map":
                    int(o),

                "retained_oxygen_substituent_maps":
                    retained,

                "support_maps":
                    _support(
                        c,
                        o,
                        *retained,
                    ),

                "mapping_scope":
                    (
                        "mapped_aromatic_c_b_source_"
                        "plus_external_o_c_formation_"
                        "with_retained_mapped_organic_"
                        "oxygen_substituent"
                    ),
            }
        )

    return out


# -----------------------------------------------------------------------------
# 2. Arylboronic acid / arylboron -> phenol
#
# Ar-B -> Ar-OH
#
# Product oxygen must be terminal hydroxyl rather than an ether oxygen.
# -----------------------------------------------------------------------------

def detect_post_v655_arylboronic_acid_to_phenol(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    for c, o in _post_v655_formed_c_o_pairs(
        facts
    ):

        if not _post_v655_aryl_c_has_boron(
            facts,
            c,
        ):
            continue

        po = _atom(
            facts,
            "product",
            o,
        )

        if po.get("element") != "O":
            continue

        if (
            int(
                po.get(
                    "total_h",
                    0,
                )
                or 0
            )
            < 1
        ):
            continue

        # Terminal phenolic oxygen:
        # no independently mapped organic substituent besides incoming aryl C.
        product_other_organic = [
            int(x["map"])
            for x in _attachments(
                facts,
                "product",
                o,
            )
            if (
                x.get("map") is not None
                and
                int(x["map"]) != int(c)
                and
                x.get("element")
                in _POST_V655_ORGANIC_O_NEIGHBORS
            )
        ]

        if product_other_organic:
            continue

        key = (
            int(c),
            int(o),
        )

        if key in seen:
            continue

        seen.add(key)

        out.append(
            {
                "event":
                    "arylboronic_acid_to_phenol",

                "role":
                    "specific_event",

                "subtype":
                    "arylboron_hydroxylation",

                "center_map":
                    int(c),

                "oxygen_map":
                    int(o),

                "support_maps":
                    _support(
                        c,
                        o,
                    ),

                "mapping_scope":
                    (
                        "mapped_aromatic_c_b_source_"
                        "to_terminal_mapped_phenolic_o_h"
                    ),
            }
        )

    return out


# -----------------------------------------------------------------------------
# 3. O-difluoromethylation
#
# Cl-CF2(H/...) donor representation -> O-CF2H product.
# -----------------------------------------------------------------------------

def detect_post_v655_o_difluoromethylation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    for c, o in _post_v655_formed_c_o_pairs(
        facts
    ):

        rc = _atom(
            facts,
            "reactant",
            c,
        )

        pc = _atom(
            facts,
            "product",
            c,
        )

        if (
            rc.get("element") != "C"
            or
            pc.get("element") != "C"
        ):
            continue

        r_att = _attachments(
            facts,
            "reactant",
            c,
        )

        p_att = _attachments(
            facts,
            "product",
            c,
        )

        r_f = sum(
            x.get("element") == "F"
            for x in r_att
        )

        p_f = sum(
            x.get("element") == "F"
            for x in p_att
        )

        r_cl = sum(
            x.get("element") == "Cl"
            for x in r_att
        )

        if not (
            r_f == 2
            and
            p_f == 2
            and
            r_cl >= 1
        ):
            continue

        if (
            int(
                pc.get(
                    "total_h",
                    0,
                )
                or 0
            )
            < 1
        ):
            continue

        if not _post_v655_reactant_oh_or_oxide(
            facts,
            o,
        ):
            continue

        key = (
            int(c),
            int(o),
        )

        if key in seen:
            continue

        seen.add(key)

        out.append(
            {
                "event":
                    "o_difluoromethylation",

                "role":
                    "specific_event",

                "subtype":
                    (
                        "chlorodifluoromethyl_"
                        "to_o_difluoromethyl"
                    ),

                "center_map":
                    int(c),

                "oxygen_map":
                    int(o),

                "support_maps":
                    _support(
                        c,
                        o,
                    ),

                "mapping_scope":
                    (
                        "mapped_c_f2_center_with_"
                        "reactant_chlorine_attachment_"
                        "and_external_o_c_formation"
                    ),
            }
        )

    return out


# -----------------------------------------------------------------------------
# 4. Nitrate ester formation
#
# ROH/O(-) + nitrate-like mapped N -> R-O-NO2
# -----------------------------------------------------------------------------

def detect_post_v655_nitrate_ester_formation(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()

    for rec in facts.get(
        "formed_bonds",
        [],
    ) or []:

        a = int(rec["a"])
        b = int(rec["b"])

        ae = _atom(
            facts,
            "product",
            a,
        ).get("element")

        be = _atom(
            facts,
            "product",
            b,
        ).get("element")

        if {
            ae,
            be,
        } != {
            "N",
            "O",
        }:
            continue

        if ae == "N":
            n, o = a, b
        else:
            n, o = b, a

        if not _post_v655_reactant_oh_or_oxide(
            facts,
            o,
        ):
            continue

        rn = _atom(
            facts,
            "reactant",
            n,
        )

        if rn.get("element") != "N":
            continue

        if (
            int(
                rn.get(
                    "charge",
                    0,
                )
                or 0
            )
            < 0
        ):
            continue

        reactant_oxygen_neighbors = [
            x
            for x in _attachments(
                facts,
                "reactant",
                n,
            )
            if x.get("element") == "O"
        ]

        if (
            len(
                reactant_oxygen_neighbors
            )
            < 2
        ):
            continue

        key = (
            int(n),
            int(o),
        )

        if key in seen:
            continue

        seen.add(key)

        out.append(
            {
                "event":
                    "nitrate_ester_formation",

                "role":
                    "specific_event",

                "subtype":
                    "alcohol_to_nitrate_ester",

                "center_map":
                    int(o),

                "nitrogen_map":
                    int(n),

                "support_maps":
                    _support(
                        o,
                        n,
                    ),

                "mapping_scope":
                    (
                        "mapped_nitrate_like_n_"
                        "plus_external_o_n_formation_"
                        "from_o_h_or_o_minus"
                    ),
            }
        )

    return out


# -----------------------------------------------------------------------------
# 5. Boekelheide rearrangement
#
# Aza-N-oxide deoxygenation accompanied by oxygen installation at the alpha
# substituent / adjacent carbon environment.
#
# This is intentionally stricter than simply promoting the existing
# aza_n_oxide_to_aza_nitrogen directional conversion.
# -----------------------------------------------------------------------------

def detect_post_v655_boekelheide_rearrangement(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:

    out: list[dict[str, Any]] = []
    seen: set[tuple[int, int, int]] = set()

    co_pairs = _post_v655_formed_c_o_pairs(
        facts
    )

    # Independent probe contract: exactly one newly formed mapped C-O bond.
    if len(co_pairs) != 1:
        return out

    c, o = co_pairs[0]

    # PATCH-RRG-05: direct aromatic-ring N-oxide oxygenation is not
    # Boekelheide alpha-sidechain oxygen transfer.
    if bool(_atom(facts, "reactant", c).get("aromatic")):
        return out

    nearby_n: set[int] = set()

    # Direct C-N relation.
    for x in _attachments(
        facts,
        "product",
        c,
    ):

        if (
            x.get("element") == "N"
            and
            x.get("map") is not None
        ):
            nearby_n.add(
                int(x["map"])
            )

        # Side-chain C attached to ring carbon adjacent to N.
        elif (
            x.get("element") == "C"
            and
            x.get("map") is not None
        ):

            cm = int(
                x["map"]
            )

            for y in _attachments(
                facts,
                "product",
                cm,
            ):

                if (
                    y.get("element") == "N"
                    and
                    y.get("map") is not None
                ):
                    nearby_n.add(
                        int(y["map"])
                    )


    for n in sorted(
        nearby_n
    ):

        rn = _atom(
            facts,
            "reactant",
            n,
        )

        pn = _atom(
            facts,
            "product",
            n,
        )

        if (
            rn.get("element") != "N"
            or
            pn.get("element") != "N"
        ):
            continue

        # True aza-N-oxide chemistry requires the N atom itself to be
        # an aromatic ring atom. This excludes exocyclic nitro N(+)-O(-),
        # which otherwise shares the same formal charge motif.
        if not (
            bool(rn.get("aromatic"))
            and
            bool(rn.get("in_ring"))
        ):
            continue

        # True aza-N-oxide chemistry requires the N atom itself to be
        # an aromatic ring atom. This excludes exocyclic nitro N(+)-O(-),
        # which otherwise shares the same formal charge motif.
        if not (
            bool(rn.get("aromatic"))
            and
            bool(rn.get("in_ring"))
        ):
            continue


        reactant_n_oxide = any(
            (
                x.get("element") == "O"
                and
                int(
                    x.get(
                        "charge",
                        0,
                    )
                    or 0
                )
                < 0
            )
            for x in _attachments(
                facts,
                "reactant",
                n,
            )
        )

        if not reactant_n_oxide:
            continue


        product_n_oxide = any(
            (
                x.get("element") == "O"
                and
                int(
                    x.get(
                        "charge",
                        0,
                    )
                    or 0
                )
                < 0
            )
            for x in _attachments(
                facts,
                "product",
                n,
            )
        )

        if product_n_oxide:
            continue


        key = (
            int(n),
            int(c),
            int(o),
        )

        if key in seen:
            continue

        seen.add(key)


        out.append(
            {
                "event":
                    "boekelheide_rearrangement",

                "role":
                    "specific_event",

                "subtype":
                    (
                        "aza_n_oxide_alpha_"
                        "oxygenation_rearrangement"
                    ),

                "center_map":
                    int(c),

                "nitrogen_map":
                    int(n),

                "oxygen_map":
                    int(o),

                "support_maps":
                    _support(
                        n,
                        c,
                        o,
                    ),

                "mapping_scope":
                    (
                        "aza_n_oxide_loss_plus_"
                        "single_mapped_c_o_formation_"
                        "at_n_adjacent_or_alpha_"
                        "carbon_environment"
                    ),
            }
        )

    return out





def detect_p6_azole_ring_atom_exchange(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """P6 scoped structural azole O->N vertex exchange (NOT a mechanism claim).

    Require one matched 5-membered aromatic O,N,N ring and N,N,N product
    with precisely four conserved ring maps and an incoming non-ring N.
    Ring matching, mapped-atom element conservation and the two original C--N
    insertion bonds are positive evidence.  Does not alter ring-lineage facts.
    """
    lineage = facts.get("ring_lineage_v2") or {}
    reactant_rings = lineage.get("reactant_rings") or []
    product_rings = lineage.get("product_rings") or []
    out: list[dict[str, Any]] = []
    seen: set[tuple[int, ...]] = set()

    def _ring_edges(ring: dict[str, Any], allowed: set[int]) -> set[tuple[int, int]]:
        return {
            tuple(sorted((int(edge[0]), int(edge[1]))))
            for edge in (ring.get("mapped_edges") or [])
            if len(edge) == 2 and int(edge[0]) in allowed and int(edge[1]) in allowed
        }

    for match in lineage.get("matches") or []:
        try:
            ri = int(match["reactant_ring_index"])
            pi = int(match["product_ring_index"])
            react = reactant_rings[ri]
            prod = product_rings[pi]
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        if int(react.get("size") or 0) != 5 or int(prod.get("size") or 0) != 5:
            continue
        if react.get("element_counts") != {"C": 2, "N": 2, "O": 1}:
            continue
        if prod.get("element_counts") != {"C": 2, "N": 3}:
            continue
        if react.get("kind") != "heterocycle" or prod.get("kind") != "heterocycle":
            continue
        r_maps = {int(x) for x in (react.get("maps") or []) if int(x) > 0}
        p_maps = {int(x) for x in (prod.get("maps") or []) if int(x) > 0}
        conserved = r_maps & p_maps
        incoming = p_maps - r_maps
        if len(r_maps) != 4 or len(p_maps) != 5 or len(conserved) != 4 or len(incoming) != 1:
            continue
        incoming_map = next(iter(incoming))
        if set(map(int, match.get("inserted_maps") or [])) != incoming:
            continue
        if match.get("lost_maps"):
            continue
        if set(map(int, match.get("conserved_maps") or [])) != conserved:
            continue
        if _ring_edges(react, conserved) != _ring_edges(prod, conserved):
            continue
        # Conserved ring atoms must be mapped to themselves, have identical
        # elements, and be aromatic on BOTH sides. The only lost ring atom
        # must be the unnumbered oxygen in the original ring.
        if any(
            _atom(facts, "reactant", m).get("element") not in {"C", "N"}
            or _atom(facts, "reactant", m).get("element") != _atom(facts, "product", m).get("element")
            or not _atom(facts, "reactant", m).get("aromatic")
            or not _atom(facts, "product", m).get("aromatic")
            for m in conserved
        ):
            continue
        if any(_atom(facts, "reactant", m).get("element") == "O" for m in r_maps):
            continue
        if _atom(facts, "reactant", incoming_map).get("element") != "N":
            continue
        if _atom(facts, "product", incoming_map).get("element") != "N":
            continue
        if not _atom(facts, "product", incoming_map).get("aromatic"):
            continue
        # Proven incoming N must not have belonged to ANY previous ring,
        # even if the previous ring is unrelated to the matched azole.
        if any(incoming_map in {int(x) for x in row.get("maps", [])} for row in reactant_rings):
            continue
        # Replaced oxygen occupies a carbon--oxygen--carbon position. Require
        # corresponding carbon--nitrogen--carbon product links on mapped atoms.
        nbrs = {
            other
            for a, b in _ring_edges(prod, p_maps)
            for other in ((b,) if a == incoming_map else (a,) if b == incoming_map else ())
        }
        if len(nbrs) != 2 or not nbrs <= conserved:
            continue
        if any(_atom(facts, "reactant", m).get("element") != "C" for m in nbrs):
            continue
        ring_key = tuple(sorted(p_maps))
        if ring_key in seen:
            continue
        seen.add(ring_key)
        out.append({
            "event": "heteroaromatic_ring_atom_exchange",
            "role": "specific_event",
            "subtype": "oxadiazole_O_to_triazole_N",
            "ring_size": 5,
            "ring_maps": sorted(p_maps),
            "conserved_ring_maps": sorted(conserved),
            "incoming_map": incoming_map,
            "replacement_neighbor_maps": sorted(nbrs),
            "support_maps": sorted(p_maps),
            "mapping_scope": "matched_five_member_azole_ring_unmapped_O_to_external_mapped_N",
        })
    return out


def run_detectors(facts: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    return {
        "substitution": detect_substitution(facts),
        "rearrangement": detect_rearrangement(facts),
        "fragmentation": detect_fragmentation(facts),
        "coupling": detect_coupling(facts),
        "ring_opening": detect_ring_opening(facts),
        "ring_formation": detect_ring_formation(facts),
        "ring_size_change": detect_ring_size_change(facts),
        "ring_rewiring": detect_ring_rewiring(facts),
        "p6_azole_ring_atom_exchange": detect_p6_azole_ring_atom_exchange(facts),
        "elimination": detect_elimination(facts),
        "condensation": detect_condensation(facts),
        "pericyclic": detect_pericyclic(facts),
        "functional_group_conversion": detect_functional_group_conversions(facts),
        "silicon_transformations": detect_silicon_transformations(facts),
        "sub05_o_silylation_product_state": detect_sub05_o_silylation_product_state(facts),
        "s1_semantic_closure": detect_s1_semantic_closure(facts),
        "addition": detect_addition_events(facts),
        "conjugate_acyl_addition": detect_conjugate_acyl_addition(facts),
        "acyl_conversion": detect_acyl_conversions(facts),
        "carbonyl_olefination": detect_carbonyl_olefination(facts),
        "oxidative_structural_conversion": detect_oxidative_structural_conversions(
            facts
        ),
        "aromatic_functionalization": detect_aromatic_functionalization(facts),
        "pi_reorganization": detect_pi_reorganizations(facts),
        "isothiocyanate_formation": detect_isothiocyanate_formation(facts),
        "s2_3a_heterocumulene_closure": detect_s23a_heterocumulene_closure(facts),
        "s2_3b_formamidine_formation": detect_s23b_formamidine_formation(facts),
        "s2_3d_alpha_hydroxyimino_amide": detect_s23d_alpha_hydroxyimino_amide(facts),
        "alpha_oximation": detect_alpha_oximation(facts),
        "alkene_sulfenyl_halogenation": detect_alkene_sulfenyl_halogenation(facts),
        "unmapped_lg_o_alkylation": detect_unmapped_lg_o_alkylation(facts),
        "unmapped_lg_n_alkylation": detect_unmapped_lg_n_cyclopropylation(facts),
        "s2_3e1_product_centered_n_methylation": detect_s23e1_product_centered_n_methylation(
            facts
        ),
        "s2_3e2_carboxylic_acid_to_aminomethylene": detect_s23e2_carboxylic_acid_to_aminomethylene(
            facts
        ),
        "s2_3e3_product_centered_urea_formation": detect_s23e3_product_centered_urea_formation(
            facts
        ),
        "s2_4_pi_system_closure": detect_s24_pi_system_closure(facts),
        "s2_5_directional_fgi": detect_s25_directional_fgi(facts),
        "post_v655_n_thp": detect_n_tetrahydropyranyl_protection(facts),
        "post_v655_disulfide": detect_post_v655_disulfide_formation(facts),
        "post_v655_n_o_acetal": detect_post_v655_n_o_acetal_formation(facts),
        "post_v655_nitro_hydroxylamine": detect_post_v655_nitro_to_hydroxylamine(facts),
        # POST_V655_FINAL_FIVE_SEMANTIC_FAMILIES
        "post_v655_arylboron_o_coupling": detect_post_v655_arylboron_o_coupling(facts),
        "post_v655_arylboronic_to_phenol": detect_post_v655_arylboronic_acid_to_phenol(facts),
        "post_v655_o_difluoromethylation": detect_post_v655_o_difluoromethylation(facts),
        "post_v655_nitrate_ester": detect_post_v655_nitrate_ester_formation(facts),
        "post_v655_boekelheide": detect_post_v655_boekelheide_rearrangement(facts),
        "n_glycosylation": detect_n_glycosylation(facts),
        "vinylic_amination": detect_vinylic_amination(facts),
        "terminal_alkyne_deprotection": detect_terminal_alkyne_deprotection(facts),
        "s_deprotection": detect_s_deprotection(facts),
        "o_heteroatom_deprotection": detect_o_heteroatom_deprotection(facts),
        "n_deprotection": detect_n_deprotection(facts),
        "n_n_bond_formation": detect_n_n_bond_formation(facts),
        "expanded_heterocumulene_addition": detect_expanded_heterocumulene_addition(
            facts
        ),
        "carbonyl_deprotection": detect_carbonyl_deprotection(facts),
        "heteroaryl_alkoxy_to_oxo": detect_heteroaryl_alkoxy_to_oxo(facts),
        "aza_michael_addition": detect_aza_michael_addition(facts),
        "imine_reduction": detect_imine_reduction(facts),
        "v67_rearrangement_extensions": detect_v67_rearrangement_extensions(facts),
        "v67_carbon_to_nitrogen_migration": detect_v67_carbon_to_nitrogen_migration(facts),
        "hofmann_rearrangement": detect_hofmann_rearrangement(facts),
        "nitrile_thioamide_and_ester": detect_nitrile_thioamide_and_ester(facts),
        "guanylation": detect_guanylation(facts),
        "chan_lam_n_arylation": detect_chan_lam_n_arylation(facts),
        "n_phosphorus_functionalization": detect_n_phosphorus_functionalization(facts),
        "o_phosphorylation": detect_o_phosphorylation(facts),
        "organozinc_c_s_coupling": detect_organozinc_c_s_coupling(facts),
        "hydrodealkoxylation": detect_hydrodealkoxylation(facts),
    }
