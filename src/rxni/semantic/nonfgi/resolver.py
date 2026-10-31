from __future__ import annotations

from collections import defaultdict
from typing import Any

from .mapped_graph import build_mapped_graph_facts
from .detectors import run_detectors, _attachments, _atom
from .rearrangement_extensions import RRG_V67_EVENTS, resolve_v67_rearrangement_primary
from .event_topology import (
    EVENT_TOPOLOGY_VERSION,
    SEMANTIC_EVENT_VERSION,
    build_semantic_events,
    whole_reaction_event_topology,
)

SYSTEM_I_VERSION = 'I_v0.1-v6.7-ring-v02'
RULEPACK_VERSION = 'i_rules_0.1-v6.7-ring-v02'
GRAPH_VERSION = "i_graph_0.3-ring-lineage-v2"

GENERIC_DEMOTIONS = {
    # H-derived concepts
    "intercomponent_bond_formation": "primitive",
    "coupling_candidate": "topology_evidence",
    "heteroatom_coupling_candidate": "topology_evidence",
    "c_c_coupling": "formed_bond_or_composite_descriptor",
    "substitution_like_topology": "topology_evidence",
    "same_center_replacement": "local_event",
    "halide_substitution_candidate": "composite_evidence",
    "carbon_hydrogen_replacement": "local_event",
    "fragment_loss_to_h": "local_event",
    "generic_bond_cleavage": "local_event",
    "generic_bond_formation": "primitive",
    "bond_order_reorganization": "primitive",
    "mapped_attachment_relocation": "topology_evidence_until_positive_migration",
    "rearrangement_like": "topology_evidence_until_positive_migration",
    "redox": "derived_orthogonal_axis_from_validated_conversion",
    "reduction": "derived_orthogonal_axis_from_validated_conversion",
    "oxidation": "derived_orthogonal_axis_from_validated_conversion",
    "functional_group_interconversion": "directional_conversion_parent_not_primary",
    "elimination_like": "mapped_topology_evidence_only",
    "condensation_like_topology": "mapped_topology_evidence_only",
    "pericyclic_compatible": "formal_topology_axis_requires_recognized_edit_pattern",
    "addition_to_pi_system": "topology_parent_not_primary_without_specific_mapped_addition",
    "addition_like_topology": "mapped_pi_consumption_evidence_only",
    # E-derived concepts
    "coupling": "broad_candidate_not_primary_by_default",
    "c_n_coupling": "formed_bond_descriptor_not_family",
    "c_o_coupling": "formed_bond_descriptor_not_family",
    "c_s_coupling": "formed_bond_descriptor_not_family",
    "ring_opening": "validated_generic_ring_topology_fallback",
    "ring_formation": "validated_generic_ring_topology_fallback",
    "heterocycle_formation": "structural_ring_family_below_specific_mechanism_or_ring_identity",
    "rearrangement": "candidate_until_positive_migration_validated",
    "fragmentation": "candidate_until_component_split_or_specific_loss",
    "bond_order_change": "local_structural_event",
    "mapping_rewire_candidate": "topology_or_mapping_qc_evidence",
    "ring_rewiring_candidate": "topology_or_mapping_qc_evidence",
}

SPECIFIC_PRIORITY = {
    # POST_V67_RING_SEMANTICS_V02
    # Specific chemistry / specific ring identity outranks generic topology.
    "epoxide_rearrangement": 106,
    "alkene_to_epoxide": 105,
    "epoxide_formation": 104,
    "epoxide_opening": 104,
    "aziridine_formation": 104,
    "aziridine_opening": 104,
    "thiirane_formation": 104,
    "thiirane_opening": 104,
    "beta_lactam_formation": 104,
    "beta_lactam_opening": 104,
    "cyclopropane_formation": 103,
    "cyclopropane_opening": 103,
    "oxetane_formation": 103,
    "oxetane_opening": 103,
    "azetidine_formation": 103,
    "azetidine_opening": 103,
    "cyclobutane_formation": 103,
    "cyclobutane_opening": 103,
    "ring_expansion": 96,
    "ring_contraction": 96,
    "ring_rewiring": 96,
    "heteroaromatic_ring_atom_exchange": 99,
    "macrocyclization": 92,
    "spirocyclization": 92,
    "heterocycle_annulation": 90,
    "carbocycle_annulation": 90,
    "intramolecular_heterocyclization": 91,
    "heterocycle_formation": 88,
    "ring_opening": 82,
    "ring_formation": 82,
    "n_glycosylation": 99,
    "n_tetrahydropyranyl_protection": 98,
    # POST_V655_FINAL_RESIDUAL_CLOSURE
    "disulfide_formation": 98,
    "n_o_acetal_formation": 98,
    "nitro_to_hydroxylamine": 98,
    # POST_V655_FINAL_FIVE_SEMANTIC_FAMILIES
    "boekelheide_rearrangement": 99,
    "arylboron_o_coupling": 98,
    "arylboronic_acid_to_phenol": 98,
    "o_difluoromethylation": 98,
    "nitrate_ester_formation": 98,
    "carbon_to_nitrogen_migration": 99,
    **{name: 99 for name in RRG_V67_EVENTS},
    "hofmann_rearrangement": 99,
    "carbonyl_deprotection": 98,
    "alpha_hydroxyimino_amide_formation": 99,
    "alpha_oximation": 98,
    "alkene_sulfenyl_halogenation": 98,
    "aza_michael_addition": 98,
    "conjugate_acyl_addition": 98,
    "chan_lam_n_arylation": 98,
    "guanylation": 97,
    "formamidine_formation": 97,
    "thiourea_formation": 97,
    "urea_formation": 97,
    "nitrile_to_thioamide": 97,
    "vinylic_amination": 96,
    "n_n_bond_formation": 96,
    "n_deprotection": 95,
    "s_deprotection": 95,
    "terminal_alkyne_deprotection": 95,
    "imine_reduction": 95,
    "n_phosphorus_functionalization": 92,
    "o_phosphorylation": 92,
    # SUB02-C product-verified C--P transformations; no mechanism assignment.
    "alkylphosphonate_formation": 96,
    "phosphorus_ylide_formation": 96,
    "nitrile_to_ester": 92,
    "organozinc_c_s_coupling": 92,
    "hydrodealkoxylation": 91,
    "intramolecular_anhydride_formation": 99,
    "anhydride_formation": 99,
    "carbonyl_olefination": 99,
    "lactamization": 98,
    "lactonization": 98,
    "weinreb_amide_to_ketone": 97,
    "aromatic_acylation": 97,
    "carbonyl_alpha_arylation": 97,
    "aryl_alkene_hydroarylation": 97,
    "vinyl_alkynyl_coupling": 97,
    "suzuki_type_cross_coupling": 97,
    "carbamate_formation": 97,
    "acyl_halide_formation": 96,
    "amidation": 96,
    "transamidation": 96,
    "transesterification": 96,
    "pericyclic_compatible": 98,
    "reductive_amination": 98,
    "carbonyl_addition": 97,
    "halocyclization": 97,
    "elimination": 96,
    "heterocumulene_addition": 96,
    "isocyanate_formation": 96,
    "isothiocyanate_formation": 96,
    "heteroarene_partial_hydrogenation": 96,
    "double_bond_migration": 95,
    "pi_heteroatom_hydroaddition": 96,
    "oxidative_esterification": 95,
    "c_h_to_aldehyde_oxidation": 95,
    "alkene_oxidative_cleavage_to_aldehyde": 97,
    "heteroaryl_hydrolysis_to_oxo": 94,
    "carbonyl_alpha_beta_dehydrogenation": 95,
    "dehydrogenation": 94,
    "alkene_hydrogenation": 94,
    "partial_alkyne_hydrogenation": 94,
    "complete_alkyne_hydrogenation": 94,
    "esterification": 94,
    "nitrile_hydration": 94,
    "aromatic_functionalization": 95,
    "cross_coupling_compatible": 95,
    "cyanation_or_nitrile_installation": 94,
    "azide_substitution": 94,
    "alkene_hydrosilylation": 94,
    "alkene_hydrohalogenation": 94,
    "c_h_halogenation": 93,
    "condensation": 93,
    "oxidative_deformylation": 93,
    "quinone_formation": 93,
    "oxidative_aromatization": 93,
    "alkene_cleavage_to_alcohol": 93,
    "etherification_or_o_alkylation": 89,
    "carbonyl_substitution": 92,
    "sulfonyl_substitution": 92,
    "sulfonamide_formation": 97,
    "sulfonate_ester_formation": 97,
    "decarboxylation": 92,
    "ester_hydrolysis_or_deesterification": 91,
    "o_deacylation_or_ester_cleavage": 90,
    "o_dealkylation": 90,
    "n_desulfonylation": 92,
    "o_desulfonylation": 92,
    "o_dephosphorylation": 92,
    "n_deformylation": 92,
    "dehalogenation": 90,
    "desulfurization": 90,
    "carbonyl_deoxygenation_or_reduction": 90,
    "n_oxide_or_nitroso_reduction": 90,
    "borylation": 88,
    # S2_1_SILICON_TRANSFORMATIONS
    "desilylative_halogenation": 96,
    "silyl_enol_ether_to_carbonyl": 95,
    "silyl_ether_to_carbonyl": 95,
    "desilylation": 94,
    "silylation": 88,
    "halogenation_or_halide_exchange": 88,
    "aryl_or_heteroaryl_amination": 90,
    "n_alkylation": 89,
    "o_arylation": 90,
    "o_alkylation": 89,
    "c_n_substitution": 86,
    "c_o_substitution": 86,
    "c_s_substitution": 86,
    # SUB02-C2: mapped C(sp3)--X + external donor C--H => product C(sp3)--C.
    "carbon_alkylation": 89,
    "substitution": 80,
    # S1 named-reaction semantic-closure additions
    "alkene_aminohydroxylation": 97,
    "aminophosphonate_multicomponent_addition": 97,
    "nitrile_to_amidine": 96,
    "nitrile_to_imidate": 96,
    "nitrile_to_orthoester": 96,
    "primary_amine_to_isocyanide": 96,
    "alpha_amino_acid_to_aldehyde": 95,
    "amine_to_diazonium": 95,
    "alkene_to_ketone_oxidation": 95,
    "allylic_oxidative_acyloxylation": 95,
    "acyl_halide_to_ketene": 95,
    "carbonyl_to_alkyne": 95,
    "organosilane_to_alcohol": 95,
    "carbon_oxygen_to_hydrogen_deoxygenation": 94,
    "alkyne_hydration_to_ketone": 94,
    "alkene_hydration": 94,
    "arene_partial_hydrogenation": 94,
    "carboxylic_acid_to_nitrile": 94,
    "amide_to_aldehyde": 94,
    "carboxylic_acid_to_aldehyde": 94,
    "n_dealkylation": 94,
    # S2_5_DIRECTIONAL_FGI
    "alkene_to_alcohol_oxygenation": 96,
    "aryl_amine_deamination": 96,
    "oxime_to_nitrile": 96,
    "enol_ether_to_carbonyl": 97,
    "rearrangement": 78,
    "fragmentation": 78,
    "multi_bond_fragmentation": 78,
}


def _flatten(detectors: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    out = []
    for detector, rows in detectors.items():
        for rec in rows:
            x = dict(rec)
            x["detector"] = detector
            out.append(x)
    return out


def _dedupe_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Merge repeated semantic records while preserving multiplicity/support.

    Batch04 intentionally collapses duplicate `ring_formation; ring_formation` style
    outputs. Distinct mapped centers/endpoints stay distinct; uncentered duplicate
    semantics are merged with a union of support maps and `multiplicity`.
    """
    merged: dict[tuple[Any, ...], dict[str, Any]] = {}
    order: list[tuple[Any, ...]] = []
    for r in records:
        if r.get("role") == "axis_event":
            key = (
                r.get("event"),
                r.get("role"),
                r.get("center_map"),
                tuple(sorted(r.get("support_maps") or [])),
            )
        else:
            key = (
                r.get("event"),
                r.get("role"),
                r.get("subtype"),
                r.get("center_map"),
                r.get("alpha_map"),
                r.get("beta_map"),
                r.get("donor_map"),
                r.get("incoming_map"),
                r.get("a"),
                r.get("b"),
                r.get("formed_bond_class"),
                r.get("directional_conversion"),
            )
        if key not in merged:
            x = dict(r)
            x["support_maps"] = sorted(set(x.get("support_maps") or []))
            x["multiplicity"] = int(x.get("multiplicity", 1))
            merged[key] = x
            order.append(key)
            continue
        x = merged[key]
        x["support_maps"] = sorted(
            set(x.get("support_maps") or []) | set(r.get("support_maps") or [])
        )
        if r.get("role") != "axis_event":
            x["multiplicity"] = int(x.get("multiplicity", 1)) + int(
                r.get("multiplicity", 1)
            )
    return [merged[k] for k in order]



def _resolve_sub01_dual_halide_cc_ownership(
    facts: dict[str, Any], chosen: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """SUB-01: primary-level granularity for exact dual-halide C-C replacement.

    Keep both specific/local detector observations unchanged; only when both
    endpoints of the *same* new intercomponent mapped C-C bond independently
    have a generic halide-to-carbon primary, replace that pair by ONE neutral,
    exact-bond-anchored primary. Never infer a coupling mechanism.

    The attachment deltas, rather than pooled topology support maps, establish
    both leaving groups and the identity of each incoming carbon. An unrelated
    chosen owner touching either endpoint blocks the replacement.
    """
    halogens = {"F", "Cl", "Br", "I"}
    deltas = {
        int(delta["center_map"]): delta
        for delta in facts.get("attachment_deltas", [])
    }
    consumed: set[int] = set()
    neutral: list[dict[str, Any]] = []

    for bond in facts.get("intercomponent_formed_bonds", []):
        a, b = sorted((int(bond["a"]), int(bond["b"])))
        if any(
            facts[side]["atoms"].get(center, {}).get("element") != "C"
            for side in ("reactant", "product")
            for center in (a, b)
        ):
            continue

        matched: list[int] = []
        for center, partner in ((a, b), (b, a)):
            delta = deltas.get(center)
            if not delta or len(delta.get("lost", [])) != 1 or len(delta.get("gained", [])) != 1:
                break
            lost, gained = delta["lost"][0], delta["gained"][0]
            halide = str(lost.get("element") or "")
            if (
                halide not in halogens
                or gained.get("element") != "C"
                or gained.get("map") is None
                or int(gained["map"]) != partner
            ):
                break

            candidates = [
                i
                for i, rec in enumerate(chosen)
                if i not in consumed
                and rec.get("event") == "substitution"
                and rec.get("role") == "specific_event"
                and rec.get("center_map") == center
                and rec.get("subtype") == f"C_{halide}_to_C"
                and partner in (rec.get("support_maps") or [])
            ]
            if len(candidates) != 1:
                break
            matched.append(candidates[0])

        if len(matched) != 2 or matched[0] == matched[1]:
            continue

        # Never demote an independently selected specific/composite owner.
        def _support(rec: dict[str, Any]) -> set[int]:
            maps = {int(x) for x in (rec.get("support_maps") or [])}
            if rec.get("center_map") is not None:
                maps.add(int(rec["center_map"]))
            return maps

        if any(
            i not in matched
            and rec.get("event") != "substitution"
            and bool(_support(rec) & {a, b})
            for i, rec in enumerate(chosen)
        ):
            continue

        consumed.update(matched)
        neutral.append(
            {
                "event": "c_c_bond_formation_unresolved",
                "role": "topology_evidence",
                "subtype": "paired_dual_halide_cc_formal",
                "a": a,
                "b": b,
                "support_maps": [a, b],
                "mapping_scope": "exact_intercomponent_cc_dual_halide_replacement",
            }
        )

    if not neutral:
        return chosen
    return [rec for i, rec in enumerate(chosen) if i not in consumed] + neutral


def _ring_primary_support(rec: dict[str, Any]) -> set[int]:
    maps = {int(x) for x in (rec.get("support_maps") or []) if x is not None}
    for key in ("ring_maps", "reactant_ring_maps", "product_ring_maps"):
        maps.update(int(x) for x in (rec.get(key) or []) if x is not None)
    if rec.get("center_map") is not None:
        maps.add(int(rec["center_map"]))
    return maps


def _resolve_v67_ring_primary_ownership(chosen: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Resolve specialized ring semantics against generic/local overlapping events.

    Detector evidence is preserved elsewhere.  This function only decides primary
    ownership after all event families have had a chance to enter ``chosen``.
    """
    if not chosen:
        return chosen

    by_event: dict[str, list[set[int]]] = defaultdict(list)
    for rec in chosen:
        by_event[str(rec.get("event"))].append(_ring_primary_support(rec))

    def overlaps_any(rec: dict[str, Any], owner_events: set[str]) -> bool:
        support = _ring_primary_support(rec)
        if not support:
            return False
        return any(
            bool(support & owner_support)
            for event in owner_events
            for owner_support in by_event.get(event, [])
        )

    small_ring = {
        "epoxide_formation", "epoxide_opening",
        "aziridine_formation", "aziridine_opening",
        "thiirane_formation", "thiirane_opening",
        "cyclopropane_formation", "cyclopropane_opening",
        "oxetane_formation", "oxetane_opening",
        "azetidine_formation", "azetidine_opening",
        "beta_lactam_formation", "beta_lactam_opening",
        "cyclobutane_formation", "cyclobutane_opening",
    }
    topology_specific = {
        "ring_expansion", "ring_contraction", "ring_rewiring",
        "macrocyclization", "spirocyclization",
        "heterocycle_annulation", "carbocycle_annulation",
        "intramolecular_heterocyclization", "heterocycle_formation",
    }
    all_ring_owners = small_ring | topology_specific | {"epoxide_rearrangement"}

    # v0.3 inserted-ring-atom closure.  When a newly formed heterocycle is proven
    # specifically by a single product-only ring atom, bond-order aromatization of
    # the atoms inside that new ring is subordinate bookkeeping rather than the
    # synthetic owner (Paal-Knorr-like closure).  Keep the aromatization record in
    # specific-event evidence; suppress only primary ownership.
    inserted_heterocycle_supports = [
        _ring_primary_support(rec)
        for rec in chosen
        if rec.get("event") == "heterocycle_formation"
        and rec.get("subtype") == "single_inserted_ring_atom_heterocycle_formation"
    ]

    weak_local = {
        "substitution", "c_n_substitution", "c_o_substitution", "c_s_substitution",
        "aryl_or_heteroaryl_amination", "n_alkylation", "o_arylation", "o_alkylation",
        "etherification_or_o_alkylation", "halogenation_or_halide_exchange",
        "o_dealkylation",
    }
    small_ring_local = weak_local | {"carbonyl_addition"}
    architecture_local = weak_local | {"carbonyl_substitution"}
    size_rewire_competitors = weak_local | {
        "rearrangement", "carbonyl_substitution", "aromatic_functionalization",
        "cross_coupling_compatible", "c_h_halogenation",
        "arene_partial_hydrogenation", "heteroarene_partial_hydrogenation",
    }

    out: list[dict[str, Any]] = []
    for rec in chosen:
        event = str(rec.get("event"))

        # A graph-proven single-inserted-atom heterocycle owns internal
        # aromatization bookkeeping on the same new ring.
        if event == "oxidative_aromatization" and inserted_heterocycle_supports:
            support = _ring_primary_support(rec)
            if support and any(support <= owner for owner in inserted_heterocycle_supports):
                continue

        # Exact alkene epoxidation is more method-informative than the generic
        # product-ring identity, while epoxide_formation remains in specific_events.
        if event == "epoxide_formation" and overlaps_any(rec, {"alkene_to_epoxide"}):
            continue

        # Payne-type epoxide migration owns both local opening/formation bookkeeping
        # and incidental O-dealkylation/generic rearrangement on the same support.
        if event in {
            "epoxide_formation", "epoxide_opening", "ring_formation", "ring_opening",
            "o_dealkylation", "rearrangement",
        } and overlaps_any(rec, {"epoxide_rearrangement"}):
            continue

        # Specific small-ring identity owns generic ring/heterocycle bookkeeping and
        # weak local substitutions/additions on the same mapped center.
        if event in (
            small_ring_local
            | {"ring_formation", "ring_opening", "heterocycle_formation",
               "heterocycle_annulation", "intramolecular_heterocyclization"}
        ) and overlaps_any(rec, small_ring):
            continue

        # Beta-lactam identity is more specific than generic lactamization.
        if event == "lactamization" and overlaps_any(rec, {"beta_lactam_formation"}):
            continue

        # Validated annulation/cyclization architecture owns broad local substitution/acyl
        # bookkeeping on the same ring support, while mechanism-level condensation,
        # pericyclic, redox and named-specific events remain independent.
        if event in architecture_local and overlaps_any(
            rec, {
                "macrocyclization", "spirocyclization", "heterocycle_annulation",
                "carbocycle_annulation", "intramolecular_heterocyclization",
                "heterocycle_formation",
            }
        ):
            continue

        # Ring-size change / validated rewiring absorbs generic open+form decomposition
        # and generic local descriptors on the same lineage.
        if event in (size_rewire_competitors | {"ring_formation", "ring_opening"}) and overlaps_any(
            rec, {"ring_expansion", "ring_contraction", "ring_rewiring"}
        ):
            continue

        # Any validated specific ring owner prevents incidental O-dealkylation from
        # becoming an independent primary interpretation on the same ring support.
        if event == "o_dealkylation" and overlaps_any(rec, all_ring_owners):
            continue

        out.append(rec)

    return out




def _resolve_p6_azole_exchange_ownership(chosen: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Remove only local O->N amination double-labels at matched azole neighbors.

    Independent cyclization at a different mapped edit edge MUST survive, even
    if the second new ring shares ring atoms with the azole replacement.
    """
    azoles = [r for r in chosen if r.get("event") == "heteroaromatic_ring_atom_exchange"
              and r.get("subtype") == "oxadiazole_O_to_triazole_N"]
    if not azoles:
        return chosen
    result: list[dict[str, Any]] = []
    for rec in chosen:
        suppress = False
        if (rec.get("event") == "aryl_or_heteroaryl_amination"
                and rec.get("subtype") == "C_O_to_N"
                and rec.get("center_map") is not None):
            center = int(rec["center_map"])
            support = {int(x) for x in (rec.get("support_maps") or [])}
            for owner in azoles:
                ring = {int(x) for x in owner["ring_maps"]}
                neighbors = {int(x) for x in owner["replacement_neighbor_maps"]}
                incoming = int(owner["incoming_map"])
                if (center in neighbors and support == {center, incoming}
                        and support <= ring):
                    suppress = True
                    break
        # P6 v0.4: the matched-ring rewiring uses product_ring_maps and
        # reactant_ring_maps, not ring_maps.  Require BOTH exact correspondences
        # to avoid absorbing independent ring rewiring elsewhere in the molecule.
        if (not suppress and rec.get("event") == "ring_rewiring"
                and rec.get("mapping_scope") == "full_ring_lineage_v2_same_size_membership_change"
                and rec.get("subtype") == "validated_same_size_ring_lineage_change"
                and rec.get("product_ring_maps") and rec.get("reactant_ring_maps")):
            product_maps = {int(x) for x in rec["product_ring_maps"]}
            reactant_maps = {int(x) for x in rec["reactant_ring_maps"]}
            for owner in azoles:
                product_azole = {int(x) for x in owner["ring_maps"]}
                conserved_core = {int(x) for x in owner["conserved_ring_maps"]}
                if product_maps == product_azole and reactant_maps == conserved_core:
                    suppress = True
                    break
        if not suppress:
            result.append(rec)
    return result

def _primary_events(
    facts: dict[str, Any], records: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    specific = [r for r in records if r.get("role") == "specific_event"]
    events = defaultdict(list)
    for r in specific:
        events[r["event"]].append(r)

    chosen: list[dict[str, Any]] = []

    # High-information mapped structural events first.
    # POST_V67_RING_SEMANTICS_V02: specialized ring events precede generic topology.
    for name in (
        "epoxide_rearrangement",
        "alkene_to_epoxide",
        "epoxide_formation",
        "epoxide_opening",
        "aziridine_formation",
        "aziridine_opening",
        "thiirane_formation",
        "thiirane_opening",
        "cyclopropane_formation",
        "cyclopropane_opening",
        "oxetane_formation",
        "oxetane_opening",
        "azetidine_formation",
        "azetidine_opening",
        "beta_lactam_formation",
        "beta_lactam_opening",
        "cyclobutane_formation",
        "cyclobutane_opening",
        "ring_expansion",
        "ring_contraction",
        "ring_rewiring",
        "heteroaromatic_ring_atom_exchange",
        "macrocyclization",
        "spirocyclization",
        "heterocycle_annulation",
        "carbocycle_annulation",
        "intramolecular_heterocyclization",
        "heterocycle_formation",
        "ring_opening",
        "ring_formation",
        "n_glycosylation",
        "n_tetrahydropyranyl_protection",
        # POST_V655_FINAL_RESIDUAL_CLOSURE
        "disulfide_formation",
        "n_o_acetal_formation",
        "nitro_to_hydroxylamine",
        # POST_V655_FINAL_FIVE_SEMANTIC_FAMILIES
        "boekelheide_rearrangement",
        "arylboron_o_coupling",
        "arylboronic_acid_to_phenol",
        "o_difluoromethylation",
        "nitrate_ester_formation",
        "carbon_to_nitrogen_migration",
        "hofmann_rearrangement",
        "carbonyl_deprotection",
        "heteroaryl_hydrolysis_to_oxo",
        "alpha_hydroxyimino_amide_formation",
        "alpha_oximation",
        "alkene_sulfenyl_halogenation",
        "aza_michael_addition",
        "conjugate_acyl_addition",
        "chan_lam_n_arylation",
        "guanylation",
        "formamidine_formation",
        "thiourea_formation",
        "urea_formation",
        "nitrile_to_thioamide",
        "nitrile_to_ester",
        "vinylic_amination",
        "n_n_bond_formation",
        "n_deprotection",
        "s_deprotection",
        "terminal_alkyne_deprotection",
        "imine_reduction",
        "n_phosphorus_functionalization",
        "o_phosphorylation",
        "alkylphosphonate_formation",
        "phosphorus_ylide_formation",
        "organozinc_c_s_coupling",
        "hydrodealkoxylation",
        "intramolecular_anhydride_formation",
        "anhydride_formation",
        "carbonyl_olefination",
        "lactamization",
        "lactonization",
        "weinreb_amide_to_ketone",
        "aromatic_acylation",
        "carbonyl_alpha_arylation",
        "aryl_alkene_hydroarylation",
        "vinyl_alkynyl_coupling",
        "suzuki_type_cross_coupling",
        "amidation",
        "carbamate_formation",
        "acyl_halide_formation",
        "transamidation",
        "transesterification",
        "pericyclic_compatible",
        "reductive_amination",
        "carbonyl_addition",
        "halocyclization",
        "elimination",
        "heterocumulene_addition",
        "isocyanate_formation",
        "isothiocyanate_formation",
        "heteroarene_partial_hydrogenation",
        "arene_partial_hydrogenation",
        "double_bond_migration",
        "pi_heteroatom_hydroaddition",
        "alkene_to_ketone_oxidation",
        "alkyne_hydration_to_ketone",
        "alkene_hydration",
        "alkene_aminohydroxylation",
        "carbonyl_alpha_beta_dehydrogenation",
        "dehydrogenation",
        "alkene_hydrogenation",
        "partial_alkyne_hydrogenation",
        "complete_alkyne_hydrogenation",
        "esterification",
        "nitrile_hydration",
        "aromatic_functionalization",
        "cross_coupling_compatible",
        "cyanation_or_nitrile_installation",
        "azide_substitution",
        "alkene_hydrosilylation",
        "alkene_hydrohalogenation",
        "c_h_halogenation",
        "condensation",
        "oxidative_deformylation",
        "quinone_formation",
        "oxidative_aromatization",
        "alkene_cleavage_to_alcohol",
        "etherification_or_o_alkylation",
        "carbonyl_substitution",
        "sulfonyl_substitution",
        "sulfonamide_formation",
        "sulfonate_ester_formation",
    ):
        chosen.extend(events.get(name, []))

    # Carbonyl addition is not cross-coupling even when the incoming carbon also loses
    # an electrophile/donor handle (e.g. an organomagnesium generated from aryl bromide).
    carbonyl_add_maps = set()
    for r in events.get("carbonyl_addition", []):
        carbonyl_add_maps.update(r.get("support_maps") or [])
    if carbonyl_add_maps:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "cross_coupling_compatible"
                and bool(set(r.get("support_maps") or []) & carbonyl_add_maps)
            )
        ]

    # S2_1B_SILYL_ENOL_ALDOL_PRIMARY_GUARD
    #
    # In Mukaiyama-type carbonyl additions, restoration of the
    # silyl-enol-ether carbonyl is subordinate structural evidence,
    # not an independent defining primary event.
    #
    # The detector record remains available in specific_events.
    if events.get("carbonyl_addition"):
        chosen = [r for r in chosen if r.get("event") != "silyl_enol_ether_to_carbonyl"]

    # Carbonyl olefination outranks aldol/Knoevenagel-like condensation and generic
    # same-center replacement at the carbonyl carbon.
    olef_centers = {
        r.get("center_map")
        for r in events.get("carbonyl_olefination", [])
        if r.get("center_map") is not None
    }
    if olef_centers:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "condensation" and r.get("center_map") in olef_centers
            )
        ]

    acetal_centers = {
        r.get("center_map")
        for r in events.get("condensation", [])
        if r.get("subtype") == "acetal_or_ketal_condensation"
        and r.get("center_map") is not None
    }
    if acetal_centers:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "carbonyl_addition"
                and r.get("center_map") in acetal_centers
            )
        ]

    # Lactonization is the specific cyclic acid->ester event; keep the orthogonal
    # directional conversion but remove generic esterification at the same acyl center.
    lactone_centers = {
        int(r["center_map"])
        for r in chosen
        if r.get("event") == "lactonization" and r.get("center_map") is not None
    }
    if lactone_centers:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "esterification"
                and r.get("center_map") in lactone_centers
            )
        ]

    # POST_V67_RING_SEMANTICS_V02 -------------------------------------------------
    # Specific ring identity/architecture owns overlapping generic ring bookkeeping.
    # Preserve all detector records in specific_events/topology; this block controls only
    # primary ownership.
    ring_specific_owners = {
        "epoxide_rearrangement",
        "epoxide_formation", "epoxide_opening",
        "aziridine_formation", "aziridine_opening",
        "thiirane_formation", "thiirane_opening",
        "cyclopropane_formation", "cyclopropane_opening",
        "oxetane_formation", "oxetane_opening",
        "azetidine_formation", "azetidine_opening",
        "beta_lactam_formation", "beta_lactam_opening",
        "cyclobutane_formation", "cyclobutane_opening",
        "ring_expansion", "ring_contraction", "ring_rewiring",
        "macrocyclization", "spirocyclization",
        "heterocycle_annulation", "carbocycle_annulation",
        "intramolecular_heterocyclization", "heterocycle_formation",
    }
    ring_owner_maps = set()
    small_ring_owner_maps = set()
    size_or_rewire_maps = set()
    for r in chosen:
        if r.get("event") not in ring_specific_owners:
            continue
        support = set(r.get("support_maps") or []) | set(r.get("ring_maps") or [])
        support |= set(r.get("reactant_ring_maps") or []) | set(r.get("product_ring_maps") or [])
        ring_owner_maps.update(support)
        if r.get("event") in {
            "epoxide_rearrangement",
            "epoxide_formation", "epoxide_opening",
            "aziridine_formation", "aziridine_opening",
            "thiirane_formation", "thiirane_opening",
            "cyclopropane_formation", "cyclopropane_opening",
            "oxetane_formation", "oxetane_opening",
            "azetidine_formation", "azetidine_opening",
            "beta_lactam_formation", "beta_lactam_opening",
            "cyclobutane_formation", "cyclobutane_opening",
        }:
            small_ring_owner_maps.update(support)
        if r.get("event") in {"epoxide_rearrangement", "ring_expansion", "ring_contraction", "ring_rewiring"}:
            size_or_rewire_maps.update(support)

    if ring_owner_maps:
        generic_ring_events = {"ring_formation", "ring_opening"}
        chosen = [
            r for r in chosen
            if not (
                r.get("event") in generic_ring_events
                and bool((set(r.get("support_maps") or []) | set(r.get("ring_maps") or [])) & ring_owner_maps)
            )
        ]

    # A local O-C cleavage on the same mapped support is subordinate to a validated
    # specific ring construction/opening/reconstruction. Keep it in specific_events,
    # but do not let o_dealkylation own the reaction.
    if ring_owner_maps:
        chosen = [
            r for r in chosen
            if not (
                r.get("event") == "o_dealkylation"
                and (
                    bool(set(r.get("support_maps") or []) & ring_owner_maps)
                    or (r.get("center_map") is not None and int(r["center_map"]) in ring_owner_maps)
                )
            )
        ]

    # Small-ring identity is more informative than broad heterocycle formation/annulation.
    if small_ring_owner_maps:
        chosen = [
            r for r in chosen
            if not (
                r.get("event") in {"heterocycle_formation", "heterocycle_annulation", "intramolecular_heterocyclization"}
                and bool((set(r.get("support_maps") or []) | set(r.get("ring_maps") or [])) & small_ring_owner_maps)
            )
        ]

    # Expansion/contraction/rewiring absorb the old open+form decomposition on the same
    # lineage and prevent incidental O-dealkylation from owning a ring reconstruction.
    if size_or_rewire_maps:
        chosen = [
            r for r in chosen
            if not (
                r.get("event") in {"ring_opening", "ring_formation", "o_dealkylation"}
                and (
                    bool(set(r.get("support_maps") or []) & size_or_rewire_maps)
                    or (r.get("center_map") is not None and int(r["center_map"]) in size_or_rewire_maps)
                )
            )
        ]

    # Batch06: generic ring_formation is topological fallback when a more specific
    # heterocycle/lactam/lactone event explains the same conserved ring support.
    specific_ring_maps = set()
    for r in chosen:
        if r.get("event") in {
            "heterocycle_annulation",
            "intramolecular_heterocyclization",
            "lactamization",
            "lactonization",
        } or (r.get("event") == "carbamate_formation" and r.get("cyclic")):
            specific_ring_maps.update(r.get("support_maps") or [])
            specific_ring_maps.update(r.get("ring_maps") or [])
    if specific_ring_maps:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "ring_formation"
                and bool(set(r.get("support_maps") or []) & specific_ring_maps)
            )
        ]

    # Batch05: broad carbonyl_substitution is a fallback descriptor only. Validated
    # acyl conversions, annulations/condensations and aromatic acyl/formyl events own
    # the synthetic interpretation at the same mapped acyl center/support.
    acyl_authoritative = {
        "intramolecular_anhydride_formation",
        "anhydride_formation",
        "lactamization",
        "lactonization",
        "amidation",
        "azide_substitution",
        "carbamate_formation",
        "acyl_halide_formation",
        "transamidation",
        "transesterification",
        "weinreb_amide_to_ketone",
        "aromatic_acylation",
        "esterification",
        "condensation",
        "ring_formation",
        "heterocycle_annulation",
        "intramolecular_heterocyclization",
        "aromatic_functionalization",
    }
    acyl_centers = set()
    acyl_support = set()
    for r in chosen:
        if r.get("event") not in acyl_authoritative:
            continue
        if r.get("center_map") is not None:
            acyl_centers.add(int(r["center_map"]))
        acyl_support.update(r.get("support_maps") or [])
        acyl_support.update(r.get("acyl_center_maps") or [])
    if acyl_centers or acyl_support:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "carbonyl_substitution"
                and (
                    (
                        r.get("center_map") is not None
                        and int(r["center_map"]) in acyl_centers
                    )
                    or (
                        r.get("center_map") is not None
                        and int(r["center_map"]) in acyl_support
                    )
                    or bool(set(r.get("support_maps") or []) & acyl_support)
                )
            )
        ]

    # Specific cleavage/reduction/substitution families retained from earlier batches.
    for name in (
        "decarboxylation",
        "ester_hydrolysis_or_deesterification",
        "o_deacylation_or_ester_cleavage",
        "o_dealkylation",
        "o_dephosphorylation",
        "o_desulfonylation",
        "n_desulfonylation",
        "n_deformylation",
        "dehalogenation",
        "desulfurization",
        "carbonyl_deoxygenation_or_reduction",
        "n_oxide_or_nitroso_reduction",
        "desilylative_halogenation",
        "silyl_enol_ether_to_carbonyl",
        "silyl_ether_to_carbonyl",
        "desilylation",
        "borylation",
        "silylation",
        "halogenation_or_halide_exchange",
        "aryl_or_heteroaryl_amination",
        "n_alkylation",
        "carbon_alkylation",
        "o_arylation",
        "o_alkylation",
        "c_n_substitution",
        "c_o_substitution",
        "c_s_substitution",
    ):
        if name == "dehalogenation" and events.get("cyanation_or_nitrile_installation"):
            continue
        chosen.extend(events.get(name, []))

    # ------------------------------------------------------------------
    # POST_V655_N_THP_OWNERSHIP
    #
    # N-THP owns only weak generic C-N bookkeeping on the same
    # local C/N support. Independent centers remain independent.
    # ------------------------------------------------------------------

    thp_owners = [
        rec for rec in chosen if rec.get("event") == "n_tetrahydropyranyl_protection"
    ]

    if thp_owners:

        weak_cn = {
            "n_alkylation",
            "c_n_substitution",
            "substitution",
        }

        owner_pairs = [
            {
                int(rec["center_map"]),
                int(rec["nitrogen_map"]),
            }
            for rec in thp_owners
            if (
                rec.get("center_map") is not None
                and rec.get("nitrogen_map") is not None
            )
        ]

        filtered = []

        for rec in chosen:

            if rec.get("event") not in weak_cn:
                filtered.append(rec)
                continue

            local = {int(x) for x in (rec.get("support_maps") or []) if x is not None}

            if rec.get("center_map") is not None:
                local.add(int(rec["center_map"]))

            if any(pair <= local for pair in owner_pairs):
                continue

            filtered.append(rec)

        chosen = filtered

    # S2_2D2_RING_OWNERSHIP_GUARD
    #
    # N-H -> N-Si remains valid structural evidence when it participates
    # in formation of a silicon-containing heterocycle.  However, when
    # the complete local N/Si support is contained in the ring support
    # of a specific heterocycle construction, the annulation owns the
    # primary synthetic interpretation.
    #
    # Examples include:
    #
    #   NH2 + Cl-Si-...-Si-Cl
    #       -> one N center gains two Si attachments and closes a ring
    #
    #   NH-...-NH + SiCl3
    #       -> two N centers attach to the same Si and close a ring
    #
    # Keep the silylation specific_event record; suppress only its
    # promotion to primary_event on locally annulation-owned support.

    n_si_ring_owners = []

    for r in chosen:

        if r.get("event") not in {
            "heterocycle_annulation",
            "intramolecular_heterocyclization",
        }:
            continue

        ring_maps = {int(x) for x in (r.get("ring_maps") or []) if x is not None}

        if ring_maps:
            n_si_ring_owners.append(ring_maps)

    if n_si_ring_owners:

        filtered = []

        for r in chosen:

            suppress = False

            if r.get("event") == "silylation" and r.get("subtype") == "N_H_to_Si":

                support = {
                    int(x) for x in (r.get("support_maps") or []) if x is not None
                }

                if r.get("center_map") is not None:

                    support.add(int(r["center_map"]))

                # Require complete local support containment.
                # Mere overlap is deliberately insufficient.
                suppress = bool(support) and any(
                    support <= ring_maps for ring_maps in n_si_ring_owners
                )

            if not suppress:
                filtered.append(r)

        chosen = filtered

    # S2_2B1_N_DESULF_CANONICAL_OWNER
    #
    # A subtype-specific N-sulfonyl deprotection and the historical
    # bare n_desulfonylation record can describe the same N-S cleavage.
    #
    # Same local support + same event hierarchy:
    # keep the more specific child as the primary owner.
    #
    # Detector evidence is untouched; this filters only chosen primary
    # candidates.

    specific_n_desulf = [
        r
        for r in chosen
        if (
            r.get("event") == "n_desulfonylation"
            and r.get("subtype") == "n_sulfonyl_deprotection"
        )
    ]

    if specific_n_desulf:

        def _local_maps_s22b1(r):
            maps = {int(x) for x in (r.get("support_maps") or []) if x is not None}

            if r.get("center_map") is not None:
                maps.add(int(r["center_map"]))

            return maps

        filtered = []

        for r in chosen:

            suppress = False

            if (
                r.get("event") == "n_desulfonylation"
                and r.get("subtype") != "n_sulfonyl_deprotection"
            ):

                r_maps = _local_maps_s22b1(r)
                r_center = r.get("center_map")

                for s in specific_n_desulf:

                    s_maps = _local_maps_s22b1(s)
                    s_center = s.get("center_map")

                    same_center = (
                        r_center is not None
                        and s_center is not None
                        and int(r_center) == int(s_center)
                    )

                    local_overlap = bool(r_maps & s_maps)

                    if same_center and local_overlap:
                        suppress = True
                        break

            if not suppress:
                filtered.append(r)

        chosen = filtered

    # S2_2_FINAL_ANNULATIVE_ESTER_O_OWNERSHIP
    #
    # Some heterocycle annulations/aromatizations transform:
    #
    #     Me-O-C(=O)-...  ->  HO-c(aromatic)-...
    #
    # The mapped oxygen loses an unmapped methyl carbon and therefore
    # fragmentation legitimately emits o_dealkylation evidence.
    #
    # However, when:
    #
    #   * the same O remains attached to the same mapped carbon;
    #   * that carbon is carbonyl-like in the reactant;
    #   * that carbon becomes aromatic in the product;
    #   * the carbon belongs to a validated heterocycle-annulation ring;
    #   * oxidative_aromatization is also present;
    #
    # the O-C cleavage is embedded in the higher-order annulative
    # aromatization and must not become an independent primary event.
    #
    # Keep o_dealkylation as specific evidence; suppress only primary
    # ownership.

    annulative_ring_maps = set()

    for r in chosen:
        if r.get("event") in {
            "heterocycle_annulation",
            "intramolecular_heterocyclization",
        }:
            annulative_ring_maps.update(
                int(x) for x in (r.get("ring_maps") or []) if x is not None
            )
            annulative_ring_maps.update(
                int(x) for x in (r.get("support_maps") or []) if x is not None
            )

    has_oxidative_aromatization = any(
        r.get("event") == "oxidative_aromatization" for r in chosen
    )

    if annulative_ring_maps and has_oxidative_aromatization:

        filtered = []

        for r in chosen:

            suppress = False

            if r.get("event") == "o_dealkylation" and r.get("center_map") is not None:

                o = int(r["center_map"])

                reactant_acyl_neighbors = []

                for att in _attachments(
                    facts,
                    "reactant",
                    o,
                ):

                    if (
                        att.get("element") != "C"
                        or att.get("map") is None
                        or not bool(att.get("carbonyl_like"))
                    ):
                        continue

                    reactant_acyl_neighbors.append(int(att["map"]))

                for c in reactant_acyl_neighbors:

                    # Same mapped O-C anchor must survive.
                    retained_in_product = any(
                        (
                            att.get("element") == "C"
                            and att.get("map") is not None
                            and int(att["map"]) == c
                        )
                        for att in _attachments(
                            facts,
                            "product",
                            o,
                        )
                    )

                    if not retained_in_product:
                        continue

                    product_c = _atom(
                        facts,
                        "product",
                        c,
                    )

                    if c in annulative_ring_maps and bool(product_c.get("aromatic")):
                        suppress = True
                        break

            if not suppress:
                filtered.append(r)

        chosen = filtered

    # Batch07: generic N-deprotection is a parent fallback. Existing specific
    # deformylation/desulfonylation labels own overlapping support.
    n_specific_maps = set()
    for r in chosen:
        if r.get("event") in {"n_deformylation", "n_desulfonylation"}:
            n_specific_maps.update(r.get("support_maps") or [])
            if r.get("center_map") is not None:
                n_specific_maps.add(int(r["center_map"]))
    if n_specific_maps:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "n_deprotection"
                and (
                    (
                        r.get("center_map") is not None
                        and int(r["center_map"]) in n_specific_maps
                    )
                    or bool(set(r.get("support_maps") or []) & n_specific_maps)
                )
            )
        ]

    # ------------------------------------------------------------------
    # S2_5_DIRECTIONAL_FGI
    #
    # Exact S1 semantic detectors that already existed in the runtime
    # but were never admitted through a primary resolver path.
    #
    # Promote only when disjoint from stronger selected ownership.
    # Weak generic substitution bookkeeping is allowed to yield to the
    # exact event on overlapping support.
    # ------------------------------------------------------------------

    s2_5_existing_exact = (
        "n_dealkylation",
        "primary_amine_to_isocyanide",
        "carbonyl_to_alkyne",
        "allylic_oxidative_acyloxylation",
        "aminophosphonate_multicomponent_addition",
    )

    s2_5_weak_generic = {
        "substitution",
        "c_n_substitution",
        "c_o_substitution",
        "c_s_substitution",
        "aryl_or_heteroaryl_amination",
        "n_alkylation",
        "o_arylation",
        "o_alkylation",
        "etherification_or_o_alkylation",
        "halogenation_or_halide_exchange",
    }

    def _s2_5_maps(
        rec: dict[str, Any],
    ) -> set[int]:

        ans = {int(x) for x in (rec.get("support_maps") or []) if x is not None}

        center = rec.get("center_map")

        if center is not None:
            ans.add(int(center))

        return ans

    for event_name in s2_5_existing_exact:

        for candidate in (
            events.get(
                event_name,
                [],
            )
            or []
        ):

            cmap = _s2_5_maps(candidate)

            if not cmap:
                continue

            stronger_overlap = False

            for owner in chosen:

                omap = _s2_5_maps(owner)

                if not (cmap & omap):
                    continue

                owner_event = owner.get("event")

                # S2_5_OWNERSHIP_REFINEMENT
                #
                # An overlapping aryl/heteroaryl amination owns the
                # shared N center.  Disjoint N-dealkylation remains an
                # independent primary event.
                if (
                    event_name == "n_dealkylation"
                    and owner_event == "aryl_or_heteroaryl_amination"
                ):

                    stronger_overlap = True
                    break

                if owner_event not in s2_5_weak_generic:

                    stronger_overlap = True
                    break

            if stronger_overlap:
                continue

            # Exact child owns overlapping generic bookkeeping.
            chosen = [
                owner
                for owner in chosen
                if not (
                    owner.get("event") in s2_5_weak_generic
                    and bool(cmap & _s2_5_maps(owner))
                )
            ]

            # Preserve multiplicity: a second exact event is retained
            # when its support is disjoint from existing exact owners.
            exact_overlap = any(
                (
                    owner.get("event") in s2_5_existing_exact
                    and bool(cmap & _s2_5_maps(owner))
                )
                for owner in chosen
            )

            if not exact_overlap:

                chosen.append(candidate)

    # Exact mapped directional conversions are authoritative FGI/redox labels. Only
    # primary-eligible conversion records participate in primary resolution; reagent-side
    # conversions stay in directional_conversions/redox_axis as evidence.
    conversion_specific = [
        r
        for r in specific
        if r.get("conversion_specific") and r.get("primary_eligible", True)
    ]

    # S2_5_DIRECTIONAL_FGI exact-child ownership.

    enol_owners = [
        r
        for r in conversion_specific
        if (r.get("directional_conversion") or r.get("subtype"))
        == "enol_ether_to_carbonyl"
    ]

    alcohol_owners = [
        r
        for r in conversion_specific
        if (r.get("directional_conversion") or r.get("subtype")) == "alkene_to_alcohol"
    ]

    if enol_owners:

        owner_maps = [_s2_5_maps(r) for r in enol_owners]

        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "alkene_to_ketone_oxidation"
                and any(bool(_s2_5_maps(r) & maps) for maps in owner_maps)
            )
        ]

    if alcohol_owners:

        owner_maps = [_s2_5_maps(r) for r in alcohol_owners]

        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "pi_heteroatom_hydroaddition"
                and any(bool(_s2_5_maps(r) & maps) for maps in owner_maps)
            )
        ]

    # S2_2A_ESTER_REDUCTION_OWNS_O_DEALKYLATION
    #
    # In mapped ester reductions the alkoxy oxygen may lose an unmapped
    # carbon substituent and gain H while the adjacent ester carbonyl
    # becomes CH2OH.  Fragmentation can therefore emit o_dealkylation,
    # but ester_to_primary_alcohol is the authoritative exact
    # directional conversion.
    #
    # Suppress only overlapping O-dealkylation primary candidates.
    # Keep the underlying detector evidence.

    def _local_maps_s22a(r):
        maps = set()

        if r.get("center_map") is not None:
            maps.add(int(r["center_map"]))

        for key, value in r.items():

            if not key.endswith("_maps"):
                continue

            if not isinstance(
                value,
                (list, tuple, set),
            ):
                continue

            for x in value:
                if x is not None:
                    try:
                        maps.add(int(x))
                    except Exception:
                        pass

        return maps

    ester_reduction_owners = [
        r for r in conversion_specific if r.get("event") == "ester_to_primary_alcohol"
    ]

    if ester_reduction_owners:

        owner_maps = [_local_maps_s22a(r) for r in ester_reduction_owners]

        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "o_dealkylation"
                and any(_local_maps_s22a(r) & maps for maps in owner_maps)
            )
        ]

    # Exact conversions replace broad reduction surrogates on the same local center.
    broad_conversion_surrogates = {
        "dehalogenation",
        "carbonyl_deoxygenation_or_reduction",
        "n_oxide_or_nitroso_reduction",
    }
    for conv in conversion_specific:
        cset = set(conv.get("support_maps") or [])
        if conv.get("center_map") is not None:
            cset.add(int(conv["center_map"]))
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") in broad_conversion_surrogates
                and (
                    (r.get("center_map") is not None and int(r["center_map"]) in cset)
                    or bool(set(r.get("support_maps") or []) & cset)
                )
            )
        ]

    # Any validated exact directional conversion outranks generic same-center
    # substitution bookkeeping. This is crucial for aldehyde -> nitrile and other
    # conversions that necessarily incorporate mapped O/N atoms.
    exact_conversion_centers = {
        int(r["center_map"])
        for r in conversion_specific
        if r.get("center_map") is not None
    }
    if exact_conversion_centers:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event")
                in {
                    "substitution",
                    "c_n_substitution",
                    "c_o_substitution",
                    "c_s_substitution",
                    "aryl_or_heteroaryl_amination",
                    "n_alkylation",
                    "o_arylation",
                    "o_alkylation",
                    "etherification_or_o_alkylation",
                    "halogenation_or_halide_exchange",
                }
                and r.get("center_map") in exact_conversion_centers
            )
        ]

    # Reductive amination explicitly outranks generic C-N substitution at the same C.
    redam_centers = {
        r.get("center_map")
        for r in chosen
        if r.get("event") == "reductive_amination" and r.get("center_map") is not None
    }
    if redam_centers:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") in {"c_n_substitution", "substitution"}
                and r.get("center_map") in redam_centers
            )
        ]

    # Azide substitution outranks generic C-N substitution at that carbon.
    az_centers = {
        r.get("center_map")
        for r in chosen
        if r.get("event") == "azide_substitution" and r.get("center_map") is not None
    }
    if az_centers:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") in {"c_n_substitution", "substitution"}
                and r.get("center_map") in az_centers
            )
        ]

    epoxide_azide_maps = set()
    for r in records:
        if (
            r.get("role") == "directional_conversion"
            and r.get("event") == "epoxide_to_azido_alcohol"
        ):
            epoxide_azide_maps.update(r.get("support_maps") or [])
            if r.get("center_map") is not None:
                epoxide_azide_maps.add(int(r["center_map"]))
    if epoxide_azide_maps:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") in {"o_dealkylation", "c_n_substitution", "substitution"}
                and (
                    r.get("center_map") in epoxide_azide_maps
                    or bool(set(r.get("support_maps") or []) & epoxide_azide_maps)
                )
            )
        ]

    # Ring-rewiring/QC vetoes isolated substitution/aromatic-functionalization calls on
    # the same mapped ring unless their center lies outside the rewiring support.
    ring_rewire_maps = set()
    for r in records:
        if r.get("event") == "ring_rewiring_candidate":
            ring_rewire_maps.update(r.get("support_maps") or [])
    if ring_rewire_maps:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event")
                in {
                    "substitution",
                    "c_n_substitution",
                    "c_o_substitution",
                    "c_s_substitution",
                    "aryl_or_heteroaryl_amination",
                    "n_alkylation",
                    "o_arylation",
                    "o_alkylation",
                    "etherification_or_o_alkylation",
                    "halogenation_or_halide_exchange",
                    "aromatic_functionalization",
                    "o_dealkylation",
                }
                and r.get("center_map") in ring_rewire_maps
            )
        ]

    # Promote exact conversions when they are not already explained by a stronger event
    # on the same mapped support. Some exact conversions are intentionally co-events.
    explained = set()
    for r in chosen:
        explained.update(r.get("support_maps") or [])
        if r.get("center_map") is not None:
            explained.add(int(r["center_map"]))

    allow_co = {
        "ester_to_carboxylic_acid",
        "nitrile_to_amide",
        "nitrile_to_carboxylic_acid",
        "alkene_to_epoxide",
        "alkene_to_diol",
        "alkene_to_alkane",
        "alkyne_to_alkene",
        "alkyne_to_alkane",
        "alkane_to_alkene",
        "methylene_to_ketone",
        "alkene_oxidative_cleavage_to_aldehyde",
    }
    for r in conversion_specific:
        support = set(r.get("support_maps") or [])
        center = r.get("center_map")
        overlap = bool(support & explained) or (
            center is not None and int(center) in explained
        )
        directional_name = r.get("directional_conversion") or r.get("event")
        if (not overlap) or directional_name in allow_co:
            chosen.append(r)
            explained.update(support)

    if not chosen and conversion_specific:
        chosen.extend(conversion_specific)

    # Acid -> ester owns the acyl transformation. Suppress donor-side O-alkylation
    # bookkeeping when that donor carbon is exactly the newly attached ester O substituent.
    ester_donor_maps = {
        int(r["ester_donor_map"])
        for r in records
        if r.get("role") == "directional_conversion"
        and r.get("event") == "carboxylic_acid_to_ester"
        and r.get("ester_donor_map") is not None
    }
    if ester_donor_maps:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event") == "etherification_or_o_alkylation"
                and r.get("center_map") in ester_donor_maps
            )
        ]

    # ------------------------------------------------------------------
    # S2_4_PRIMARY_OWNERSHIP
    #
    # These two S2.4 semantic records are valid local structural evidence,
    # but they are secondary interpretations when their support overlaps an
    # already-selected primary event.
    #
    # Broad audit evidence:
    #
    #   alkene_to_ketone_oxidation:
    #       31 overlapping prior primaries
    #        6 disjoint prior primaries
    #
    #   arene_partial_hydrogenation:
    #        3 overlapping prior primaries
    #
    # Therefore:
    #   * overlapping support -> retain as specific evidence only
    #   * disjoint support    -> retain independent primary multiplicity
    #
    # Do not suppress the detector record itself.
    # ------------------------------------------------------------------

    s2_4_secondary_primary_events = {
        "alkene_to_ketone_oxidation",
        "arene_partial_hydrogenation",
    }

    def _s2_4_primary_support(
        rec: dict[str, Any],
    ) -> set[int]:

        support = {int(x) for x in (rec.get("support_maps") or []) if x is not None}

        center = rec.get("center_map")

        if center is not None:
            support.add(int(center))

        return support

    # Snapshot support from every already-selected non-secondary
    # interpretation.  Secondary S2.4 records must not suppress one
    # another; this also preserves genuine same-event multicenter
    # multiplicity.
    s2_4_prior_owner_supports = [
        _s2_4_primary_support(r)
        for r in chosen
        if (
            r.get("event") not in s2_4_secondary_primary_events
            and _s2_4_primary_support(r)
        )
    ]

    if s2_4_prior_owner_supports:

        filtered = []

        for r in chosen:

            if r.get("event") not in s2_4_secondary_primary_events:

                filtered.append(r)

                continue

            support = _s2_4_primary_support(r)

            overlaps_prior_owner = any(
                bool(support & owner_support)
                for owner_support in s2_4_prior_owner_supports
            )

            if not overlaps_prior_owner:

                # Independent center: retain.
                filtered.append(r)

        chosen = filtered

    # Strong composite events suppress generic substitution bookkeeping on their maps.
    strong_events = {
        "condensation",
        "elimination",
        "pericyclic_compatible",
        "aromatic_functionalization",
        "cross_coupling_compatible",
        "suzuki_type_cross_coupling",
        "carbonyl_alpha_arylation",
        "aryl_alkene_hydroarylation",
        "vinyl_alkynyl_coupling",
        "ring_opening",
        "ring_formation",
        "ring_expansion",
        "ring_contraction",
        "ring_rewiring",
        "epoxide_rearrangement",
        "epoxide_formation",
        "epoxide_opening",
        "aziridine_formation",
        "aziridine_opening",
        "cyclopropane_formation",
        "cyclopropane_opening",
        "oxetane_formation",
        "oxetane_opening",
        "beta_lactam_formation",
        "beta_lactam_opening",
        "heterocycle_annulation",
        "carbocycle_annulation",
        "intramolecular_heterocyclization",
        "heterocycle_formation",
        "macrocyclization",
        "spirocyclization",
        "reductive_amination",
        "carbonyl_addition",
        "carbonyl_olefination",
        "intramolecular_anhydride_formation",
        "anhydride_formation",
        "lactamization",
        "lactonization",
        "amidation",
        "carbamate_formation",
        "acyl_halide_formation",
        "transamidation",
        "transesterification",
        "weinreb_amide_to_ketone",
        "aromatic_acylation",
        "heterocumulene_addition",
        "isocyanate_formation",
        "isothiocyanate_formation",
        "thiourea_formation",
        "urea_formation",
        "alpha_hydroxyimino_amide_formation",
        "alpha_oximation",
        "alkene_sulfenyl_halogenation",
        "n_glycosylation",
        "vinylic_amination",
        "n_n_bond_formation",
        "n_deprotection",
        "carbonyl_deprotection",
        "heteroaryl_hydrolysis_to_oxo",
        "aza_michael_addition",
        "conjugate_acyl_addition",
        "imine_reduction",
        "carbon_to_nitrogen_migration",
        "hofmann_rearrangement",
        "nitrile_to_thioamide",
        "nitrile_to_ester",
        "guanylation",
        "formamidine_formation",
        "chan_lam_n_arylation",
        "n_phosphorus_functionalization",
        "o_phosphorylation",
        "alkylphosphonate_formation",
        "phosphorus_ylide_formation",
        "organozinc_c_s_coupling",
        "hydrodealkoxylation",
        "halocyclization",
        "azide_substitution",
        "heteroarene_partial_hydrogenation",
        "arene_partial_hydrogenation",
        "double_bond_migration",
        "pi_heteroatom_hydroaddition",
        "alkene_to_ketone_oxidation",
        "alkyne_hydration_to_ketone",
        "alkene_hydration",
        "alkene_aminohydroxylation",
    }
    strong_maps = set()
    for r in chosen:
        if r.get("event") in strong_events:
            strong_maps.update(r.get("support_maps") or [])
            if r.get("center_map") is not None:
                strong_maps.add(int(r["center_map"]))
    if strong_maps:
        chosen = [
            r
            for r in chosen
            if not (
                r.get("event")
                in {
                    "substitution",
                    "c_n_substitution",
                    "c_o_substitution",
                    "c_s_substitution",
                    "aryl_or_heteroaryl_amination",
                    "n_alkylation",
                    "o_arylation",
                    "o_alkylation",
                    "etherification_or_o_alkylation",
                    "halogenation_or_halide_exchange",
                }
                and r.get("center_map") in strong_maps
            )
        ]

    # Generic substitution only if no specific interpretation exists at the center.
    specific_centers = {
        r.get("center_map") for r in chosen if r.get("center_map") is not None
    }
    explained_maps = set()
    for r in chosen:
        explained_maps.update(r.get("support_maps") or [])
    for r in records:
        if r.get("event") in {"mapping_rewire_candidate", "ring_rewiring_candidate"}:
            explained_maps.update(r.get("support_maps") or [])
    if not events.get("cross_coupling_compatible"):
        for r in events.get("substitution", []):
            center = r.get("center_map")
            if center not in specific_centers and center not in explained_maps:
                chosen.append(r)

    # Positive-evidence rearrangement/fragmentation remain available.
    chosen.extend(events.get("rearrangement", []))
    chosen.extend(events.get("fragmentation", []))
    chosen.extend(events.get("multi_bond_fragmentation", []))

    # multicomponent_incorporation is deliberately NOT a primary event in batch04.

    # If nothing survives but there is C-C intercomponent formation, retain a neutral
    # unresolved structural event rather than force a coupling family.
    if not chosen:
        cc = [r for r in records if r.get("event") == "c_c_intercomponent_formation"]
        if cc:
            chosen.extend({**r, "event": "c_c_bond_formation_unresolved"} for r in cc)

    # S2_1C_MUKAIYAMA_PRIMARY_OWNERSHIP
    #
    # A silyl-enol-ether -> carbonyl change is subordinate when its
    # original enol partner is exactly the incoming carbon participating
    # in a mapped carbonyl-addition event.
    #
    # Example: Mukaiyama aldol.
    #
    # Keep the silicon record in specific_events, but do not expose it as
    # a second defining primary event. Independent silicon conversions at
    # unrelated centers are unaffected.

    carbonyl_add_incoming_maps = {
        int(r["incoming_map"])
        for r in events.get(
            "carbonyl_addition",
            [],
        )
        if r.get("incoming_map") is not None
    }

    if carbonyl_add_incoming_maps:

        filtered = []

        for r in chosen:

            if r.get("event") == "silyl_enol_ether_to_carbonyl":
                enol_partner_maps = {
                    int(x)
                    for x in (
                        r.get(
                            "enol_partner_maps",
                            [],
                        )
                        or []
                    )
                }

                if enol_partner_maps & carbonyl_add_incoming_maps:
                    continue

            filtered.append(r)

        chosen = filtered

    # POST_V67_RRG: validated positive-evidence children; retain orthogonal events.
    for event_name in RRG_V67_EVENTS:
        chosen.extend(events.get(event_name, []))
    chosen = resolve_v67_rearrangement_primary(chosen)
    chosen = _resolve_sub01_dual_halide_cc_ownership(facts, chosen)
    chosen = _resolve_v67_ring_primary_ownership(chosen)
    chosen = _resolve_p6_azole_exchange_ownership(chosen)

    # P1_STAGE2_ACYL_OWNERSHIP:
    # The oxygen-side C(acyl)-O bond loss is valid specific/raw evidence,
    # but not an independent O-deacylation *primary* if the identical
    # carbonyl carbon is independently validated as undergoing an acyl
    # substitution (e.g. acid esterification or lactone aminolysis).
    # Require the exact broken O-C(acyl) mapped bond, not mere support
    # overlap or another acyl event elsewhere in the reaction.
    acyl_exchange_centers = {
        int(rec["center_map"])
        for rec in chosen
        if (
            rec.get("event") in {
                "esterification", "amidation", "transamidation", "transesterification",
            }
            or (
                rec.get("event") == "acyl_heteroatom_relocation"
                and rec.get("subtype") == "acyl_o_to_n"
            )
        )
        and rec.get("center_map") is not None
    }
    if acyl_exchange_centers:
        filtered = []
        for rec in chosen:
            suppress = False
            if (
                rec.get("event") == "o_deacylation_or_ester_cleavage"
                and rec.get("center_map") is not None
            ):
                oxygen = int(rec["center_map"])
                # A real departing ester/lactone O must have lost the
                # specific mapped acyl-carbon attachment.
                for att in _attachments(facts, "reactant", oxygen):
                    if att.get("element") != "C" or att.get("map") is None:
                        continue
                    acyl = int(att["map"])
                    if acyl not in acyl_exchange_centers:
                        continue
                    if not bool(att.get("carbonyl_like")):
                        continue
                    retained = any(
                        x.get("element") == "C"
                        and x.get("map") is not None
                        and int(x["map"]) == acyl
                        for x in _attachments(facts, "product", oxygen)
                    )
                    if not retained:
                        suppress = True
                        break
            if not suppress:
                filtered.append(rec)
        chosen = filtered

    chosen = _dedupe_records(chosen)
    chosen.sort(key=lambda r: SPECIFIC_PRIORITY.get(r.get("event"), 50), reverse=True)
    return chosen


def _component_provenance_records(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    if len(facts.get("intercomponent_formed_bonds", [])) >= 2:
        out.append(
            {
                "event": "multicomponent_incorporation",
                "role": "component_provenance",
                "support_maps": sorted(
                    {
                        m
                        for b in facts["intercomponent_formed_bonds"]
                        for m in (int(b["a"]), int(b["b"]))
                    }
                ),
                "n_intercomponent_formed_bonds": len(
                    facts["intercomponent_formed_bonds"]
                ),
                "mapping_scope": "mapped_intercomponent_endpoints",
            }
        )
    return out


def classify_system_i(reaction: str) -> dict[str, Any]:
    facts = build_mapped_graph_facts(reaction)
    detector_outputs = run_detectors(facts)
    records = _dedupe_records(_flatten(detector_outputs))

    # Structural oxidation/reduction events may carry a redox_direction even when
    # they are not simple directional FG conversions. Derive the axis from that
    # validated event, never from raw local oxidation-state arithmetic.
    derived_axes = []
    for r in records:
        direction = r.get("redox_direction")
        if r.get("role") == "specific_event" and direction in {
            "oxidation",
            "reduction",
        }:
            derived_axes.append(
                {
                    "event": direction,
                    "role": "axis_event",
                    "subtype": r.get("event"),
                    "center_map": r.get("center_map"),
                    "support_maps": r.get("support_maps") or [],
                    "mapping_scope": "derived_from_validated_specific_event",
                }
            )
    records = _dedupe_records(records + derived_axes)

    by_role = defaultdict(list)
    for r in records:
        by_role[r.get("role", "unknown")].append(r)

    component_provenance = _dedupe_records(_component_provenance_records(facts))
    primary = [] if not facts["mapping"]["valid"] else _primary_events(facts, records)

    semantic_events = (
        build_semantic_events(facts, primary) if facts["mapping"]["valid"] else []
    )
    event_topology = whole_reaction_event_topology(facts)

    mapping_notes = []
    if facts["mapping"]["element_mismatches"]:
        mapping_notes.append("element_map_mismatch")
    if (
        facts["mapping"]["duplicate_maps_reactant"]
        or facts["mapping"]["duplicate_maps_product"]
    ):
        mapping_notes.append("duplicate_atom_maps")
    if (
        facts["mapping"]["parse_failures_reactant"]
        or facts["mapping"]["parse_failures_product"]
    ):
        mapping_notes.append("rdkit_parse_failure")

    return {
        "system": SYSTEM_I_VERSION,
        "versions": {
            "graph": GRAPH_VERSION,
            "rules": RULEPACK_VERSION,
            "event_topology": EVENT_TOPOLOGY_VERSION,
            "semantic_events": SEMANTIC_EVENT_VERSION,
        },
        "reaction": reaction,
        "mapping_qc": {
            "valid": facts["mapping"]["valid"],
            "notes": mapping_notes,
            "n_common_maps": len(facts["mapping"]["common_maps"]),
            "common_maps": facts["mapping"]["common_maps"],
            "reactant_only_maps": facts["mapping"]["reactant_only_maps"],
            "product_only_maps": facts["mapping"]["product_only_maps"],
        },
        "primary_events": primary,
        "semantic_events": semantic_events,
        "semantic_event_count": len(semantic_events),
        "multi_semantic_event": len(semantic_events) > 1,
        "event_topology": event_topology,
        "specific_events": by_role["specific_event"],
        "local_events": by_role["local_event"],
        "primitives": by_role["primitive"],
        "topology_evidence": by_role["topology_evidence"],
        "directional_conversions": by_role["directional_conversion"],
        "redox_axis": by_role["axis_event"],
        "component_provenance": component_provenance,
        "formed_bond_classes": facts["formed_bond_classes"],
        "graph_summary": {
            "formed_bonds": facts["formed_bonds"],
            "broken_bonds": facts["broken_bonds"],
            "bond_order_changes": facts["bond_order_changes"],
            "intercomponent_formed_bonds": facts["intercomponent_formed_bonds"],
            "component_splits": facts["component_splits"],
            "component_merges": facts["component_merges"],
            "common_cycle_rank_reactant": facts["common_cycle_rank_reactant"],
            "common_cycle_rank_product": facts["common_cycle_rank_product"],
            "common_cycle_rank_delta": facts["common_cycle_rank_delta"],
            "validated_ring_openings": facts["ring_openings"],
            "validated_ring_formations": facts["ring_formations"],
            "ring_lineage_v2": facts.get("ring_lineage_v2", {}),
        },
        "detector_outputs": detector_outputs,
        "generic_label_policy": GENERIC_DEMOTIONS,
    }
