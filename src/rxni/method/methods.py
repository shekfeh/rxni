from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any, Iterable

from .reagents import (
    extract_reagents_from_reaction,
    reagent_families,
    reagent_label,
    resolve_reagent_token,
)
from .reagent_roles import classify_reagent_role, role_profile, partition_condition_components


def _def(display_name: str, parent_event: str, *queries: str) -> dict[str, Any]:
    return {
        "display_name": display_name,
        "parent_event": parent_event,
        "knowledge_queries": list(queries) or [display_name],
    }


METHOD_DEFINITIONS: dict[str, dict[str, Any]] = {
    # Esterification
    "esterification.mitsunobu": _def("Mitsunobu esterification", "esterification", "Mitsunobu", "Mitsunobu esterification"),
    "esterification.steglich": _def("Steglich esterification", "esterification", "Steglich esterification", "DCC DMAP esterification"),
    "esterification.carbodiimide": _def("Carbodiimide-mediated esterification", "esterification", "carbodiimide esterification", "DCC esterification", "EDC esterification"),
    "esterification.thionyl_chloride": _def("SOCl2-assisted esterification", "esterification", "thionyl chloride esterification", "acid chloride esterification"),
    "esterification.acyl_chloride": _def("Acyl-chloride esterification", "esterification", "acyl chloride alcohol esterification"),
    "esterification.anhydride": _def("Anhydride acylation to ester", "esterification", "anhydride esterification"),
    "esterification.acid_catalyzed": _def("Acid-catalyzed direct esterification", "esterification", "Fischer esterification", "acid catalyzed esterification"),
    "esterification.transesterification": _def("Transesterification", "esterification", "transesterification"),

    # Amidation / acyl substitution
    "amidation.carbodiimide": _def("Carbodiimide-mediated amidation", "amidation", "carbodiimide amidation", "DCC peptide coupling", "EDC amidation"),
    "amidation.carbodiimide_additive": _def("Carbodiimide amidation with coupling additive", "amidation", "EDC HOBt amidation", "DCC HOBt peptide coupling", "Oxyma coupling"),
    "amidation.uronium": _def("Uronium-mediated amidation", "amidation", "HATU amidation", "HBTU amidation", "TBTU amidation"),
    "amidation.phosphonium": _def("Phosphonium-mediated amidation", "amidation", "PyBOP amidation", "phosphonium peptide coupling"),
    "amidation.acid_chloride": _def("Acid-chloride amidation", "amidation", "acid chloride amidation"),
    "amidation.anhydride": _def("Anhydride-mediated amidation", "amidation", "mixed anhydride amidation", "anhydride amidation"),

    # Hydrolysis / deprotection
    "hydrolysis.basic": _def("Basic hydrolysis / saponification", "hydrolysis_deprotection", "saponification", "basic ester hydrolysis"),
    "hydrolysis.acidic": _def("Acidic hydrolysis", "hydrolysis_deprotection", "acid hydrolysis"),
    "deprotection.hydrogenolysis": _def("Hydrogenolytic deprotection", "hydrolysis_deprotection", "hydrogenolysis deprotection", "Cbz hydrogenolysis"),
    "deprotection.acidolysis": _def("Acid-mediated deprotection", "hydrolysis_deprotection", "acid deprotection", "Boc deprotection TFA"),
    "deprotection.desilylation_fluoride": _def("Fluoride-mediated desilylation", "silyl_deprotection", "TBAF desilylation", "fluoride silyl ether deprotection"),
    "deprotection.o_dealkylation_boron_halide": _def("Boron-halide O-dealkylation", "o_dealkylation", "BBr3 demethylation", "boron tribromide ether cleavage"),

    # Alkyl substitution
    "substitution.base_assisted": _def("Base-assisted alkyl substitution", "substitution", "SN2 alkylation", "base assisted alkylation"),
    "substitution.phase_transfer": _def("Phase-transfer substitution", "substitution", "phase transfer alkylation"),
    "o_alkylation.mitsunobu": _def("Mitsunobu O-alkylation", "o_alkylation", "Mitsunobu etherification"),

    # Sulfonyl substitution
    "sulfonyl_substitution.base_assisted": _def("Base-assisted sulfonyl substitution", "sulfonyl_substitution", "sulfonyl chloride amine base", "sulfonamide formation"),
    "sulfonyl_substitution.nucleophilic_catalysis": _def("DMAP-assisted sulfonyl substitution", "sulfonyl_substitution", "DMAP sulfonylation"),

    # Aryl / heteroaryl substitution and coupling
    "aryl_coupling.suzuki": _def("Suzuki-type cross-coupling", "aryl_coupling", "Suzuki coupling", "Suzuki Miyaura"),
    "aryl_coupling.buchwald_hartwig": _def("Buchwald-Hartwig amination", "aryl_coupling", "Buchwald Hartwig amination"),
    "aryl_coupling.ullmann": _def("Copper-mediated aryl substitution", "aryl_coupling", "Ullmann coupling", "copper aryl amination"),
    "aryl_substitution.snar": _def("SNAr-compatible aromatic substitution", "aryl_coupling", "nucleophilic aromatic substitution", "SNAr"),
    "aryl_coupling.pd_catalyzed_generic": _def("Generic Pd-catalyzed aryl coupling", "aryl_coupling", "palladium catalyzed aryl coupling"),
    "aryl_coupling.nickel_catalyzed_generic": _def("Generic Ni-catalyzed aryl coupling", "aryl_coupling", "nickel catalyzed aryl coupling"),

    # Reduction / reductive chemistry
    "reduction.borohydride": _def("Borohydride reduction", "reduction", "sodium borohydride reduction", "borohydride carbonyl reduction"),
    "reduction.aluminum_hydride": _def("Aluminum-hydride reduction", "reduction", "LiAlH4 reduction"),
    "reduction.dibal": _def("DIBAL-type reduction", "reduction", "DIBAL reduction"),
    "hydrogenation.catalytic": _def("Catalytic hydrogenation", "reduction", "hydrogenation palladium", "H2 Pd hydrogenation"),
    "reductive_amination.borohydride": _def("Borohydride reductive amination", "reductive_amination", "sodium triacetoxyborohydride reductive amination", "cyanoborohydride reductive amination"),

    # Oxidation / special observed families
    "oxidation.peracid": _def("Peracid oxidation", "oxidation", "mCPBA oxidation", "peracid sulfoxidation"),
    "halogenation.n_halosuccinimide": _def("N-halosuccinimide halogenation", "halogenation", "NBS bromination", "NCS chlorination", "NIS iodination"),
    "cyanation.cyanide_source": _def("Cyanide-source nitrile installation", "cyanation", "cyanation", "nitrile installation cyanide"),
}


def method_definition(method_id: str) -> dict[str, Any]:
    d = dict(METHOD_DEFINITIONS.get(method_id, {}))
    d.setdefault("display_name", method_id)
    d.setdefault("parent_event", method_id.split(".", 1)[0])
    d.setdefault("knowledge_queries", [d["display_name"]])
    return d


def _mol_has(smiles: str, smarts: str) -> bool:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        pat = Chem.MolFromSmarts(smarts)
        return bool(mol is not None and pat is not None and mol.HasSubstructMatch(pat))
    except Exception:
        return False


def _left_components(reaction: str) -> list[str]:
    left = str(reaction).split(">", 1)[0]
    return [x for x in left.split(".") if x]


def _substrate_flags(reaction: str) -> dict[str, bool]:
    comps = _left_components(reaction)
    return {
        "carboxylic_acid": any(_mol_has(x, "[CX3](=O)[OX2H1]") for x in comps),
        "alcohol": any(_mol_has(x, "[OX2H1][#6]") for x in comps),
        "amine": any(_mol_has(x, "[N;H1,H2;!$(NC=O)]") for x in comps),
        "acyl_chloride": any(_mol_has(x, "[CX3](=O)Cl") for x in comps),
        "anhydride": any(_mol_has(x, "[CX3](=O)O[CX3](=O)") for x in comps),
        "ester": any(_mol_has(x, "[CX3](=O)O[#6]") for x in comps),
        "aryl_halide": any(_mol_has(x, "[c][F,Cl,Br,I]") for x in comps),
        "alkyl_halide": any(_mol_has(x, "[CX4][Cl,Br,I]") for x in comps),
        "sulfonyl_chloride": any(_mol_has(x, "S(=O)(=O)Cl") for x in comps),
        "organoboron": any(_mol_has(x, "[#6][B]") or _mol_has(x, "[B][#6]") for x in comps),
        "carbonyl": any(_mol_has(x, "[CX3]=O") for x in comps),
    }


def _pred(
    method_id: str,
    confidence: float,
    defining: list[str],
    supporting: list[str] | None = None,
    discriminator: str | None = None,
) -> dict[str, Any]:
    d = method_definition(method_id)
    # `label` remains as a compatibility alias for pre-v6.5.2 callers.
    legacy = {
        "esterification.mitsunobu": "mitsunobu_like",
        "esterification.steglich": "steglich_like",
    }.get(method_id, method_id)
    return {
        "label": legacy,
        "method_id": method_id,
        "display_name": d["display_name"],
        "parent_event": d["parent_event"],
        "confidence": float(confidence),
        "defining_evidence": defining,
        "supporting_evidence": supporting or [],
        "discriminator": discriminator,
    }


def _families(reagents: Iterable[str]) -> set[str]:
    out: set[str] = set()
    for r in reagents:
        out.update(reagent_families(r))
    return out


def _has_quaternary_ammonium(reagents: Iterable[str]) -> bool:
    return any(_mol_has(r, "[N+;X4]") for r in reagents)


def _event_in(event: str, values: set[str]) -> bool:
    return event in values


HYDROLYSIS_EVENTS = {
    "ester_to_carboxylic_acid", "o_dealkylation", "n_deprotection",
    "carbonyl_deprotection", "heteroaryl_hydrolysis_to_oxo",
    "o_deacylation_or_ester_cleavage", "silyl_deprotection",
}
SUBSTITUTION_EVENTS = {
    "n_alkylation", "substitution", "c_n_substitution",
    "etherification_or_o_alkylation", "o_alkylation",
}
ARYL_EVENTS = {
    "aryl_or_heteroaryl_amination", "o_arylation", "suzuki_type_cross_coupling",
    "cross_coupling_compatible", "carbonyl_alpha_arylation", "aromatic_functionalization",
}
REDUCTION_EVENTS = {
    "aldehyde_to_primary_alcohol", "ester_to_primary_alcohol",
    "ketone_to_secondary_alcohol", "alkene_hydrogenation",
    "aryl_halide_to_hydrogen", "reductive_amination",
}


def infer_methods(
    primary_event: str | None,
    reaction: str,
    reagents: Iterable[str] | None = None,
    *,
    solvents: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    """Infer method hypotheses only from event-compatible condition signatures.

    Reagents cannot define the structural transformation.  The primary event is
    always an upstream prerequisite/gate.
    """
    event = (primary_event or "").lower()
    rs = set(reagents if reagents is not None else extract_reagents_from_reaction(reaction))
    fam = _families(rs)
    flags = _substrate_flags(reaction)
    labels = {reagent_label(x) for x in rs}
    out: list[dict[str, Any]] = []

    pph3 = resolve_reagent_token("PPh3")
    dmap = resolve_reagent_token("DMAP")
    socl2 = resolve_reagent_token("SOCl2")

    # Esterification
    if event in {"esterification", "carboxylic_acid_to_ester"}:
        if "phosphine" in fam and "azodicarboxylate" in fam:
            out.append(_pred("esterification.mitsunobu", .98, ["phosphine", "azodicarboxylate"], discriminator="phosphine + azodicarboxylate within an esterification transformation"))
        if "carbodiimide" in fam and dmap in rs:
            out.append(_pred("esterification.steglich", .97, ["carbodiimide", "DMAP"], discriminator="carbodiimide + DMAP"))
        elif "carbodiimide" in fam:
            out.append(_pred("esterification.carbodiimide", .92, ["carbodiimide"], discriminator="carbodiimide without full Steglich signature"))
        if socl2 in rs:
            out.append(_pred("esterification.thionyl_chloride", .82, ["SOCl2"]))
        if flags["acyl_chloride"]:
            out.append(_pred("esterification.acyl_chloride", .96, ["acyl chloride reactant"]))
        if flags["anhydride"]:
            out.append(_pred("esterification.anhydride", .94, ["carboxylic anhydride reactant"]))
        if flags["ester"] and flags["alcohol"] and not flags["carboxylic_acid"] and not out:
            out.append(_pred("esterification.transesterification", .82, ["ester reactant", "alcohol reactant"]))
        if "strong_acid" in fam and flags["carboxylic_acid"] and flags["alcohol"] and not out:
            out.append(_pred("esterification.acid_catalyzed", .88, ["carboxylic acid", "alcohol", "strong acid"]))

    # Amidation / transamidation / carbamate context
    if event in {"amidation", "transamidation", "carbamate_formation"}:
        if "carbodiimide" in fam and "coupling_additive" in fam:
            out.append(_pred("amidation.carbodiimide_additive", .97, ["carbodiimide", "coupling additive"]))
        elif "carbodiimide" in fam:
            out.append(_pred("amidation.carbodiimide", .93, ["carbodiimide"]))
        if "uronium_coupling_reagent" in fam:
            out.append(_pred("amidation.uronium", .96, ["uronium coupling reagent"]))
        if "phosphonium_coupling_reagent" in fam:
            out.append(_pred("amidation.phosphonium", .96, ["phosphonium coupling reagent"]))
        if flags["acyl_chloride"]:
            out.append(_pred("amidation.acid_chloride", .96, ["acyl chloride reactant"]))
        if flags["anhydride"]:
            out.append(_pred("amidation.anhydride", .90, ["anhydride reactant"]))

    # Hydrolysis / deprotection
    if event in HYDROLYSIS_EVENTS:
        if event == "silyl_deprotection" and "fluoride_source" in fam:
            out.append(_pred("deprotection.desilylation_fluoride", .98, ["fluoride source", "silyl deprotection event"]))
        if event == "o_dealkylation" and "lewis_acid_dealkylation" in fam:
            out.append(_pred("deprotection.o_dealkylation_boron_halide", .97, ["BBr3/BCl3-like Lewis acid", "O-dealkylation event"]))
        if "hydroxide_base" in fam:
            out.append(_pred("hydrolysis.basic", .92, ["hydroxide base"]))
        if "strong_acid" in fam:
            if event in {"n_deprotection", "carbonyl_deprotection", "o_dealkylation"}:
                out.append(_pred("deprotection.acidolysis", .83, ["strong acid", event]))
            else:
                out.append(_pred("hydrolysis.acidic", .85, ["strong acid"]))
        if "hydrogen" in fam and fam & {"palladium_catalyst", "platinum_catalyst", "nickel_catalyst"}:
            out.append(_pred("deprotection.hydrogenolysis", .93, ["H2", "metal hydrogenation catalyst"]))

    # N/O alkylation and substitution
    if event in SUBSTITUTION_EVENTS:
        if "phosphine" in fam and "azodicarboxylate" in fam and event in {"o_alkylation", "etherification_or_o_alkylation"}:
            out.append(_pred("o_alkylation.mitsunobu", .94, ["phosphine", "azodicarboxylate", event]))
        if fam & {"hydroxide_base", "carbonate_base", "tertiary_amine_base"}:
            out.append(_pred("substitution.base_assisted", .74, ["base", event], discriminator="generic base-assisted substitution; not a named-method assignment"))
        if _has_quaternary_ammonium(rs) and fam & {"hydroxide_base", "carbonate_base"}:
            out.append(_pred("substitution.phase_transfer", .80, ["quaternary ammonium", "inorganic base"]))

    # Sulfonyl substitution
    if event == "sulfonyl_substitution":
        if fam & {"tertiary_amine_base", "carbonate_base", "hydroxide_base"}:
            out.append(_pred("sulfonyl_substitution.base_assisted", .87, ["base", "sulfonyl substitution event"]))
        if dmap in rs:
            out.append(_pred("sulfonyl_substitution.nucleophilic_catalysis", .88, ["DMAP", "sulfonyl substitution event"]))

    # Aryl / heteroaryl substitution and coupling
    if event in ARYL_EVENTS:
        metal = fam & {"palladium_catalyst", "copper_catalyst", "nickel_catalyst"}
        if event == "suzuki_type_cross_coupling" or ("palladium_catalyst" in fam and ("organoboron_or_boron_reagent" in fam or flags["organoboron"])):
            out.append(_pred("aryl_coupling.suzuki", .96 if event == "suzuki_type_cross_coupling" else .88, ["organoboron", "Pd-compatible coupling context"]))
        if event == "aryl_or_heteroaryl_amination" and "palladium_catalyst" in fam and "phosphine" in fam:
            out.append(_pred("aryl_coupling.buchwald_hartwig", .92, ["Pd", "phosphine", "C-N arylation event"]))
        if event in {"aryl_or_heteroaryl_amination", "o_arylation"} and "copper_catalyst" in fam:
            out.append(_pred("aryl_coupling.ullmann", .86, ["Cu", event]))
        if event in {"aryl_or_heteroaryl_amination", "o_arylation", "aromatic_functionalization"} and flags["aryl_halide"] and not metal:
            if fam & {"carbonate_base", "hydroxide_base", "tertiary_amine_base"}:
                out.append(_pred("aryl_substitution.snar", .68, ["aryl halide", "nucleophile/base", "no detected transition-metal catalyst"], discriminator="SNAr-compatible, not definitive"))
        if event in {"cross_coupling_compatible", "carbonyl_alpha_arylation", "aromatic_functionalization"}:
            if "palladium_catalyst" in fam and not any(x["method_id"] in {"aryl_coupling.suzuki", "aryl_coupling.buchwald_hartwig"} for x in out):
                out.append(_pred("aryl_coupling.pd_catalyzed_generic", .72, ["Pd catalyst", event], discriminator="generic Pd-catalyzed coupling evidence; donor-specific method unresolved"))
            if "nickel_catalyst" in fam:
                out.append(_pred("aryl_coupling.nickel_catalyzed_generic", .72, ["Ni catalyst", event], discriminator="generic Ni-catalyzed coupling evidence"))

    # Reduction / reductive amination
    if event in REDUCTION_EVENTS:
        if event == "reductive_amination" and "borohydride" in fam:
            out.append(_pred("reductive_amination.borohydride", .91, ["borohydride family", "reductive amination event"]))
        elif "borohydride" in fam:
            out.append(_pred("reduction.borohydride", .90, ["borohydride family"]))
        if "aluminum_hydride" in fam:
            if "DIBAL-H" in labels:
                out.append(_pred("reduction.dibal", .91, ["DIBAL-H-like reagent"]))
            else:
                out.append(_pred("reduction.aluminum_hydride", .91, ["aluminum hydride family"]))
        if "hydrogen" in fam and fam & {"palladium_catalyst", "platinum_catalyst", "nickel_catalyst"}:
            out.append(_pred("hydrogenation.catalytic", .95, ["H2", "Pd/Pt/Ni catalyst"]))

    if event in {"sulfide_to_sulfoxide", "aldehyde_to_carboxylic_acid"} and "peracid_oxidant" in fam:
        out.append(_pred("oxidation.peracid", .92, ["peracid oxidant", event]))

    if event == "halogenation_or_halide_exchange" and "n_halosuccinimide" in fam:
        out.append(_pred("halogenation.n_halosuccinimide", .90, ["N-halosuccinimide", event]))

    if event == "cyanation_or_nitrile_installation" and "cyanide_source" in fam:
        out.append(_pred("cyanation.cyanide_source", .94, ["cyanide source", event]))

    # Deterministic deduplication by method id, keeping strongest evidence.
    best: dict[str, dict[str, Any]] = {}
    for p in out:
        mid = p["method_id"]
        if mid not in best or p["confidence"] > best[mid]["confidence"]:
            best[mid] = p
    return sorted(best.values(), key=lambda x: (-x["confidence"], x["method_id"]))


def build_method_landscape(
    member_rows: Iterable[dict[str, Any]],
    *,
    primary_event: str | None,
    top_reagents: int = 20,
    top_signatures: int = 20,
    landscape_type: str = "global",
    sort_by_weighted_support: bool = False,
) -> dict[str, Any]:
    method_reactions: Counter[str] = Counter()
    method_occurrences: Counter[str] = Counter()
    method_weight: Counter[str] = Counter()
    method_examples = defaultdict(list)
    method_meta: dict[str, dict[str, Any]] = {}
    agent_counts = Counter(); sig_counts = Counter(); reagents_per_rxn = []
    n = 0; n_with = 0; recognized = 0; occurrence_total = 0

    for row in member_rows:
        n += 1
        rxn = str(row.get("reaction") or row.get("forward_reaction_smiles") or "")
        occ = int(row.get("membership_occurrence_count") or 1)
        occurrence_total += occ
        agents = extract_reagents_from_reaction(rxn)
        condition_roles = partition_condition_components(agents)
        active_components = condition_roles.active_components
        reagents_per_rxn.append(agents)
        if agents:
            n_with += 1
        agent_counts.update(set(agents))
        if agents:
            sig_counts[tuple(sorted(agents))] += 1
        preds = infer_methods(primary_event, rxn, active_components, solvents=condition_roles.solvents)
        if preds:
            recognized += 1
        chem_score = float(row.get("chemistry_score") or row.get("score") or 1.0)
        for pred in preds:
            mid = pred["method_id"]
            method_reactions[mid] += 1
            method_occurrences[mid] += occ
            method_weight[mid] += chem_score * pred["confidence"]
            method_meta.setdefault(mid, pred)
            if len(method_examples[mid]) < 5:
                method_examples[mid].append({
                    "reaction_key": row.get("reaction_key"),
                    "reaction_id": row.get("reaction_id"),
                    "template_hash": row.get("template_hash"),
                    "chemistry_score": chem_score,
                    "reagents": [reagent_label(x) for x in condition_roles.reagents],
                    "solvents": [reagent_label(x) for x in condition_roles.solvents],
                    "catalysts": [reagent_label(x) for x in condition_roles.catalysts],
                    "bases": [reagent_label(x) for x in condition_roles.bases],
                    "additives": [reagent_label(x) for x in condition_roles.additives],
                })

    mids = list(method_reactions)
    if sort_by_weighted_support:
        mids.sort(key=lambda mid: (-method_weight[mid], -method_reactions[mid], mid))
    else:
        mids.sort(key=lambda mid: (-method_reactions[mid], -method_weight[mid], mid))

    methods = []
    for mid in mids:
        count = method_reactions[mid]
        meta = method_meta[mid]
        methods.append({
            **meta,
            "support_reactions": count,
            "support_fraction": count / n if n else 0.0,
            "support_occurrences": method_occurrences[mid],
            "weighted_support": method_weight[mid],
            "representative_members": method_examples[mid],
        })

    roles = role_profile(reagents_per_rxn, n_reactions=n)
    return {
        "schema_version": "method_landscape_v1",
        "schema_revision": "v653_conditioned_landscapes",
        "landscape_type": landscape_type,
        "parent_event": primary_event,
        "evidence_pool": {
            "candidate_reactions": n,
            "membership_occurrences": occurrence_total,
            "with_agent_field": n_with,
            "agent_coverage_fraction": n_with / n if n else 0.0,
            "recognized_reactions": recognized,
            "recognized_fraction": recognized / n if n else 0.0,
            "unresolved_reactions": n - recognized,
        },
        "methods": methods,
        "agent_roles": roles["roles"],
        "reagent_families": roles.get("families", []),
        "top_reagents": [
            {
                "reagent": k,
                "label": reagent_label(k),
                "role": classify_reagent_role(k).primary_role,
                "families": classify_reagent_role(k).families,
                "count": v,
                "fraction": v / n if n else 0.0,
            }
            for k, v in agent_counts.most_common(top_reagents)
        ],
        "top_agent_signatures": [
            {
                "reagents": list(k),
                "labels": [reagent_label(x) for x in k],
                "count": v,
                "fraction": v / n if n else 0.0,
            }
            for k, v in sig_counts.most_common(top_signatures)
        ],
    }


def build_conditioned_method_landscape(
    member_rows: Iterable[dict[str, Any]],
    *,
    primary_event: str | None,
    query_reaction: str,
    condition_evidence: Any,
) -> dict[str, Any]:
    """Build a query-conditioned landscape while preserving the global prior.

    Member rows should carry a condition-aware chemistry_score.  Direct query
    method signatures are reported separately and can add a strongly supported
    method even when it is rare in the retrieved empirical neighborhood.
    """
    land = build_method_landscape(
        member_rows,
        primary_event=primary_event,
        landscape_type="query_conditioned",
        sort_by_weighted_support=True,
    )
    qreagents = list(getattr(condition_evidence, "active_components", []) or [])
    if not qreagents:
        qreagents = list(getattr(condition_evidence, "reagents", []) or [])
    qsolvents = list(getattr(condition_evidence, "solvents", []) or [])
    direct = infer_methods(primary_event, query_reaction, qreagents, solvents=qsolvents)
    direct_by = {x["method_id"]: x for x in direct}

    existing = {m["method_id"] for m in land["methods"]}
    for m in land["methods"]:
        q = direct_by.get(m["method_id"])
        m["query_condition_support"] = 0.0 if q is None else q["confidence"]
        m["query_condition_evidence"] = [] if q is None else q["defining_evidence"]
        m["conditioned_rank_score"] = float(m["weighted_support"]) * (1.0 + m["query_condition_support"])

    for mid, q in direct_by.items():
        if mid in existing:
            continue
        land["methods"].append({
            **q,
            "support_reactions": 0,
            "support_fraction": 0.0,
            "support_occurrences": 0,
            "weighted_support": 0.0,
            "representative_members": [],
            "query_condition_support": q["confidence"],
            "query_condition_evidence": q["defining_evidence"],
            "conditioned_rank_score": q["confidence"],
            "empirical_support_status": "query_signature_only_no_retrieved_member_support",
        })

    land["methods"].sort(
        key=lambda m: (
            -float(m.get("conditioned_rank_score", 0.0)),
            -float(m.get("query_condition_support", 0.0)),
            -float(m.get("weighted_support", 0.0)),
            m["method_id"],
        )
    )
    land["query_condition_predictions"] = direct
    land["query_condition_evidence"] = (
        condition_evidence.as_dict() if hasattr(condition_evidence, "as_dict") else dict(condition_evidence or {})
    )
    land["interpretation"] = "conditioned_support_is_ranking_evidence_not_a_calibrated_probability"
    return land
