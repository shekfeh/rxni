from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from typing import Any, Iterable

from .reagents import reagent_families, reagent_label, resolve_reagent_token


def _canon_many(values):
    return {resolve_reagent_token(v) for v in values}


ROLE_REGISTRY = {
    "solvent": _canon_many([
        "THF", "C1CCOC1", "ClCCl", "CN(C)C=O", "CC#N", "Cc1ccccc1",
        "CCOC(C)=O", "CC(C)=O", "CO", "CCO", "O", "CCOCC",
        "O=C1OCCO1", "C1COCCO1", "CC(=O)N(C)C",
    ]),
    "phosphine": _canon_many(["PPh3"]),
    "azo_reagent": _canon_many(["DIAD", "DEAD"]),
    "carbodiimide": _canon_many(["DCC", "DIC", "EDC"]),
    "coupling_additive": _canon_many(["HOBt", "HOAt", "Oxyma"]),
    "acyl_activation_reagent": _canon_many(["SOCl2", "O=C(Cl)C(=O)Cl"]),
    "nucleophilic_catalyst": _canon_many(["DMAP"]),
    "base": _canon_many([
        "Et3N", "DIPEA", "pyridine", "NMM", "CC(C)(C)[O-]", "[OH-]",
    ]),
    "acid": _canon_many(["H2SO4", "TFA", "HCl"]),
    "spectator_or_counterion": _canon_many([
        "[Na+]", "[K+]", "[Li+]", "[Cs+]", "[Cl-]", "[Br-]", "[F-]", "[I-]",
    ]),
}

FAMILY_TO_ROLES = {
    "palladium_catalyst": ["metal_catalyst"],
    "platinum_catalyst": ["metal_catalyst"],
    "nickel_catalyst": ["metal_catalyst"],
    "copper_catalyst": ["metal_catalyst"],
    "phosphine": ["ligand_or_reagent", "phosphine"],
    "azodicarboxylate": ["azo_reagent"],
    "carbodiimide": ["coupling_reagent", "carbodiimide"],
    "coupling_additive": ["coupling_additive"],
    "uronium_coupling_reagent": ["coupling_reagent"],
    "phosphonium_coupling_reagent": ["coupling_reagent"],
    "tertiary_amine_base": ["base"],
    "hydroxide_base": ["base"],
    "carbonate_base": ["base"],
    "strong_acid": ["acid"],
    "lewis_acid_dealkylation": ["lewis_acid"],
    "fluoride_source": ["fluoride_source"],
    "borohydride": ["reducing_agent"],
    "aluminum_hydride": ["reducing_agent"],
    "hydrogen": ["reducing_agent"],
    "peracid_oxidant": ["oxidant"],
    "n_halosuccinimide": ["halogenating_agent"],
    "cyanide_source": ["cyanide_source"],
    "quaternary_ammonium": ["phase_transfer_catalyst"],
    "alkali_counterion": ["spectator_or_counterion"],
    "halide_counterion": ["spectator_or_counterion"],
}

ROLE_PRIORITY = (
    "metal_catalyst", "nucleophilic_catalyst", "phase_transfer_catalyst",
    "coupling_reagent", "coupling_additive", "reducing_agent", "oxidant",
    "halogenating_agent", "cyanide_source", "fluoride_source", "lewis_acid",
    "phosphine", "azo_reagent", "carbodiimide", "acyl_activation_reagent",
    "base", "acid", "ligand_or_reagent", "solvent", "spectator_or_counterion",
)


@dataclass(frozen=True)
class ReagentRole:
    smiles: str
    label: str
    primary_role: str
    all_roles: list[str]
    families: list[str]

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ConditionComponentRoles:
    """Role-resolved view of a reaction middle-field component pool.

    The partition is deliberately conservative.  It does not claim that a
    component can have only one chemical function; ``component_roles`` retains
    the full multi-role classification.  The buckets are a practical view for
    condition matching/reporting so that obvious solvents do not dilute reagent
    similarity and bases/catalysts/additives can be compared independently.
    """

    components: list[str]
    reagents: list[str]
    solvents: list[str]
    catalysts: list[str]
    bases: list[str]
    additives: list[str]
    acids: list[str]
    spectators: list[str]
    unresolved: list[str]
    component_roles: list[dict[str, Any]]

    @property
    def active_components(self) -> list[str]:
        out: list[str] = []
        for vals in (self.reagents, self.catalysts, self.bases, self.additives, self.acids):
            for value in vals:
                if value not in out:
                    out.append(value)
        return out

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["active_components"] = self.active_components
        d["labels"] = {
            name: [reagent_label(x) for x in getattr(self, name)]
            for name in (
                "reagents", "solvents", "catalysts", "bases", "additives",
                "acids", "spectators", "unresolved",
            )
        }
        return d


def classify_reagent_role(smiles: str) -> ReagentRole:
    c = resolve_reagent_token(smiles)
    families = list(reagent_families(smiles))
    roles = {r for r, vals in ROLE_REGISTRY.items() if c in vals}
    for fam in families:
        roles.update(FAMILY_TO_ROLES.get(fam, []))
    ordered = [r for r in ROLE_PRIORITY if r in roles]
    ordered += sorted(roles - set(ordered))
    primary = ordered[0] if ordered else "unresolved"
    return ReagentRole(c, reagent_label(c), primary, ordered or ["unresolved"], families)


def partition_condition_components(components: Iterable[str]) -> ConditionComponentRoles:
    canon: list[str] = []
    for raw in components:
        c = resolve_reagent_token(str(raw))
        if c and c not in canon:
            canon.append(c)

    buckets: dict[str, list[str]] = {
        "reagents": [], "solvents": [], "catalysts": [], "bases": [],
        "additives": [], "acids": [], "spectators": [], "unresolved": [],
    }
    details: list[dict[str, Any]] = []

    for c in canon:
        rr = classify_reagent_role(c)
        details.append(rr.as_dict())
        p = rr.primary_role
        if p == "solvent":
            bucket = "solvents"
        elif p in {"metal_catalyst", "nucleophilic_catalyst", "phase_transfer_catalyst"}:
            bucket = "catalysts"
        elif p == "base":
            bucket = "bases"
        elif p == "coupling_additive":
            bucket = "additives"
        elif p in {"acid", "lewis_acid"}:
            bucket = "acids"
        elif p == "spectator_or_counterion":
            bucket = "spectators"
        elif p == "unresolved":
            bucket = "unresolved"
        else:
            bucket = "reagents"
        buckets[bucket].append(c)

    return ConditionComponentRoles(
        components=canon,
        reagents=buckets["reagents"],
        solvents=buckets["solvents"],
        catalysts=buckets["catalysts"],
        bases=buckets["bases"],
        additives=buckets["additives"],
        acids=buckets["acids"],
        spectators=buckets["spectators"],
        unresolved=buckets["unresolved"],
        component_roles=details,
    )


def role_profile(reagent_lists: Iterable[Iterable[str]], *, n_reactions: int | None = None) -> dict[str, Any]:
    role_counts: Counter[str] = Counter()
    reagent_counts: Counter[tuple[str, str]] = Counter()
    family_counts: Counter[str] = Counter()
    observed = 0
    for reagents in reagent_lists:
        seen_roles = set(); seen_reagents = set(); seen_families = set()
        for r in reagents:
            rr = classify_reagent_role(r)
            seen_roles.add(rr.primary_role)
            seen_reagents.add((rr.smiles, rr.label))
            seen_families.update(rr.families)
        if seen_reagents:
            observed += 1
        role_counts.update(seen_roles)
        reagent_counts.update(seen_reagents)
        family_counts.update(seen_families)
    denom = n_reactions if n_reactions is not None else observed
    return {
        "n_reactions": denom,
        "n_with_reagents": observed,
        "roles": [
            {"role": k, "count": v, "fraction": v / denom if denom else 0.0}
            for k, v in role_counts.most_common()
        ],
        "families": [
            {"family": k, "count": v, "fraction": v / denom if denom else 0.0}
            for k, v in family_counts.most_common()
        ],
        "reagents": [
            {
                "reagent": s, "label": lab, "count": v,
                "fraction": v / denom if denom else 0.0,
                "role": classify_reagent_role(s).primary_role,
                "families": classify_reagent_role(s).families,
            }
            for (s, lab), v in reagent_counts.most_common()
        ],
    }
