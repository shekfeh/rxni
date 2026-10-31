from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from functools import lru_cache
from typing import Any, Iterable

# Query aliases are conveniences, not chemistry truth.  Complex reagents whose
# exact salt/solvate encoding varies are intentionally allowed to remain as
# stable textual tokens; family recognition below provides the chemistry layer.
REAGENT_ALIASES: dict[str, str] = {
    "pph3": "c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    "triphenylphosphine": "c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    "diad": "CC(C)OC(=O)N=NC(=O)OC(C)C",
    "dead": "CCOC(=O)N=NC(=O)OCC",
    "thf": "C1CCOC1",
    "dcc": "C1CCC(CC1)N=C=NC2CCCCC2",
    "dmap": "CN(C)c1ccncc1",
    "thionyl_chloride": "O=S(Cl)Cl",
    "socl2": "O=S(Cl)Cl",
    "edc": "CCN=C=NCCCN(C)C",
    "edci": "CCN=C=NCCCN(C)C",
    "edic": "CCN=C=NCCCN(C)C",
    "dic": "CC(C)N=C=NC(C)C",
    "hobt": "HOBt",
    "hoat": "HOAt",
    "oxyma": "Oxyma",
    "hatu": "HATU",
    "hbtu": "HBTU",
    "tbtu": "TBTU",
    "pybop": "PyBOP",
    "et3n": "CCN(CC)CC",
    "triethylamine": "CCN(CC)CC",
    "dipea": "CCN(C(C)C)C(C)C",
    "hunigs_base": "CCN(C(C)C)C(C)C",
    "nmm": "CN1CCOCC1",
    "pyridine": "c1ccncc1",
    "dmf": "CN(C)C=O",
    "dcm": "ClCCl",
    "mecn": "CC#N",
    "toluene": "Cc1ccccc1",
    "methanol": "CO",
    "meoh": "CO",
    "ethanol": "CCO",
    "etoh": "CCO",
    "water": "O",
    "h2so4": "O=S(=O)(O)O",
    "tfa": "O=C(O)C(F)(F)F",
    "hcl": "Cl",
    "bbr3": "B(Br)(Br)Br",
    "bcl3": "B(Cl)(Cl)Cl",
    "tbaf": "CCCC[N+](CCCC)(CCCC)CCCC.[F-]",
    "naf": "[Na+].[F-]",
    "kf": "[K+].[F-]",
    "csf": "[Cs+].[F-]",
    "naoh": "[Na+].[OH-]",
    "koh": "[K+].[OH-]",
    "lioh": "[Li+].[OH-]",
    "k2co3": "[K+].[K+].[O-]C(=O)[O-]",
    "cs2co3": "[Cs+].[Cs+].[O-]C(=O)[O-]",
    "nabh4": "[Na+].[BH4-]",
    "sodium_borohydride": "[Na+].[BH4-]",
    "nabh3cn": "[Na+].[BH3-]C#N",
    "sodium_cyanoborohydride": "[Na+].[BH3-]C#N",
    "nabh(oac)3": "NaBH(OAc)3",
    "stabh": "NaBH(OAc)3",
    "lialh4": "[Li+].[AlH4-]",
    "lah": "[Li+].[AlH4-]",
    "dibal": "DIBAL-H",
    "dibal-h": "DIBAL-H",
    "h2": "[H][H]",
    "pd": "[Pd]",
    "pd/c": "[Pd]",
    "pt": "[Pt]",
    "ni": "[Ni]",
    "cu": "[Cu]",
    "cui": "[Cu]I",
    "mcpba": "O=C(OO)c1cccc(Cl)c1",
    "nbs": "O=C1CCC(=O)N1Br",
    "ncs": "O=C1CCC(=O)N1Cl",
    "nis": "O=C1CCC(=O)N1I",
    "nacn": "[Na+].[C-]#N",
    "kcn": "[K+].[C-]#N",
    "cucn": "[Cu]C#N",
    "zn(cn)2": "N#C[Zn]C#N",
}

TEXTUAL_REAGENT_TOKENS = {
    "HOBt", "HOAt", "Oxyma", "HATU", "HBTU", "TBTU", "PyBOP",
    "DIBAL-H", "NaBH(OAc)3",
}

DISPLAY_ALIASES: dict[str, str] = {
    "c1ccc(P(c2ccccc2)c2ccccc2)cc1": "PPh3",
    "CC(C)OC(=O)N=NC(=O)OC(C)C": "DIAD",
    "CCOC(=O)N=NC(=O)OCC": "DEAD",
    "C1CCOC1": "THF",
    "C1CCC(CC1)N=C=NC2CCCCC2": "DCC",
    "CN(C)c1ccncc1": "DMAP",
    "O=S(Cl)Cl": "SOCl2",
    "CCN=C=NCCCN(C)C": "EDC",
    "CC(C)N=C=NC(C)C": "DIC",
    "CCN(CC)CC": "Et3N",
    "CCN(C(C)C)C(C)C": "DIPEA",
    "CN1CCOCC1": "NMM",
    "c1ccncc1": "pyridine",
    "CN(C)C=O": "DMF",
    "ClCCl": "DCM",
    "CC#N": "MeCN",
    "Cc1ccccc1": "toluene",
    "CO": "MeOH",
    "CCO": "EtOH",
    "O": "water",
    "O=S(=O)(O)O": "H2SO4",
    "O=C(O)C(F)(F)F": "TFA",
    "B(Br)(Br)Br": "BBr3",
    "B(Cl)(Cl)Cl": "BCl3",
    "[H][H]": "H2",
    "[Pd]": "Pd",
    "[Pt]": "Pt",
    "[Ni]": "Ni",
    "[Cu]": "Cu",
    "O=C(OO)c1cccc(Cl)c1": "mCPBA",
    "O=C1CCC(=O)N1Br": "NBS",
    "O=C1CCC(=O)N1Cl": "NCS",
    "O=C1CCC(=O)N1I": "NIS",
    "HATU": "HATU",
    "HBTU": "HBTU",
    "TBTU": "TBTU",
    "PyBOP": "PyBOP",
    "DIBAL-H": "DIBAL-H",
    "NaBH(OAc)3": "NaBH(OAc)3",
}

# Exact aliases supplement structural family recognition.  They intentionally
# describe reagent *families*, not named-reaction conclusions.
FAMILY_ALIASES: dict[str, set[str]] = {
    "phosphine": {"pph3", "triphenylphosphine"},
    "azodicarboxylate": {"diad", "dead"},
    "carbodiimide": {"dcc", "edc", "edci", "edic", "dic"},
    "coupling_additive": {"hobt", "hoat", "oxyma"},
    "uronium_coupling_reagent": {"hatu", "hbtu", "tbtu"},
    "phosphonium_coupling_reagent": {"pybop"},
    "tertiary_amine_base": {"et3n", "triethylamine", "dipea", "hunigs_base", "nmm", "pyridine"},
    "strong_acid": {"h2so4", "tfa", "hcl"},
    "lewis_acid_dealkylation": {"bbr3", "bcl3"},
    "fluoride_source": {"tbaf", "naf", "kf", "csf"},
    "hydroxide_base": {"naoh", "koh", "lioh"},
    "carbonate_base": {"k2co3", "cs2co3"},
    "borohydride": {"nabh4", "sodium_borohydride", "nabh3cn", "sodium_cyanoborohydride", "nabh(oac)3", "stabh"},
    "aluminum_hydride": {"lialh4", "lah", "dibal", "dibal-h"},
    "hydrogen": {"h2"},
    "palladium_catalyst": {"pd", "pd/c"},
    "platinum_catalyst": {"pt"},
    "nickel_catalyst": {"ni"},
    "copper_catalyst": {"cu", "cui"},
    "peracid_oxidant": {"mcpba"},
    "n_halosuccinimide": {"nbs", "ncs", "nis"},
    "cyanide_source": {"nacn", "kcn", "cucn", "zn(cn)2"},
}


def _rdkit_canonical(smiles: str) -> str:
    if smiles in TEXTUAL_REAGENT_TOKENS:
        return smiles
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        pass
    return smiles.strip()


@lru_cache(maxsize=65536)
def resolve_reagent_token(token: str) -> str:
    token = token.strip()
    if not token:
        return ""
    raw = REAGENT_ALIASES.get(token.lower(), token)
    return _rdkit_canonical(raw)


def resolve_reagent_components(token: str) -> list[str]:
    """Resolve a known salt/reagent alias into observed molecular components.

    This is important for matching CLI aliases such as TBAF/NaBH4 against
    reaction SMILES where ionic components are stored separately.  Arbitrary
    user-supplied dotted SMILES are still kept intact by parse_reagent_args.
    """
    text=str(token).strip()
    raw=REAGENT_ALIASES.get(text.lower())
    if raw is not None and "." in raw:
        out=[]
        for part in raw.split("."):
            c=_rdkit_canonical(part)
            if c and c not in out:
                out.append(c)
        return out
    c=resolve_reagent_token(text)
    return [c] if c else []


def parse_reagent_args(values: Iterable[str] | None) -> list[str]:
    """Normalize repeatable CLI/query reagent values conservatively."""
    out: list[str] = []
    for value in values or []:
        if value is None:
            continue
        chunks: list[str] = []
        for part in str(value).replace(";", ",").split(","):
            part = part.strip()
            if not part:
                continue
            dotted = part.split(".")
            if len(dotted) > 1 and all(x.lower() in REAGENT_ALIASES for x in dotted):
                chunks.extend(dotted)
            else:
                chunks.append(part)
        for chunk in chunks:
            for c in resolve_reagent_components(chunk):
                if c and c not in out:
                    out.append(c)
    return out


def extract_agent_field(reaction: str) -> str:
    parts = str(reaction).split(">")
    if len(parts) != 3:
        return ""
    return parts[1].strip()


def extract_reagents_from_reaction(reaction: str) -> list[str]:
    field = extract_agent_field(reaction)
    if not field:
        return []
    out: list[str] = []
    for token in field.split("."):
        c = resolve_reagent_token(token)
        if c and c not in out:
            out.append(c)
    return out


def reagent_label(canonical_smiles: str) -> str:
    c = resolve_reagent_token(canonical_smiles)
    for raw, label in DISPLAY_ALIASES.items():
        if resolve_reagent_token(raw) == c:
            return label
    return canonical_smiles


def _has_smarts(smiles: str, smarts: str) -> bool:
    if smiles in TEXTUAL_REAGENT_TOKENS:
        return False
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        pat = Chem.MolFromSmarts(smarts)
        return bool(mol is not None and pat is not None and mol.HasSubstructMatch(pat))
    except Exception:
        return False


@lru_cache(maxsize=65536)
def reagent_families(smiles_or_alias: str) -> tuple[str, ...]:
    """Return chemistry families without inferring a reaction method."""
    raw = str(smiles_or_alias).strip()
    low = raw.lower()
    c = resolve_reagent_token(raw)
    fam: set[str] = set()
    for family, aliases in FAMILY_ALIASES.items():
        if low in aliases or any(resolve_reagent_token(a) == c for a in aliases):
            fam.add(family)

    # Structural fallbacks recognize stored reaction components independent of
    # spelling/alias choice.
    if _has_smarts(c, "[N]=[C]=[N]"):
        fam.add("carbodiimide")
    if _has_smarts(c, "[P;X3]([#6])([#6])[#6]"):
        fam.add("phosphine")
    if _has_smarts(c, "N=NC(=O)O") and _has_smarts(c, "OC(=O)N=N"):
        fam.add("azodicarboxylate")
    if c == "[F-]" or _has_smarts(c, "[F-]"):
        fam.add("fluoride_source")
    if c == "[OH-]" or _has_smarts(c, "[OH-]"):
        fam.add("hydroxide_base")
    if _has_smarts(c, "[O-]C(=O)[O-]"):
        fam.add("carbonate_base")
    if "[BH4-]" in raw or "[BH3-]" in raw or "NaBH(OAc)3" in raw:
        fam.add("borohydride")
    if "[AlH4-]" in raw or raw == "DIBAL-H":
        fam.add("aluminum_hydride")
    if c in {"[Pd]", "[Pd+2]", "[Pd+]"} or "[Pd" in raw:
        fam.add("palladium_catalyst")
    if c in {"[Pt]", "[Pt+2]"} or "[Pt" in raw:
        fam.add("platinum_catalyst")
    if c in {"[Ni]", "[Ni+2]"} or "[Ni" in raw:
        fam.add("nickel_catalyst")
    if c in {"[Cu]", "[Cu+]", "[Cu+2]"} or "[Cu" in raw:
        fam.add("copper_catalyst")
    if _has_smarts(c, "[#6][B]") or _has_smarts(c, "[B][#6]"):
        fam.add("organoboron_or_boron_reagent")
    if _has_smarts(c, "[N+;X4]"):
        fam.add("quaternary_ammonium")
    if _has_smarts(c, "[C-]#N") or c in {"N#C[Cu]", "N#C[Zn]C#N"}:
        fam.add("cyanide_source")
    if _has_smarts(c, "C(=O)OO"):
        fam.add("peracid_oxidant")
    if _has_smarts(c, "[N]([Cl,Br,I])C(=O)"):
        fam.add("n_halosuccinimide")
    if c == "[H][H]":
        fam.add("hydrogen")
    if c in {"[Na+]", "[K+]", "[Li+]", "[Cs+]"}:
        fam.add("alkali_counterion")
    if c in {"[Cl-]", "[Br-]", "[I-]", "[F-]"}:
        fam.add("halide_counterion")
    return tuple(sorted(fam))


def reconstruct_ionic_reagent_evidence(reagents: Iterable[str]) -> list[str]:
    """Conservative set-level reconstruction for split ionic components.

    These are evidence signatures, not fabricated molecular structures or
    stoichiometric assignments.  They let `[K+]` + carbonate carry more meaning
    than either component alone while preserving the observed components.
    """
    vals=list(dict.fromkeys(resolve_reagent_token(x) for x in reagents if x))
    fam=_families_for_set(vals)
    sig=[]
    has_alkali=bool(fam & {"alkali_counterion"})
    if has_alkali and "hydroxide_base" in fam:
        sig.append("alkali_hydroxide_like")
    if has_alkali and "carbonate_base" in fam:
        sig.append("alkali_carbonate_like")
    if has_alkali and "borohydride" in fam:
        sig.append("alkali_borohydride_like")
    if "fluoride_source" in fam and any(_has_smarts(x, "[N+;X4]") for x in vals):
        sig.append("quaternary_ammonium_fluoride_like")
    return sig


def _families_for_set(reagents: Iterable[str]) -> set[str]:
    out:set[str]=set()
    for r in reagents:
        out.update(reagent_families(r))
    return out


@dataclass(frozen=True)
class ReagentRecord:
    reagent_id: str
    canonical_smiles: str
    canonical_name: str
    reagent_families: list[str]
    informativeness: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def describe_reagent(smiles_or_alias: str) -> ReagentRecord:
    c = resolve_reagent_token(smiles_or_alias)
    families = list(reagent_families(smiles_or_alias))
    high = {
        "azodicarboxylate", "carbodiimide", "uronium_coupling_reagent",
        "phosphonium_coupling_reagent", "borohydride", "aluminum_hydride",
        "palladium_catalyst", "copper_catalyst", "cyanide_source",
        "lewis_acid_dealkylation", "peracid_oxidant",
    }
    low_families = {"alkali_counterion", "halide_counterion"}
    if set(families) & high:
        info = "high"
    elif families and set(families) <= low_families:
        info = "low"
    elif families:
        info = "moderate"
    else:
        info = "unknown"
    return ReagentRecord(
        reagent_id=c,
        canonical_smiles=c,
        canonical_name=reagent_label(c),
        reagent_families=families,
        informativeness=info,
    )


@dataclass(frozen=True)
class ReagentMatch:
    score: float
    query_coverage: float
    jaccard: float
    exact_query_coverage: float
    exact_jaccard: float
    n_query: int
    n_candidate: int
    n_overlap: int
    n_family_overlap: int
    matched: list[str]
    family_matches: list[dict[str, Any]]
    query_reagents: list[str]
    candidate_reagents: list[str]
    match_mode: str = "exact+family"

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["matched_labels"] = [reagent_label(x) for x in self.matched]
        d["family_match_labels"] = [
            f"{reagent_label(x['query'])}~{reagent_label(x['candidate'])} [{x['family']}]"
            for x in self.family_matches
        ]
        d["query_labels"] = [reagent_label(x) for x in self.query_reagents]
        d["candidate_labels"] = [reagent_label(x) for x in self.candidate_reagents]
        return d


# Family equivalence is deliberately weaker than exact molecular identity.
# Only method-informative families are included.  Counterions and generic
# spectators never create family-level reagent similarity.
FAMILY_EQUIVALENCE_WEIGHTS: dict[str, float] = {
    "azodicarboxylate": 0.92,
    "carbodiimide": 0.88,
    "coupling_additive": 0.85,
    "uronium_coupling_reagent": 0.90,
    "phosphonium_coupling_reagent": 0.90,
    "phosphine": 0.85,
    "tertiary_amine_base": 0.70,
    "strong_acid": 0.65,
    "lewis_acid_dealkylation": 0.85,
    "fluoride_source": 0.82,
    "hydroxide_base": 0.75,
    "carbonate_base": 0.75,
    "borohydride": 0.86,
    "aluminum_hydride": 0.82,
    "hydrogen": 0.80,
    "palladium_catalyst": 0.82,
    "platinum_catalyst": 0.82,
    "nickel_catalyst": 0.82,
    "copper_catalyst": 0.82,
    "peracid_oxidant": 0.86,
    "n_halosuccinimide": 0.86,
    "cyanide_source": 0.86,
    "organoboron_or_boron_reagent": 0.80,
}


def _reagent_pair_similarity(query: str, candidate: str) -> tuple[float, str | None]:
    if query == candidate:
        return 1.0, None
    qf = set(reagent_families(query))
    cf = set(reagent_families(candidate))
    shared = qf & cf & set(FAMILY_EQUIVALENCE_WEIGHTS)
    if not shared:
        return 0.0, None
    family = max(shared, key=lambda f: (FAMILY_EQUIVALENCE_WEIGHTS[f], f))
    return float(FAMILY_EQUIVALENCE_WEIGHTS[family]), family


def _soft_one_to_one_matches(query: list[str], candidate: list[str]) -> list[tuple[int, int, float, str | None]]:
    """Deterministic one-to-one exact/family matching.

    Candidate condition lists are normally short.  A descending greedy matching
    is sufficient here and, importantly, prevents one candidate component from
    satisfying several query components at once.
    """
    pairs: list[tuple[float, int, int, str | None]] = []
    for qi, q in enumerate(query):
        for ci, c in enumerate(candidate):
            score, family = _reagent_pair_similarity(q, c)
            if score > 0.0:
                pairs.append((score, qi, ci, family))
    pairs.sort(key=lambda x: (-x[0], x[1], x[2], x[3] or ""))
    used_q: set[int] = set()
    used_c: set[int] = set()
    out: list[tuple[int, int, float, str | None]] = []
    for score, qi, ci, family in pairs:
        if qi in used_q or ci in used_c:
            continue
        used_q.add(qi); used_c.add(ci)
        out.append((qi, ci, float(score), family))
    return out


def reagent_match(query_reagents: Iterable[str], candidate_reagents: Iterable[str]) -> ReagentMatch:
    q = list(dict.fromkeys(resolve_reagent_token(x) for x in query_reagents if x))
    c = list(dict.fromkeys(resolve_reagent_token(x) for x in candidate_reagents if x))
    qs, cs = set(q), set(c)
    exact = qs & cs
    exact_union = qs | cs
    exact_coverage = len(exact) / len(qs) if qs else 0.0
    exact_jaccard = len(exact) / len(exact_union) if exact_union else 0.0

    matched_pairs = _soft_one_to_one_matches(q, c)
    credit = sum(x[2] for x in matched_pairs)
    coverage = credit / len(q) if q else 0.0
    soft_union = len(q) + len(c) - credit
    jaccard = credit / soft_union if soft_union > 0 else 0.0
    score = 0.80 * coverage + 0.20 * jaccard if q else 0.0

    family_matches: list[dict[str, Any]] = []
    for qi, ci, pair_score, family in matched_pairs:
        if family is None:
            continue
        family_matches.append({
            "query": q[qi],
            "candidate": c[ci],
            "family": family,
            "similarity": pair_score,
        })

    return ReagentMatch(
        score=float(score),
        query_coverage=float(coverage),
        jaccard=float(jaccard),
        exact_query_coverage=float(exact_coverage),
        exact_jaccard=float(exact_jaccard),
        n_query=len(q),
        n_candidate=len(c),
        n_overlap=len(exact),
        n_family_overlap=len(family_matches),
        matched=sorted(exact),
        family_matches=family_matches,
        query_reagents=q,
        candidate_reagents=c,
    )


def infer_method_labels(primary_event: str | None, candidate_reagents: Iterable[str], reaction: str = "") -> list[dict[str, Any]]:
    """Backward-compatible wrapper around the v6.5.3 method classifier."""
    from .methods import infer_methods
    return infer_methods(primary_event, reaction, candidate_reagents)


def profile_reagents(reactions: Iterable[str], *, primary_event: str | None = None, top_n: int = 25) -> dict[str, Any]:
    from .methods import build_method_landscape
    rows = [{"reaction": r, "membership_occurrence_count": 1} for r in reactions]
    land = build_method_landscape(rows, primary_event=primary_event, top_reagents=top_n, top_signatures=top_n)
    return {
        "n_reactions": land["evidence_pool"]["candidate_reactions"],
        "n_with_agent_field": land["evidence_pool"]["with_agent_field"],
        "agent_coverage_fraction": land["evidence_pool"]["agent_coverage_fraction"],
        "top_reagents": land["top_reagents"],
        "top_agent_signatures": land["top_agent_signatures"],
        "method_labels": [
            {
                "label": m["method_id"], "method_id": m["method_id"],
                "display_name": m["display_name"], "confidence": m["confidence"],
                "evidence": m["defining_evidence"], "count": m["support_reactions"],
                "fraction": m["support_fraction"],
            }
            for m in land["methods"]
        ],
        "method_landscape": land,
    }
