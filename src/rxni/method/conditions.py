from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Iterable

from .reagents import (
    extract_reagents_from_reaction,
    parse_reagent_args,
    reagent_label,
    resolve_reagent_token,
    resolve_reagent_components,
    reconstruct_ionic_reagent_evidence,
    reagent_match,
)
from .reagent_roles import partition_condition_components


def _dedupe(values: Iterable[str]) -> list[str]:
    out: list[str] = []
    for value in values:
        if value and value not in out:
            out.append(value)
    return out


def _canon_many(values: Iterable[str] | None) -> list[str]:
    return _dedupe(resolve_reagent_token(str(x)) for x in (values or []) if str(x).strip())


def _canon_reagent_many(values: Iterable[str] | None) -> list[str]:
    out: list[str] = []
    for value in values or []:
        if not str(value).strip():
            continue
        out.extend(resolve_reagent_components(str(value)))
    return _dedupe(out)


def _without(values: Iterable[str], excluded: Iterable[str]) -> list[str]:
    ex = set(excluded)
    return [x for x in values if x not in ex]


@dataclass(frozen=True)
class QueryConditionEvidence:
    """Canonical query-side condition evidence for Runtime-v6.5.3 PATCH 1.

    Structural interpretation never consumes this object.  It is constructed only
    after Mapping Correction + semantic resolution and is used by condition-aware
    ranking, method-landscape conditioning, and knowledge retrieval.

    PATCH 1 resolves obvious roles from mixed patent middle fields so solvents,
    catalysts, bases and additives no longer dilute core reagent similarity.
    """

    schema_version: str = "query_condition_evidence_v1_1"
    components: list[str] = field(default_factory=list)
    reagents: list[str] = field(default_factory=list)
    solvents: list[str] = field(default_factory=list)
    catalysts: list[str] = field(default_factory=list)
    bases: list[str] = field(default_factory=list)
    additives: list[str] = field(default_factory=list)
    acids: list[str] = field(default_factory=list)
    spectators: list[str] = field(default_factory=list)
    unresolved: list[str] = field(default_factory=list)
    reconstructed_reagent_signatures: list[str] = field(default_factory=list)
    temperature_c: float | None = None
    reaction_time_h: float | None = None
    atmosphere: str | None = None
    pressure: float | str | None = None
    sources: dict[str, Any] = field(default_factory=dict)

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
        for name in (
            "components", "reagents", "solvents", "catalysts", "bases",
            "additives", "acids", "spectators", "unresolved",
        ):
            d[f"{name[:-1] if name.endswith('s') else name}_labels" if name != "components" else "component_labels"] = [
                reagent_label(x) for x in getattr(self, name)
            ]
        # Stable historical field names.
        d["reagent_labels"] = [reagent_label(x) for x in self.reagents]
        d["solvent_labels"] = [reagent_label(x) for x in self.solvents]
        d["catalyst_labels"] = [reagent_label(x) for x in self.catalysts]
        d["base_labels"] = [reagent_label(x) for x in self.bases]
        d["additive_labels"] = [reagent_label(x) for x in self.additives]
        d["acid_labels"] = [reagent_label(x) for x in self.acids]
        return d

    @property
    def has_evidence(self) -> bool:
        return bool(
            self.active_components
            or self.solvents
            or self.temperature_c is not None
            or self.reaction_time_h is not None
            or self.atmosphere
            or self.pressure is not None
        )


@dataclass(frozen=True)
class ConditionMatch:
    score: float
    reagent_score: float | None
    solvent_score: float | None
    catalyst_score: float | None
    base_score: float | None
    additive_score: float | None
    acid_score: float | None
    dimensions_used: list[str]
    unavailable_dimensions: list[str]
    reagent_match: dict[str, Any] | None = None
    solvent_match: dict[str, Any] | None = None
    catalyst_match: dict[str, Any] | None = None
    base_match: dict[str, Any] | None = None
    additive_match: dict[str, Any] | None = None
    acid_match: dict[str, Any] | None = None
    candidate_roles: dict[str, Any] | None = None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


_DIMENSION_WEIGHTS = {
    "reagents": 0.40,
    "solvents": 0.10,
    "catalysts": 0.15,
    "bases": 0.15,
    "additives": 0.10,
    "acids": 0.10,
}


def condition_match(query: QueryConditionEvidence, candidate_components: Iterable[str]) -> ConditionMatch:
    """Role-aware comparison of query evidence with a V5 middle-field pool.

    Exact molecular identity receives full credit.  Within a role bucket,
    chemically informative reagent-family equivalents receive partial credit
    (for example DIAD~DEAD as azodicarboxylates).  Obvious solvents are compared
    only against query solvents and do not dilute reagent Jaccard similarity.
    """
    roles = partition_condition_components(candidate_components)
    parts: list[tuple[float, float]] = []
    used: list[str] = []
    matches: dict[str, dict[str, Any] | None] = {
        "reagents": None, "solvents": None, "catalysts": None,
        "bases": None, "additives": None, "acids": None,
    }

    for name in ("reagents", "solvents", "catalysts", "bases", "additives", "acids"):
        qvals = list(getattr(query, name) or [])
        if not qvals:
            continue
        cvals = list(getattr(roles, name) or [])
        match = reagent_match(qvals, cvals)
        matches[name] = match.as_dict()
        parts.append((_DIMENSION_WEIGHTS[name], match.score))
        used.append(name)

    denom = sum(w for w, _ in parts)
    score = sum(w * x for w, x in parts) / denom if denom else 0.0
    unavailable: list[str] = []
    if query.temperature_c is not None:
        unavailable.append("temperature_c")
    if query.reaction_time_h is not None:
        unavailable.append("reaction_time_h")
    if query.atmosphere:
        unavailable.append("atmosphere")
    if query.pressure is not None:
        unavailable.append("pressure")

    return ConditionMatch(
        score=float(score),
        reagent_score=None if matches["reagents"] is None else float(matches["reagents"]["score"]),
        solvent_score=None if matches["solvents"] is None else float(matches["solvents"]["score"]),
        catalyst_score=None if matches["catalysts"] is None else float(matches["catalysts"]["score"]),
        base_score=None if matches["bases"] is None else float(matches["bases"]["score"]),
        additive_score=None if matches["additives"] is None else float(matches["additives"]["score"]),
        acid_score=None if matches["acids"] is None else float(matches["acids"]["score"]),
        dimensions_used=used,
        unavailable_dimensions=unavailable,
        reagent_match=matches["reagents"],
        solvent_match=matches["solvents"],
        catalyst_match=matches["catalysts"],
        base_match=matches["bases"],
        additive_match=matches["additives"],
        acid_match=matches["acids"],
        candidate_roles=roles.as_dict(),
    )


def build_query_condition_evidence(
    reaction: str,
    *,
    cli_reagents: Iterable[str] | None = None,
    solvents: Iterable[str] | None = None,
    catalysts: Iterable[str] | None = None,
    bases: Iterable[str] | None = None,
    additives: Iterable[str] | None = None,
    acids: Iterable[str] | None = None,
    temperature_c: float | None = None,
    reaction_time_h: float | None = None,
    atmosphere: str | None = None,
    pressure: float | str | None = None,
    orderly_record: Any | None = None,
) -> QueryConditionEvidence:
    """Merge and role-resolve condition evidence without feeding it upstream.

    Provenance is retained per source.  ORDerly fields are accepted either as an
    OrderlyRecord instance or a dictionary with matching field names.
    """

    embedded = extract_reagents_from_reaction(reaction)
    cli = parse_reagent_args(cli_reagents)

    orderly_agents: list[str] = []
    orderly_solvents: list[str] = []
    orderly_temperature = None
    orderly_time = None
    orderly_atmosphere = None
    orderly_pressure = None
    if orderly_record is not None:
        getter = orderly_record.get if isinstance(orderly_record, dict) else lambda k, default=None: getattr(orderly_record, k, default)
        orderly_agents = _canon_reagent_many(getter("agents", []) or [])
        orderly_solvents = _canon_many(getter("solvents", []) or [])
        orderly_temperature = getter("temperature_c", None)
        orderly_time = getter("reaction_time_h", None)
        orderly_atmosphere = getter("atmosphere", None)
        orderly_pressure = getter("pressure", None)

    explicit_solvents = _canon_many(solvents)
    explicit_catalysts = _canon_reagent_many(catalysts)
    explicit_bases = _canon_reagent_many(bases)
    explicit_additives = _canon_reagent_many(additives)
    explicit_acids = _canon_reagent_many(acids)

    merged_components = _dedupe([
        *embedded, *cli, *orderly_agents,
        *explicit_catalysts, *explicit_bases, *explicit_additives, *explicit_acids,
    ])
    resolved = partition_condition_components(merged_components)

    # Explicit role arguments override a generic classification bucket while the
    # original molecular component remains present in ``components``/provenance.
    forced = set(explicit_catalysts + explicit_bases + explicit_additives + explicit_acids)
    core_reagents = _without(resolved.reagents, forced)
    role_catalysts = _dedupe([*resolved.catalysts, *explicit_catalysts])
    role_bases = _dedupe([*resolved.bases, *explicit_bases])
    role_additives = _dedupe([*resolved.additives, *explicit_additives])
    role_acids = _dedupe([*resolved.acids, *explicit_acids])
    merged_solvents = _dedupe([*resolved.solvents, *explicit_solvents, *orderly_solvents])

    temp = temperature_c if temperature_c is not None else orderly_temperature
    rtime = reaction_time_h if reaction_time_h is not None else orderly_time
    atm = atmosphere if atmosphere is not None else orderly_atmosphere
    pres = pressure if pressure is not None else orderly_pressure

    return QueryConditionEvidence(
        components=merged_components,
        reagents=core_reagents,
        solvents=merged_solvents,
        catalysts=role_catalysts,
        bases=role_bases,
        additives=role_additives,
        acids=role_acids,
        spectators=list(resolved.spectators),
        unresolved=list(resolved.unresolved),
        reconstructed_reagent_signatures=reconstruct_ionic_reagent_evidence(merged_components),
        temperature_c=None if temp is None else float(temp),
        reaction_time_h=None if rtime is None else float(rtime),
        atmosphere=None if atm in (None, "") else str(atm),
        pressure=pres,
        sources={
            "reaction_agent_field": embedded,
            "cli_reagents": cli,
            "explicit_solvents": explicit_solvents,
            "explicit_catalysts": explicit_catalysts,
            "explicit_bases": explicit_bases,
            "explicit_additives": explicit_additives,
            "explicit_acids": explicit_acids,
            "orderly_agents": orderly_agents,
            "orderly_solvents": orderly_solvents,
            "role_resolution": resolved.as_dict(),
            "temperature": (
                "explicit" if temperature_c is not None
                else "orderly" if orderly_temperature is not None
                else None
            ),
            "reaction_time": (
                "explicit" if reaction_time_h is not None
                else "orderly" if orderly_time is not None
                else None
            ),
            "atmosphere": (
                "explicit" if atmosphere is not None
                else "orderly" if orderly_atmosphere is not None
                else None
            ),
            "pressure": (
                "explicit" if pressure is not None
                else "orderly" if orderly_pressure is not None
                else None
            ),
        },
    )
