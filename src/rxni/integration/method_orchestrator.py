"""Runtime-v6.7-local post-retrieval condition/method/ChemKN orchestration.

M7E2 intentionally does not implement structural candidate generation.
Its input is a structurally retrieved member envelope produced by the frozen
the preserved Retrieval-V2 layer.

Ownership:
    runtime_v6_7.retrieval -> local copy of frozen Retrieval-V2 algorithms
    runtime_v6_7.method    -> local condition evidence, compatibility, method landscapes
    runtime_v6_7.chemkn    -> local exact/related knowledge enrichment
    this module            -> orchestration only

No Runtime-v6.5 import is permitted.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Iterable

from rxni.method.conditions import (
    build_query_condition_evidence,
    condition_match,
)
from rxni.method.reagent_roles import (
    partition_condition_components,
)
from rxni.method.methods import (
    infer_methods,
    build_method_landscape,
    build_conditioned_method_landscape,
)
from rxni.chemkn.bridge import enrich_method_landscape


# --- mechanically extracted historical orchestration helpers ---
def _score_value(scores: dict[str, Any], mode: str) -> float:
    aliases = {
        "local": ("local", "local_score"),
        "multi": ("multi", "multi_score"),
        "semantic": ("semantic", "semantic_score"),
        "sem_multi": ("sem_multi", "sem_multi_score"),
        "reagent_conditioned": ("reagent_conditioned", "conditioned"),
        "conditioned": ("conditioned", "reagent_conditioned"),
    }
    if mode == "exact_edit":
        return 1.0
    for k in aliases.get(mode, (mode,)):
        if k in scores and scores[k] is not None:
            return float(scores[k])
    raise KeyError(f"member_scores() did not return a score for mode={mode!r}; keys={sorted(scores)}")

def _semantic_primary_event(record: dict[str, Any]) -> str | None:
    sem = record.get("semantics") or {}
    events = sem.get("primary_events") or []
    if events and isinstance(events[0], dict):
        return events[0].get("event") or events[0].get("subtype")
    prof = record.get("semantic_profile") or {}
    for key in prof:
        if str(key).startswith("P:"):
            return str(key).split(":", 1)[1]
    return None

def _condition_relation(query_evidence, cm, candidate_components: list[str]) -> str:
    if not query_evidence.has_evidence:
        return "not_conditioned"
    if not candidate_components:
        return "no_candidate_conditions"
    matches = [
        cm.reagent_match, cm.solvent_match, cm.catalyst_match,
        cm.base_match, cm.additive_match, cm.acid_match,
    ]
    matches = [m for m in matches if m is not None]
    if not matches:
        return "different_conditions"
    exact = all(float(m.get("exact_query_coverage", 0.0)) >= 0.999999 for m in matches)
    if exact:
        return "exact_condition_match"
    full_soft = all(float(m.get("query_coverage", 0.0)) >= 0.85 for m in matches)
    any_family = any(bool(m.get("family_matches")) for m in matches)
    if full_soft and any_family:
        return "family_equivalent_condition_match"
    if cm.score > 0.0:
        return "partial_condition_match"
    return "different_conditions"

def _candidate_category(query_primary: str | None, candidate_primary: str | None, condition_relation: str) -> tuple[str, str]:
    if query_primary and candidate_primary == query_primary:
        chemistry_relation = "same_primary_event"
    elif candidate_primary:
        chemistry_relation = "related_different_primary_event"
    else:
        chemistry_relation = "candidate_semantics_unresolved"
    return chemistry_relation, f"{chemistry_relation}__{condition_relation}"
# --- end mechanically extracted helpers ---


CONDITION_SCORE_KEYS = (
    "condition_score",
    "conditioned",
    "reagent_conditioned",
    "reagent_score",
    "reagent_query_coverage",
    "reagent_exact_query_coverage",
    "reagent_jaccard",
    "reagent_exact_jaccard",
)


def _candidate_components(member: dict[str, Any]) -> list[str]:
    metadata = member.get("metadata") or {}
    values = metadata.get("reagents") or member.get("reagents") or []
    return [str(x) for x in values if str(x)]


def _candidate_primary_event(
    member: dict[str, Any],
) -> str | None:
    metadata = member.get("metadata") or {}
    value = metadata.get("candidate_primary_event")
    if value:
        return str(value)

    # Optional fallback for richer native retrieval envelopes.
    semantics = (
        metadata.get("semantics")
        or member.get("semantics")
        or {}
    )
    try:
        return _semantic_primary_event(semantics)
    except Exception:
        return None


def _chemistry_relation(
    query_primary_event: str | None,
    candidate_primary_event: str | None,
) -> str:
    if (
        query_primary_event
        and candidate_primary_event
        and query_primary_event == candidate_primary_event
    ):
        return "same_primary_event"
    if candidate_primary_event:
        return "different_primary_event"
    return "unknown_primary_event"


def condition_member(
    member: dict[str, Any],
    *,
    primary_event: str | None,
    condition_evidence: Any,
    mode: str = "sem_multi",
    chemistry_weight: float = 0.80,
    reagent_weight: float = 0.20,
) -> dict[str, Any]:
    """Apply query-condition evidence to one already-retrieved member.

    Structural scores and candidate identity are copied, never recomputed.
    """
    out = deepcopy(member)
    scores = dict(out.get("scores") or {})
    metadata = dict(out.get("metadata") or {})

    candidate_components = _candidate_components(out)
    candidate_roles = partition_condition_components(
        candidate_components
    )
    metadata["condition_roles"] = candidate_roles.as_dict()

    base = float(_score_value(scores, mode))

    if condition_evidence.has_evidence:
        cm = condition_match(
            condition_evidence,
            candidate_components,
        )
        cm_dict = cm.as_dict()

        conditioned = (
            float(chemistry_weight) * base
            + float(reagent_weight) * float(cm.score)
        )

        scores["condition_score"] = float(cm.score)
        scores["conditioned"] = float(conditioned)
        scores["reagent_conditioned"] = float(conditioned)

        reagent_match = cm_dict.get("reagent_match")
        metadata["condition_match"] = cm_dict
        metadata["reagent_match"] = reagent_match

        if isinstance(reagent_match, dict):
            scores["reagent_score"] = float(
                reagent_match.get("score") or 0.0
            )
            scores["reagent_query_coverage"] = float(
                reagent_match.get("query_coverage") or 0.0
            )
            scores["reagent_exact_query_coverage"] = float(
                reagent_match.get(
                    "exact_query_coverage"
                ) or 0.0
            )
            scores["reagent_jaccard"] = float(
                reagent_match.get("jaccard") or 0.0
            )
            scores["reagent_exact_jaccard"] = float(
                reagent_match.get("exact_jaccard") or 0.0
            )

        condition_relation = _condition_relation(
            condition_evidence,
            cm,
            candidate_components,
        )
    else:
        condition_relation = "no_query_condition_evidence"

    candidate_event = _candidate_primary_event(out)
    chemistry_relation, category = _candidate_category(
        primary_event,
        candidate_event,
        condition_relation,
    )

    metadata["chemistry_relation"] = chemistry_relation
    metadata["condition_relation"] = condition_relation
    metadata["candidate_category"] = category

    predictions = infer_methods(
        primary_event,
        str(out.get("reaction") or ""),
        candidate_roles.active_components,
        solvents=candidate_roles.solvents,
    )
    metadata["method_predictions"] = predictions
    metadata["method_labels"] = deepcopy(predictions)

    out["scores"] = scores
    out["metadata"] = metadata
    return out


def condition_members(
    members: Iterable[dict[str, Any]],
    *,
    primary_event: str | None,
    condition_evidence: Any,
    mode: str = "sem_multi",
    chemistry_weight: float = 0.80,
    reagent_weight: float = 0.20,
    sort_conditioned: bool = True,
) -> list[dict[str, Any]]:
    out = [
        condition_member(
            member,
            primary_event=primary_event,
            condition_evidence=condition_evidence,
            mode=mode,
            chemistry_weight=chemistry_weight,
            reagent_weight=reagent_weight,
        )
        for member in members
    ]

    if sort_conditioned and condition_evidence.has_evidence:
        out.sort(
            key=lambda row: float(
                (row.get("scores") or {}).get(
                    "conditioned",
                    _score_value(row.get("scores") or {}, mode),
                )
            ),
            reverse=True,
        )

    return out


def _landscape_row(
    member: dict[str, Any],
    *,
    score: float,
) -> dict[str, Any]:
    metadata = member.get("metadata") or {}
    return {
        "reaction_key": member.get("reaction_key"),
        "reaction": member.get("reaction"),
        "reaction_id": metadata.get("reaction_id"),
        "template_hash": member.get("template_hash"),
        "membership_occurrence_count": (
            member.get("membership_occurrence_count") or 1
        ),
        "chemistry_score": float(score),
    }


def build_landscapes(
    members: Iterable[dict[str, Any]],
    *,
    primary_event: str | None,
    query_reaction: str,
    condition_evidence: Any,
    mode: str = "sem_multi",
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    members = list(members)
    if not members:
        return None, None

    global_members = sorted(
        members,
        key=lambda member: (
            -float(
                _score_value(
                    member.get("scores") or {},
                    mode,
                )
            ),
            str(member.get("reaction_key") or ""),
        ),
    )

    global_rows = [
        _landscape_row(
            member,
            score=float(
                _score_value(
                    member.get("scores") or {},
                    mode,
                )
            ),
        )
        for member in global_members
    ]

    global_landscape = build_method_landscape(
        global_rows,
        primary_event=primary_event,
        landscape_type="global",
    )

    conditioned_landscape = None
    if condition_evidence.has_evidence:
        conditioned_rows = [
            _landscape_row(
                member,
                score=float(
                    _score_value(
                        member.get("scores") or {},
                        "conditioned",
                    )
                ),
            )
            for member in members
        ]

        conditioned_landscape = (
            build_conditioned_method_landscape(
                conditioned_rows,
                primary_event=primary_event,
                query_reaction=query_reaction,
                condition_evidence=condition_evidence,
            )
        )

    return global_landscape, conditioned_landscape


def orchestrate_post_retrieval(
    members: Iterable[dict[str, Any]],
    *,
    reaction: str,
    primary_event: str | None,
    mode: str = "sem_multi",
    chemistry_weight: float = 0.80,
    reagent_weight: float = 0.20,
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
    chemkn_registry: Any | None = None,
    chemkn_link_store: Any | None = None,
    chemkn_include_candidates: bool = True,
    chemkn_max_candidates: int = 6,
    exact_method_threshold: float = 0.85,
) -> dict[str, Any]:
    """v6.7-local orchestration preserving the frozen retrieval-envelope contract."""

    evidence = build_query_condition_evidence(
        reaction,
        cli_reagents=cli_reagents,
        solvents=solvents,
        catalysts=catalysts,
        bases=bases,
        additives=additives,
        acids=acids,
        temperature_c=temperature_c,
        reaction_time_h=reaction_time_h,
        atmosphere=atmosphere,
        pressure=pressure,
        orderly_record=orderly_record,
    )

    query_method_predictions = infer_methods(
        primary_event,
        reaction,
        evidence.active_components,
        solvents=evidence.solvents,
    )

    input_members = [deepcopy(x) for x in members]
    conditioned_members = condition_members(
        input_members,
        primary_event=primary_event,
        condition_evidence=evidence,
        mode=mode,
        chemistry_weight=chemistry_weight,
        reagent_weight=reagent_weight,
        sort_conditioned=True,
    )

    global_landscape, conditioned_landscape = build_landscapes(
        conditioned_members,
        primary_event=primary_event,
        query_reaction=reaction,
        condition_evidence=evidence,
        mode=mode,
    )

    result = {
        "query_condition_evidence": evidence.as_dict(),
        "query_reagents": list(evidence.active_components),
        "query_method_predictions": query_method_predictions,
        "reagent_weight": (
            float(reagent_weight)
            if evidence.has_evidence
            else 0.0
        ),
        "condition_weight": (
            float(reagent_weight)
            if evidence.has_evidence
            else 0.0
        ),
        "ranking_mode": (
            "conditioned"
            if evidence.has_evidence
            else mode
        ),
        "global_method_landscape": global_landscape,
        "conditioned_method_landscape": conditioned_landscape,
        "method_profile": global_landscape,
        "members": conditioned_members,
        "mode": mode,
    }

    if chemkn_registry is not None:
        if global_landscape is not None:
            result["global_method_landscape"] = enrich_method_landscape(
                global_landscape,
                chemkn_registry,
                retrieval_result=result,
                link_store=chemkn_link_store,
                include_candidates=chemkn_include_candidates,
                max_candidates=chemkn_max_candidates,
                exact_method_threshold=exact_method_threshold,
            )
            result["method_profile"] = result[
                "global_method_landscape"
            ]

        if conditioned_landscape is not None:
            result["conditioned_method_landscape"] = enrich_method_landscape(
                conditioned_landscape,
                chemkn_registry,
                retrieval_result=result,
                link_store=chemkn_link_store,
                include_candidates=chemkn_include_candidates,
                max_candidates=chemkn_max_candidates,
                exact_method_threshold=exact_method_threshold,
            )

    return result
