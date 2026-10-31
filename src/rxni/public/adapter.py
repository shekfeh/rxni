"""Adapter from the native Runtime-v6.6 retrieval result to the public API."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from .models import (
    AnalysisOptions,
    AnalysisRequest,
    AnalysisResult,
    CandidateCategory,
    CandidateHit,
    ChemKNSection,
    ConditionsSection,
    DiagnosticsSection,
    MethodLandscapeView,
    ProvenanceSection,
    ReactionSection,
    RetrievalSection,
    RuntimeInfo,
    SemanticsSection,
    StructureSection,
    SummarySection,
    TemplateHit,
)


def _as_dict(obj: Any) -> dict[str, Any]:
    if obj is None:
        return {}
    if isinstance(obj, dict):
        return dict(obj)
    if isinstance(obj, Mapping):
        return dict(obj)
    if hasattr(obj, "as_dict"):
        value = obj.as_dict()
        if isinstance(value, Mapping):
            return dict(value)
    if hasattr(obj, "__dict__"):
        try:
            return dict(vars(obj))
        except Exception:
            pass
    return {}


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def _dict_list(value: Any) -> list[dict[str, Any]]:
    return [_as_dict(x) for x in _as_list(value) if _as_dict(x)]


def _float_dict(value: Any) -> dict[str, float]:
    out: dict[str, float] = {}
    for k, v in _as_dict(value).items():
        try:
            out[str(k)] = float(v)
        except Exception:
            continue
    return out


def _int_dict(value: Any) -> dict[str, int]:
    out: dict[str, int] = {}
    for k, v in _as_dict(value).items():
        try:
            out[str(k)] = int(v)
        except Exception:
            continue
    return out


def _optional_nonnegative_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        out = int(value)
    except Exception:
        return None
    return max(0, out)


def _extract_query_record(raw: dict[str, Any]) -> dict[str, Any]:
    return _as_dict(raw.get("query_record"))


def _template_score(row: dict[str, Any], ranking_mode: str) -> float | None:
    for key in (
        "score",
        "template_score",
        "aggregate_score",
        "final_score",
        ranking_mode,
    ):
        value = row.get(key)
        if value is not None:
            try:
                return float(value)
            except Exception:
                pass

    for key in (
        ranking_mode,
        "sem_multi",
        "multi",
        "semantic",
        "local",
        "exact_edit",
    ):
        value = _as_dict(row.get("scores")).get(key)
        if value is not None:
            try:
                return float(value)
            except Exception:
                pass

    return None


def _normalize_predictions(value: Any) -> list[dict[str, Any]]:
    rows = []
    for x in _as_list(value):
        d = _as_dict(x)
        if d:
            rows.append(d)
        elif x is not None:
            rows.append({"value": x})
    return rows


def _normalize_categories(value: Any) -> list[CandidateCategory]:
    out: list[CandidateCategory] = []
    for item in _as_list(value):
        d = _as_dict(item)
        if not d:
            continue
        try:
            count = int(d.get("count", 0))
        except Exception:
            count = 0
        fraction = d.get("fraction")
        try:
            fraction = float(fraction) if fraction is not None else None
        except Exception:
            fraction = None
        out.append(
            CandidateCategory(
                category=str(d.get("category") or "unclassified"),
                count=count,
                fraction=fraction,
            )
        )
    return out


def _normalize_candidates(
    value: Any,
    *,
    include: bool,
) -> list[CandidateHit]:
    if not include:
        return []

    out: list[CandidateHit] = []
    for rank, item in enumerate(_as_list(value), 1):
        d = _as_dict(item)
        md = _as_dict(d.get("metadata"))

        template_number = d.get("template_number")
        if template_number is None:
            template_number = d.get("retro_template_number")

        channels = [str(x) for x in _as_list(d.get("generation_channels"))]

        out.append(
            CandidateHit(
                rank=rank,
                reaction_key=(
                    str(d["reaction_key"])
                    if d.get("reaction_key") is not None
                    else None
                ),
                reaction=(
                    str(d.get("reaction") or d.get("forward_reaction_smiles"))
                    if (d.get("reaction") or d.get("forward_reaction_smiles"))
                    else None
                ),
                template_hash=(
                    str(d.get("template_hash") or d.get("retro_template_hash"))
                    if (d.get("template_hash") or d.get("retro_template_hash"))
                    else None
                ),
                template_number=template_number,
                exact_edit=(
                    bool(d["exact_edit"]) if d.get("exact_edit") is not None else None
                ),
                generation_channels=channels,
                scores=_float_dict(d.get("scores")),
                candidate_category=(
                    str(md["candidate_category"])
                    if md.get("candidate_category") is not None
                    else None
                ),
                chemistry_relation=(
                    str(md["chemistry_relation"])
                    if md.get("chemistry_relation") is not None
                    else None
                ),
                condition_relation=(
                    str(md["condition_relation"])
                    if md.get("condition_relation") is not None
                    else None
                ),
                metadata=md,
            )
        )
    return out


def _normalize_templates(
    value: Any,
    *,
    include: bool,
    ranking_mode: str,
) -> list[TemplateHit]:
    if not include:
        return []

    out: list[TemplateHit] = []
    for rank, item in enumerate(_as_list(value), 1):
        d = _as_dict(item)
        known = {
            "template_hash",
            "retro_template_hash",
            "template_number",
            "retro_template_number",
            "score",
            "template_score",
            "aggregate_score",
            "final_score",
            "member_count",
            "n_members",
            "support",
            "unique_member_count",
            "scores",
        }
        metadata = {k: v for k, v in d.items() if k not in known}

        member_count = None
        for key in ("member_count", "n_members", "support", "unique_member_count"):
            if d.get(key) is not None:
                try:
                    member_count = int(d[key])
                except Exception:
                    pass
                break

        out.append(
            TemplateHit(
                rank=rank,
                template_hash=(
                    str(d.get("template_hash") or d.get("retro_template_hash"))
                    if (d.get("template_hash") or d.get("retro_template_hash"))
                    else None
                ),
                template_number=(
                    d.get("template_number")
                    if d.get("template_number") is not None
                    else d.get("retro_template_number")
                ),
                score=_template_score(d, ranking_mode),
                member_count=member_count,
                scores=_float_dict(d.get("scores")),
                metadata=metadata,
            )
        )
    return out


def _landscape_view(value: Any, landscape_type: str) -> MethodLandscapeView:
    d = _as_dict(value)
    return MethodLandscapeView(
        available=bool(d),
        landscape_type=landscape_type,
        data=d,
    )


def _dedup_dicts(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    import json
    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        try:
            key = json.dumps(row, sort_keys=True, default=str)
        except Exception:
            key = repr(row)
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def _semantic_diagnostics(sem: dict[str, Any]) -> dict[str, Any]:
    tokens = (
        "diagnostic", "validat", "hypoth", "veto", "reject", "ambigu",
        "mapping_control", "mapping_correction", "recovery", "abst",
    )
    out: dict[str, Any] = {}
    for key, value in sem.items():
        kl = str(key).lower()
        if any(tok in kl for tok in tokens):
            out[str(key)] = value
    return out


def _extract_chemkn(raw: dict[str, Any], scope: str = "all") -> ChemKNSection:
    """Conservatively surface already-enriched ChemKN cards if present.

    No compatibility inference is performed here.  If the runtime/service did
    not enrich the method landscapes, the section remains explicitly
    ``not_requested``.
    """

    exact: list[dict[str, Any]] = []
    related: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {}

    def visit(obj: Any) -> None:
        if isinstance(obj, Mapping):
            for key, value in obj.items():
                key_l = str(key).lower()
                if key_l in {
                    "exact_cards", "exact_chemkn_cards", "chemkn_exact"
                }:
                    exact.extend(_dict_list(value))
                elif key_l in {
                    "related_cards", "candidate_cards",
                    "related_chemkn_cards", "chemkn_related",
                }:
                    related.extend(_dict_list(value))
                elif "chemkn" in key_l and isinstance(value, Mapping):
                    metadata[str(key)] = dict(value)
                else:
                    visit(value)
        elif isinstance(obj, list):
            for item in obj:
                visit(item)

    visit(raw.get("global_method_landscape"))
    visit(raw.get("conditioned_method_landscape"))

    exact = _dedup_dicts(exact)
    related = _dedup_dicts(related)
    if scope == "exact":
        related = []
    elif scope == "related":
        exact = []
    elif scope == "none":
        exact = []
        related = []

    status = str(raw.get("_chemkn_status") or (
        "available" if (exact or related or metadata) else "not_requested"
    ))
    return ChemKNSection(
        status=status,
        exact_cards=exact,
        related_cards=related,
        metadata=metadata,
    )


def adapt_runtime_result(
    *,
    raw: dict[str, Any],
    request: AnalysisRequest,
    semantic_policy: str,
    library_id: str,
    service_version: str,
    retrieval_manifest: dict[str, Any] | None = None,
    retrieval_representation_semantic_policy: str | None = None,
) -> AnalysisResult:
    qrec = _extract_query_record(raw)
    ranking_mode = str(raw.get("ranking_mode") or request.options.mode.value)
    mode = str(raw.get("mode") or request.options.mode.value)

    runtime_version = str(raw.get("runtime_version") or "6.6")

    representation_version = None
    manifest = dict(retrieval_manifest or {})
    for key in (
        "representation_version",
        "retrieval_representation_version",
        "schema_version",
        "build_id",
    ):
        if manifest.get(key):
            representation_version = str(manifest[key])
            break

    consistent = None
    if retrieval_representation_semantic_policy is not None:
        consistent = (
            str(retrieval_representation_semantic_policy)
            == str(semantic_policy)
        )

    semantic_profile = _as_dict(qrec.get("semantic_profile"))
    canonical_sem = _as_dict(qrec.get("semantics"))

    reaction_centers = _dict_list(
        qrec.get("reaction_centers") or qrec.get("centers")
    )

    def sem_rows(key: str) -> list[dict[str, Any]]:
        return _dict_list(
            canonical_sem.get(key)
            or qrec.get(key)
            or raw.get(key)
        )

    directional = sem_rows("directional_conversions")
    specific = sem_rows("specific_events")
    primary_rows = sem_rows("primary_events")
    semantic_events = sem_rows("semantic_events")
    local_events = sem_rows("local_events")
    primitives = sem_rows("primitives")
    topology_evidence = sem_rows("topology_evidence")
    redox_axis = sem_rows("redox_axis")
    formed_bond_classes = [
        str(x) for x in _as_list(
            canonical_sem.get("formed_bond_classes")
            or qrec.get("formed_bond_classes")
        ) if x is not None
    ]

    warnings = []
    for value in (
        qrec.get("warnings"),
        raw.get("warnings"),
        raw.get("_cli_warnings"),
        raw.get("_service_warnings"),
    ):
        for item in _as_list(value):
            if item is not None:
                warnings.append(str(item))

    categories = _normalize_categories(raw.get("candidate_categories"))

    templates = _normalize_templates(
        raw.get("templates"),
        include=request.options.include_templates,
        ranking_mode=ranking_mode,
    )
    candidates = _normalize_candidates(
        raw.get("members"),
        include=request.options.include_candidates,
    )

    normalized_evidence = _as_dict(raw.get("query_condition_evidence"))

    global_landscape = _landscape_view(
        raw.get("global_method_landscape"),
        "global",
    )
    conditioned_landscape = _landscape_view(
        raw.get("conditioned_method_landscape"),
        "conditioned",
    )

    timings = _float_dict(raw.get("timings"))

    raw_runtime_result = (
        raw if request.options.include_raw_runtime_result else None
    )

    diagnostics = DiagnosticsSection(
        timings=timings if request.options.include_diagnostics else {},
        missing_provider_count=(
            _optional_nonnegative_int(raw.get("missing_provider"))
            if request.options.include_diagnostics
            else None
        ),
        bulk_provider=(
            _as_dict(raw.get("bulk_provider"))
            if request.options.include_diagnostics
            else {}
        ),
        cli_or_service_warnings=warnings,
        raw_runtime_result=raw_runtime_result,
    )

    try:
        selected = int(raw["selected_candidate_keys"])
    except Exception:
        selected = None

    try:
        members = int(raw["candidate_member_count"])
    except Exception:
        members = None

    try:
        elapsed = float(raw["elapsed_seconds"])
    except Exception:
        elapsed = None

    return AnalysisResult(
        runtime=RuntimeInfo(
            runtime_version=runtime_version,
            semantic_policy=semantic_policy,
            retrieval_representation_version=representation_version,
            retrieval_representation_semantic_policy=(
                retrieval_representation_semantic_policy
            ),
            semantic_representation_consistent=consistent,
        ),
        reaction=ReactionSection(
            input_reaction=request.reaction,
            query_reaction_key=(
                str(raw["query_reaction_key"])
                if raw.get("query_reaction_key") is not None
                else None
            ),
        ),
        summary=SummarySection(
            primary_event=(
                str(raw["primary_event"])
                if raw.get("primary_event") is not None
                else None
            ),
            mode=mode,
            ranking_mode=ranking_mode,
            selected_candidate_count=selected,
            scored_member_count=members,
            returned_template_count=len(templates),
            elapsed_seconds=elapsed,
        ),
        structure=StructureSection(
            structural_keys=_as_dict(raw.get("structural_keys")),
            mapping_qc=_as_dict(raw.get("mapping_qc") or canonical_sem.get("mapping_qc")),
            mapping_correction=_as_dict(raw.get("mapping_correction")),
            reaction_centers=reaction_centers,
        ),
        semantics=SemanticsSection(
            primary_event=(
                str(raw["primary_event"])
                if raw.get("primary_event") is not None
                else None
            ),
            semantic_profile=semantic_profile,
            primary_events=primary_rows,
            semantic_events=semantic_events,
            specific_events=specific,
            local_events=local_events,
            directional_conversions=directional,
            primitives=primitives,
            topology_evidence=topology_evidence,
            redox_axis=redox_axis,
            formed_bond_classes=formed_bond_classes,
            event_topology=_as_dict(canonical_sem.get("event_topology")),
            graph_summary=_as_dict(canonical_sem.get("graph_summary")),
            component_provenance=(
                canonical_sem.get("component_provenance")
                if isinstance(canonical_sem.get("component_provenance"), (dict, list))
                else {}
            ),
            detector_outputs=_as_dict(canonical_sem.get("detector_outputs")),
            diagnostic_fields=_semantic_diagnostics(canonical_sem),
            semantic_event_count=(
                int(canonical_sem["semantic_event_count"])
                if canonical_sem.get("semantic_event_count") is not None
                else None
            ),
            multi_semantic_event=(
                bool(canonical_sem["multi_semantic_event"])
                if canonical_sem.get("multi_semantic_event") is not None
                else None
            ),
            canonical_semantics=canonical_sem,
            warnings=warnings,
        ),
        retrieval=RetrievalSection(
            mode=mode,
            ranking_mode=ranking_mode,
            candidate_policy=str(
                raw.get("candidate_policy")
                or request.options.candidate_policy.value
            ),
            candidate_budget=(
                int(raw["candidate_budget"])
                if raw.get("candidate_budget") is not None
                else request.options.candidate_budget
            ),
            generation_counts=_int_dict(raw.get("generation_counts")),
            selected_channel_counts=_int_dict(
                raw.get("selected_channel_counts")
            ),
            selected_candidate_count=selected,
            scored_member_count=members,
            categories=categories,
            templates=templates,
            candidates=candidates,
        ),
        conditions=ConditionsSection(
            supplied=request.conditions,
            has_evidence=bool(
                normalized_evidence.get(
                    "has_evidence",
                    request.conditions.has_evidence,
                )
            ),
            normalized_evidence=normalized_evidence,
            query_method_predictions=_normalize_predictions(
                raw.get("query_method_predictions")
            ),
            reagent_weight=float(raw.get("reagent_weight") or 0.0),
            global_method_landscape=global_landscape,
            conditioned_method_landscape=conditioned_landscape,
        ),
        chemkn=_extract_chemkn(raw, request.options.chemkn_scope.value),
        diagnostics=diagnostics,
        provenance=ProvenanceSection(
            service_version=service_version,
            library_id=library_id,
            retrieval_manifest=manifest,
        ),
    )
