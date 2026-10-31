"""RXNI application service for reaction intelligence.

This is the single application-facing orchestration boundary.  CLI, HTTP,
notebooks, and future UI code should call this service rather than importing
the native retriever directly.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Mapping

from rxni.integration.retriever import V66Retriever
from rxni.retrieval.libraries import resolve_library, LibraryFormatError

from .adapter import adapt_runtime_result
from .models import (
    AnalysisOptions,
    AnalysisRequest,
    AnalysisResult,
    ConditionInput,
    MappingControl,
    model_json_schema_compat,
    model_validate_compat,
)


SERVICE_VERSION = "reaction-intelligence-service-v0.2"


class ServiceConfigurationError(RuntimeError):
    pass


def _enrich_chemkn(raw: dict[str, Any], options: AnalysisOptions) -> None:
    if not options.include_chemkn or options.chemkn_scope.value == "none":
        raw["_chemkn_status"] = "not_requested"
        return

    try:
        from rxni.chemkn.knowledge_registry import KnowledgeRegistry
        from rxni.chemkn.method_links import MethodLinkStore
        from rxni.chemkn.bridge import enrich_method_landscape
    except Exception as exc:
        raw["_chemkn_status"] = "unavailable"
        raw.setdefault("_service_warnings", []).append(
            f"ChemKN modules unavailable: {type(exc).__name__}: {exc}"
        )
        return

    cards_path = os.environ.get("RXNI_KNOWLEDGE_CARDS")
    method_links = os.environ.get("RXNI_METHOD_LINKS")
    if not cards_path:
        raw["_chemkn_status"] = "unavailable"
        raw.setdefault("_service_warnings", []).append(
            "ChemKN requested without configured provider: set RXNI_KNOWLEDGE_CARDS "
            "to a normalized, locally authorized JSONL collection."
        )
        return
    if not Path(cards_path).expanduser().is_file():
        raw["_chemkn_status"] = "unavailable"
        raw.setdefault("_service_warnings", []).append(
            "Configured RXNI_KNOWLEDGE_CARDS file does not exist."
        )
        return
    if method_links and not Path(method_links).expanduser().is_file():
        raw["_chemkn_status"] = "unavailable"
        raw.setdefault("_service_warnings", []).append(
            "Configured RXNI_METHOD_LINKS file does not exist."
        )
        return
    try:
        registry = KnowledgeRegistry.from_jsonl(Path(cards_path).expanduser())
        links = MethodLinkStore.load(method_links) if method_links else None
        include_related = options.chemkn_scope.value in {"all", "related"}
        limit = max(0, int(options.chemkn_card_limit))
        if raw.get("global_method_landscape"):
            raw["global_method_landscape"] = enrich_method_landscape(
                raw["global_method_landscape"],
                registry,
                retrieval_result=raw,
                link_store=links,
                include_candidates=include_related,
                max_candidates=limit,
            )
            raw["method_profile"] = raw["global_method_landscape"]
        if raw.get("conditioned_method_landscape"):
            raw["conditioned_method_landscape"] = enrich_method_landscape(
                raw["conditioned_method_landscape"],
                registry,
                retrieval_result=raw,
                link_store=links,
                include_candidates=include_related,
                max_candidates=limit,
            )
        raw["_chemkn_status"] = "available"
    except Exception as exc:
        raw["_chemkn_status"] = "unavailable"
        raw.setdefault("_service_warnings", []).append(
            f"ChemKN enrichment failed: {type(exc).__name__}: {exc}"
        )


class ReactionIntelligenceService:
    """RXNI application service preserving the Retrieval-V2 contract.

    ``mapping_control`` is process-level because the current semantic policy
    seam is controlled by ``RUNTIME_V67_MAPPING_CONTROL``.  A web process
    should therefore use one semantic policy for its lifetime rather than
    switching policies per request.
    """

    def __init__(
        self,
        *,
        library: str | Path,
        mapping_control: str | MappingControl = MappingControl.CONVERSION_VALIDATED,
        library_id: str | None = None,
        retrieval_representation_semantic_policy: str | None = None,
        retriever: Any | None = None,
    ) -> None:
        self.library = Path(library)
        self.mapping_control = (
            mapping_control.value
            if isinstance(mapping_control, MappingControl)
            else str(mapping_control)
        )

        if self.mapping_control not in {
            MappingControl.CONVERSION_VALIDATED.value,
            MappingControl.STRICT.value,
        }:
            raise ServiceConfigurationError(
                f"Unsupported mapping_control={self.mapping_control!r}"
            )

        # Semantic mapping policy remains process-global.
        os.environ["RUNTIME_V67_MAPPING_CONTROL"] = self.mapping_control

        self.library_id = library_id or self.library.name
        self.retrieval_representation_semantic_policy = (
            retrieval_representation_semantic_policy
        )
        self.retrieval_manifest = self._load_retrieval_manifest()

        if retriever is not None:
            self._retriever = retriever
        else:
            if not self.library.exists():
                raise ServiceConfigurationError(
                    f"Template library does not exist: {self.library}"
                )
            self._retriever = V66Retriever(str(self.library))

    def _load_retrieval_manifest(self) -> dict[str, Any]:
        try:
            path = resolve_library(self.library).retrieval_manifest
        except (LibraryFormatError, FileNotFoundError):
            # Preserve the original optional-manifest behavior when a custom
            # retriever is injected or the library is unavailable.
            return {}
        if not path.exists():
            return {}

        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return {}

        if not isinstance(payload, dict):
            return {}

        # Public provenance intentionally excludes server-local paths and
        # arbitrary wide build payloads.
        preferred = (
            "build_id",
            "schema_version",
            "representation_version",
            "retrieval_representation_version",
            "semantic_authority",
            "semantic_policy",
            "scoring_contract",
            "created_at",
            "created_utc",
        )
        out = {k: payload[k] for k in preferred if k in payload}

        # Preserve a compact runtime section when present.
        if isinstance(payload.get("runtime"), dict):
            out["runtime"] = dict(payload["runtime"])

        return out

    def analyze(
        self,
        reaction: str,
        *,
        conditions: ConditionInput | Mapping[str, Any] | None = None,
        options: AnalysisOptions | Mapping[str, Any] | None = None,
    ) -> AnalysisResult:
        condition_model = self._coerce_conditions(conditions)
        option_model = self._coerce_options(options)

        request = AnalysisRequest(
            reaction=reaction,
            conditions=condition_model,
            options=option_model,
        )
        return self.analyze_request(request)

    def analyze_request(
        self,
        request: AnalysisRequest | Mapping[str, Any],
    ) -> AnalysisResult:
        if not isinstance(request, AnalysisRequest):
            request = model_validate_compat(AnalysisRequest, request)  # type: ignore[assignment]

        conditions = request.conditions
        options = request.options

        raw = self._retriever.retrieve(
            request.reaction,
            mode=options.mode.value,
            top_k=options.top_k,
            exclude_template=options.exclude_template,
            max_members=options.candidate_budget,
            per_channel_limit=options.per_channel_limit,
            candidate_policy=options.candidate_policy.value,
            member_limit=options.member_limit,
            reagents=list(conditions.reagents),
            solvents=list(conditions.solvents),
            catalysts=list(conditions.catalysts),
            bases=list(conditions.bases),
            additives=list(conditions.additives),
            acids=list(conditions.acids),
            temperature_c=conditions.temperature_c,
            reaction_time_h=conditions.reaction_time_h,
            atmosphere=conditions.atmosphere,
            pressure=conditions.pressure,
            reagent_weight=options.reagent_weight,
            profile_methods=(options.profile_methods or options.include_chemkn),
        )

        if not isinstance(raw, dict):
            if hasattr(raw, "as_dict"):
                raw = raw.as_dict()
            else:
                raise TypeError(
                    "V66Retriever.retrieve returned an unsupported object "
                    f"of type {type(raw)!r}"
                )

        _enrich_chemkn(raw, options)

        return adapt_runtime_result(
            raw=raw,
            request=request,
            semantic_policy=self.mapping_control,
            library_id=self.library_id,
            service_version=SERVICE_VERSION,
            retrieval_manifest=self.retrieval_manifest,
            retrieval_representation_semantic_policy=(
                self.retrieval_representation_semantic_policy
            ),
        )

    @staticmethod
    def _coerce_conditions(
        value: ConditionInput | Mapping[str, Any] | None,
    ) -> ConditionInput:
        if value is None:
            return ConditionInput()
        if isinstance(value, ConditionInput):
            return value
        return model_validate_compat(ConditionInput, value)  # type: ignore[return-value]

    @staticmethod
    def _coerce_options(
        value: AnalysisOptions | Mapping[str, Any] | None,
    ) -> AnalysisOptions:
        if value is None:
            return AnalysisOptions()
        if isinstance(value, AnalysisOptions):
            return value
        return model_validate_compat(AnalysisOptions, value)  # type: ignore[return-value]

    @staticmethod
    def request_json_schema() -> dict[str, Any]:
        return model_json_schema_compat(AnalysisRequest)

    @staticmethod
    def result_json_schema() -> dict[str, Any]:
        return model_json_schema_compat(AnalysisResult)
