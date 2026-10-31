"""Public typed schema for RXNI.

The models in this module are intentionally decoupled from internal retrieval,
semantic, ChemKN, and storage dataclasses.  They form the stable boundary used
by the service layer, CLI, HTTP API, and future UI.
"""

from __future__ import annotations

from enum import Enum
from typing import Any, Mapping

from pydantic import BaseModel, Field

try:  # Pydantic v2
    from pydantic import ConfigDict
    _PYDANTIC_V2 = hasattr(BaseModel, "model_dump")
except Exception:  # pragma: no cover
    ConfigDict = None
    _PYDANTIC_V2 = False


SCHEMA_VERSION = "reaction-intelligence-api-v1"


class PublicModel(BaseModel):
    """Base class with a closed public schema."""

    if _PYDANTIC_V2:
        model_config = ConfigDict(extra="forbid")  # type: ignore[misc]
    else:  # pragma: no cover - Pydantic v1 compatibility
        class Config:
            extra = "forbid"


class MappingControl(str, Enum):
    CONVERSION_VALIDATED = "conversion_validated"
    STRICT = "strict"


class RetrievalMode(str, Enum):
    LOCAL = "local"
    MULTI = "multi"
    SEMANTIC = "semantic"
    SEM_MULTI = "sem_multi"
    EXACT_EDIT = "exact_edit"


class CandidatePolicy(str, Enum):
    SPECIFICITY = "specificity"
    LEGACY = "legacy"


class ChemKNScope(str, Enum):
    ALL = "all"
    EXACT = "exact"
    RELATED = "related"
    NONE = "none"


class ConditionInput(PublicModel):
    reagents: list[str] = Field(default_factory=list)
    solvents: list[str] = Field(default_factory=list)
    catalysts: list[str] = Field(default_factory=list)
    bases: list[str] = Field(default_factory=list)
    additives: list[str] = Field(default_factory=list)
    acids: list[str] = Field(default_factory=list)
    temperature_c: float | None = None
    reaction_time_h: float | None = None
    atmosphere: str | None = None
    pressure: str | None = None

    @property
    def has_evidence(self) -> bool:
        return bool(
            self.reagents
            or self.solvents
            or self.catalysts
            or self.bases
            or self.additives
            or self.acids
            or self.temperature_c is not None
            or self.reaction_time_h is not None
            or self.atmosphere
            or self.pressure
        )


class AnalysisOptions(PublicModel):
    mode: RetrievalMode = RetrievalMode.SEM_MULTI
    top_k: int = Field(default=20, ge=1)
    candidate_budget: int = Field(default=3000, ge=1)
    per_channel_limit: int = Field(default=5000, ge=1)
    candidate_policy: CandidatePolicy = CandidatePolicy.SPECIFICITY
    member_limit: int | None = Field(default=None, ge=0)
    exclude_template: str | None = None
    reagent_weight: float = Field(default=0.20, ge=0.0, le=1.0)
    profile_methods: bool = False
    include_candidates: bool = True
    include_templates: bool = True
    include_diagnostics: bool = False
    include_raw_runtime_result: bool = False
    include_chemkn: bool = False
    chemkn_scope: ChemKNScope = ChemKNScope.ALL
    chemkn_card_limit: int = Field(default=8, ge=0)


class AnalysisRequest(PublicModel):
    reaction: str = Field(min_length=1)
    conditions: ConditionInput = Field(default_factory=ConditionInput)
    options: AnalysisOptions = Field(default_factory=AnalysisOptions)


class RuntimeInfo(PublicModel):
    runtime_version: str
    schema_version: str = SCHEMA_VERSION
    semantic_policy: str
    retrieval_representation_version: str | None = None
    retrieval_representation_semantic_policy: str | None = None
    semantic_representation_consistent: bool | None = None


class ReactionSection(PublicModel):
    input_reaction: str
    query_reaction_key: str | None = None


class SummarySection(PublicModel):
    primary_event: str | None = None
    mode: str
    ranking_mode: str
    selected_candidate_count: int | None = None
    scored_member_count: int | None = None
    returned_template_count: int = 0
    elapsed_seconds: float | None = None


class StructureSection(PublicModel):
    structural_keys: dict[str, Any] = Field(default_factory=dict)
    mapping_qc: dict[str, Any] = Field(default_factory=dict)
    mapping_correction: dict[str, Any] = Field(default_factory=dict)
    reaction_centers: list[dict[str, Any]] = Field(default_factory=list)


class SemanticsSection(PublicModel):
    primary_event: str | None = None
    semantic_profile: dict[str, Any] = Field(default_factory=dict)
    primary_events: list[dict[str, Any]] = Field(default_factory=list)
    semantic_events: list[dict[str, Any]] = Field(default_factory=list)
    specific_events: list[dict[str, Any]] = Field(default_factory=list)
    local_events: list[dict[str, Any]] = Field(default_factory=list)
    directional_conversions: list[dict[str, Any]] = Field(default_factory=list)
    primitives: list[dict[str, Any]] = Field(default_factory=list)
    topology_evidence: list[dict[str, Any]] = Field(default_factory=list)
    redox_axis: list[dict[str, Any]] = Field(default_factory=list)
    formed_bond_classes: list[str] = Field(default_factory=list)
    event_topology: dict[str, Any] = Field(default_factory=dict)
    graph_summary: dict[str, Any] = Field(default_factory=dict)
    component_provenance: dict[str, Any] | list[Any] = Field(default_factory=dict)
    detector_outputs: dict[str, Any] = Field(default_factory=dict)
    diagnostic_fields: dict[str, Any] = Field(default_factory=dict)
    semantic_event_count: int | None = None
    multi_semantic_event: bool | None = None
    canonical_semantics: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)


class CandidateCategory(PublicModel):
    category: str
    count: int
    fraction: float | None = None


class CandidateHit(PublicModel):
    rank: int
    reaction_key: str | None = None
    reaction: str | None = None
    template_hash: str | None = None
    template_number: int | str | None = None
    exact_edit: bool | None = None
    generation_channels: list[str] = Field(default_factory=list)
    scores: dict[str, float] = Field(default_factory=dict)
    candidate_category: str | None = None
    chemistry_relation: str | None = None
    condition_relation: str | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)


class TemplateHit(PublicModel):
    rank: int
    template_hash: str | None = None
    template_number: int | str | None = None
    score: float | None = None
    member_count: int | None = None
    scores: dict[str, float] = Field(default_factory=dict)
    metadata: dict[str, Any] = Field(default_factory=dict)


class RetrievalSection(PublicModel):
    mode: str
    ranking_mode: str
    candidate_policy: str
    candidate_budget: int | None = None
    generation_counts: dict[str, int] = Field(default_factory=dict)
    selected_channel_counts: dict[str, int] = Field(default_factory=dict)
    selected_candidate_count: int | None = None
    scored_member_count: int | None = None
    categories: list[CandidateCategory] = Field(default_factory=list)
    templates: list[TemplateHit] = Field(default_factory=list)
    candidates: list[CandidateHit] = Field(default_factory=list)


class MethodLandscapeView(PublicModel):
    available: bool
    landscape_type: str
    data: dict[str, Any] = Field(default_factory=dict)


class ConditionsSection(PublicModel):
    supplied: ConditionInput
    has_evidence: bool
    normalized_evidence: dict[str, Any] = Field(default_factory=dict)
    query_method_predictions: list[dict[str, Any]] = Field(default_factory=list)
    reagent_weight: float = 0.0
    global_method_landscape: MethodLandscapeView
    conditioned_method_landscape: MethodLandscapeView


class ChemKNSection(PublicModel):
    """Stable ChemKN envelope.

    RXNI retrieval does not itself guarantee that textbook/reaction-card
    enrichment has been requested.  The public schema therefore represents the
    state explicitly instead of overloading method-landscape fields.
    """

    status: str = "not_requested"
    exact_cards: list[dict[str, Any]] = Field(default_factory=list)
    related_cards: list[dict[str, Any]] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)


class DiagnosticsSection(PublicModel):
    timings: dict[str, float] = Field(default_factory=dict)
    missing_provider_count: int | None = Field(default=None, ge=0)
    bulk_provider: dict[str, Any] = Field(default_factory=dict)
    cli_or_service_warnings: list[str] = Field(default_factory=list)
    raw_runtime_result: dict[str, Any] | None = None


class ProvenanceSection(PublicModel):
    service_version: str
    library_id: str
    retrieval_manifest: dict[str, Any] = Field(default_factory=dict)


class AnalysisResult(PublicModel):
    schema_version: str = SCHEMA_VERSION
    runtime: RuntimeInfo
    reaction: ReactionSection
    summary: SummarySection
    structure: StructureSection
    semantics: SemanticsSection
    retrieval: RetrievalSection
    conditions: ConditionsSection
    chemkn: ChemKNSection
    diagnostics: DiagnosticsSection
    provenance: ProvenanceSection


def model_dump_compat(model: BaseModel, **kwargs: Any) -> dict[str, Any]:
    if hasattr(model, "model_dump"):
        return model.model_dump(**kwargs)  # type: ignore[attr-defined]
    return model.dict(**kwargs)  # pragma: no cover


def model_validate_compat(model_cls: type[BaseModel], value: Any) -> BaseModel:
    if hasattr(model_cls, "model_validate"):
        return model_cls.model_validate(value)  # type: ignore[attr-defined]
    return model_cls.parse_obj(value)  # pragma: no cover


def model_json_schema_compat(model_cls: type[BaseModel]) -> dict[str, Any]:
    if hasattr(model_cls, "model_json_schema"):
        return model_cls.model_json_schema()  # type: ignore[attr-defined]
    return model_cls.schema()  # pragma: no cover
