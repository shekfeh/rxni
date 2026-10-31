from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Iterable


@dataclass
class CompatibilitySignature:
    """Generic optional compatibility overlay, independent of knowledge source."""
    required_events: list[str] = field(default_factory=list)
    required_event_mode: str = "none"          # none | soft | hard
    event_alternatives: list[dict[str, Any]] = field(default_factory=list)
    event_alternatives_mode: str = "none"      # none | soft | hard
    direction: str | None = None
    allowed_donor_elements: list[str] = field(default_factory=list)
    allowed_electrophile_classes: list[str] = field(default_factory=list)
    required_motifs: list[str] = field(default_factory=list)
    preferred_motifs: list[str] = field(default_factory=list)
    forbidden_motifs: list[str] = field(default_factory=list)
    specificity: str | None = None
    knowledge_role: str | None = None
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class KnowledgeCard:
    """Common runtime-facing card model.

    `source_payload` always retains the original record so the adapter never
    destroys source-specific information.
    """
    source: str
    source_id: str
    canonical_name: str
    source_kind: str                         # e.g. reaction_card | literature_record | curated_rule
    transformation: str | None = None
    reaction_family: str | None = None
    aliases: list[str] = field(default_factory=list)
    retrieval_tags: list[str] = field(default_factory=list)

    compatibility: CompatibilitySignature = field(default_factory=CompatibilitySignature)

    reagents: list[str] = field(default_factory=list)
    catalysts: list[str] = field(default_factory=list)
    solvents: list[str] = field(default_factory=list)

    mechanism: list[str] = field(default_factory=list)
    scope: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    competing_reactions: list[str] = field(default_factory=list)

    recognition: dict[str, Any] = field(default_factory=dict)
    ontology: dict[str, Any] = field(default_factory=dict)
    condition_regimes: list[dict[str, Any]] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)
    curation: dict[str, Any] = field(default_factory=dict)
    extensions: dict[str, Any] = field(default_factory=dict)

    source_payload: dict[str, Any] = field(default_factory=dict)

    @property
    def stable_ref(self) -> str:
        return f"{self.source}:{self.source_id}"

    def as_dict(self, *, include_raw: bool = False) -> dict[str, Any]:
        d = asdict(self)
        d["stable_ref"] = self.stable_ref
        if not include_raw:
            d.pop("source_payload", None)
        return d


@dataclass
class RuntimeCompatibilityQuery:
    """Minimal runtime evidence needed to evaluate formal card compatibility."""
    events: set[str] = field(default_factory=set)
    direction: str | None = None
    donor_elements: set[str] = field(default_factory=set)
    electrophile_classes: set[str] = field(default_factory=set)
    motifs: set[str] = field(default_factory=set)

    @classmethod
    def from_values(
        cls,
        *,
        events: Iterable[str] = (),
        direction: str | None = None,
        donor_elements: Iterable[str] = (),
        electrophile_classes: Iterable[str] = (),
        motifs: Iterable[str] = (),
    ) -> "RuntimeCompatibilityQuery":
        return cls(
            events={str(x) for x in events if str(x)},
            direction=None if direction is None else str(direction),
            donor_elements={str(x) for x in donor_elements if str(x)},
            electrophile_classes={str(x) for x in electrophile_classes if str(x)},
            motifs={str(x) for x in motifs if str(x)},
        )


@dataclass
class CompatibilityDecision:
    compatible: bool
    score: float
    hard_failures: list[str] = field(default_factory=list)
    matched: list[str] = field(default_factory=list)
    missing_soft: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)
