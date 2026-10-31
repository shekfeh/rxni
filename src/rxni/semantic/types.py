"""Canonical semantic datatypes.

This module contains data only.  It performs no detection, validation,
resolution, redox inference, or retrieval scoring.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any, Mapping


class SemanticRole(str, Enum):
    PRIMARY = "primary_event"
    SPECIFIC = "specific_event"
    DIRECTIONAL = "directional_conversion"
    AXIS = "axis_event"
    LOCAL = "local_event"
    PRIMITIVE = "primitive"
    TOPOLOGY = "topology_evidence"
    UNKNOWN = "unknown"


class ValidationStatus(str, Enum):
    ACCEPTED = "accepted"
    REJECTED = "rejected"
    NOT_RUN = "not_run"
    AMBIGUOUS = "ambiguous"


@dataclass(frozen=True)
class SemanticRecord:
    """One semantic claim or structural-semantic record.

    Known schema fields are explicit.  Historical or detector-specific fields
    are retained losslessly in ``extra`` during migration.
    """

    event: str
    role: str = "unknown"
    subtype: str | None = None
    center_map: int | None = None
    support_maps: tuple[int, ...] = ()
    source_functional_group: str | None = None
    target_functional_group: str | None = None
    conversion_family: str | None = None
    redox_direction: str | None = None
    primary_eligible: bool | None = None
    mapping_scope: str | None = None
    multiplicity: int = 1
    extra: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, record: Mapping[str, Any]) -> "SemanticRecord":
        known = {
            "event", "role", "subtype", "center_map", "support_maps",
            "source_functional_group", "target_functional_group",
            "conversion_family", "redox_direction", "primary_eligible",
            "mapping_scope", "multiplicity",
        }
        extra = {k: v for k, v in record.items() if k not in known}
        support = tuple(sorted({int(x) for x in (record.get("support_maps") or [])}))
        center = record.get("center_map")
        return cls(
            event=str(record.get("event") or ""),
            role=str(record.get("role") or "unknown"),
            subtype=(None if record.get("subtype") is None else str(record.get("subtype"))),
            center_map=(None if center is None else int(center)),
            support_maps=support,
            source_functional_group=record.get("source_functional_group"),
            target_functional_group=record.get("target_functional_group"),
            conversion_family=record.get("conversion_family"),
            redox_direction=record.get("redox_direction"),
            primary_eligible=record.get("primary_eligible"),
            mapping_scope=record.get("mapping_scope"),
            multiplicity=int(record.get("multiplicity", 1) or 1),
            extra=extra,
        )

    def as_dict(self) -> dict[str, Any]:
        out = {
            "event": self.event,
            "role": self.role,
            "subtype": self.subtype,
            "center_map": self.center_map,
            "support_maps": list(self.support_maps),
            "source_functional_group": self.source_functional_group,
            "target_functional_group": self.target_functional_group,
            "conversion_family": self.conversion_family,
            "redox_direction": self.redox_direction,
            "primary_eligible": self.primary_eligible,
            "mapping_scope": self.mapping_scope,
            "multiplicity": self.multiplicity,
        }
        out.update(dict(self.extra))
        return {k: v for k, v in out.items() if v is not None}


@dataclass(frozen=True)
class ValidationCheck:
    name: str
    passed: bool | None
    detail: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class SemanticValidation:
    label: str
    status: ValidationStatus
    validator_class: str | None = None
    checks: tuple[ValidationCheck, ...] = ()
    veto_reasons: tuple[str, ...] = ()
    contract_version: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "status": self.status.value,
            "validator_class": self.validator_class,
            "checks": [x.as_dict() for x in self.checks],
            "veto_reasons": list(self.veto_reasons),
            "contract_version": self.contract_version,
        }


@dataclass
class SemanticResult:
    """Canonical public semantic result for Runtime v6.6+."""

    reaction: str
    engine_version: str
    rules_version: str
    mapping_qc: dict[str, Any] = field(default_factory=dict)

    primary_events: list[SemanticRecord] = field(default_factory=list)
    specific_events: list[SemanticRecord] = field(default_factory=list)
    directional_conversions: list[SemanticRecord] = field(default_factory=list)
    redox_axis: list[SemanticRecord] = field(default_factory=list)

    semantic_events: list[dict[str, Any]] = field(default_factory=list)
    local_events: list[SemanticRecord] = field(default_factory=list)
    primitives: list[SemanticRecord] = field(default_factory=list)
    topology_evidence: list[SemanticRecord] = field(default_factory=list)
    event_topology: dict[str, Any] = field(default_factory=dict)
    component_provenance: list[dict[str, Any]] = field(default_factory=list)
    formed_bond_classes: list[str] = field(default_factory=list)
    graph_summary: dict[str, Any] = field(default_factory=dict)

    # Raw detector material is kept for diagnostics but is not semantic truth.
    raw_hypotheses: dict[str, Any] = field(default_factory=dict)
    validations: list[SemanticValidation] = field(default_factory=list)
    rejected_hypotheses: list[dict[str, Any]] = field(default_factory=list)

    semantic_profile: dict[str, float] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        def rows(xs):
            return [x.as_dict() if hasattr(x, "as_dict") else x for x in xs]

        return {
            "reaction": self.reaction,
            "engine_version": self.engine_version,
            "rules_version": self.rules_version,
            "mapping_qc": self.mapping_qc,
            "primary_events": rows(self.primary_events),
            "specific_events": rows(self.specific_events),
            "directional_conversions": rows(self.directional_conversions),
            "redox_axis": rows(self.redox_axis),
            "semantic_events": self.semantic_events,
            "local_events": rows(self.local_events),
            "primitives": rows(self.primitives),
            "topology_evidence": rows(self.topology_evidence),
            "event_topology": self.event_topology,
            "component_provenance": self.component_provenance,
            "formed_bond_classes": self.formed_bond_classes,
            "graph_summary": self.graph_summary,
            "raw_hypotheses": self.raw_hypotheses,
            "validations": rows(self.validations),
            "rejected_hypotheses": self.rejected_hypotheses,
            "semantic_profile": self.semantic_profile,
            "diagnostics": self.diagnostics,
        }
