"""Native Runtime-v6.6 reconstruction of the frozen v6.5.6 semantic reference.

This module contains no Runtime-v6.5 imports.

It adapts the mechanically migrated v6.5.6 non-FGI rule engine into the
Runtime-v6.6 SemanticResult contract.  It intentionally does NOT yet replace
``reference_v656`` as the audit/reference backend.

M6C2D1 must first prove:
  1. reference scientific-field parity;
  2. current Stage-E canonical parity when this native reference is supplied
     to the existing v6.6 enforcement layer.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from .nonfgi.resolver import classify_system_i
from .types import SemanticRecord, SemanticResult
from .version import SEMANTIC_ENGINE_VERSION, SEMANTIC_RULESET_VERSION


def _records(rows: list[dict[str, Any]] | None) -> list[SemanticRecord]:
    return [SemanticRecord.from_mapping(x) for x in (rows or [])]


def classify_native_reference(reaction: str) -> SemanticResult:
    """Return the migrated v6.5.6 rule-engine result as a v6.6 SemanticResult."""
    sem = deepcopy(classify_system_i(reaction))
    versions = sem.get("versions") or {}

    return SemanticResult(
        reaction=reaction,
        engine_version=SEMANTIC_ENGINE_VERSION,
        rules_version=SEMANTIC_RULESET_VERSION,
        mapping_qc=deepcopy(sem.get("mapping_qc") or {}),
        primary_events=_records(sem.get("primary_events")),
        specific_events=_records(sem.get("specific_events")),
        directional_conversions=_records(sem.get("directional_conversions")),
        redox_axis=_records(sem.get("redox_axis")),
        semantic_events=deepcopy(sem.get("semantic_events") or []),
        local_events=_records(sem.get("local_events")),
        primitives=_records(sem.get("primitives")),
        topology_evidence=_records(sem.get("topology_evidence")),
        event_topology=deepcopy(sem.get("event_topology") or {}),
        component_provenance=deepcopy(
            sem.get("component_provenance") or []
        ),
        formed_bond_classes=list(
            sem.get("formed_bond_classes") or []
        ),
        graph_summary=deepcopy(sem.get("graph_summary") or {}),
        raw_hypotheses=deepcopy(sem.get("detector_outputs") or {}),
        diagnostics={
            "backend": "native__reference",
            "reference_semantic_system": sem.get("system"),
            "reference_rule_version": versions.get("rules"),
            "reference_graph_version": versions.get("graph"),
            "mapping_recovery": {
                "source": "native_v6_6_fact_boundary",
                "legacy_mapping_correction_invoked": False,
            },
        },
    )

# M6C2D2 compatibility name used by the preserved Stage-E pipeline.
classify_reference = classify_native_reference
