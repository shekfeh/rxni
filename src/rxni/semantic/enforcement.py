"""Canonical validated-FGI enforcement for Runtime v6.6.

This module is the boundary between raw native directional hypotheses and
canonical semantic truth.

Policy:
- ACCEPTED hypotheses enter ``directional_conversions``.
- REJECTED hypotheses are diagnostic only.
- AMBIGUOUS hypotheses abstain and are diagnostic only.
- NOT_RUN is a hard canonicalization error because Stage E is enabled only
  after the 82/82 proof ledger is complete.
- Conversion-derived reference primary/specific rows survive only when their
  underlying directional conversion has an accepted native instance.
- No generic primary is promoted as fallback after an exact conversion is
  rejected.  Rejection therefore produces conservative abstention rather than
  semantic invention.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
from typing import Any, Iterable

from .detectors import directional_hypotheses, run_directional_hypothesis_detectors
from .structure import build_mapped_graph_facts
from .types import SemanticRecord, SemanticResult, ValidationStatus
from .validation import ValidationEngine


class CanonicalValidationError(RuntimeError):
    """Raised when canonical enforcement encounters an unvalidated hypothesis."""


def _support(rec: SemanticRecord) -> tuple[int, ...]:
    return tuple(sorted(int(x) for x in rec.support_maps))


def directional_instance_key(rec: SemanticRecord) -> tuple[str, int | None, tuple[int, ...]]:
    return (str(rec.event), rec.center_map, _support(rec))


def _conversion_child_name(rec: SemanticRecord) -> str | None:
    if not bool(rec.extra.get("conversion_specific")):
        return None
    child = rec.extra.get("directional_conversion")
    if child:
        return str(child)
    if rec.subtype:
        return str(rec.subtype)
    return str(rec.event) if rec.event else None


def conversion_instance_key(
    rec: SemanticRecord,
) -> tuple[str, int | None, tuple[int, ...]] | None:
    child = _conversion_child_name(rec)
    if not child:
        return None
    return (child, rec.center_map, _support(rec))


def validate_directional_hypotheses(
    facts: dict[str, Any],
    *,
    engine: ValidationEngine | None = None,
) -> tuple[
    list[SemanticRecord],
    list[Any],
    list[dict[str, Any]],
    dict[str, Any],
]:
    """Validate every native directional hypothesis instance.

    Returns:
      accepted records,
      SemanticValidation rows,
      non-accepted diagnostic rows,
      detector output dictionary.
    """
    engine = engine or ValidationEngine()
    detector_outputs = run_directional_hypothesis_detectors(facts)
    raw = []
    for rows in detector_outputs.values():
        raw.extend(
            SemanticRecord.from_mapping(r)
            for r in rows
            if r.get("role") == "directional_conversion"
        )

    accepted: list[SemanticRecord] = []
    validations = []
    rejected: list[dict[str, Any]] = []

    for rec in raw:
        val = engine.validate_directional(facts, rec)
        validations.append(val)

        if val.status == ValidationStatus.ACCEPTED:
            accepted.append(rec)
            continue

        if val.status == ValidationStatus.NOT_RUN:
            raise CanonicalValidationError(
                f"Canonical enforcement encountered NOT_RUN validator: "
                f"event={rec.event!r} center_map={rec.center_map!r}"
            )

        rejected.append(
            {
                **rec.as_dict(),
                "validation_status": val.status.value,
                "validator_class": val.validator_class,
                "veto_reasons": list(val.veto_reasons),
                "checks": [x.as_dict() for x in val.checks],
                "contract_version": val.contract_version,
            }
        )

    return accepted, validations, rejected, detector_outputs


def _accepted_indexes(
    accepted: Iterable[SemanticRecord],
) -> tuple[
    set[tuple[str, int | None, tuple[int, ...]]],
    set[tuple[str, int | None]],
]:
    exact = set()
    event_center = set()
    for rec in accepted:
        exact.add(directional_instance_key(rec))
        event_center.add((str(rec.event), rec.center_map))
    return exact, event_center


def _conversion_row_is_accepted(
    rec: SemanticRecord,
    *,
    accepted_exact: set[tuple[str, int | None, tuple[int, ...]]],
    accepted_event_center: set[tuple[str, int | None]],
) -> bool:
    key = conversion_instance_key(rec)
    if key is None:
        return True

    # Exact instance identity is authoritative.  Event+center fallback is only
    # used when a historical specific/primary row has reduced support metadata.
    if key in accepted_exact:
        return True
    child, center, support = key
    if not support and (child, center) in accepted_event_center:
        return True
    return False


def filter_conversion_derived_rows(
    rows: Iterable[SemanticRecord],
    accepted: Iterable[SemanticRecord],
) -> list[SemanticRecord]:
    """Remove conversion-derived rows whose exact native child was not accepted."""
    exact, event_center = _accepted_indexes(accepted)
    return [
        rec
        for rec in rows
        if _conversion_row_is_accepted(
            rec,
            accepted_exact=exact,
            accepted_event_center=event_center,
        )
    ]


def _semantic_event_matches_primary(
    event_row: dict[str, Any],
    primary: SemanticRecord,
) -> bool:
    if str(event_row.get("event") or "") != str(primary.event or ""):
        return False
    if event_row.get("center_map") != primary.center_map:
        return False

    ev_sub = str(event_row.get("subtype") or "")
    p_sub = str(primary.subtype or "")
    if ev_sub and p_sub and ev_sub != p_sub:
        return False

    ev_support = tuple(sorted(int(x) for x in (event_row.get("support_maps") or [])))
    p_support = _support(primary)
    if ev_support and p_support and ev_support != p_support:
        return False
    return True


def filter_semantic_events_to_primary(
    semantic_events: Iterable[dict[str, Any]],
    primary: Iterable[SemanticRecord],
) -> list[dict[str, Any]]:
    prim = list(primary)
    return [
        deepcopy(row)
        for row in semantic_events
        if any(_semantic_event_matches_primary(row, p) for p in prim)
    ]


def enforce_validated_directional_semantics(
    reference: SemanticResult,
    reaction: str,
) -> SemanticResult:
    """Canonical Stage-E semantics with strict or conversion-validated mapping policy."""
    import os

    # Runtime-v6.7 owns the active semantic policy. Fall back to the old
    # variable for compatibility with scripts that have not yet migrated.
    mapping_control = (
        os.environ.get(
            "RUNTIME_V67_MAPPING_CONTROL",
            os.environ.get(
                "RUNTIME_V66_MAPPING_CONTROL",
                "conversion_validated",
            ),
        )
        .strip()
        .lower()
    )
    allowed_controls = {"strict", "conversion_validated"}
    if mapping_control not in allowed_controls:
        raise CanonicalValidationError(
            f"Unsupported RUNTIME_V67_MAPPING_CONTROL={mapping_control!r}; "
            f"expected one of {sorted(allowed_controls)}"
        )

    facts = build_mapped_graph_facts(reaction)
    strict_accepted, validations, rejected, detector_outputs = (
        validate_directional_hypotheses(facts)
    )

    accepted = list(strict_accepted)
    effective_rejected = list(rejected)
    recoveries = []

    def _check_map(diagnostic):
        return {
            str(row.get("name")): row.get("passed")
            for row in (diagnostic.get("checks") or [])
        }

    def _required_true(checks, names):
        return all(checks.get(name) is True for name in names)

    def _conversion_validated_mapping_recovery(diagnostic):
        validator = str(diagnostic.get("validator_class") or "")
        vetoes = {
            str(x)
            for x in (diagnostic.get("veto_reasons") or [])
        }
        checks = _check_map(diagnostic)

        if not vetoes:
            return False, ""

        if validator == "VC2_ACYL_CENTER":
            if vetoes != {"GRAPH_EDIT_COMPATIBLE"}:
                return False, ""
            required = (
                "mapped_center_present",
                "mapped_center_is_carbon",
                "SOURCE_PRESENT",
                "SOURCE_DESTROYED",
                "TARGET_ABSENT_BEFORE",
                "TARGET_CREATED",
            )
            if _required_true(checks, required):
                return True, "vc2_graph_only_equivalent_atom_mapping"
            return False, ""

        if validator == "VC4_ATTACHMENT_REPLACEMENT":
            allowed = {
                "SOURCE_PRESENT",
                "GRAPH_EDIT_COMPATIBLE",
            }
            if not vetoes.issubset(allowed):
                return False, ""
            required = (
                "mapped_center_present",
                "SOURCE_DESTROYED",
                "TARGET_ABSENT_BEFORE",
                "TARGET_CREATED",
                "event_specific_guard",
                "mapped_carbon_scaffold_preserved",
            )
            if _required_true(checks, required):
                return True, "vc4_boundary_mapping_with_valid_attachment_replacement"
            return False, ""

        if validator == "VC5_HYDROGEN_REPLACEMENT":
            allowed = {
                "SOURCE_PRESENT",
                "TARGET_ABSENT_BEFORE",
                "GRAPH_EDIT_COMPATIBLE",
            }
            if not vetoes.issubset(allowed):
                return False, ""
            required = (
                "mapped_center_present",
                "SOURCE_DESTROYED",
                "TARGET_CREATED",
                "event_specific_guard",
                "mapped_scaffold_preserved",
            )
            if _required_true(checks, required):
                return True, "vc5_boundary_mapping_with_valid_h_replacement"
            return False, ""

        if validator == "VC6_SPECIAL_MULTICENTER":
            event = str(diagnostic.get("event") or "")

            if event == "benzylic_halide_to_aldehyde":
                required_vetoes = {
                    "SOURCE_PRESENT",
                    "SOURCE_DESTROYED",
                    "GRAPH_EDIT_COMPATIBLE",
                }
                if vetoes != required_vetoes:
                    return False, ""
                required = (
                    "mapped_center_present",
                    "TARGET_ABSENT_BEFORE",
                    "TARGET_CREATED",
                    "benzylic_aromatic_attachment_preserved",
                    "mapped_carbon_scaffold_preserved",
                )
                if _required_true(checks, required):
                    return True, "vc6_benzylic_halide_boundary_mapping_to_aldehyde"
                return False, ""

            allowed = {
                "SOURCE_PRESENT",
                "GRAPH_EDIT_COMPATIBLE",
            }
            if not vetoes.issubset(allowed):
                return False, ""
            required = (
                "mapped_center_present",
                "SOURCE_DESTROYED",
                "TARGET_ABSENT_BEFORE",
                "TARGET_CREATED",
                "mapped_nitrogen_gained",
                "hydrogen_gain",
                "mapped_carbon_scaffold_preserved",
            )
            if _required_true(checks, required):
                return True, "vc6_boundary_mapping_with_valid_multicenter_conversion"
            return False, ""

        return False, ""

    if mapping_control == "conversion_validated":
        reference_keys = {
            key
            for rec in list(reference.primary_events) + list(reference.specific_events)
            for key in [conversion_instance_key(rec)]
            if key is not None
        }

        raw_native = {}
        for rows in detector_outputs.values():
            for row in rows:
                if row.get("role") != "directional_conversion":
                    continue
                rec = SemanticRecord.from_mapping(row)
                raw_native[directional_instance_key(rec)] = rec

        rejected_by_key = {}
        for row in rejected:
            rec = SemanticRecord.from_mapping(row)
            rejected_by_key[directional_instance_key(rec)] = row

        accepted_keys = {directional_instance_key(rec) for rec in accepted}
        recovered_keys = set()

        for key, diagnostic in rejected_by_key.items():
            if key not in reference_keys or key not in raw_native:
                continue

            eligible, reason = _conversion_validated_mapping_recovery(diagnostic)
            if not eligible:
                continue

            if key not in accepted_keys:
                accepted.append(raw_native[key])
                accepted_keys.add(key)

            recovered_keys.add(key)
            recoveries.append(
                {
                    "directional_instance_key": [
                        key[0],
                        key[1],
                        list(key[2]),
                    ],
                    "event": key[0],
                    "center_map": key[1],
                    "support_maps": list(key[2]),
                    "strict_validation_status":
                        diagnostic.get("validation_status"),
                    "strict_validator_class":
                        diagnostic.get("validator_class"),
                    "strict_veto_reasons":
                        list(diagnostic.get("veto_reasons") or []),
                    "recovery_policy":
                        "conversion_validated_mapping_tolerance",
                    "recovery_reason": reason,
                }
            )

        if recovered_keys:
            kept = []
            for row in rejected:
                rec = SemanticRecord.from_mapping(row)
                if directional_instance_key(rec) not in recovered_keys:
                    kept.append(row)
            effective_rejected = kept

    deduped = []
    seen = set()
    for rec in accepted:
        key = directional_instance_key(rec)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(rec)
    accepted = deduped

    primary = filter_conversion_derived_rows(
        reference.primary_events,
        accepted,
    )
    specific = filter_conversion_derived_rows(
        reference.specific_events,
        accepted,
    )

    result = deepcopy(reference)
    result.primary_events = primary
    result.specific_events = specific
    result.directional_conversions = list(accepted)
    result.redox_axis = []
    result.semantic_events = filter_semantic_events_to_primary(
        reference.semantic_events,
        primary,
    )

    result.validations = list(validations)
    result.rejected_hypotheses = effective_rejected

    result.raw_hypotheses = {
        "reference": deepcopy(reference.raw_hypotheses),
        "native_directional": deepcopy(detector_outputs),
    }

    counts = Counter(v.status.value for v in validations)
    result.diagnostics.update(
        {
            "backend": "validated_native_fgi_overlay",
            "canonical_fgi_enforcement": True,
            "mapping_control": mapping_control,
            "mapping_control_version": "mapping_control_v3",
            "mapping_control_default": "conversion_validated",
            "mapping_control_recovery_count": len(recoveries),
            "mapping_control_recoveries": recoveries,
            "strict_accepted_directional_instances":
                len(strict_accepted),
            "effective_accepted_directional_instances":
                len(accepted),
            "directional_hypotheses_total": len(validations),
            "directional_validation_counts":
                dict(sorted(counts.items())),
            "accepted_directional_instances": len(accepted),
            "nonaccepted_directional_instances":
                len(effective_rejected),
            "conversion_primary_rows_before": sum(
                bool(r.extra.get("conversion_specific"))
                for r in reference.primary_events
            ),
            "conversion_primary_rows_after": sum(
                bool(r.extra.get("conversion_specific"))
                for r in primary
            ),
            "conversion_specific_rows_before": sum(
                bool(r.extra.get("conversion_specific"))
                for r in reference.specific_events
            ),
            "conversion_specific_rows_after": sum(
                bool(r.extra.get("conversion_specific"))
                for r in specific
            ),
            "canonical_abstention_policy": (
                "STRICT keeps only native ACCEPTED directionals. "
                "CONVERSION_VALIDATED additionally permits validator-specific "
                "mapping tolerance only when chemistry-defining positive checks "
                "prove the source-to-target conversion despite bounded mapping "
                "ambiguity or boundary loss."
            ),
        }
    )
    return result
