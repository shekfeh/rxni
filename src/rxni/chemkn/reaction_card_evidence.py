from __future__ import annotations
from typing import Any


def reaction_card_summary(card) -> dict[str, Any]:
    """Compact runtime-facing view of a normalized reaction card."""
    if card.source_kind != "reaction_card":
        raise ValueError("reaction_card_summary expects a reaction_card source_kind")

    regimes = []
    for r in card.condition_regimes:
        regimes.append({
            "regime_id": r.get("regime_id"),
            "label": r.get("label"),
            "regime_class": r.get("regime_class"),
            "confidence": r.get("confidence"),
            "scope": r.get("regime_scope"),
            "generality": r.get("generality"),
            "usage": r.get("usage"),
            "applicability": (r.get("applicability") or {}).get("summary"),
            "generalized_conditions": (r.get("generalized_conditions") or {}).get("summary"),
            "evidence_refs": r.get("evidence_refs") or [],
            "warnings": r.get("warnings") or [],
        })

    manifest = (card.extensions or {}).get("manifest") or {}
    return {
        "reaction_id": card.source_id,
        "canonical_name": card.canonical_name,
        "reaction_family": card.reaction_family,
        "transformation": card.transformation,
        "recognition": card.recognition,
        "condition_regimes": regimes,
        "evidence_confidence": (card.evidence or {}).get("overall_confidence"),
        "curation": {
            "status": (card.curation or {}).get("status"),
            "confidence": (card.curation or {}).get("confidence"),
            "needs_review": (card.curation or {}).get("needs_review"),
            "review_reasons": (card.curation or {}).get("review_reasons") or [],
            "unresolved_issues": (card.curation or {}).get("unresolved_issues") or [],
            "strict_ledger_status": manifest.get("strict_ledger_status"),
        },
        "master_reaction_name_id": (card.extensions or {}).get("master_reaction_name_id"),
    }
