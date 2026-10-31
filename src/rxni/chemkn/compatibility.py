from __future__ import annotations
from .knowledge_schema import (
    CompatibilityDecision,
    CompatibilitySignature,
    RuntimeCompatibilityQuery,
)


_ANY_DIRECTIONS = {"", "any", "none", None}


def evaluate_compatibility(
    signature: CompatibilitySignature,
    query: RuntimeCompatibilityQuery,
) -> CompatibilityDecision:
    """Evaluate optional source-neutral compatibility.

    This deliberately treats hard and soft evidence differently. A hard miss
    excludes a card; a soft miss lowers support but does not exclude it.
    """
    hard: list[str] = []
    matched: list[str] = []
    soft_missing: list[str] = []
    score = 1.0

    req = set(signature.required_events)
    overlap = req & query.events
    if req:
        if overlap:
            matched.append(f"required_event:{sorted(overlap)}")
        elif signature.required_event_mode == "hard":
            hard.append(f"required_event_any:{sorted(req)}")
        elif signature.required_event_mode == "soft":
            score -= 0.20
            soft_missing.append(f"required_event_any:{sorted(req)}")

    if signature.event_alternatives:
        alt_hit = False
        for alt in signature.event_alternatives:
            events = set(alt.get("events_any") or [])
            if events & query.events:
                alt_hit = True
                matched.append(f"event_alternative:{alt.get('name') or sorted(events)}")
                break
        if not alt_hit:
            if signature.event_alternatives_mode == "hard":
                hard.append("event_alternatives")
            elif signature.event_alternatives_mode == "soft":
                score -= 0.15
                soft_missing.append("event_alternatives")

    if signature.direction not in _ANY_DIRECTIONS and query.direction:
        if signature.direction == query.direction:
            matched.append(f"direction:{query.direction}")
        else:
            hard.append(f"direction:{signature.direction}")

    if signature.allowed_donor_elements and query.donor_elements:
        allowed = set(signature.allowed_donor_elements)
        if allowed & query.donor_elements:
            matched.append(f"donor_element:{sorted(allowed & query.donor_elements)}")
        else:
            hard.append(f"allowed_donor_elements:{sorted(allowed)}")

    if signature.allowed_electrophile_classes and query.electrophile_classes:
        allowed = set(signature.allowed_electrophile_classes)
        if allowed & query.electrophile_classes:
            matched.append(f"electrophile:{sorted(allowed & query.electrophile_classes)}")
        else:
            hard.append(f"allowed_electrophile_classes:{sorted(allowed)}")

    motifs = query.motifs
    if signature.required_motifs:
        missing = set(signature.required_motifs) - motifs
        if missing:
            hard.append(f"required_motifs:{sorted(missing)}")
        else:
            matched.append("required_motifs")

    if signature.forbidden_motifs:
        forbidden = set(signature.forbidden_motifs) & motifs
        if forbidden:
            hard.append(f"forbidden_motifs:{sorted(forbidden)}")

    if signature.preferred_motifs:
        pref = set(signature.preferred_motifs)
        if pref & motifs:
            matched.append(f"preferred_motifs:{sorted(pref & motifs)}")
        else:
            score -= 0.05
            soft_missing.append("preferred_motifs")

    if hard:
        score = 0.0

    return CompatibilityDecision(
        compatible=not hard,
        score=max(0.0, min(1.0, score)),
        hard_failures=hard,
        matched=matched,
        missing_soft=soft_missing,
    )
