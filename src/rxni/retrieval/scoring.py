"""Runtime-v6.6 retrieval scoring.

The scientific equations are a mechanically migrated, hash-audited copy of the
frozen Runtime-v6.5 scorer. M5C changes namespace ownership only.
"""

from __future__ import annotations

from typing import Any

from ._core.algorithms.scoring import member_scores as _member_scores_native


def member_scores(
    query_record: dict[str, Any],
    candidate_record: dict[str, Any],
) -> dict[str, float]:
    return dict(_member_scores_native(query_record, candidate_record))


def template_score(
    scores: list[float],
    *,
    max_weight: float = 0.75,
    mean_weight: float = 0.25,
    top_n: int = 3,
) -> float:
    if not scores:
        return 0.0
    vals = sorted((float(x) for x in scores), reverse=True)
    best = vals[0]
    mean_top = sum(vals[:top_n]) / len(vals[:top_n])
    return max_weight * best + mean_weight * mean_top
