"""Frozen System-I retrieval scoring.

Initial behavior is preserved from the validated Tool-93 / relaxed-retrieval
experiments. Do not alter weights or similarity definitions without a new
retrieval experiment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class RetrievalScoringConfig:
    local_r_weight: float = 0.30
    local_p_weight: float = 0.30
    local_edit_weight: float = 0.40

    struct_weight: float = 0.75
    semantic_weight: float = 0.25

    local_combined_weight: float = 0.55
    multi_combined_weight: float = 0.20
    semantic_combined_weight: float = 0.25

    template_max_weight: float = 0.75
    template_mean_weight: float = 0.25
    template_top_n: int = 3

    def validate(self) -> None:
        checks = (
            (
                self.local_r_weight
                + self.local_p_weight
                + self.local_edit_weight,
                "local weights",
            ),
            (
                self.struct_weight + self.semantic_weight,
                "structural/semantic weights",
            ),
            (
                self.local_combined_weight
                + self.multi_combined_weight
                + self.semantic_combined_weight,
                "combined weights",
            ),
            (
                self.template_max_weight
                + self.template_mean_weight,
                "template weights",
            ),
        )

        for value, name in checks:
            if abs(value - 1.0) > 1e-8:
                raise ValueError(f"{name} must sum to 1; got {value}")


def sparse_dict(
    payload: dict[str, Any] | None,
) -> dict[int, float]:
    if not payload:
        return {}

    idx = payload.get("indices") or []
    cnt = payload.get("counts") or []

    return {
        int(i): float(c)
        for i, c in zip(idx, cnt)
        if float(c) != 0.0
    }


def count_tanimoto(
    a: dict[str, Any] | None,
    b: dict[str, Any] | None,
) -> float:
    """Min/max Tanimoto for non-negative sparse count vectors."""

    da = sparse_dict(a)
    db = sparse_dict(b)

    if not da and not db:
        return 1.0

    keys = set(da) | set(db)

    numerator = sum(
        min(da.get(k, 0.0), db.get(k, 0.0))
        for k in keys
    )

    denominator = sum(
        max(da.get(k, 0.0), db.get(k, 0.0))
        for k in keys
    )

    return numerator / denominator if denominator else 0.0


def weighted_jaccard(
    a: dict[str, float],
    b: dict[str, float],
) -> float:

    if not a and not b:
        return 1.0

    keys = set(a) | set(b)

    numerator = sum(
        min(float(a.get(k, 0.0)), float(b.get(k, 0.0)))
        for k in keys
    )

    denominator = sum(
        max(float(a.get(k, 0.0)), float(b.get(k, 0.0)))
        for k in keys
    )

    return numerator / denominator if denominator else 0.0


def center_pair_score(
    a: dict[str, Any],
    b: dict[str, Any],
    *,
    r_weight: float = 0.30,
    p_weight: float = 0.30,
    edit_weight: float = 0.40,
) -> float:

    ac = a.get("channels") or {}
    bc = b.get("channels") or {}

    return (
        r_weight
        * count_tanimoto(
            ac.get("reactant_local"),
            bc.get("reactant_local"),
        )
        + p_weight
        * count_tanimoto(
            ac.get("product_local"),
            bc.get("product_local"),
        )
        + edit_weight
        * count_tanimoto(
            ac.get("edit_local"),
            bc.get("edit_local"),
        )
    )


def local_score(
    a: dict[str, Any],
    b: dict[str, Any],
    *,
    r_weight: float = 0.30,
    p_weight: float = 0.30,
    edit_weight: float = 0.40,
) -> float:
    """Symmetric best-center average."""

    ca = (
        (a.get("localized_fingerprints") or {})
        .get("centers")
        or []
    )

    cb = (
        (b.get("localized_fingerprints") or {})
        .get("centers")
        or []
    )

    if not ca or not cb:
        return 0.0

    matrix = [
        [
            center_pair_score(
                x,
                y,
                r_weight=r_weight,
                p_weight=p_weight,
                edit_weight=edit_weight,
            )
            for y in cb
        ]
        for x in ca
    ]

    q_to_c = (
        sum(max(row) for row in matrix)
        / len(ca)
    )

    c_to_q = (
        sum(
            max(matrix[i][j] for i in range(len(ca)))
            for j in range(len(cb))
        )
        / len(cb)
    )

    return 0.5 * (q_to_c + c_to_q)


def multi_score(
    a: dict[str, Any],
    b: dict[str, Any],
) -> float:

    ac = (
        (a.get("multichannel_fingerprint") or {})
        .get("channels")
        or {}
    )

    bc = (
        (b.get("multichannel_fingerprint") or {})
        .get("channels")
        or {}
    )

    weights = {
        "reactant_global": 0.20,
        "product_global": 0.20,
        "reactant_local": 0.20,
        "product_local": 0.20,
        "edit": 0.20,
    }

    return sum(
        weight * count_tanimoto(ac.get(key), bc.get(key))
        for key, weight in weights.items()
    )


def semantic_score(
    a: dict[str, Any],
    b: dict[str, Any],
) -> float:

    return weighted_jaccard(
        a.get("semantic_profile") or {},
        b.get("semantic_profile") or {},
    )


def member_scores(
    a: dict[str, Any],
    b: dict[str, Any],
    config: RetrievalScoringConfig | None = None,
) -> dict[str, float]:

    cfg = config or RetrievalScoringConfig()
    cfg.validate()

    local = local_score(
        a,
        b,
        r_weight=cfg.local_r_weight,
        p_weight=cfg.local_p_weight,
        edit_weight=cfg.local_edit_weight,
    )

    multi = multi_score(a, b)
    semantic = semantic_score(a, b)

    return {
        "local": local,
        "multi": multi,
        "semantic": semantic,

        "sem_local": (
            cfg.struct_weight * local
            + cfg.semantic_weight * semantic
        ),

        "sem_multi": (
            cfg.struct_weight * multi
            + cfg.semantic_weight * semantic
        ),

        "sem_local_multi": (
            cfg.local_combined_weight * local
            + cfg.multi_combined_weight * multi
            + cfg.semantic_combined_weight * semantic
        ),
    }


def template_score(
    scores: list[float],
    config: RetrievalScoringConfig | None = None,
) -> float:

    cfg = config or RetrievalScoringConfig()
    cfg.validate()

    values = sorted(
        (float(x) for x in scores),
        reverse=True,
    )

    if not values:
        return 0.0

    top = values[: max(1, cfg.template_top_n)]

    return (
        cfg.template_max_weight * values[0]
        + cfg.template_mean_weight
        * (sum(top) / len(top))
    )
