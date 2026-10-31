"""Runtime-v6.7 live query retrieval-record construction.

The retrieval representation algorithms remain the frozen Runtime-v6.6
generation. The semantic classification used to build the live *query* record
is intentionally Runtime-v6.7 so current detector/resolver changes propagate
into the CLI and query semantic profile.

Stored Library-v5 / Retrieval-V2 candidate records remain v6.6 generation until
an explicit representation rebuild is performed.

This module does not reproduce a historical semantic-profile function.
Canonical ``SemanticResult.semantic_profile`` from Runtime-v6.7 is used
directly.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from rxni.semantic import classify_reaction

from ._core.algorithms.retrieval_representation import (
    extract_reaction_centers,
    generate_localized_fingerprints,
    generate_multichannel_fingerprint,
)
from .version import LIVE_RECORD_VERSION, RETRIEVAL_API_VERSION


def build_retrieval_record(
    reaction: str,
    *,
    record_id: str | None = None,
    radius: int = 2,
    fp_size: int = 2048,
) -> dict[str, Any]:
    canonical = classify_reaction(
        reaction,
        backend="validated",
    )
    sem = canonical.as_dict()

    centers = extract_reaction_centers(
        reaction,
        classification=sem,
    )
    localized = generate_localized_fingerprints(
        reaction,
        classification=sem,
        radius=radius,
        fp_size=fp_size,
    )
    multi = generate_multichannel_fingerprint(
        reaction,
        classification=sem,
        radius=radius,
        fp_size=fp_size,
    )

    return {
        "record_id": "" if record_id is None else str(record_id),
        "reaction": reaction,
        "semantic_profile": deepcopy(canonical.semantic_profile),
        "semantics": sem,
        "reaction_centers": centers,
        "localized_fingerprints": localized,
        "multichannel_fingerprint": multi,
        "metadata": {
            "retrieval_api_version": RETRIEVAL_API_VERSION,
            "retrieval_record_version": LIVE_RECORD_VERSION,
            "semantic_engine_version": canonical.engine_version,
            "semantic_rules_version": canonical.rules_version,
            "semantic_backend": "validated",
            "semantic_runtime": "runtime_v6_7",
            "retrieval_generation": "runtime_v6_6_retrieval_v2",
            "structural_representation_implementation":
                "rxni.retrieval._core",
            "m5c_native_retrieval": True,
        },
    }
