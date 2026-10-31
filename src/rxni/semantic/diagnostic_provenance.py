"""Canonical diagnostic provenance for the self-contained Runtime-v6.6 semantic API.

M6C3 removes retriever/reference-bridge provenance from the active semantic
package.  This module changes diagnostics only; scientific semantic fields are
not modified.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


LEGACY_DIAGNOSTIC_KEYS = {
    "reference_semantic_profile",
    "legacy_reference_backend",
    "legacy_retriever",
}


def normalize_diagnostic_provenance(
    diagnostics: dict[str, Any] | None,
    *,
    backend: str,
) -> dict[str, Any]:
    out = deepcopy(diagnostics or {})

    for key in LEGACY_DIAGNOSTIC_KEYS:
        out.pop(key, None)

    # Normalize the provenance-bearing keys created during staged migration.
    out["semantic_runtime"] = "runtime_v6_6"
    out["semantic_backend"] = str(backend)
    out["legacy_runtime_v6_5_dependency"] = False
    out["mapped_fact_source"] = (
        "runtime_v6_6.semantic.structure.mapped_graph"
    )
    out["nonfgi_rule_engine"] = (
        "runtime_v6_6.semantic.nonfgi."
    )

    mapping_recovery = dict(out.get("mapping_recovery") or {})
    mapping_recovery.update(
        {
            "source": "native_v6_6_fact_boundary",
            "legacy_mapping_correction_invoked": False,
        }
    )
    out["mapping_recovery"] = mapping_recovery

    # Keep enforcement/validation diagnostics already emitted by Stage E.
    return out
