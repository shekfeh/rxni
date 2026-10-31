"""Conservative structural semantic-extension hypotheses.

These rules are retained for parity while they are decomposed into chemical
domains. The public v6.6 name is architectural, not historical.
"""

from __future__ import annotations

from typing import Any


def detect_structural_semantic_extensions(
    facts: dict[str, Any],
) -> list[dict[str, Any]]:
    from . import _generated_fgi_rules as r
    return r.detect_structural_semantic_extensions(facts)
