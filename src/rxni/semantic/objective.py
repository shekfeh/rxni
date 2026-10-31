"""Native objective semantic inputs for Runtime-v6.6.

This module owns objective fields that are direct functions of the mapped
reaction graph and therefore should not be inherited from the v6.5 semantic
reference path.

M6C1B scope:
    formed_bond_classes only.

Primary/specific semantic labels are intentionally out of scope.
"""

from __future__ import annotations

from typing import Any

from .structure.mapped_graph import build_mapped_graph_facts


def native_formed_bond_classes(mapped_reaction: str) -> list[Any]:
    """Return formed-bond classes from native Runtime-v6.6 mapped facts."""
    facts = build_mapped_graph_facts(mapped_reaction)
    return list(facts.get("formed_bond_classes") or [])
