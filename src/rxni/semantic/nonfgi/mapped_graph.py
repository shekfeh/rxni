"""Runtime-v6.7 mapped-fact boundary for migrated non-FGI rules.

Runtime-v6.7 now owns the mapped-graph implementation locally.  The implementation
preserves the frozen v6.6 base-fact contract and includes the additive v6.7
full-ring lineage extension used by repaired/specialized ring detectors.
"""

from rxni.semantic.structure.mapped_graph import (
    build_mapped_graph_facts as _build_mapped_graph_facts_v67,
    canonical_bond_class,
)


def build_mapped_graph_facts(reaction: str):
    """Return v6.7-owned graph facts preserving the v6.6 base contract."""
    return _build_mapped_graph_facts_v67(reaction)


__all__ = [
    "build_mapped_graph_facts",
    "canonical_bond_class",
]
