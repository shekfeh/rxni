"""Native structural-facts public boundary.

The implementation is vendored into `structure/mapped_graph.py` from the frozen
v6.5.6 mapped-graph module by the Milestone-2 installer. It is copied, not
imported at runtime, so v6.6 owns its structural implementation after install.
"""

from __future__ import annotations

from typing import Any

try:
    from .mapped_graph import build_mapped_graph_facts as _build
except ImportError as exc:  # installer not yet run
    _build = None
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def build_structural_facts(reaction: str) -> dict[str, Any]:
    if _build is None:
        raise RuntimeError(
            "Native v6.6 mapped_graph.py is not installed. "
            "Run tools/install_v660_milestone2.py from the repository root."
        ) from _IMPORT_ERROR
    return _build(reaction)


# Stable public alias for callers that want the historical function spelling.
build_mapped_graph_facts = build_structural_facts
