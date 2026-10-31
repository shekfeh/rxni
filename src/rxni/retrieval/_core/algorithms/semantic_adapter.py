"""Runtime-v6.7 semantic adapter for the frozen structural retrieval code.

The structural representation implementation remains the validated v6.6
generation, but semantic classification must resolve through Runtime-v6.7 so
current detector/resolver changes are used for live query construction.
"""

from __future__ import annotations

from typing import Any

from rxni.semantic import classify_reaction


def classify_for_retrieval(mapped_reaction: str) -> dict[str, Any]:
    return classify_reaction(
        mapped_reaction,
        backend="validated",
    ).as_dict()
