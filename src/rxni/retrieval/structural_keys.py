"""Runtime-v6.6 exact structural candidate keys.

M4 proved structural-key parity. M5C mechanically migrates the frozen algorithm
into the v6.6-owned native compatibility closure.
"""

from __future__ import annotations

from typing import Any

from ._core.algorithms.query_structural_keys import (
    derive_query_structural_keys as _derive_native,
)


def derive_structural_keys(
    reaction: str,
    *,
    store_readable: bool = False,
) -> dict[str, Any]:
    return _derive_native(
        reaction,
        store_readable=store_readable,
    )
