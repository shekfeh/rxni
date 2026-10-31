"""Shared helpers for the Template Library v5 builder."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable


def reaction_key(reaction: str) -> str:
    """Stable key for the exact mapped reaction string.

    This intentionally matches the current reranker's exact-reaction deduplication
    semantics: no remapping, canonicalization, or chemistry normalization is added.
    """
    return hashlib.sha256(reaction.encode("utf-8")).hexdigest()


def u64_to_sql(value: Any) -> str | None:
    if value is None:
        return None
    try:
        # Zero padding makes SQLite TEXT ordering agree with unsigned numeric order.
        return f"{int(value):020d}"
    except Exception:
        return None


def sql_to_u64(value: Any) -> int | None:
    if value in (None, ""):
        return None
    return int(value)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_json(path: str | Path, payload: Any) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def chunked(values: Iterable[Any], n: int):
    buf = []
    for value in values:
        buf.append(value)
        if len(buf) >= n:
            yield buf
            buf = []
    if buf:
        yield buf
