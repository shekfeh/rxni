"""Small shared utilities for builder_v5."""

from __future__ import annotations

import hashlib
import json
import os
import numbers
from pathlib import Path
from typing import Any


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_json(path: str | Path, payload: Any) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, p)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def reaction_key(reaction: str) -> str:
    return sha256_text(reaction)


def require_u64(value: Any, field: str) -> int | None:
    """Return exact uint64 or fail loudly.

    Floats are intentionally forbidden. A hash that touched float64 may have
    silently lost low bits; accepting it would reintroduce the historical V4
    corruption mode this builder was created to remove.
    """
    if value is None:
        return None
    if isinstance(value, bool):
        raise TypeError(f"{field}: bool is not a valid uint64 hash")
    if isinstance(value, numbers.Real) and not isinstance(value, numbers.Integral):
        raise TypeError(f"{field}: floating-point hash forbidden: {value!r}")
    iv = int(value)
    if iv < 0 or iv > 0xFFFFFFFFFFFFFFFF:
        raise OverflowError(f"{field}: outside uint64 range: {iv}")
    return iv


def u64_to_sql(value: Any, field: str = "hash") -> str | None:
    iv = require_u64(value, field)
    return None if iv is None else str(iv)


def sql_to_u64(value: Any, field: str = "hash") -> int | None:
    if value is None:
        return None
    if isinstance(value, float):
        raise TypeError(f"{field}: float recovered from SQL: {value!r}")
    return require_u64(int(str(value)), field)


def file_sha256(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()
