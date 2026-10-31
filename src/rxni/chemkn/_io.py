from __future__ import annotations
import json
from pathlib import Path
from typing import Any, Iterator


def read_jsonl(path: str | Path) -> Iterator[dict[str, Any]]:
    p = Path(path)
    with p.open(encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def index_jsonl(path: str | Path, key: str) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for row in read_jsonl(path):
        value = row.get(key)
        if value is not None:
            out[str(value)] = row
    return out


def str_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(x) for x in value if x is not None and str(x).strip()]
    s = str(value).strip()
    return [s] if s else []
