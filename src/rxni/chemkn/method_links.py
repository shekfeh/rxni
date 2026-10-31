from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class MethodLinkStore:
    """Explicit method_id -> knowledge source crosswalk.

    No fuzzy named-reaction inference is performed here. That is deliberate:
    benchmarkable method IDs must have auditable links.
    """

    def __init__(self, rows: list[dict[str, Any]]):
        self.rows = rows
        self.by_method: dict[str, list[dict[str, Any]]] = {}
        for row in rows:
            mid = str(row["method_id"])
            self.by_method.setdefault(mid, []).extend(row.get("links") or [])

    @classmethod
    def load(cls, path: str | Path):
        rows = []
        with Path(path).open(encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    rows.append(json.loads(line))
        return cls(rows)

    def links_for(self, method_id: str) -> list[dict[str, Any]]:
        return list(self.by_method.get(str(method_id), []))


def resolve_method_links(method_id: str, registry, link_store: MethodLinkStore):
    out = []
    for link in link_store.links_for(method_id):
        source = link.get("source")
        source_id = link.get("source_id")
        card = registry.get(source, source_id) if source and source_id else None
        out.append({
            "method_id": method_id,
            "link": link,
            "resolved": card is not None,
            "card": None if card is None else card.as_dict(),
        })
    return out
