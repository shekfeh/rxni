"""Parquet + SQLite provider for precomputed Library-v5 System-I records.

The class implements the Runtime-v6.5 RetrievalRecordProvider protocol and can
be passed directly to system_i_reranker.rerank_edit_candidates().
"""

from __future__ import annotations

import sqlite3
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import pyarrow.parquet as pq

from rxni.retrieval._core.foundation.common import chunked, reaction_key

from .retrieval_provider import OnTheFlyRecordProvider
from .stored_record_codec import retrieval_record_from_storage_row


@dataclass
class LibraryRecordProvider:
    representations_path: str | Path
    locator_db_path: str | Path
    fallback_provider: OnTheFlyRecordProvider | None = field(
        default_factory=OnTheFlyRecordProvider
    )

    cache: dict[str, dict[str, Any]] = field(default_factory=dict)
    cache_hits: int = 0
    cache_misses: int = 0
    read_seconds: float = 0.0
    row_groups_read: int = 0
    prefetched_records: int = 0
    missing_stored_records: int = 0

    def __post_init__(self) -> None:
        self.representations_path = Path(self.representations_path)
        self.locator_db_path = Path(self.locator_db_path)
        if not self.representations_path.exists():
            raise FileNotFoundError(self.representations_path)
        if not self.locator_db_path.exists():
            raise FileNotFoundError(self.locator_db_path)
        self._pf = pq.ParquetFile(self.representations_path)

    def _lookup_locations(self, keys: list[str]) -> dict[str, tuple[int, int]]:
        if not keys:
            return {}
        out: dict[str, tuple[int, int]] = {}
        conn = sqlite3.connect(self.locator_db_path)
        try:
            for part in chunked(keys, 800):
                placeholders = ",".join("?" for _ in part)
                query = (
                    "SELECT reaction_key, row_group_id, row_in_group "
                    f"FROM locator WHERE reaction_key IN ({placeholders})"
                )
                for key, rg, ri in conn.execute(query, part):
                    out[str(key)] = (int(rg), int(ri))
        finally:
            conn.close()
        return out

    def _load_locations(
        self,
        requested: dict[str, tuple[str, int, int]],
    ) -> None:
        # requested: reaction_key -> (reaction_string, row_group, row_in_group)
        by_group: dict[int, list[tuple[str, str, int]]] = defaultdict(list)
        for key, (reaction, rg, ri) in requested.items():
            by_group[rg].append((key, reaction, ri))

        for rg in sorted(by_group):
            t0 = time.perf_counter()
            table = self._pf.read_row_group(rg)
            self.read_seconds += time.perf_counter() - t0
            self.row_groups_read += 1

            for key, requested_reaction, row_in_group in by_group[rg]:
                row = table.slice(row_in_group, 1).to_pylist()[0]
                if str(row.get("reaction_key")) != key:
                    raise RuntimeError(
                        f"Locator mismatch for {key}: row group {rg}, row {row_in_group}"
                    )
                record = retrieval_record_from_storage_row(row)
                stored_reaction = str(record.get("reaction") or "")
                if stored_reaction != requested_reaction:
                    raise RuntimeError(
                        "Stored reaction string does not match reaction_key request"
                    )
                self.cache[requested_reaction] = record
                self.prefetched_records += 1

    def prefetch(self, requests: Iterable[tuple[str, str]]) -> None:
        unique: dict[str, str] = {}
        for reaction, _record_id in requests:
            if reaction in self.cache:
                continue
            unique[reaction_key(reaction)] = reaction

        locations = self._lookup_locations(list(unique))
        requested: dict[str, tuple[str, int, int]] = {}
        for key, reaction in unique.items():
            loc = locations.get(key)
            if loc is None:
                self.missing_stored_records += 1
                continue
            requested[key] = (reaction, loc[0], loc[1])

        self._load_locations(requested)

    def get(self, reaction: str, *, record_id: str = "") -> dict[str, Any]:
        if reaction in self.cache:
            self.cache_hits += 1
            return self.cache[reaction]

        self.cache_misses += 1
        key = reaction_key(reaction)
        loc = self._lookup_locations([key]).get(key)
        if loc is not None:
            self._load_locations({key: (reaction, loc[0], loc[1])})
            return self.cache[reaction]

        self.missing_stored_records += 1
        if self.fallback_provider is None:
            raise KeyError(f"No stored System-I representation for reaction {key}")

        record = self.fallback_provider.get(reaction, record_id=record_id)
        self.cache[reaction] = record
        return record

    def snapshot(self) -> dict[str, Any]:
        fallback = self.fallback_provider.snapshot() if self.fallback_provider is not None else {}
        return {
            "provider": "library_v5",
            "cache_size": len(self.cache),
            "cache_hits": self.cache_hits,
            "cache_misses": self.cache_misses,
            # Common field expected by system_i_reranker diagnostics.
            "build_seconds": float(fallback.get("build_seconds", 0.0)),
            "read_seconds": self.read_seconds,
            "row_groups_read": self.row_groups_read,
            "prefetched_records": self.prefetched_records,
            "missing_stored_records": self.missing_stored_records,
            "fallback_cache_size": int(fallback.get("cache_size", 0)),
            "fallback_cache_misses": int(fallback.get("cache_misses", 0)),
        }
