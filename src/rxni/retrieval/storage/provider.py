"""Build-scoped stored provider for retrieval_representations_v2.parquet."""

from __future__ import annotations

import json
import os
import sqlite3
from collections import OrderedDict
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from .codec import retrieval_record_from_storage_row_v2
from .custom_profiles import load_custom_profile, validate_schema, normalize_row
from ..version import STORED_REPRESENTATION_GENERATION


class V2StoredRecordProvider:
    representation_generation = STORED_REPRESENTATION_GENERATION

    def __init__(
        self,
        representations_path: str | Path,
        locator_path: str | Path,
        manifest_path: str | Path,
        *,
        row_group_cache_size: int = 4,
        schema_profile: str | Path | None = None,
    ):
        self.representations_path = Path(representations_path)
        self.locator_path = Path(locator_path)
        self.manifest_path = Path(manifest_path)
        self.row_group_cache_size = max(1, int(row_group_cache_size))

        if not self.representations_path.exists():
            raise FileNotFoundError(self.representations_path)
        if not self.locator_path.exists():
            raise FileNotFoundError(self.locator_path)
        if not self.manifest_path.exists():
            raise FileNotFoundError(self.manifest_path)

        self.manifest = json.loads(self.manifest_path.read_text())
        if (
            self.manifest.get("representation_generation")
            != STORED_REPRESENTATION_GENERATION
        ):
            raise ValueError(
                "Manifest representation generation mismatch"
            )
        self.build_id = str(self.manifest["build_id"])
        self._pf = pq.ParquetFile(self.representations_path)
        # Explicit per-provider profile, or opt-in process configuration.
        # Default keeps all previous canonical and historical reads unchanged.
        selected = schema_profile if schema_profile is not None else os.environ.get("RXNI_V2_SCHEMA_PROFILE")
        self._custom_profile = load_custom_profile(selected) if selected else None
        if self._custom_profile is not None:
            validate_schema(self._custom_profile, self._pf.schema_arrow)
            if self.manifest.get("physical_schema_profile") != self._custom_profile.profile_id:
                raise ValueError("Custom schema profile not declared in library manifest")
            if self.manifest.get("custom_schema_profile_sha256") != self._custom_profile.sha256:
                raise ValueError("Custom schema profile SHA256 not attested by manifest")
        self._cache: OrderedDict[int, list[dict[str, Any]]] = OrderedDict()

    def _lookup(self, reaction_key: str) -> tuple[int, int]:
        conn = sqlite3.connect(self.locator_path)
        try:
            row = conn.execute(
                """
                SELECT row_group_id, row_in_group
                FROM representation_locator_v2
                WHERE reaction_key=?
                  AND build_id=?
                  AND representation_generation=?
                """,
                (
                    str(reaction_key),
                    self.build_id,
                    STORED_REPRESENTATION_GENERATION,
                ),
            ).fetchone()
        finally:
            conn.close()
        if row is None:
            raise KeyError(reaction_key)
        return int(row[0]), int(row[1])

    def _row_group(self, row_group_id: int) -> list[dict[str, Any]]:
        if row_group_id in self._cache:
            rows = self._cache.pop(row_group_id)
            self._cache[row_group_id] = rows
            return rows

        rows = self._pf.read_row_group(row_group_id).to_pylist()
        if self._custom_profile is not None:
            rows = [normalize_row(self._custom_profile, row) for row in rows]
        self._cache[row_group_id] = rows
        while len(self._cache) > self.row_group_cache_size:
            self._cache.popitem(last=False)
        return rows

    def get(
        self,
        reaction: str | None = None,
        *,
        record_id: str | None = None,
    ) -> dict[str, Any]:
        if not record_id:
            raise ValueError(
                "V2StoredRecordProvider requires record_id=reaction_key"
            )
        rg, ri = self._lookup(str(record_id))
        rows = self._row_group(rg)
        if not (0 <= ri < len(rows)):
            raise IndexError(
                f"locator row_in_group={ri} outside row group {rg}"
            )
        row = rows[ri]
        actual = str(row.get("reaction_key") or row.get("record_id") or "")
        if actual and actual != str(record_id):
            raise RuntimeError(
                "Locator/store key mismatch: "
                f"requested={record_id} stored={actual}"
            )
        return retrieval_record_from_storage_row_v2(
            row,
            expected_build_id=self.build_id,
        )

    def prefetch(self, reaction_keys) -> None:
        """Warm relevant row groups without changing scientific behavior."""
        keys = [str(k) for k in reaction_keys]
        if not keys:
            return
        conn = sqlite3.connect(self.locator_path)
        try:
            qmarks = ",".join("?" for _ in keys)
            rows = conn.execute(
                f"""
                SELECT DISTINCT row_group_id
                FROM representation_locator_v2
                WHERE build_id=?
                  AND representation_generation=?
                  AND reaction_key IN ({qmarks})
                """,
                [
                    self.build_id,
                    STORED_REPRESENTATION_GENERATION,
                    *keys,
                ],
            ).fetchall()
        finally:
            conn.close()
        for (rg,) in rows:
            self._row_group(int(rg))

    def snapshot(self) -> dict[str, Any]:
        return {
            **({"physical_schema_profile": self._custom_profile.profile_id}
               if self._custom_profile is not None else {}),
            "provider": type(self).__name__,
            "representation_generation": self.representation_generation,
            "build_id": self.build_id,
            "representations_path": str(self.representations_path),
            "locator_path": str(self.locator_path),
        }
