"""Versioned *source* schema adapters for RXNI retrieval V2 representation builds.

This module does not define or change the on-disk Retrieval-V2 ABI. Source
profiles only specify how to obtain (opaque key, forward reaction SMILES) from
input Parquet/CSV rows. They do NOT make a source library retrieval-ready:
external structural/metadata indices must use the same keys.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

SCHEMA_VERSION = 1
FIELDS = ("reaction_key", "reaction_smiles", "valid", "feature_error")


class SourceSchemaError(ValueError):
    """Invalid source profile or incompatible source data."""


@dataclass(frozen=True)
class SourceProfile:
    source_format: str
    key_column: str
    reaction_column: str
    valid_column: str | None = None
    feature_error_column: str | None = None
    profile_id: str = "external"

    def columns(self) -> list[str]:
        return list(dict.fromkeys(c for c in (self.key_column, self.reaction_column,
                     self.valid_column, self.feature_error_column) if c))

    def as_dict(self) -> dict[str, Any]:
        return {"schema_version": SCHEMA_VERSION, "profile_id": self.profile_id,
                "source_format": self.source_format,
                "columns": {"reaction_key": self.key_column,
                            "reaction_smiles": self.reaction_column,
                            "valid": self.valid_column,
                            "feature_error": self.feature_error_column}}


def load_profile(path: str | Path) -> SourceProfile:
    path = Path(path)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceSchemaError(f"Cannot load source profile {path}: {exc}") from exc
    if not isinstance(raw, dict) or raw.get("schema_version") != SCHEMA_VERSION:
        raise SourceSchemaError("Source schema profile requires schema_version=1")
    allowed = {"schema_version", "profile_id", "source_format", "columns"}
    if set(raw) - allowed:
        raise SourceSchemaError(f"Unknown profile fields: {sorted(set(raw)-allowed)}")
    fmt = raw.get("source_format")
    if fmt not in {"parquet", "csv"}:
        raise SourceSchemaError("source_format must be 'parquet' or 'csv'")
    cols = raw.get("columns")
    if not isinstance(cols, dict) or set(cols)-set(FIELDS):
        raise SourceSchemaError("columns must contain recognized canonical fields only")
    for key in ("reaction_key", "reaction_smiles"):
        if not isinstance(cols.get(key), str) or not cols[key].strip():
            raise SourceSchemaError(f"Required column mapping missing: {key}")
    for key in ("valid", "feature_error"):
        if cols.get(key) is not None and (not isinstance(cols[key], str) or not cols[key].strip()):
            raise SourceSchemaError(f"Invalid optional column mapping: {key}")
    for k, value in cols.items():
        if value is not None and ("\x00" in value or "\n" in value or "\r" in value):
            raise SourceSchemaError(f"Unsafe column name: {k}")
    if cols["reaction_key"] == cols["reaction_smiles"]:
        raise SourceSchemaError("reaction_key and reaction_smiles must use different columns")
    profile_id = raw.get("profile_id", "external")
    if not isinstance(profile_id, str) or not profile_id or len(profile_id) > 128:
        raise SourceSchemaError("profile_id must be a nonempty string of <=128 characters")
    return SourceProfile(fmt, cols["reaction_key"], cols["reaction_smiles"],
                         cols.get("valid"), cols.get("feature_error"), profile_id)


def _check_suffix(path: Path, fmt: str) -> None:
    expected = ".parquet" if fmt == "parquet" else ".csv"
    if path.suffix.lower() != expected:
        raise SourceSchemaError(f"Profile format {fmt!r} requires {expected} input: {path}")


def source_columns(path: Path, profile: SourceProfile) -> list[str]:
    _check_suffix(path, profile.source_format)
    if not path.is_file():
        raise FileNotFoundError(path)
    if profile.source_format == "parquet":
        import pyarrow.parquet as pq
        schema = pq.read_schema(path)
        names = schema.names
        # Reject missing columns before invoking a worker pool.
        for col in (profile.key_column, profile.reaction_column):
            field = schema.field(col) if col in names else None
            if field is None:
                continue
            import pyarrow as pa
            if not (pa.types.is_string(field.type) or pa.types.is_large_string(field.type)):
                raise SourceSchemaError(f"Expected a UTF-8 string field for {col!r}, got {field.type}")
    else:
        with path.open("r", newline="", encoding="utf-8-sig") as stream:
            names = csv.DictReader(stream).fieldnames or []
    missing = set(profile.columns()) - set(names)
    if missing:
        raise SourceSchemaError(f"Missing source columns: {sorted(missing)}; available: {names}")
    return names


def _accepted(value: Any, column: str) -> bool:
    if isinstance(value, bool):
        return value
    if value in (1, 0):
        return bool(value)
    if isinstance(value, str):
        match = value.strip().lower()
        if match in ("true", "1"):
            return True
        if match in ("false", "0"):
            return False
    raise SourceSchemaError(f"Unknown validity flag {value!r} in column {column}")


def _iter_rows(path: Path, profile: SourceProfile) -> Iterator[dict[str, Any]]:
    if profile.source_format == "parquet":
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(columns=profile.columns(), batch_size=8192):
            yield from batch.to_pylist()
    else:
        with path.open("r", newline="", encoding="utf-8-sig") as stream:
            yield from csv.DictReader(stream)


def iter_source(path: str | Path, profile: SourceProfile, limit: int | None = None):
    """Yield deterministic (ordinal,key,reaction), preserving opaque source keys.

    Invalid source rows are skipped as in the historical builder. Malformed
    validity markers fail closed rather than being silently treated as true.
    """
    path = Path(path).expanduser().resolve()
    source_columns(path, profile)
    emitted = 0
    for row in _iter_rows(path, profile):
        if profile.valid_column and not _accepted(row.get(profile.valid_column), profile.valid_column):
            continue
        if profile.feature_error_column and str(row.get(profile.feature_error_column) or "").strip():
            continue
        key, reaction = row.get(profile.key_column), row.get(profile.reaction_column)
        if key is None or reaction is None or not str(key).strip() or not str(reaction).strip():
            continue
        # This API requires the forward, three-component reaction form. A
        # standalone SMILES is not safely interchangeable with a reaction.
        reaction = str(reaction)
        if reaction.count(">") != 2:
            raise SourceSchemaError(f"Expected reactants>agents>products reaction; key={key!r}")
        yield emitted, str(key), reaction
        emitted += 1
        if limit is not None and emitted >= limit:
            return


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def describe(path: Path, profile: SourceProfile) -> dict[str, Any]:
    available = source_columns(path, profile)
    return {"status": "SOURCE_SCHEMA_COMPATIBLE", "source": str(path.resolve()),
            "profile": profile.as_dict(), "available_columns": available,
            "notice": "Source schema accepted; full retrieval readiness requires matching structural index, metadata lookup, and V2 materialization."}


def main():
    p = argparse.ArgumentParser(description="Inspect RXNI builder input schemas; does not modify source files")
    p.add_argument("action", choices=["inspect"])
    p.add_argument("--source-core", required=True)
    p.add_argument("--input-schema", required=True)
    args = p.parse_args()
    print(json.dumps(describe(Path(args.source_core), load_profile(args.input_schema)), indent=2))


if __name__ == "__main__":
    main()
