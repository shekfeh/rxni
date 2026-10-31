"""Explicit, fail-closed physical-schema mappings for RXNI V2 read adapters.

No scientific reclassification, JSON rewriting, inferred aliases, or codec changes.
The canonical 10-column writer is immutable. Additional physical columns must
be declared as provenance/annotations and are not retrieval features.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from .schema_profiles import CANONICAL_COLUMNS, sha256

CANONICAL_NAMES = tuple(name for name, _ in CANONICAL_COLUMNS)
_IDENTIFIER = re.compile(r"^[a-zA-Z][a-zA-Z0-9_.-]{0,127}$")
_STRING_TYPES = {"string", "large_string"}


@dataclass(frozen=True)
class CustomV2Profile:
    profile_id: str
    columns: dict[str, str]
    extra_columns: dict[str, dict[str, str]]
    sha256: str


def load_custom_profile(path: str | Path) -> CustomV2Profile:
    p = Path(path).expanduser().resolve()
    if not p.is_file() or p.stat().st_size > 65536:
        raise ValueError("Schema profile file missing or too large")
    data = json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or set(data) != {"schema_profile_version", "profile_id", "columns", "extra_columns"}:
        raise ValueError("Unknown schema profile keys or missing required keys")
    if data["schema_profile_version"] != 1 or isinstance(data["schema_profile_version"], bool):
        raise ValueError("Unsupported schema_profile_version")
    pid = data["profile_id"]
    if not isinstance(pid, str) or not _IDENTIFIER.fullmatch(pid) or not pid.startswith("rxni_custom_"):
        raise ValueError("Custom profile ID must begin with rxni_custom_")
    mapping = data["columns"]
    if not isinstance(mapping, dict) or set(mapping) != set(CANONICAL_NAMES):
        raise ValueError("All ten canonical fields require an explicit source-column mapping")
    if any(not isinstance(x, str) or not _IDENTIFIER.fullmatch(x) for x in mapping.values()):
        raise ValueError("Invalid mapped source column")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Source column mapping must be one-to-one")
    extra = data["extra_columns"]
    if not isinstance(extra, dict) or len(extra) > 256:
        raise ValueError("extra_columns must be a bounded mapping")
    for name, spec in extra.items():
        if not isinstance(name, str) or not _IDENTIFIER.fullmatch(name) or name in mapping.values():
            raise ValueError(f"Invalid or ambiguous extra column {name!r}")
        if not isinstance(spec, dict) or set(spec) != {"type", "role"}:
            raise ValueError(f"Invalid extra column declaration {name!r}")
        if not isinstance(spec["type"], str) or not spec["type"] or len(spec["type"]) > 100:
            raise ValueError(f"Invalid extra column type {name!r}")
        if spec["role"] not in ("provenance", "annotation"):
            raise ValueError("Extra fields may only have provenance/annotation roles")
    return CustomV2Profile(pid, dict(mapping), dict(extra), sha256(p))


def validate_schema(profile: CustomV2Profile, schema: Any) -> dict[str, Any]:
    found = [(field.name, str(field.type)) for field in schema]
    names = [name for name, _ in found]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate Parquet source columns")
    observed = dict(found)
    allowed = set(profile.columns.values()) | set(profile.extra_columns)
    unexpected = set(observed) - allowed
    missing = allowed - set(observed)
    if unexpected or missing:
        raise ValueError(f"Physical schema has unexpected={sorted(unexpected)} missing={sorted(missing)}")
    for canonical, source in profile.columns.items():
        if observed[source] not in _STRING_TYPES:
            raise ValueError(f"Canonical field {canonical} must map to a string column, got {observed[source]}")
    for source, spec in profile.extra_columns.items():
        if observed[source] != spec["type"]:
            raise ValueError(f"Extra column {source} type mismatch: {observed[source]} != {spec['type']}")
    return {"profile": profile.profile_id, "physical_columns": len(found), "extra_columns": len(profile.extra_columns)}


def normalize_row(profile: CustomV2Profile, row: Mapping[str, Any]) -> dict[str, Any]:
    """Project only explicitly mapped canonical fields; do not coerce values."""
    result: dict[str, Any] = {}
    for canonical, source in profile.columns.items():
        if source not in row:
            raise ValueError(f"Missing source field {source}")
        value = row[source]
        if value is not None and not isinstance(value, str):
            raise TypeError(f"Expected string field {source}, found {type(value).__name__}")
        result[canonical] = value
    for required in ("record_id", "reaction_key", "reaction", "representation_generation", "representation_build_id"):
        if not result[required]:
            raise ValueError(f"Empty required field {required}")
    return result


def validate_custom_parquet(path: str | Path, profile_path: str | Path, manifest: str | Path | None = None) -> dict[str, Any]:
    import pyarrow.parquet as pq

    profile = load_custom_profile(profile_path)
    path = Path(path).expanduser().resolve()
    pf = pq.ParquetFile(path)
    summary = validate_schema(profile, pf.schema_arrow)
    requested = list(profile.columns.values())
    total = 0
    ids = set()
    for batch in pf.iter_batches(columns=requested, batch_size=8192):
        for row in batch.to_pylist():
            normalized = normalize_row(profile, row)
            if normalized["representation_generation"] != "retrieval_representations_v2":
                raise ValueError("Incorrect stored representation generation")
            ids.add(normalized["representation_build_id"])
            total += 1
    if total < 1 or total != pf.metadata.num_rows or len(ids) != 1:
        raise ValueError("Empty/incorrect row count or inconsistent build IDs")
    result = {**summary, "status": "CUSTOM_SCHEMA_READ_VALIDATED", "rows": total,
              "build_id": next(iter(ids)), "parquet_sha256": sha256(path),
              "profile_sha256": profile.sha256, "writer_is_canonical_10": True}
    if manifest is not None:
        mf = json.loads(Path(manifest).read_text(encoding="utf-8"))
        if mf.get("physical_schema_profile") != profile.profile_id:
            raise ValueError("Manifest physical_schema_profile differs from supplied custom profile")
        if mf.get("custom_schema_profile_sha256") != profile.sha256:
            raise ValueError("Manifest does not attest exact custom profile SHA256")
        for k, expected in (("build_id", result["build_id"]),
                            ("representations_sha256", result["parquet_sha256"]),
                            ("rows", total), ("row_groups", pf.num_row_groups),
                            ("representation_generation", "retrieval_representations_v2")):
            if mf.get(k) != expected:
                raise ValueError(f"Manifest metadata mismatch: {k}")
        result["manifest_checked"] = True
    else:
        result["manifest_checked"] = False
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate explicitly mapped RXNI V2 custom Parquet layout")
    parser.add_argument("parquet", type=Path)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    print(json.dumps(validate_custom_parquet(args.parquet, args.profile, args.manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
