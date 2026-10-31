"""Strict RXNI Retrieval-V2 Parquet physical-schema profiles.

The public writer uses an explicit, compact 10-field schema. A historical
V6.7 writer inherited four extra V1 fields (14 total). Reading the historical
layout is supported only after validating that its extra fields are null.

No scientific normalization, codec rewrite, or transformation is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

CANONICAL_SCHEMA_PROFILE = "rxni_retrieval_v2_compact10_v1"
LEGACY_SCHEMA_PROFILE = "rxni_retrieval_v2_inherited14_v1"

CANONICAL_COLUMNS = (
    ("record_id", "string"),
    ("reaction", "string"),
    ("semantic_profile_json", "string"),
    ("semantics_json", "string"),
    ("reaction_centers_json", "string"),
    ("localized_fingerprints_json", "string"),
    ("multichannel_fingerprint_json", "string"),
    ("reaction_key", "string"),
    ("representation_generation", "string"),
    ("representation_build_id", "string"),
)
LEGACY_COLUMNS = (
    ("reaction_key", "string"),
    ("forward_edit_signature_hash", "uint64"),
    ("record_id", "string"),
    ("reaction", "string"),
    ("semantic_profile_json", "string"),
    ("semantics_json", "string"),
    ("reaction_centers_json", "string"),
    ("localized_fingerprints_json", "string"),
    ("multichannel_fingerprint_json", "string"),
    ("library_schema_version", "string"),
    ("reaction_key_version", "string"),
    ("system_i_representation_version", "string"),
    ("representation_generation", "string"),
    ("representation_build_id", "string"),
)
LEGACY_ONLY_COLUMNS = (
    "forward_edit_signature_hash", "library_schema_version",
    "reaction_key_version", "system_i_representation_version",
)


def classify_physical_schema(schema: Any) -> str:
    """Match exact field names, Arrow types, and order; reject unknown layouts."""
    found = tuple((field.name, str(field.type)) for field in schema)
    if found == CANONICAL_COLUMNS:
        return CANONICAL_SCHEMA_PROFILE
    if found == LEGACY_COLUMNS:
        return LEGACY_SCHEMA_PROFILE
    raise ValueError(f"Unsupported V2 physical schema: {found!r}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_parquet(path: str | Path, manifest: str | Path | None = None) -> dict[str, Any]:
    """Read-only validation: schema, legacy-null guard, generation and manifest.

    This is a storage compatibility check, NOT chemical correctness or ranking
    validation. It does not silently project or drop legacy data.
    """
    import pyarrow.parquet as pq

    path = Path(path).expanduser().resolve()
    pf = pq.ParquetFile(path)
    profile = classify_physical_schema(pf.schema_arrow)
    extra_nonnull = {name: 0 for name in LEGACY_ONLY_COLUMNS}
    cols = ["record_id", "reaction_key", "representation_generation", "representation_build_id"]
    if profile == LEGACY_SCHEMA_PROFILE:
        cols.extend(LEGACY_ONLY_COLUMNS)
    count = 0
    ids = set()
    generations = set()
    for batch in pf.iter_batches(columns=cols, batch_size=8192):
        # These are small metadata-only columns, not expensive semantics JSON.
        values = batch.to_pydict()
        n = batch.num_rows
        count += n
        for field in ("reaction_key", "record_id", "representation_generation", "representation_build_id"):
            if any(not isinstance(v, str) or not v for v in values[field]):
                raise ValueError(f"Empty or invalid required V2 field: {field}")
        ids.update(values["representation_build_id"])
        generations.update(values["representation_generation"])
        for name in extra_nonnull:
            if name in values:
                extra_nonnull[name] += sum(v is not None for v in values[name])
    if profile == LEGACY_SCHEMA_PROFILE and any(extra_nonnull.values()):
        raise ValueError(f"Non-null legacy-only fields require explicit migration audit: {extra_nonnull}")
    if len(ids) != 1 or len(generations) != 1 or generations != {"retrieval_representations_v2"}:
        raise ValueError(f"Inconsistent representation generation or build IDs: {ids}, {generations}")
    if count != pf.metadata.num_rows:
        raise ValueError("Parquet row-count mismatch")
    result: dict[str, Any] = {
        "status": "SCHEMA_PROFILE_VALIDATED",
        "profile": profile,
        "canonical_for_new_writes": profile == CANONICAL_SCHEMA_PROFILE,
        "rows": count,
        "row_groups": pf.num_row_groups,
        "build_id": next(iter(ids)),
        "generation": next(iter(generations)),
        "legacy_only_nonnull_counts": extra_nonnull if profile == LEGACY_SCHEMA_PROFILE else {},
        "sha256": sha256(path),
    }
    if manifest is not None:
        source = json.loads(Path(manifest).read_text(encoding="utf-8"))
        for field, value in (
            ("representations_sha256", result["sha256"]),
            ("rows", count),
            ("row_groups", pf.num_row_groups),
            ("build_id", result["build_id"]),
            ("representation_generation", result["generation"]),
        ):
            if source.get(field) != value:
                raise ValueError(f"Manifest mismatch for {field}: {source.get(field)!r} != {value!r}")
        declared = source.get("physical_schema_profile")
        if declared is not None and declared != profile:
            raise ValueError(f"Manifest declares profile {declared!r} but file is {profile!r}")
        result["manifest_checked"] = True
    else:
        result["manifest_checked"] = False
    return result


def main():
    ap = argparse.ArgumentParser(description="Validate known RXNI V2 Parquet physical-schema profiles (read-only)")
    ap.add_argument("parquet", type=Path)
    ap.add_argument("--manifest", type=Path)
    args = ap.parse_args()
    print(json.dumps(validate_parquet(args.parquet, args.manifest), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
