"""Compact persisted codec for Runtime-v6.6 retrieval representation v2.

The physical V1 structural codec has been mechanically migrated into the v6.6
retrieval namespace. M5C therefore has no Runtime-v6.5 import edge here.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from .._core.algorithms.stored_record_codec import (
    retrieval_record_from_storage_row as _decode_native,
    retrieval_record_to_storage_row as _encode_native,
)

from ..version import STORED_REPRESENTATION_GENERATION


V2_EXTRA_COLUMNS = (
    "representation_generation",
    "representation_build_id",
)


def retrieval_record_to_storage_row_v2(
    record: dict[str, Any],
    *,
    build_id: str,
) -> dict[str, Any]:
    row = dict(_encode_native(record))
    row["representation_generation"] = STORED_REPRESENTATION_GENERATION
    row["representation_build_id"] = str(build_id)
    return row


def retrieval_record_from_storage_row_v2(
    row: dict[str, Any],
    *,
    expected_build_id: str | None = None,
) -> dict[str, Any]:
    generation = str(row.get("representation_generation") or "")
    if generation != STORED_REPRESENTATION_GENERATION:
        raise ValueError(
            "Representation generation mismatch: "
            f"expected={STORED_REPRESENTATION_GENERATION!r} "
            f"actual={generation!r}"
        )

    build_id = str(row.get("representation_build_id") or "")
    if expected_build_id is not None and build_id != str(expected_build_id):
        raise ValueError(
            "Representation build mismatch: "
            f"expected={expected_build_id!r} actual={build_id!r}"
        )

    base = dict(row)
    base.pop("representation_generation", None)
    base.pop("representation_build_id", None)

    rec = deepcopy(_decode_native(base))
    meta = dict(rec.get("metadata") or {})
    meta.update(
        {
            "representation_generation": generation,
            "representation_build_id": build_id,
            "stored_representation": True,
        }
    )
    rec["metadata"] = meta
    return rec


# Physical V2 storage order is explicit and independent of the source library.
# All ten persisted fields are strings in the validated JSON-serialized codec.
# Do not reorder them without a new storage-format contract.
V2_STORAGE_COLUMNS = (
    "record_id",
    "reaction",
    "semantic_profile_json",
    "semantics_json",
    "reaction_centers_json",
    "localized_fingerprints_json",
    "multichannel_fingerprint_json",
    "reaction_key",
    *V2_EXTRA_COLUMNS,
)


def retrieval_storage_schema():
    """Return the existing V2 Parquet schema without a legacy V1 input file.

    The schema matches the insertion order of the canonical stored-record
    encoder plus the reaction locator and generation columns. Arrow is
    imported lazily so the codec can be examined without a storage backend.
    """
    import pyarrow as pa

    return pa.schema([pa.field(name, pa.string()) for name in V2_STORAGE_COLUMNS])
