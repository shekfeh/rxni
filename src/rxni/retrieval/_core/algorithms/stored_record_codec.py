"""Canonical storage codec for System-I retrieval records.

The scorer consumes one canonical in-memory representation regardless of
whether the representation was generated on demand or loaded from a future
precomputed library.

The decoder intentionally accepts either:

    nested Python/Arrow objects

or:

    JSON-serialized nested fields

so the future library builder is not forced into one physical Parquet layout.
"""

from __future__ import annotations

import copy
import json

from typing import Any, Mapping


REQUIRED_RECORD_FIELDS = (
    "record_id",
    "reaction",
    "semantic_profile",
    "reaction_centers",
    "localized_fingerprints",
    "multichannel_fingerprint",
)


def _first(
    row: Mapping[str, Any],
    *names: str,
) -> Any:

    for name in names:
        if name in row:
            value = row[name]

            if value is not None:
                return value

    return None


def _decode_nested(
    value: Any,
    *,
    default: Any,
) -> Any:

    if value is None:
        return copy.deepcopy(default)

    if isinstance(
        value,
        (dict, list),
    ):
        return copy.deepcopy(value)

    if isinstance(value, bytes):
        value = value.decode("utf-8")

    if isinstance(value, str):
        text = value.strip()

        if not text:
            return copy.deepcopy(default)

        return json.loads(text)

    raise TypeError(
        "Unsupported stored nested value: "
        f"{type(value)!r}"
    )


def validate_sparse_count_fp(
    payload: Any,
    *,
    label: str = "fingerprint",
) -> None:

    if payload is None:
        return

    if not isinstance(payload, Mapping):
        raise TypeError(
            f"{label} must be a mapping"
        )

    indices = payload.get("indices") or []
    counts = payload.get("counts") or []

    if len(indices) != len(counts):
        raise ValueError(
            f"{label}: indices/counts length mismatch "
            f"{len(indices)} != {len(counts)}"
        )

    for i in indices:
        int(i)

    for c in counts:
        if float(c) < 0:
            raise ValueError(
                f"{label}: negative count {c}"
            )


def validate_retrieval_record(
    record: Mapping[str, Any],
) -> None:

    missing = [
        key
        for key in REQUIRED_RECORD_FIELDS
        if key not in record
    ]

    if missing:
        raise ValueError(
            "Retrieval record missing fields: "
            f"{missing}"
        )

    if not isinstance(
        record["semantic_profile"],
        Mapping,
    ):
        raise TypeError(
            "semantic_profile must be a mapping"
        )

    if not isinstance(
        record["reaction_centers"],
        Mapping,
    ):
        raise TypeError(
            "reaction_centers must be a mapping"
        )

    if not isinstance(
        record["localized_fingerprints"],
        Mapping,
    ):
        raise TypeError(
            "localized_fingerprints must be a mapping"
        )

    if not isinstance(
        record["multichannel_fingerprint"],
        Mapping,
    ):
        raise TypeError(
            "multichannel_fingerprint must be a mapping"
        )

    channels = (
        record[
            "multichannel_fingerprint"
        ].get("channels")
        or {}
    )

    for name, fp in channels.items():
        validate_sparse_count_fp(
            fp,
            label=f"multichannel.{name}",
        )

    centers = (
        record[
            "localized_fingerprints"
        ].get("centers")
        or []
    )

    for ci, center in enumerate(centers):
        cchannels = (
            center.get("channels")
            or {}
        )

        for name, fp in cchannels.items():
            validate_sparse_count_fp(
                fp,
                label=(
                    f"center[{ci}].{name}"
                ),
            )


def retrieval_record_from_storage_row(
    row: Mapping[str, Any],
) -> dict[str, Any]:
    """Decode a stored row into the canonical Tool-92/93 record."""

    record_id = str(
        _first(
            row,
            "record_id",
            "reaction_id",
            "forward_reaction_hash",
        )
        or ""
    )

    reaction = str(
        _first(
            row,
            "reaction",
            "forward_reaction_smiles",
            "mapped_reaction",
        )
        or ""
    )

    semantic_profile = _decode_nested(
        _first(
            row,
            "semantic_profile",
            "semantic_profile_json",
        ),
        default={},
    )

    reaction_centers = _decode_nested(
        _first(
            row,
            "reaction_centers",
            "reaction_centers_json",
        ),
        default={},
    )

    localized_fingerprints = (
        _decode_nested(
            _first(
                row,
                "localized_fingerprints",
                "localized_fingerprints_json",
            ),
            default={},
        )
    )

    multichannel_fingerprint = (
        _decode_nested(
            _first(
                row,
                "multichannel_fingerprint",
                "multichannel_fingerprint_json",
            ),
            default={},
        )
    )

    semantics = _decode_nested(
        _first(
            row,
            "semantics",
            "semantics_json",
        ),
        default={},
    )

    record = {
        "record_id":
            record_id,

        "reaction":
            reaction,

        "semantic_profile":
            semantic_profile,

        "semantics":
            semantics,

        "reaction_centers":
            reaction_centers,

        "localized_fingerprints":
            localized_fingerprints,

        "multichannel_fingerprint":
            multichannel_fingerprint,
    }

    validate_retrieval_record(
        record
    )

    return record


def retrieval_record_to_storage_row(
    record: Mapping[str, Any],
) -> dict[str, Any]:
    """Reference serialization used for parity tests.

    The future library builder may use Arrow structs/lists instead.
    The decoder accepts both.
    """

    validate_retrieval_record(
        record
    )

    def dump(value: Any) -> str:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
        )

    return {
        "record_id":
            str(
                record.get(
                    "record_id",
                    "",
                )
            ),

        "reaction":
            str(
                record.get(
                    "reaction",
                    "",
                )
            ),

        "semantic_profile_json":
            dump(
                record[
                    "semantic_profile"
                ]
            ),

        "semantics_json":
            dump(
                record.get(
                    "semantics",
                    {},
                )
            ),

        "reaction_centers_json":
            dump(
                record[
                    "reaction_centers"
                ]
            ),

        "localized_fingerprints_json":
            dump(
                record[
                    "localized_fingerprints"
                ]
            ),

        "multichannel_fingerprint_json":
            dump(
                record[
                    "multichannel_fingerprint"
                ]
            ),
    }
