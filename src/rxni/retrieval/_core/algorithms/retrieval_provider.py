"""Retrieval-record providers for System-I.

The reranker depends only on RetrievalRecordProvider.

Current:
    OnTheFlyRecordProvider

Validation/future-storage simulation:
    StoredDictRecordProvider

Future:
    LibraryRecordProvider
"""

from __future__ import annotations

import time

from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Protocol, runtime_checkable


from .retrieval_record import (
    build_retrieval_record,
)

from .stored_record_codec import (
    retrieval_record_from_storage_row,
)


@runtime_checkable
class RetrievalRecordProvider(
    Protocol
):
    def get(
        self,
        reaction: str,
        *,
        record_id: str = "",
    ) -> dict[str, Any]:
        ...

    def prefetch(
        self,
        requests: Iterable[
            tuple[str, str]
        ],
    ) -> None:
        ...

    def snapshot(
        self,
    ) -> dict[str, Any]:
        ...


@dataclass
class OnTheFlyRecordProvider:
    """Generate System-I representations lazily."""

    radius: int = 2
    fp_size: int = 2048

    cache: dict[
        str,
        dict[str, Any],
    ] = field(
        default_factory=dict
    )

    cache_hits: int = 0
    cache_misses: int = 0
    build_seconds: float = 0.0


    def get(
        self,
        reaction: str,
        *,
        record_id: str = "",
    ) -> dict[str, Any]:

        if reaction in self.cache:
            self.cache_hits += 1
            return self.cache[reaction]

        self.cache_misses += 1

        t0 = time.perf_counter()

        record = build_retrieval_record(
            reaction,
            radius=self.radius,
            fp_size=self.fp_size,
            record_id=record_id,
        )

        self.build_seconds += (
            time.perf_counter()
            - t0
        )

        self.cache[reaction] = record

        return record


    def prefetch(
        self,
        requests: Iterable[
            tuple[str, str]
        ],
    ) -> None:
        # Deliberately no-op.
        #
        # We do not eagerly generate thousands of fingerprints here.
        # A future LibraryRecordProvider will use this hook for a
        # batched Parquet/index read.
        return None


    def snapshot(
        self,
    ) -> dict[str, Any]:

        return {
            "provider":
                "on_the_fly",

            "cache_size":
                len(self.cache),

            "cache_hits":
                self.cache_hits,

            "cache_misses":
                self.cache_misses,

            "build_seconds":
                self.build_seconds,
        }


@dataclass
class StoredDictRecordProvider:
    """Provider over already stored/serialized retrieval records.

    This acts as the test double for the future library-backed provider.
    """

    rows_by_reaction: Mapping[
        str,
        Mapping[str, Any],
    ]

    cache: dict[
        str,
        dict[str, Any],
    ] = field(
        default_factory=dict
    )

    cache_hits: int = 0
    cache_misses: int = 0
    read_seconds: float = 0.0


    def get(
        self,
        reaction: str,
        *,
        record_id: str = "",
    ) -> dict[str, Any]:

        if reaction in self.cache:
            self.cache_hits += 1
            return self.cache[reaction]

        self.cache_misses += 1

        if reaction not in self.rows_by_reaction:
            raise KeyError(
                "No stored retrieval record "
                f"for reaction: {reaction[:100]}"
            )

        t0 = time.perf_counter()

        record = (
            retrieval_record_from_storage_row(
                self.rows_by_reaction[
                    reaction
                ]
            )
        )

        self.read_seconds += (
            time.perf_counter()
            - t0
        )

        if (
            record_id
            and not record.get(
                "record_id"
            )
        ):
            record["record_id"] = (
                record_id
            )

        self.cache[reaction] = record

        return record


    def prefetch(
        self,
        requests: Iterable[
            tuple[str, str]
        ],
    ) -> None:

        # Decode requested rows now.
        #
        # The future Parquet-backed provider will replace this with
        # a real batched read.
        for reaction, record_id in requests:

            if reaction in self.cache:
                continue

            self.get(
                reaction,
                record_id=record_id,
            )


    def snapshot(
        self,
    ) -> dict[str, Any]:

        return {
            "provider":
                "stored_dict",

            "cache_size":
                len(self.cache),

            "cache_hits":
                self.cache_hits,

            "cache_misses":
                self.cache_misses,

            # Keep the common metric so reranker diagnostics remain stable.
            "build_seconds":
                0.0,

            "read_seconds":
                self.read_seconds,
        }
