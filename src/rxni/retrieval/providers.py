"""Runtime-v6.6 retrieval-record providers.

No module in this file imports Runtime-v6.5.  The old V1 representation remains
readable through a mechanically migrated provider/codec implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

from .record import build_retrieval_record
from .version import LIVE_RECORD_VERSION, V5_COMPAT_REPRESENTATION


class RetrievalRecordProvider(Protocol):
    representation_generation: str

    def get(
        self,
        reaction: str,
        *,
        record_id: str | None = None,
    ) -> dict[str, Any]:
        ...


@dataclass
class LiveV66RecordProvider:
    representation_generation: str = LIVE_RECORD_VERSION

    def get(
        self,
        reaction: str,
        *,
        record_id: str | None = None,
    ) -> dict[str, Any]:
        return build_retrieval_record(
            reaction,
            record_id=record_id,
        )


class V5CompatibilityProvider:
    """Read-only adapter around immutable V1 representations."""

    representation_generation = V5_COMPAT_REPRESENTATION

    def __init__(self, library: str | Path):
        from ._core.algorithms.library_record_provider import (
            LibraryRecordProvider,
        )

        library = Path(library)
        self.library = library
        self._provider = LibraryRecordProvider(
            library / "system_i_representations_v1.parquet",
            library / "representation_locator_v1.sqlite",
            fallback_provider=None,
        )

    def get(
        self,
        reaction: str,
        *,
        record_id: str | None = None,
    ) -> dict[str, Any]:
        return self._provider.get(
            reaction,
            record_id=record_id,
        )
