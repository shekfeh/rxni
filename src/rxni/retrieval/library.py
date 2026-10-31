"""Library-v5 data foundation contract for Runtime-v6.6 retrieval.

M4 decision:
KEEP immutable core/membership/provenance/structural-index artifacts.
REBUILD only stored retrieval representations under a new generation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class LibraryLayout:
    root: Path

    @classmethod
    def from_path(cls, path: str | Path) -> "LibraryLayout":
        return cls(Path(path))

    @property
    def reaction_core(self) -> Path:
        return self.root / "reaction_core_v5.parquet"

    @property
    def membership(self) -> Path:
        return self.root / "reaction_template_membership_v5.parquet"

    @property
    def provenance(self) -> Path:
        return self.root / "reaction_provenance_v5.parquet"

    @property
    def template_library(self) -> Path:
        return self.root / "template_library_v5.parquet"

    @property
    def structural_index(self) -> Path:
        return self.root / "structural_index_v1.sqlite"

    @property
    def v5_representations(self) -> Path:
        return self.root / "system_i_representations_v1.parquet"

    @property
    def v5_locator(self) -> Path:
        return self.root / "representation_locator_v1.sqlite"

    @property
    def v66_representations(self) -> Path:
        return self.root / "retrieval_representations_v2.parquet"

    @property
    def v66_locator(self) -> Path:
        return self.root / "representation_locator_v2.sqlite"

    def validate_foundation(self) -> None:
        required = [
            self.reaction_core,
            self.membership,
            self.provenance,
            self.structural_index,
            self.v5_representations,
            self.v5_locator,
        ]
        missing = [str(p) for p in required if not p.exists()]
        if missing:
            raise FileNotFoundError(
                "Library foundation is incomplete:\n"
                + "\n".join(missing)
            )
