"""Machine-readable directional-FGI contract registry.

Contracts describe what a validator must prove.  They are not detectors and do
not themselves inspect chemistry.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from ..version import FGI_CONTRACT_VERSION


@dataclass(frozen=True)
class FGIContract:
    label: str
    validator_class: str
    source_states: tuple[str, ...]
    target_states: tuple[str, ...]
    anchor: str
    mapping_contract: str
    one_sided_mapping_policy: str
    expected_graph_edit_contract: str
    redox_consequence: str
    status: str
    priority: str
    event_specific_guard: str = ""

    @classmethod
    def from_mapping(cls, label: str, row: Mapping[str, Any]) -> "FGIContract":
        split = lambda x: tuple(v for v in str(x or "").split(";") if v)
        return cls(
            label=label,
            validator_class=str(row.get("validator_class") or ""),
            source_states=split(row.get("source_state")),
            target_states=split(row.get("target_state")),
            anchor=str(row.get("anchor") or ""),
            mapping_contract=str(row.get("mapping_contract") or ""),
            one_sided_mapping_policy=str(row.get("one_sided_mapping_policy") or ""),
            expected_graph_edit_contract=str(row.get("expected_graph_edit_contract") or ""),
            redox_consequence=str(row.get("redox_consequence") or ""),
            status=str(row.get("status") or ""),
            priority=str(row.get("priority") or ""),
            event_specific_guard=str(row.get("event_specific_guard") or ""),
        )


class ContractRegistry:
    def __init__(self, contracts: Iterable[FGIContract]):
        rows = list(contracts)
        self._contracts = {x.label: x for x in rows}
        if len(self._contracts) != len(rows):
            raise ValueError("Duplicate FGI contract label")

    @property
    def version(self) -> str:
        return FGI_CONTRACT_VERSION

    def get(self, label: str) -> FGIContract | None:
        return self._contracts.get(label)

    def require(self, label: str) -> FGIContract:
        out = self.get(label)
        if out is None:
            raise KeyError(f"No FGI contract registered for {label!r}")
        return out

    def labels(self) -> tuple[str, ...]:
        return tuple(sorted(self._contracts))

    def __len__(self) -> int:
        return len(self._contracts)


def load_default_registry() -> ContractRegistry:
    path = Path(__file__).with_name("fgi_contracts_v1.json")
    raw = json.loads(path.read_text(encoding="utf-8"))
    return ContractRegistry(
        FGIContract.from_mapping(label, row)
        for label, row in raw.items()
    )
