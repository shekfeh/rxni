"""Detector interface.

Detectors PROPOSE hypotheses.  They do not validate, resolve ownership, derive
redox, or build retrieval profiles.
"""

from __future__ import annotations

from typing import Protocol

from ..structure.facts import StructuralFacts
from ..types import SemanticRecord


class SemanticDetector(Protocol):
    name: str

    def detect(self, facts: StructuralFacts) -> list[SemanticRecord]:
        ...
