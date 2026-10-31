"""Opt-in registry for normalized knowledge cards.

RXNI v0.1 does not auto-discover or bundle publisher-specific knowledge bases.
External adapters should normalize licensed sources to the documented JSONL format.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from ._io import read_jsonl
from .compatibility import evaluate_compatibility
from .knowledge_schema import KnowledgeCard, CompatibilitySignature, RuntimeCompatibilityQuery

_TOKEN = re.compile(r"[A-Za-z0-9_+\-]+")
_CARD_FIELDS = set(KnowledgeCard.__dataclass_fields__)
_COMPAT_FIELDS = set(CompatibilitySignature.__dataclass_fields__)
_REQUIRED = ("source", "source_id", "canonical_name", "source_kind")


def _tokens(text: str) -> set[str]:
    return {x.lower() for x in _TOKEN.findall(text or "") if len(x) > 1}


def normalize_card(row: dict, *, line: int | None = None) -> KnowledgeCard:
    """Validate a normalized record without executing arbitrary provider code."""
    if not isinstance(row, dict):
        raise ValueError(f"knowledge record {line}: expected JSON object")
    missing = [k for k in _REQUIRED if not isinstance(row.get(k), str) or not row[k].strip()]
    if missing:
        raise ValueError(f"knowledge record {line}: missing/empty required fields {missing}")
    unknown = set(row) - _CARD_FIELDS - {"stable_ref"}
    if unknown:
        raise ValueError(f"knowledge record {line}: unknown fields {sorted(unknown)}")
    record = {k: v for k, v in row.items() if k != "stable_ref"}
    signature = record.get("compatibility")
    if signature is None:
        record["compatibility"] = CompatibilitySignature()
    elif isinstance(signature, dict):
        unknown_sig = set(signature) - _COMPAT_FIELDS
        if unknown_sig:
            raise ValueError(f"knowledge record {line}: unknown compatibility fields {sorted(unknown_sig)}")
        record["compatibility"] = CompatibilitySignature(**signature)
    elif not isinstance(signature, CompatibilitySignature):
        raise ValueError(f"knowledge record {line}: compatibility must be an object")
    card = KnowledgeCard(**record)
    if row.get("stable_ref") and row["stable_ref"] != card.stable_ref:
        raise ValueError(f"knowledge record {line}: inconsistent stable_ref")
    return card


class KnowledgeRegistry:
    def __init__(self, cards: Iterable[KnowledgeCard] = ()):
        self.cards = list(cards)
        self.by_ref = {}
        self.by_source_id = defaultdict(dict)
        for card in self.cards:
            if card.stable_ref in self.by_ref:
                raise ValueError(f"Duplicate knowledge card {card.stable_ref}")
            self.by_ref[card.stable_ref] = card
            self.by_source_id[card.source][card.source_id] = card

    @classmethod
    def from_jsonl(cls, path: str | Path) -> "KnowledgeRegistry":
        return cls(normalize_card(row, line=i) for i, row in enumerate(read_jsonl(path), 1))

    @classmethod
    def from_paths(cls, *, knowledge_cards: str | Path | None = None) -> "KnowledgeRegistry":
        return cls.from_jsonl(knowledge_cards) if knowledge_cards is not None else cls()

    def get(self, source: str, source_id: str) -> KnowledgeCard | None:
        return self.by_source_id.get(source, {}).get(str(source_id))

    def search(self, text: str, *, top_k: int = 20, sources: set[str] | None = None):
        q = _tokens(text)
        scored = []
        for card in self.cards:
            if sources is not None and card.source not in sources:
                continue
            hay = " ".join([
                card.canonical_name, card.transformation or "", card.reaction_family or "",
                " ".join(card.aliases), " ".join(card.retrieval_tags)
            ])
            overlap = q & _tokens(hay)
            if overlap:
                scored.append((len(overlap) / max(1, len(q)), card))
        scored.sort(key=lambda x: (-x[0], x[1].stable_ref))
        return [{"score": score, "card": card.as_dict()} for score, card in scored[:top_k]]

    def compatible_cards(
        self, query: RuntimeCompatibilityQuery, *, sources: set[str] | None = None,
        top_k: int = 100,
    ):
        out = []
        for card in self.cards:
            if sources is not None and card.source not in sources:
                continue
            decision = evaluate_compatibility(card.compatibility, query)
            if decision.compatible:
                out.append({
                    "compatibility_score": decision.score,
                    "decision": decision.as_dict(),
                    "card": card.as_dict(),
                })
        out.sort(key=lambda x: (-x["compatibility_score"], x["card"]["stable_ref"]))
        return out[:top_k]
