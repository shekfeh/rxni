from __future__ import annotations

from typing import Any

from .compatibility import evaluate_compatibility
from .method_links import MethodLinkStore, resolve_method_links
from .reaction_card_evidence import reaction_card_summary
from .runtime_compatibility import from_retrieval_result
from rxni.method.methods import method_definition


def _gate(card_dict: dict[str, Any], registry, runtime_query):
    source = card_dict.get("source"); sid = card_dict.get("source_id")
    card = registry.get(source, sid) if source and sid else None
    if card is None:
        return {"compatible": False, "score": 0.0, "hard_failures": ["unresolved_card"], "matched": [], "missing_soft": []}
    return evaluate_compatibility(card.compatibility, runtime_query).as_dict()


def _candidate_cards(registry, method_id: str, runtime_query, *, top_k: int = 8):
    d = method_definition(method_id); seen = set(); out = []
    for q in d.get("knowledge_queries") or [d.get("display_name", method_id)]:
        for hit in registry.search(q, top_k=top_k * 2):
            card = hit["card"]; ref = card["stable_ref"]
            if ref in seen:
                continue
            seen.add(ref)
            gate = _gate(card, registry, runtime_query)
            if not gate["compatible"]:
                continue
            out.append({
                "link_status": "candidate",
                "link_basis": "knowledge_query+compatibility_gate",
                "relation_type": "transformation_family",
                "query": q,
                "search_score": hit["score"],
                "compatibility": gate,
                "card": card,
            })
    return sorted(
        out,
        key=lambda x: (-(x["search_score"] * x["compatibility"]["score"]), x["card"]["stable_ref"]),
    )[:top_k]


def _general_transformation_cards(registry, parent_event: str | None, runtime_query, *, top_k: int = 12):
    if not parent_event:
        return []
    out = []
    seen = set()
    for hit in registry.search(str(parent_event).replace("_", " "), top_k=top_k * 3):
        card = hit["card"]
        ref = card["stable_ref"]
        if ref in seen:
            continue
        seen.add(ref)
        gate = _gate(card, registry, runtime_query)
        if not gate["compatible"]:
            continue
        out.append({
            "relation_type": "transformation_family",
            "link_status": "general",
            "link_basis": "primary_transformation_query+compatibility_gate",
            "search_score": hit["score"],
            "compatibility": gate,
            "card": card,
        })
    return sorted(
        out,
        key=lambda x: (-(x["search_score"] * x["compatibility"]["score"]), x["card"]["stable_ref"]),
    )[:top_k]


def _unique_items(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    seen = set()
    for item in items:
        card = item.get("card") or {}
        key = card.get("stable_ref") or (card.get("source"), card.get("source_id"))
        if key in seen:
            continue
        seen.add(key)
        out.append(item)
    return out


def enrich_method_landscape(
    landscape: dict[str, Any],
    registry,
    *,
    retrieval_result: dict[str, Any] | None = None,
    link_store: MethodLinkStore | None = None,
    include_candidates: bool = True,
    max_candidates: int = 6,
    exact_method_threshold: float = 0.85,
) -> dict[str, Any]:
    """Attach relation-aware ChemKN evidence.

    A reviewed method link establishes that a knowledge card is an authoritative
    card for a runtime method_id.  It does *not* establish that the current query
    used that method.  Therefore an `exact_method` query relation is emitted only
    for a query-conditioned landscape with strong direct method evidence.

    Global-landscape method cards remain `transformation_family` evidence.
    """
    runtime_query = from_retrieval_result(
        retrieval_result or {"primary_event": landscape.get("parent_event")}
    )
    landscape_type = str(landscape.get("landscape_type") or "global")
    methods = []
    exact_all: list[dict[str, Any]] = []
    related_all: list[dict[str, Any]] = []

    for method in landscape.get("methods") or []:
        m = dict(method); mid = m["method_id"]
        query_support = float(m.get("query_condition_support") or 0.0)
        query_exact = landscape_type == "query_conditioned" and query_support >= exact_method_threshold
        resolved = resolve_method_links(mid, registry, link_store) if link_store is not None else []
        exact = []; related_explicit = []; rejected = []
        for r in resolved:
            if not r.get("resolved") or not r.get("card"):
                continue
            card = r["card"]; gate = _gate(card, registry, runtime_query)
            declared_relation = (r.get("link") or {}).get("relation_type") or "exact_method"
            if declared_relation == "exact_method":
                effective_relation = "exact_method" if query_exact else "transformation_family"
            else:
                effective_relation = declared_relation
            item = {
                "link_status": "explicit",
                "declared_relation_type": declared_relation,
                "relation_type": effective_relation,
                "query_exact_method_supported": query_exact,
                "query_condition_support": query_support,
                "link": r["link"],
                "compatibility": gate,
                "card": card,
            }
            if card.get("source_kind") == "reaction_card":
                source_card = registry.get(card["source"], card["source_id"])
                item["reaction_card_summary"] = reaction_card_summary(source_card)
            if not gate["compatible"]:
                rejected.append(item)
            elif effective_relation == "exact_method":
                exact.append(item)
            else:
                related_explicit.append(item)

        candidates = [] if exact or not include_candidates else _candidate_cards(
            registry, mid, runtime_query, top_k=max_candidates
        )
        compatible_explicit = exact + related_explicit
        m["chemkn_evidence"] = {
            "explicit_links": compatible_explicit,
            "exact_links": exact,
            "related_explicit_links": related_explicit,
            "rejected_explicit_links": rejected,
            "candidate_links": candidates,
        }
        m["chemkn_exact"] = exact
        m["chemkn_related"] = related_explicit + candidates
        exact_all.extend(exact)
        related_all.extend(related_explicit)
        related_all.extend(candidates)
        methods.append(m)

    general = _general_transformation_cards(
        registry, landscape.get("parent_event"), runtime_query, top_k=max(8, max_candidates)
    )
    related_all.extend(general)

    out = dict(landscape)
    out["methods"] = methods
    out["chemkn_exact"] = _unique_items(exact_all)
    out["chemkn_related"] = _unique_items(related_all)
    out["chemkn_general_transformation"] = general
    out["knowledge_sources"] = {k: len(v) for k, v in registry.by_source_id.items()}
    out["runtime_compatibility_events"] = sorted(runtime_query.events)
    out["exact_method_threshold"] = float(exact_method_threshold)
    out["relation_types"] = [
        "exact_method", "transformation_family", "mechanistically_related", "contextual"
    ]
    return out
