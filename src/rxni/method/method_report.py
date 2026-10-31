from __future__ import annotations
from typing import Any


def _pct(x):
    return f"{100.0 * float(x):5.1f}%"


def _card_ref(item: dict[str, Any]) -> tuple[Any, Any, Any]:
    c = item.get("card") or {}
    return (c.get("stable_ref"), c.get("source"), c.get("source_id"))


def _unique(items):
    out=[]; seen=set()
    for x in items or []:
        k=_card_ref(x)
        if k in seen:
            continue
        seen.add(k); out.append(x)
    return out


def _render_chemkn(lines: list[str], land: dict[str, Any], *, scope: str = "all", limit: int = 8):
    if scope == "none":
        return
    exact = _unique(land.get("chemkn_exact") or [])
    general = _unique(land.get("chemkn_related") or [])
    exact_refs={_card_ref(x) for x in exact}
    general=[x for x in general if _card_ref(x) not in exact_refs]

    if scope in {"all", "exact"}:
        lines.append("STRICTLY RELATED / EXACT-METHOD CHEMKN CARDS")
        lines.append("-" * 88)
        if not exact:
            lines.append("(none: no method card is supported as exact for this query)")
        for x in exact[:max(0,int(limit))]:
            c=x.get("card") or {}
            lines.append(
                f"{c.get('source')}:{c.get('source_id')} — {c.get('canonical_name') or c.get('title') or '-'}"
            )
            lines.append(
                f"    relation=exact_method | query_support={float(x.get('query_condition_support') or 0.0):.2f}"
            )
        lines.append("")

    if scope in {"all", "related"}:
        lines.append("GENERAL / RELATED CHEMKN KNOWLEDGE")
        lines.append("-" * 88)
        if not general:
            lines.append("(none)")
        for x in general[:max(0,int(limit))]:
            c=x.get("card") or {}
            lines.append(
                f"{c.get('source')}:{c.get('source_id')} — {c.get('canonical_name') or c.get('title') or '-'}"
            )
            lines.append(
                f"    relation={x.get('relation_type','transformation_family')} | basis={x.get('link_basis') or x.get('link_status') or '-'}"
            )
        lines.append("")


def _render_landscape(lines: list[str], land: dict[str, Any], title: str):
    ev = land.get("evidence_pool") or {}
    lines.append(title)
    lines.append("-" * 88)
    lines.append(f"Transformation : {land.get('parent_event') or 'unknown'}")
    lines.append(f"Evidence pool  : {ev.get('candidate_reactions',0)} reactions | {ev.get('membership_occurrences',0)} source occurrences")
    lines.append(f"Agent coverage : {ev.get('with_agent_field',0)}/{ev.get('candidate_reactions',0)} ({_pct(ev.get('agent_coverage_fraction',0))})")
    lines.append(f"Recognized     : {ev.get('recognized_reactions',0)}/{ev.get('candidate_reactions',0)} ({_pct(ev.get('recognized_fraction',0))})")
    lines.append(f"{'method_id':38s} {'n':>5s} {'pool':>8s} {'conf':>6s} {'q':>6s}  defining evidence")
    methods = land.get("methods") or []
    if not methods:
        lines.append("(no high-confidence method family recognized)")
    for m in methods:
        evtxt = ", ".join(m.get("defining_evidence") or [])
        lines.append(
            f"{m.get('method_id','')[:38]:38s} {m.get('support_reactions',0):5d} "
            f"{_pct(m.get('support_fraction',0)):>8s} {m.get('confidence',0):6.2f} "
            f"{m.get('query_condition_support',0):6.2f}  {evtxt}"
        )
    lines.append("")


def render_method_landscape(
    result: dict[str, Any], *, chemkn_scope: str = "all", chemkn_card_limit: int = 8
) -> str:
    global_land = result.get("global_method_landscape") or result.get("method_profile") or result.get("method_landscape") or {}
    conditioned = result.get("conditioned_method_landscape") or {}
    lines = ["REACTION METHOD LANDSCAPES", "=" * 88]
    _render_landscape(lines, global_land, "GLOBAL METHOD LANDSCAPE")
    if conditioned:
        _render_landscape(lines, conditioned, "QUERY-CONDITIONED METHOD LANDSCAPE")
        q = result.get("query_condition_evidence") or {}
        lines.append("QUERY CONDITION EVIDENCE")
        lines.append("-" * 88)
        lines.append("reagents : " + ", ".join(q.get("reagent_labels") or []))
        lines.append("solvents : " + ", ".join(q.get("solvent_labels") or []))
        lines.append("catalysts: " + ", ".join(q.get("catalyst_labels") or []))
        lines.append("bases    : " + ", ".join(q.get("base_labels") or []))
        lines.append("additives: " + ", ".join(q.get("additive_labels") or []))
        if q.get("temperature_c") is not None:
            lines.append(f"temperature: {q['temperature_c']} C")
        if q.get("reaction_time_h") is not None:
            lines.append(f"time       : {q['reaction_time_h']} h")
        lines.append("")

    # For a conditioned query the strict/general knowledge distinction should be
    # reported from the conditioned landscape; otherwise use the global one.
    knowledge_land = conditioned or global_land
    _render_chemkn(lines, knowledge_land, scope=chemkn_scope, limit=chemkn_card_limit)

    lines.append("TOP GLOBAL METHOD/CONDITION COMPONENTS")
    lines.append("-" * 88)
    for r in (global_land.get("top_reagents") or [])[:15]:
        lines.append(f"{r.get('label','')[:26]:26s} {r.get('role','')[:22]:22s} {r.get('count',0):6d}  {_pct(r.get('fraction',0))}")
    lines.append("")
    t = result.get("timings") or {}
    if t:
        lines.append(f"Runtime: {t.get('total', result.get('elapsed_seconds',0)):.3f} s | candidates={result.get('candidate_member_count',0)}")
    return "\n".join(lines)
