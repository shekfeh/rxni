from __future__ import annotations

import json
from typing import Any, Iterable

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from rxni.public import AnalysisResult, model_dump_compat


def _short(value: Any, n: int = 96) -> str:
    if value is None:
        return "-"
    if isinstance(value, (dict, list, tuple)):
        text = json.dumps(value, ensure_ascii=False, default=str)
    else:
        text = str(value)
    text = " ".join(text.split())
    return text if len(text) <= n else text[: n - 1] + "…"


def _score(candidate: Any, ranking_mode: str) -> str:
    scores = getattr(candidate, "scores", {}) or {}
    for key in (ranking_mode, "conditioned", "sem_multi", "multi", "semantic", "local", "exact_edit"):
        if key in scores:
            try:
                return f"{float(scores[key]):.5f}"
            except Exception:
                return str(scores[key])
    return "-"


def render_summary(console: Console, result: AnalysisResult) -> None:
    s = result.summary
    r = result.runtime
    grid = Table.grid(padding=(0, 2))
    grid.add_column(style="bold")
    grid.add_column()
    rows = [
        ("Runtime", r.runtime_version),
        ("Schema", result.schema_version),
        ("Semantic policy", r.semantic_policy),
        ("Primary event", s.primary_event or "ABSTAIN"),
        ("Mode / ranking", f"{s.mode} / {s.ranking_mode}"),
        ("Candidates", str(s.selected_candidate_count if s.selected_candidate_count is not None else "-")),
        ("Scored members", str(s.scored_member_count if s.scored_member_count is not None else "-")),
        ("Templates", str(s.returned_template_count)),
        ("Elapsed", f"{s.elapsed_seconds:.3f} s" if s.elapsed_seconds is not None else "-"),
    ]
    if r.retrieval_representation_semantic_policy:
        rows.append(("Retrieval semantic policy", r.retrieval_representation_semantic_policy))
        rows.append(("Semantic lineage aligned", str(r.semantic_representation_consistent)))
    for k, v in rows:
        grid.add_row(k, v)
    console.print(Panel(grid, title="Reaction Intelligence", border_style="cyan"))


def render_semantics(console: Console, result: AnalysisResult, *, deep: bool = False) -> None:
    sem = result.semantics
    table = Table(title="Semantic interpretation", box=box.SIMPLE_HEAVY)
    table.add_column("Layer", style="bold")
    table.add_column("Registered events")

    def names(rows: Iterable[dict[str, Any]]) -> str:
        out = []
        for row in rows:
            ev = row.get("event") or row.get("subtype") or row.get("name") or row.get("label")
            st = row.get("subtype")
            if ev:
                label = str(ev)
                if st and str(st) != label:
                    label += f" [{st}]"
                out.append(label)
        return "; ".join(out) if out else "-"

    table.add_row("Primary", names(sem.primary_events) or (sem.primary_event or "-"))
    table.add_row("Specific", names(sem.specific_events))
    table.add_row("Local", names(sem.local_events))
    table.add_row("Directional / FGI", names(sem.directional_conversions))
    table.add_row("Semantic events", names(sem.semantic_events))
    table.add_row("Primitives", names(sem.primitives))
    table.add_row("Topology evidence", names(sem.topology_evidence))
    table.add_row("Redox axis", names(sem.redox_axis))
    table.add_row("Formed bond classes", ", ".join(sem.formed_bond_classes) if sem.formed_bond_classes else "-")
    console.print(table)

    if sem.semantic_profile:
        p = Table(title="Semantic retrieval profile", box=box.SIMPLE)
        p.add_column("Axis")
        p.add_column("Token")
        p.add_column("Weight", justify="right")
        axis_map = {
            "P": "primary",
            "PS": "primary subtype",
            "D": "directional",
            "R": "redox",
            "B": "bond class",
            "S": "specific",
            "L": "local",
            "PRIM": "primitive",
            "TOPO": "topology",
        }
        for token, weight in sorted(sem.semantic_profile.items()):
            prefix = str(token).split(":", 1)[0]
            p.add_row(axis_map.get(prefix, prefix), str(token), str(weight))
        console.print(p)

    if deep:
        for title, value in (
            ("Event topology", sem.event_topology),
            ("Graph summary", sem.graph_summary),
            ("Component provenance", sem.component_provenance),
            ("Validation / diagnostic fields", sem.diagnostic_fields),
        ):
            if value not in ({}, [], None):
                console.print(Panel.fit(json.dumps(value, indent=2, ensure_ascii=False, default=str), title=title))


def render_retrieval(
    console: Console,
    result: AnalysisResult,
    *,
    template_limit: int = 10,
    candidate_limit: int = 0,
) -> None:
    ret = result.retrieval
    if ret.categories:
        cats = Table(title="Candidate categories", box=box.SIMPLE)
        cats.add_column("Category")
        cats.add_column("Count", justify="right")
        cats.add_column("Fraction", justify="right")
        for row in ret.categories:
            cats.add_row(row.category, str(row.count), f"{row.fraction:.3f}" if row.fraction is not None else "-")
        console.print(cats)

    if ret.templates and template_limit:
        t = Table(title="Top templates", box=box.SIMPLE_HEAVY)
        t.add_column("#", justify="right")
        t.add_column("Score", justify="right")
        t.add_column("Members", justify="right")
        t.add_column("Template")
        for row in ret.templates[:template_limit]:
            t.add_row(
                str(row.rank),
                f"{row.score:.5f}" if row.score is not None else "-",
                str(row.member_count if row.member_count is not None else "-"),
                _short(row.template_hash, 34),
            )
        console.print(t)

    if ret.candidates and candidate_limit:
        c = Table(title="Top candidates", box=box.SIMPLE_HEAVY)
        c.add_column("#", justify="right")
        c.add_column("Score", justify="right")
        c.add_column("Category")
        c.add_column("Template")
        c.add_column("Reaction")
        for row in ret.candidates[:candidate_limit]:
            c.add_row(
                str(row.rank),
                _score(row, ret.ranking_mode),
                _short(row.candidate_category, 30),
                _short(row.template_hash, 22),
                _short(row.reaction, 82),
            )
        console.print(c)


def _method_rows(view: Any) -> list[dict[str, Any]]:
    if not getattr(view, "available", False):
        return []
    data = getattr(view, "data", {}) or {}
    rows = data.get("methods") or []
    return [dict(x) for x in rows if isinstance(x, dict)]


def render_methods(console: Console, result: AnalysisResult, *, limit: int = 12) -> None:
    for title, view in (
        ("Global method landscape", result.conditions.global_method_landscape),
        ("Query-conditioned method landscape", result.conditions.conditioned_method_landscape),
    ):
        rows = _method_rows(view)
        if not rows:
            console.print(Panel("No method landscape available.", title=title))
            continue
        table = Table(title=title, box=box.SIMPLE_HEAVY)
        table.add_column("#", justify="right")
        table.add_column("Method")
        table.add_column("Confidence", justify="right")
        table.add_column("Support", justify="right")
        table.add_column("Fraction", justify="right")
        table.add_column("Discriminator")
        for i, row in enumerate(rows[:limit], 1):
            table.add_row(
                str(i),
                str(row.get("display_name") or row.get("method_id") or row.get("name") or "-"),
                _short(row.get("confidence") or row.get("score"), 12),
                str(row.get("support_reactions") or row.get("support_count") or row.get("count") or "-"),
                _short(row.get("support_fraction") or row.get("fraction"), 12),
                _short(row.get("discriminator"), 60),
            )
        console.print(table)


def _card_parts(item: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    card = item.get("card") if isinstance(item.get("card"), dict) else item
    title = (
        card.get("title")
        or card.get("name")
        or card.get("reaction_name")
        or card.get("card_id")
        or card.get("source_id")
        or "ChemKN card"
    )
    return str(title), card


def render_chemkn(console: Console, result: AnalysisResult, *, limit: int = 8) -> None:
    ck = result.chemkn
    console.print(Panel(f"status = {ck.status}", title="ChemKN"))
    for title, rows in (("Exact method", ck.exact_cards), ("Related / family", ck.related_cards)):
        if not rows:
            continue
        console.print(Text(title, style="bold magenta"))
        for i, item in enumerate(rows[:limit], 1):
            name, card = _card_parts(item)
            relation = item.get("relation_type") if isinstance(item, dict) else None
            body = {
                "relation": relation,
                "query_condition_support": item.get("query_condition_support") if isinstance(item, dict) else None,
                "link_basis": item.get("link_basis") if isinstance(item, dict) else None,
                "card": card,
            }
            console.print(Panel(json.dumps(body, indent=2, ensure_ascii=False, default=str), title=f"{i}. {name}"))


def render_conditions(console: Console, result: AnalysisResult) -> None:
    cond = result.conditions
    table = Table(title="Query condition evidence", box=box.SIMPLE)
    table.add_column("Field", style="bold")
    table.add_column("Value")
    supplied = model_dump_compat(cond.supplied)
    for key, value in supplied.items():
        if value not in (None, [], ""):
            table.add_row(key, _short(value, 120))
    if len(table.rows) == 0:
        table.add_row("evidence", "none supplied")
    table.add_row("has_evidence", str(cond.has_evidence))
    table.add_row("reagent_weight", str(cond.reagent_weight))
    console.print(table)
    if cond.query_method_predictions:
        console.print(Panel(json.dumps(cond.query_method_predictions, indent=2, default=str), title="Query method predictions"))


def render_diagnostics(console: Console, result: AnalysisResult, *, raw_semantics: bool = False) -> None:
    d = result.diagnostics
    if d.timings:
        table = Table(title="Timings", box=box.SIMPLE)
        table.add_column("Stage")
        table.add_column("Seconds", justify="right")
        for k, v in sorted(d.timings.items()):
            table.add_row(k, f"{v:.5f}")
        console.print(table)
    if d.bulk_provider:
        console.print(Panel(json.dumps(d.bulk_provider, indent=2, default=str), title="Provider stats"))
    if d.cli_or_service_warnings:
        console.print(Panel("\n".join(f"• {x}" for x in d.cli_or_service_warnings), title="Warnings", border_style="yellow"))
    if raw_semantics:
        console.print(Panel(json.dumps(result.semantics.canonical_semantics, indent=2, ensure_ascii=False, default=str), title="Canonical semantics JSON"))
        if result.semantics.detector_outputs:
            console.print(Panel(json.dumps(result.semantics.detector_outputs, indent=2, ensure_ascii=False, default=str), title="Detector outputs"))


def render_analysis(
    console: Console,
    result: AnalysisResult,
    *,
    full: bool = False,
    template_limit: int = 10,
    candidate_limit: int = 0,
    methods: bool = False,
    chemkn: bool = False,
) -> None:
    render_summary(console, result)
    render_semantics(console, result, deep=full)
    render_retrieval(console, result, template_limit=template_limit, candidate_limit=candidate_limit)
    if full:
        render_conditions(console, result)
    if methods:
        render_methods(console, result)
    if chemkn:
        render_chemkn(console, result)
    if full:
        render_diagnostics(console, result)


def render_json(console: Console, result: AnalysisResult) -> None:
    console.print_json(json.dumps(model_dump_compat(result), ensure_ascii=False, default=str))
