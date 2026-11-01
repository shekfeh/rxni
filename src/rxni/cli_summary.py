#!/usr/bin/env python3
"""RXNI compact analysis and retrieval reporting CLI.

Use `python -m rxni.cli_summary --help` for available commands.
This is a presentation layer over the frozen scientific engines."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import pydoc
import shutil
import sys
from pathlib import Path
from typing import Any, Sequence

from rxni.integration.retriever import V66Retriever

CLI_VERSION = "rxni-summary/0.1.1"

try:
    from rxni.chemkn.knowledge_registry import KnowledgeRegistry
    from rxni.chemkn.method_links import MethodLinkStore
    from rxni.chemkn.bridge import enrich_method_landscape
except Exception:  # optional card layer
    KnowledgeRegistry = MethodLinkStore = enrich_method_landscape = None


def _dict(x: Any) -> dict[str, Any]:
    if isinstance(x, dict):
        return x
    if hasattr(x, "as_dict"):
        try:
            y = x.as_dict()
            if isinstance(y, dict):
                return y
        except Exception:
            pass
    return {}


def _get(x: Any, *names: str, default=None):
    d = _dict(x)
    for name in names:
        if d.get(name) is not None:
            return d[name]
    return default


def _short(x: Any, n: int = 80) -> str:
    if x is None:
        return "-"
    if isinstance(x, (dict, list, tuple)):
        s = json.dumps(x, ensure_ascii=False, default=str)
    else:
        s = str(x)
    s = " ".join(s.split())
    return s if len(s) <= n else s[: n - 1] + "…"


def _f(x: Any, digits: int = 4) -> str:
    try:
        return f"{float(x):.{digits}f}"
    except Exception:
        return "-"


def _section(title: str) -> str:
    return f"{title}\n{'-' * len(title)}"


def _kv(label: str, value: Any) -> str:
    return f"{label:<26} {value}"


def _existing(path: str | None) -> str | None:
    return path if path and Path(path).exists() else None


def _many(values: Sequence[str] | None) -> list[str]:
    out: list[str] = []
    for raw in values or []:
        out.extend(x.strip() for x in str(raw).split(",") if x.strip())
    return out


def _library(args) -> str:
    if args.library:
        return args.library
    for env in ("REACTION_INTELLIGENCE_LIBRARY", "V5"):
        if os.environ.get(env):
            return os.environ[env]
    raise SystemExit(
        "Library not specified. Use --library PATH or set "
        "REACTION_INTELLIGENCE_LIBRARY / V5."
    )


def _reaction(args) -> str:
    if args.reaction:
        return args.reaction.strip()
    p = Path(args.reaction_file)
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            return line
    raise SystemExit(f"No reaction found in {p}")


def _write_json(path: str, result: dict[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )


def _add_runtime_args(p: argparse.ArgumentParser) -> None:
    rxn = p.add_mutually_exclusive_group(required=True)
    rxn.add_argument("--reaction", help="Reaction SMILES.")
    rxn.add_argument(
        "--reaction-file", help="Read the first reaction from a text file."
    )

    p.add_argument(
        "--library",
        help="Library-v5 path; optional if REACTION_INTELLIGENCE_LIBRARY or V5 is set.",
    )
    p.add_argument(
        "--mapping-control",
        choices=["conversion_validated", "strict"],
        help="Override RUNTIME_V67_MAPPING_CONTROL for this invocation.",
    )
    p.add_argument(
        "--mode",
        default="sem_multi",
        choices=["local", "multi", "semantic", "sem_multi", "exact_edit"],
    )
    p.add_argument("--top-k", type=int, default=20)
    p.add_argument("--exclude-template")

    g = p.add_argument_group("retrieval")
    g.add_argument("--candidate-budget", type=int, default=3000)
    g.add_argument("--per-channel-limit", type=int, default=5000)
    g.add_argument(
        "--candidate-policy", choices=["specificity", "legacy"], default="specificity"
    )
    g.add_argument(
        "--member-limit",
        type=int,
        help="Members serialized to result JSON; 0 means all.",
    )

    g = p.add_argument_group("query conditions")
    g.add_argument(
        "--reagent", "--reagents", dest="reagents", action="append", default=[]
    )
    g.add_argument("--solvent", action="append", default=[])
    g.add_argument("--catalyst", action="append", default=[])
    g.add_argument("--base", action="append", default=[])
    g.add_argument("--additive", action="append", default=[])
    g.add_argument("--acid", action="append", default=[])
    g.add_argument("--temperature-c", type=float)
    g.add_argument("--reaction-time-h", type=float)
    g.add_argument("--atmosphere")
    g.add_argument("--pressure")
    g.add_argument("--reagent-weight", type=float, default=0.20)

    g = p.add_argument_group("methods / ChemKN")
    g.add_argument("--profile-methods", action="store_true")
    g.add_argument("--method-display-limit", type=int, default=8)
    g.add_argument(
        "--chemkn-cards", choices=["all", "exact", "related", "none"], default="all"
    )
    g.add_argument("--chemkn-card-limit", type=int, default=8)
    g.add_argument("--knowledge-cards", help="Explicit normalized JSONL knowledge collection")
    g.add_argument("--method-links")

    g = p.add_argument_group("reporting")
    g.add_argument("--format", choices=["text", "json", "csv", "tsv"], default="text")
    g.add_argument(
        "--summary", choices=["compact", "standard", "full", "none"], default="standard"
    )
    g.add_argument("--output", help="Write complete result JSON.")
    g.add_argument("--report-output", help="Write the rendered report.")
    g.add_argument("--show-templates", action="store_true")
    g.add_argument("--template-display-limit", type=int, default=10)
    g.add_argument("--show-candidates", action="store_true")
    g.add_argument("--candidate-format", choices=["table", "detail"], default="table")
    g.add_argument("--candidate-display-limit", type=int, default=20)
    g.add_argument("--full-reactions", action="store_true")
    g.add_argument("--pager", choices=["auto", "never", "always"], default="auto")


def _enrich(result: dict[str, Any], args) -> dict[str, Any]:
    if (
        not args.profile_methods
        or args.chemkn_cards == "none"
        or not result.get("method_profile")
    ):
        return result
    if KnowledgeRegistry is None:
        result.setdefault("_cli_warnings", []).append(
            "Native ChemKN enrichment modules unavailable."
        )
        return result

    cards = _existing(args.knowledge_cards or os.environ.get("RXNI_KNOWLEDGE_CARDS"))
    if not cards:
        return result
    reg = KnowledgeRegistry.from_jsonl(cards)
    links = (
        MethodLinkStore.load(args.method_links)
        if _existing(args.method_links)
        else None
    )
    related = args.chemkn_cards in {"all", "related"}
    limit = max(0, int(args.chemkn_card_limit))

    if result.get("global_method_landscape"):
        result["global_method_landscape"] = enrich_method_landscape(
            result["global_method_landscape"],
            reg,
            retrieval_result=result,
            link_store=links,
            include_candidates=related,
            max_candidates=limit,
        )
        result["method_profile"] = result["global_method_landscape"]
    if result.get("conditioned_method_landscape"):
        result["conditioned_method_landscape"] = enrich_method_landscape(
            result["conditioned_method_landscape"],
            reg,
            retrieval_result=result,
            link_store=links,
            include_candidates=related,
            max_candidates=limit,
        )
    return result


def _score(row: dict[str, Any], mode: str) -> Any:
    scores = row.get("scores")
    if isinstance(scores, dict):
        for key in (mode, "conditioned", "sem_multi", "multi", "semantic", "local"):
            if scores.get(key) is not None:
                return scores[key]
    for key in ("score", "template_score", "aggregate_score", "final_score"):
        if row.get(key) is not None:
            return row[key]
    return None


def _template_table(result: dict[str, Any], limit: int) -> str:
    rows = result.get("templates") or []
    if not rows:
        return "No ranked templates."
    mode = str(result.get("ranking_mode") or result.get("mode") or "sem_multi")
    rows = rows if limit == 0 else rows[:limit]
    out = [f"{'#':>3} {'score':>8} {'members':>7} {'template':<22} event"]
    for i, raw in enumerate(rows, 1):
        row = _dict(raw)
        out.append(
            f"{i:>3} {_f(_score(row, mode), 5):>8} "
            f"{str(_get(row, 'member_count', 'n_members', 'support', default='-')):>7} "
            f"{_short(_get(row, 'template_hash', 'retro_template_hash', default='-'), 22):<22} "
            f"{_short(_get(row, 'primary_event', 'dominant_primary_event', 'event', default='-'), 40)}"
        )
    return "\n".join(out)


def _candidate_table(result: dict[str, Any], args) -> str:
    rows = result.get("members") or []
    if not rows:
        return "No serialized candidates."
    rows = (
        rows
        if args.candidate_display_limit == 0
        else rows[: args.candidate_display_limit]
    )
    mode = str(result.get("ranking_mode") or result.get("mode") or args.mode)

    if args.candidate_format == "detail":
        blocks = []
        for i, raw in enumerate(rows, 1):
            row, md = _dict(raw), _dict(_dict(raw).get("metadata"))
            reaction = _get(row, "reaction", "forward_reaction_smiles", default="-")
            if not args.full_reactions:
                reaction = _short(reaction, 150)
            blocks.append(
                "\n".join(
                    [
                        f"[{i}] {row.get('reaction_key', '-')}",
                        _kv("Score", _f(_score(row, mode), 6)),
                        _kv(
                            "Template",
                            _get(
                                row, "template_hash", "retro_template_hash", default="-"
                            ),
                        ),
                        _kv("Category", md.get("candidate_category", "-")),
                        _kv("Chemistry relation", md.get("chemistry_relation", "-")),
                        _kv("Condition relation", md.get("condition_relation", "-")),
                        _kv("Reaction", reaction),
                    ]
                )
            )
        return "\n\n".join(blocks)

    out = [f"{'#':>3} {'score':>8} {'category':<22} {'template':<18} reaction"]
    for i, raw in enumerate(rows, 1):
        row, md = _dict(raw), _dict(_dict(raw).get("metadata"))
        reaction = _get(row, "reaction", "forward_reaction_smiles", default="-")
        if not args.full_reactions:
            reaction = _short(reaction, 85)
        out.append(
            f"{i:>3} {_f(_score(row, mode), 5):>8} "
            f"{_short(md.get('candidate_category', '-'), 22):<22} "
            f"{_short(_get(row, 'template_hash', 'retro_template_hash', default='-'), 18):<18} "
            f"{reaction}"
        )
    return "\n".join(out)


def _method_rows(x: Any) -> list[dict[str, Any]]:
    if isinstance(x, list):
        return [_dict(y) for y in x]
    d = _dict(x)
    for key in ("methods", "ranked_methods", "entries", "rows", "landscape"):
        if isinstance(d.get(key), list):
            return [_dict(y) for y in d[key]]
    return [d] if d else []


def _method_table(x: Any, limit: int) -> str:
    rows = _method_rows(x)
    if not rows:
        return "No method landscape available."
    rows = rows if limit == 0 else rows[:limit]
    out = [f"{'#':>3} {'score':>8} {'support':>7} method"]
    for i, row in enumerate(rows, 1):
        out.append(
            f"{i:>3} {_f(_get(row, 'score', 'weight', 'fraction', 'conditioned_score'), 5):>8} "
            f"{str(_get(row, 'support', 'support_count', 'count', 'n', default='-')):>7} "
            f"{_short(_get(row, 'method_id', 'method', 'label', 'name', default='-'), 80)}"
        )
    return "\n".join(out)


def _render_text(args, result: dict[str, Any]) -> str:
    if args.summary == "none":
        return json.dumps(result, indent=2, ensure_ascii=False, default=str)

    compact = args.summary == "compact"
    full = args.summary == "full"
    mapping = os.environ.get(
        "RUNTIME_V67_MAPPING_CONTROL",
        args.mapping_control or "conversion_validated",
    )
    lines = [
        _section("REACTION INTELLIGENCE"),
        _kv("Runtime", result.get("runtime_version", "6.7")),
        _kv("Mapping policy", mapping),
        _kv("Primary event", result.get("primary_event", "-")),
        _kv(
            "Mode / ranking",
            f"{result.get('mode', args.mode)} / {result.get('ranking_mode', '-') }",
        ),
        _kv("Candidates", result.get("selected_candidate_keys", "-")),
        _kv("Scored members", result.get("candidate_member_count", "-")),
        _kv("Templates", len(result.get("templates") or [])),
        _kv("Elapsed", f"{_f(result.get('elapsed_seconds'), 3)} s"),
    ]
    if compact:
        return "\n".join(lines)

    qrec = _dict(result.get("query_record"))
    profile = qrec.get("semantic_profile")
    lines += ["", _section("REACTION / SEMANTICS")]
    lines.append(_kv("Query reaction key", result.get("query_reaction_key", "-")))
    if isinstance(profile, dict) and profile:
        keys = sorted(
            profile,
            key=lambda k: (
                -float(profile[k]) if isinstance(profile[k], (int, float)) else 0,
                str(k),
            ),
        )
        lines.append(_kv("Semantic profile", ", ".join(map(str, keys[:12]))))
    if result.get("mapping_qc") not in (None, {}, []):
        lines.append(_kv("Mapping QC", _short(result.get("mapping_qc"), 110)))
    if result.get("mapping_correction") not in (None, {}, []):
        lines.append(
            _kv("Mapping correction", _short(result.get("mapping_correction"), 110))
        )

    lines += ["", _section("RETRIEVAL")]
    lines.append(
        _kv("Candidate policy", result.get("candidate_policy", args.candidate_policy))
    )
    lines.append(
        _kv("Candidate budget", result.get("candidate_budget", args.candidate_budget))
    )
    if result.get("selected_channel_counts") or result.get("generation_counts"):
        lines.append(
            _kv(
                "Channels",
                _short(
                    result.get("selected_channel_counts")
                    or result.get("generation_counts"),
                    110,
                ),
            )
        )
    if result.get("candidate_categories"):
        cats = [
            f"{_dict(x).get('category','-')}:{_dict(x).get('count','-')}"
            for x in result["candidate_categories"][:8]
        ]
        lines.append(_kv("Candidate categories", "; ".join(cats)))

    if args.show_templates or full:
        lines += [
            "",
            _section("TOP TEMPLATES"),
            _template_table(result, args.template_display_limit),
        ]

    evidence = _dict(result.get("query_condition_evidence"))
    if evidence and evidence.get("has_evidence"):
        lines += ["", _section("QUERY CONDITIONS")]
        for label, key in (
            ("Active components", "active_components"),
            ("Solvents", "solvents"),
            ("Catalysts", "catalysts"),
            ("Bases", "bases"),
            ("Temperature", "temperature_c"),
            ("Time", "reaction_time_h"),
        ):
            if evidence.get(key) not in (None, "", [], {}):
                lines.append(_kv(label, _short(evidence[key], 100)))
    if result.get("query_method_predictions"):
        lines.append(
            _kv("Method predictions", _short(result["query_method_predictions"], 120))
        )

    if args.profile_methods:
        lines += [
            "",
            _section("GLOBAL METHOD LANDSCAPE"),
            _method_table(
                result.get("global_method_landscape"), args.method_display_limit
            ),
        ]
        if result.get("conditioned_method_landscape"):
            lines += [
                "",
                _section("QUERY-CONDITIONED METHOD LANDSCAPE"),
                _method_table(
                    result.get("conditioned_method_landscape"),
                    args.method_display_limit,
                ),
            ]

    if args.show_candidates:
        lines += ["", _section("TOP CANDIDATES"), _candidate_table(result, args)]

    if full and result.get("timings"):
        lines += ["", _section("TIMINGS")]
        lines += [_kv(k, f"{_f(v, 4)} s") for k, v in result["timings"].items()]
    if full and result.get("_cli_warnings"):
        lines += ["", _section("CLI WARNINGS")]
        lines += [f"- {x}" for x in result["_cli_warnings"]]
    return "\n".join(lines)


def _render_delimited(args, result: dict[str, Any]) -> str:
    rows = result.get("members") or []
    rows = (
        rows
        if args.candidate_display_limit == 0
        else rows[: args.candidate_display_limit]
    )
    mode = str(result.get("ranking_mode") or result.get("mode") or args.mode)
    buf = io.StringIO()
    w = csv.writer(
        buf, delimiter="," if args.format == "csv" else "\t", lineterminator="\n"
    )
    w.writerow(
        ["rank", "reaction_key", "score", "template_hash", "category", "reaction"]
    )
    for i, raw in enumerate(rows, 1):
        row, md = _dict(raw), _dict(_dict(raw).get("metadata"))
        w.writerow(
            [
                i,
                row.get("reaction_key"),
                _score(row, mode),
                _get(row, "template_hash", "retro_template_hash"),
                md.get("candidate_category"),
                _get(row, "reaction", "forward_reaction_smiles"),
            ]
        )
    return buf.getvalue()


def _render(args, result: dict[str, Any]) -> str:
    if args.format == "json":
        return json.dumps(result, indent=2, ensure_ascii=False, default=str)
    if args.format in {"csv", "tsv"}:
        return _render_delimited(args, result)
    return _render_text(args, result)


def _page_or_print(text: str, mode: str) -> None:
    if mode == "always":
        pydoc.pager(text)
        return
    if mode == "auto" and sys.stdout.isatty():
        height = shutil.get_terminal_size((120, 35)).lines
        if text.count("\n") + 1 > max(20, height - 3):
            pydoc.pager(text)
            return
    print(text)


def _run(args) -> dict[str, Any]:
    if args.mapping_control:
        os.environ["RUNTIME_V67_MAPPING_CONTROL"] = args.mapping_control

    member_limit = args.member_limit
    if member_limit is None and (args.show_candidates or args.format in {"csv", "tsv"}):
        member_limit = args.candidate_display_limit

    rt = V66Retriever(_library(args))
    result = rt.retrieve(
        _reaction(args),
        mode=args.mode,
        top_k=args.top_k,
        exclude_template=args.exclude_template,
        max_members=args.candidate_budget,
        per_channel_limit=args.per_channel_limit,
        candidate_policy=args.candidate_policy,
        member_limit=member_limit,
        reagents=_many(args.reagents),
        solvents=_many(args.solvent),
        catalysts=_many(args.catalyst),
        bases=_many(args.base),
        additives=_many(args.additive),
        acids=_many(args.acid),
        temperature_c=args.temperature_c,
        reaction_time_h=args.reaction_time_h,
        atmosphere=args.atmosphere,
        pressure=args.pressure,
        reagent_weight=args.reagent_weight,
        profile_methods=args.profile_methods,
    )
    if not isinstance(result, dict):
        result = result.as_dict()
    result.setdefault(
        "_cli",
        {
            "cli_version": CLI_VERSION,
            "command": args.cmd,
            "mapping_control": os.environ.get(
                "RUNTIME_V67_MAPPING_CONTROL", "runtime-default"
            ),
        },
    )
    return _enrich(result, args)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="rxni")
    ap.add_argument("--version", action="version", version=CLI_VERSION)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name, help_text in (
        ("retrieve", "RXNI structural and semantic retrieval."),
        ("analyze", "Retrieval with human-facing summary."),
        ("conditions", "Method/condition-focused view."),
        ("run", "Backward-friendly analysis alias."),
    ):
        p = sp.add_parser(name, help=help_text)
        _add_runtime_args(p)
    return ap


def main(argv: Sequence[str] | None = None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)

    if args.cmd == "conditions":
        args.profile_methods = True
        if args.summary == "standard":
            args.summary = "full"
    if args.cmd in {"analyze", "run"} and args.format == "text":
        args.show_templates = True

    result = _run(args)
    if args.output:
        _write_json(args.output, result)

    rendered = _render(args, result)
    if args.report_output:
        p = Path(args.report_output)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(
            rendered if rendered.endswith("\n") else rendered + "\n", encoding="utf-8"
        )

    if args.format in {"csv", "tsv"}:
        print(rendered, end="" if rendered.endswith("\n") else "\n")
    elif args.format == "json":
        print(rendered)
    else:
        _page_or_print(rendered, args.pager)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
