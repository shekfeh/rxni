#!/usr/bin/env python3
"""RXNI analysis and review command-line interface.

Provides analyze, retrieve, conditions, and run operations.
Use `python -m rxni.cli_analysis --help` for the actual option contract.
The scientific semantic and retrieval engines are unchanged by this interface."""

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

CLI_VERSION = "rxni-analysis/0.1.1"

try:
    from rxni.chemkn.knowledge_registry import KnowledgeRegistry
    from rxni.chemkn.method_links import MethodLinkStore
    from rxni.chemkn.bridge import enrich_method_landscape
except Exception:  # optional card layer
    KnowledgeRegistry = MethodLinkStore = enrich_method_landscape = None

try:
    from rxni.method.method_report import (
        render_method_landscape as _native_render_method_landscape,
    )
except Exception:  # optional human-facing native report helper
    _native_render_method_landscape = None


REVIEW_CHOICES = (
    "events",
    "fgi",
    "axes",
    "graph",
    "validation",
    "detectors",
    "methods",
    "chemkn",
    "retrieval",
    "semantics-json",
    "all",
    "audit",
)


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

    g = p.add_argument_group("deep chemistry review")
    g.add_argument(
        "--review",
        action="append",
        choices=REVIEW_CHOICES,
        default=[],
        help=(
            "Repeatable deep-review section. Examples: --review events "
            "--review fgi --review axes --review chemkn. "
            "'all' shows all curated review sections; 'audit' also dumps "
            "raw detector/semantic diagnostics."
        ),
    )
    g.add_argument(
        "--comprehensive",
        action="store_true",
        help="Shortcut for --review all while keeping the normal summary first.",
    )
    g.add_argument(
        "--audit",
        action="store_true",
        help=(
            "Maximum semantic audit: all review sections plus raw detector "
            "outputs and canonical semantics JSON."
        ),
    )
    g.add_argument(
        "--show-all-events",
        action="store_true",
        help="Shortcut for --review events.",
    )
    g.add_argument(
        "--show-fgi",
        action="store_true",
        help="Shortcut for --review fgi.",
    )
    g.add_argument(
        "--show-axes",
        action="store_true",
        help="Shortcut for --review axes.",
    )
    g.add_argument(
        "--show-validation",
        action="store_true",
        help="Shortcut for --review validation.",
    )
    g.add_argument(
        "--show-detectors",
        action="store_true",
        help="Shortcut for --review detectors.",
    )
    g.add_argument(
        "--show-method-details",
        action="store_true",
        help="Shortcut for --review methods.",
    )
    g.add_argument(
        "--show-chemkn",
        action="store_true",
        help="Shortcut for --review chemkn; automatically enables method profiling.",
    )
    g.add_argument(
        "--method-review-limit",
        type=int,
        default=0,
        help="Detailed method rows shown in fallback renderer; 0 means all.",
    )
    g.add_argument(
        "--raw-review-json",
        action="store_true",
        help="Include raw canonical semantics JSON after curated review sections.",
    )


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
            f"{i:>3} {_f(_get(row, 'score', 'weighted_support', 'weight', 'support_fraction', 'fraction', 'conditioned_score'), 5):>8} "
            f"{str(_get(row, 'support_reactions', 'support', 'support_count', 'count', 'n', default='-')):>7} "
            f"{_short(_get(row, 'method_id', 'method', 'label', 'name', default='-'), 80)}"
        )
    return "\n".join(out)


def _review_set(args) -> set[str]:
    review = set(getattr(args, "review", []) or [])
    if getattr(args, "comprehensive", False):
        review.add("all")
    if getattr(args, "audit", False):
        review.add("audit")
    if getattr(args, "show_all_events", False):
        review.add("events")
    if getattr(args, "show_fgi", False):
        review.add("fgi")
    if getattr(args, "show_axes", False):
        review.add("axes")
    if getattr(args, "show_validation", False):
        review.add("validation")
    if getattr(args, "show_detectors", False):
        review.add("detectors")
    if getattr(args, "show_method_details", False):
        review.add("methods")
    if getattr(args, "show_chemkn", False):
        review.add("chemkn")
    if getattr(args, "raw_review_json", False):
        review.add("semantics-json")

    if "audit" in review:
        review.update(
            {
                "events",
                "fgi",
                "axes",
                "graph",
                "validation",
                "detectors",
                "methods",
                "chemkn",
                "retrieval",
                "semantics-json",
            }
        )
    elif "all" in review:
        review.update(
            {
                "events",
                "fgi",
                "axes",
                "graph",
                "validation",
                "methods",
                "chemkn",
                "retrieval",
            }
        )
    return review


def _semantics(result: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    qrec = _dict(result.get("query_record"))
    sem = _dict(qrec.get("semantics"))
    return qrec, sem


def _json_text(x: Any) -> str:
    return json.dumps(x, indent=2, ensure_ascii=False, default=str)


def _event_name(row: dict[str, Any]) -> str:
    return str(
        _get(
            row,
            "event",
            "subtype",
            "name",
            "label",
            "type",
            "directional_conversion",
            default="-",
        )
    )


def _event_bucket_table(title: str, value: Any) -> str:
    rows = (
        value if isinstance(value, list) else ([] if value in (None, {}) else [value])
    )
    if not rows:
        return f"{title}: none"

    out = [
        title,
        f"{'#':>3} {'event':<34} {'subtype':<28} {'role/status':<22} {'center':<8} support/source",
    ]
    for i, raw in enumerate(rows, 1):
        row = _dict(raw)
        if not row:
            out.append(f"{i:>3} {_short(raw, 120)}")
            continue

        role_status = _get(
            row,
            "role",
            "status",
            "validation_status",
            "resolution_status",
            default="-",
        )
        center = _get(row, "center_map", "center", "map", default="-")
        support = _get(
            row,
            "support_maps",
            "support_edit_ids",
            "support",
            "source",
            "detector",
            default="-",
        )
        out.append(
            f"{i:>3} "
            f"{_short(_event_name(row), 34):<34} "
            f"{_short(row.get('subtype', '-'), 28):<28} "
            f"{_short(role_status, 22):<22} "
            f"{_short(center, 8):<8} "
            f"{_short(support, 80)}"
        )
    return "\n".join(out)


def _render_event_inventory(result: dict[str, Any]) -> str:
    qrec, sem = _semantics(result)
    if not sem:
        return "Canonical semantics object is unavailable in query_record."

    buckets = (
        ("PRIMARY EVENTS", "primary_events"),
        ("SEMANTIC EVENTS", "semantic_events"),
        ("SPECIFIC EVENTS", "specific_events"),
        ("LOCAL EVENTS", "local_events"),
        ("DIRECTIONAL CONVERSIONS / FGI", "directional_conversions"),
        ("PRIMITIVES", "primitives"),
        ("TOPOLOGY EVIDENCE", "topology_evidence"),
        ("REDOX AXIS", "redox_axis"),
    )
    blocks = [
        _kv("Semantic event count", sem.get("semantic_event_count", "-")),
        _kv("Multi semantic event", sem.get("multi_semantic_event", "-")),
    ]
    for title, key in buckets:
        blocks += ["", _event_bucket_table(title, sem.get(key))]
    return "\n".join(blocks)


def _profile_axis_name(token: str) -> str:
    prefix = token.split(":", 1)[0] if ":" in token else token
    return {
        "PS": "primary subtype",
        "P": "primary event",
        "D": "directional conversion",
        "R": "redox",
        "B": "formed bond class",
        "S": "specific event",
        "L": "local event",
        "PRIM": "primitive",
        "TOPO": "topology",
        "FAMILY": "family",
        "SUBFAMILY": "subfamily",
        "FGI": "FGI",
        "REDOX": "redox",
    }.get(prefix, prefix)


def _render_axes(result: dict[str, Any]) -> str:
    qrec, sem = _semantics(result)
    profile = qrec.get("semantic_profile") or {}
    lines = []

    lines.append("RETRIEVAL SEMANTIC PROFILE")
    if isinstance(profile, dict) and profile:
        lines.append(f"{'axis':<22} {'weight':>9} token")
        for token, weight in sorted(
            profile.items(),
            key=lambda kv: (
                _profile_axis_name(str(kv[0])),
                str(kv[0]),
            ),
        ):
            lines.append(
                f"{_profile_axis_name(str(token)):<22} " f"{_f(weight, 4):>9} {token}"
            )
    else:
        lines.append("none")

    explicit_keys = (
        "family",
        "families",
        "subfamily",
        "subfamilies",
        "fgi",
        "redox",
        "redox_axis",
        "formed_bond_classes",
        "topology_evidence",
        "event_topology",
        "component_provenance",
        "semantic_event_count",
        "multi_semantic_event",
    )
    lines += ["", "CANONICAL SEMANTIC DIMENSIONS"]
    found = False
    for key in explicit_keys:
        if key in sem and sem.get(key) not in (None, "", [], {}):
            found = True
            lines.append(f"\n[{key}]")
            lines.append(_json_text(sem.get(key)))
    if not found:
        lines.append("No explicit dimension fields found beyond semantic_profile.")
    return "\n".join(lines)


def _render_graph_review(result: dict[str, Any]) -> str:
    qrec, sem = _semantics(result)
    blocks = []
    for title, value in (
        ("STRUCTURAL KEYS", result.get("structural_keys")),
        ("MAPPING QC", result.get("mapping_qc") or sem.get("mapping_qc")),
        ("MAPPING CORRECTION", result.get("mapping_correction")),
        ("GRAPH SUMMARY", sem.get("graph_summary")),
        ("EVENT TOPOLOGY", sem.get("event_topology")),
        ("REACTION CENTERS", qrec.get("reaction_centers")),
        ("COMPONENT PROVENANCE", sem.get("component_provenance")),
        ("FORMED BOND CLASSES", sem.get("formed_bond_classes")),
    ):
        if value not in (None, "", [], {}):
            blocks += [_section(title), _json_text(value), ""]
    return (
        "\n".join(blocks).rstrip()
        if blocks
        else "No graph/topology review fields found."
    )


def _top_matching_fields(
    obj: dict[str, Any],
    tokens: tuple[str, ...],
    *,
    skip: tuple[str, ...] = (),
) -> list[tuple[str, Any]]:
    out = []
    for key, value in obj.items():
        key_l = str(key).lower()
        if any(s in key_l for s in skip):
            continue
        if any(t in key_l for t in tokens):
            out.append((str(key), value))
    return out


def _render_fgi_review(result: dict[str, Any]) -> str:
    _, sem = _semantics(result)
    blocks = [
        _event_bucket_table(
            "CANONICAL ACCEPTED DIRECTIONAL CONVERSIONS",
            sem.get("directional_conversions"),
        )
    ]

    # Stage-E/Mapping-Control audit fields evolve; show every relevant
    # top-level field rather than hard-coding one transient schema.
    fields = _top_matching_fields(
        sem,
        (
            "direction",
            "fgi",
            "hypoth",
            "validat",
            "veto",
            "reject",
            "ambigu",
            "mapping_control",
            "mapping_correction",
            "recovery",
        ),
        skip=("directional_conversions",),
    )
    if fields:
        blocks.append("")
        blocks.append("FGI / VALIDATION DIAGNOSTICS")
        for key, value in fields:
            blocks += [f"\n[{key}]", _json_text(value)]

    if result.get("mapping_correction") not in (None, {}, []):
        blocks += [
            "",
            "RETRIEVAL MAPPING-CORRECTION TRACE",
            _json_text(result.get("mapping_correction")),
        ]
    return "\n".join(blocks)


def _walk_matching(
    obj: Any,
    *,
    tokens: tuple[str, ...],
    path: str = "",
    depth: int = 0,
    max_depth: int = 5,
    skip_keys: tuple[str, ...] = ("detector_outputs",),
) -> list[tuple[str, Any]]:
    if depth > max_depth:
        return []
    out: list[tuple[str, Any]] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            key_s = str(key)
            key_l = key_s.lower()
            child = f"{path}.{key_s}" if path else key_s
            if any(x in key_l for x in skip_keys):
                continue
            if any(t in key_l for t in tokens):
                out.append((child, value))
            if isinstance(value, (dict, list)):
                out.extend(
                    _walk_matching(
                        value,
                        tokens=tokens,
                        path=child,
                        depth=depth + 1,
                        max_depth=max_depth,
                        skip_keys=skip_keys,
                    )
                )
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            child = f"{path}[{i}]"
            if isinstance(value, (dict, list)):
                out.extend(
                    _walk_matching(
                        value,
                        tokens=tokens,
                        path=child,
                        depth=depth + 1,
                        max_depth=max_depth,
                        skip_keys=skip_keys,
                    )
                )
    return out


def _render_validation_review(result: dict[str, Any]) -> str:
    _, sem = _semantics(result)
    matches = _walk_matching(
        sem,
        tokens=(
            "diagnostic",
            "validat",
            "status",
            "veto",
            "reject",
            "ambigu",
            "hypoth",
            "mapping_control",
            "recovery",
        ),
    )
    if result.get("mapping_qc") not in (None, {}, []):
        matches.insert(0, ("retrieval.mapping_qc", result.get("mapping_qc")))
    if result.get("mapping_correction") not in (None, {}, []):
        matches.insert(
            1, ("retrieval.mapping_correction", result.get("mapping_correction"))
        )

    if not matches:
        return "No explicit validation/diagnostic fields found."

    seen = set()
    lines = []
    for path, value in matches:
        fingerprint = (path, json.dumps(value, sort_keys=True, default=str))
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        lines += [f"[{path}]", _json_text(value), ""]
    return "\n".join(lines).rstrip()


def _render_detector_review(result: dict[str, Any]) -> str:
    _, sem = _semantics(result)
    outputs = sem.get("detector_outputs")
    if outputs in (None, {}, []):
        return "No detector_outputs field found."
    return _json_text(outputs)


def _render_retrieval_review(result: dict[str, Any]) -> str:
    keys = (
        "generation_counts",
        "generation_available_counts",
        "selected_channel_counts",
        "candidate_categories",
        "bulk_provider",
        "missing_provider",
        "serialized_member_limit",
        "structural_keys",
    )
    out = {}
    for key in keys:
        if result.get(key) not in (None, {}, []):
            out[key] = result.get(key)
    return _json_text(out) if out else "No retrieval diagnostic fields found."


def _method_detail_fallback(landscape: Any, limit: int = 0) -> str:
    rows = _method_rows(landscape)
    if not rows:
        return "No method landscape available."
    if limit:
        rows = rows[:limit]
    blocks = []
    for i, row in enumerate(rows, 1):
        blocks.append(
            "\n".join(
                [
                    f"[{i}] {_get(row, 'method_id', 'display_name', 'name', default='-')}",
                    _kv("Display name", _get(row, "display_name", "name", default="-")),
                    _kv("Parent event", row.get("parent_event", "-")),
                    _kv("Confidence", _get(row, "confidence", "score", default="-")),
                    _kv(
                        "Support reactions",
                        _get(row, "support_reactions", "support_count", default="-"),
                    ),
                    _kv(
                        "Support fraction",
                        _get(row, "support_fraction", "fraction", default="-"),
                    ),
                    _kv("Support occurrences", row.get("support_occurrences", "-")),
                    _kv("Weighted support", row.get("weighted_support", "-")),
                    _kv("Defining evidence", _short(row.get("defining_evidence"), 140)),
                    _kv(
                        "Supporting evidence",
                        _short(row.get("supporting_evidence"), 140),
                    ),
                    _kv("Discriminator", _short(row.get("discriminator"), 140)),
                ]
            )
        )
    return "\n\n".join(blocks)


def _render_method_review(result: dict[str, Any], args) -> str:
    if _native_render_method_landscape is not None:
        try:
            return _native_render_method_landscape(
                result,
                chemkn_scope=args.chemkn_cards,
                chemkn_card_limit=args.chemkn_card_limit,
            )
        except Exception as exc:
            fallback_warning = (
                f"[native method renderer failed: {type(exc).__name__}: {exc}]"
            )
        else:
            fallback_warning = ""
    else:
        fallback_warning = "[native method renderer unavailable]"

    blocks = [fallback_warning]
    blocks += [
        "GLOBAL METHOD LANDSCAPE",
        _method_detail_fallback(
            result.get("global_method_landscape"),
            args.method_review_limit,
        ),
    ]
    if result.get("conditioned_method_landscape"):
        blocks += [
            "",
            "QUERY-CONDITIONED METHOD LANDSCAPE",
            _method_detail_fallback(
                result.get("conditioned_method_landscape"),
                args.method_review_limit,
            ),
        ]
    return "\n".join(x for x in blocks if x)


def _chemkn_items(landscape: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    d = _dict(landscape)
    exact = [_dict(x) for x in (d.get("chemkn_exact") or [])]
    related = [_dict(x) for x in (d.get("chemkn_related") or [])]
    return exact, related


def _dedup_chemkn(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    seen = set()
    for item in items:
        card = _dict(item.get("card"))
        key = (
            item.get("relation_type"),
            card.get("card_id") or card.get("source_id") or item.get("card_id"),
            json.dumps(card or item, sort_keys=True, default=str),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(item)
    return out


def _render_chemkn_review(result: dict[str, Any], args) -> str:
    global_exact, global_related = _chemkn_items(result.get("global_method_landscape"))
    cond_exact, cond_related = _chemkn_items(result.get("conditioned_method_landscape"))
    exact = _dedup_chemkn(global_exact + cond_exact)
    related = _dedup_chemkn(global_related + cond_related)

    blocks = []
    for title, rows in (
        ("EXACT-METHOD CHEMKN CARDS", exact),
        ("RELATED / TRANSFORMATION-FAMILY CHEMKN CARDS", related),
    ):
        blocks.append(title)
        if not rows:
            blocks.append("none")
            blocks.append("")
            continue
        for i, item in enumerate(rows, 1):
            card = _dict(item.get("card"))
            display = (
                card.get("title")
                or card.get("name")
                or card.get("card_id")
                or card.get("source_id")
                or item.get("card_id")
                or "-"
            )
            blocks += [
                f"\n[{i}] {display}",
                _kv("Relation", item.get("relation_type", "-")),
                _kv("Declared relation", item.get("declared_relation_type", "-")),
                _kv("Query support", item.get("query_condition_support", "-")),
                _kv("Link basis", _short(item.get("link_basis"), 120)),
                _kv(
                    "Source",
                    _short(
                        card.get("source_id")
                        or card.get("source")
                        or card.get("chapter"),
                        120,
                    ),
                ),
                _kv("Card ID", card.get("card_id", "-")),
                _json_text(card or item),
            ]
        blocks.append("")

    if not exact and not related:
        blocks += [
            "No enriched ChemKN cards found. "
            "Check that --profile-methods is active and that Carey/March/"
            "reaction-card resources are present or passed explicitly."
        ]
    return "\n".join(blocks).rstrip()


def _render_deep_reviews(args, result: dict[str, Any]) -> str:
    review = _review_set(args)
    if not review:
        return ""

    qrec, sem = _semantics(result)
    blocks = []

    def add(title: str, body: str):
        blocks.extend(["", _section(title), body])

    if "events" in review:
        add("COMPLETE REGISTERED EVENT INVENTORY", _render_event_inventory(result))
    if "fgi" in review:
        add("FGI / DIRECTIONAL-CONVERSION REVIEW", _render_fgi_review(result))
    if "axes" in review:
        add("SEMANTIC AXES / DIMENSIONS", _render_axes(result))
    if "graph" in review:
        add("GRAPH / TOPOLOGY / CENTER REVIEW", _render_graph_review(result))
    if "validation" in review:
        add(
            "VALIDATION / ABSTENTION / MAPPING REVIEW",
            _render_validation_review(result),
        )
    if "retrieval" in review:
        add("RETRIEVAL DIAGNOSTICS", _render_retrieval_review(result))
    if "methods" in review:
        add("DETAILED METHOD LANDSCAPE", _render_method_review(result, args))
    if "chemkn" in review:
        add("CHEMKN KNOWLEDGE CARDS", _render_chemkn_review(result, args))
    if "detectors" in review:
        add("RAW DETECTOR OUTPUTS", _render_detector_review(result))
    if "semantics-json" in review:
        add("RAW CANONICAL SEMANTICS JSON", _json_text(sem or qrec))
    return "\n".join(blocks)


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

    deep = _render_deep_reviews(args, result)
    if deep:
        lines.append(deep)
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
    ap = argparse.ArgumentParser(
        prog="rxni",
        description="RXNI reaction retrieval, semantic analysis and chemistry reviews.",
        epilog=(
            "Examples: --review events --review fgi --review axes; "
            "--show-chemkn; --comprehensive; --audit. "
            "Default output remains concise."
        ),
    )
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

    requested_reviews = _review_set(args)
    if requested_reviews.intersection({"methods", "chemkn"}):
        args.profile_methods = True

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
