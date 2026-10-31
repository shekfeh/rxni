#!/usr/bin/env python3
"""RXNI batch orchestration for user-selected single-reaction scripts.

Usage: `python -m rxni.batch_cli --input reactions.tsv --output-dir results --cli-script /path/to/cli_analysis.py --cli-args "analyze --format json"`.
Pass --review-script only when a compatible independent review CLI exists.
No reaction inference logic is implemented in this orchestration module."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------


def _guess_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    if path.suffix.lower() == ".csv":
        return ","
    sample = path.read_text(encoding="utf-8", errors="replace")[:8192]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error:
        return "\t"


def _iter_table(
    path: Path,
    reaction_column: str,
    case_id_column: str | None,
) -> Iterable[dict[str, Any]]:
    delim = _guess_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames:
            raise SystemExit(f"No header found in {path}")
        if reaction_column not in reader.fieldnames:
            raise SystemExit(
                f"Missing reaction column {reaction_column!r}; "
                f"available={reader.fieldnames}"
            )
        for i, row in enumerate(reader):
            reaction = str(row.get(reaction_column) or "").strip()
            if not reaction:
                continue
            case_id = (
                str(row.get(case_id_column) or "").strip()
                if case_id_column and case_id_column in row
                else ""
            )
            yield {
                "row_index": i,
                "case_id": case_id or f"row_{i:07d}",
                "reaction": reaction,
                "input": row,
            }


def _iter_jsonl(
    path: Path,
    reaction_column: str,
    case_id_column: str | None,
) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            obj = json.loads(line)
            if not isinstance(obj, dict):
                raise SystemExit(f"JSONL row {i+1} is not an object")
            reaction = str(obj.get(reaction_column) or "").strip()
            if not reaction:
                continue
            case_id = (
                str(obj.get(case_id_column) or "").strip() if case_id_column else ""
            )
            yield {
                "row_index": i,
                "case_id": case_id or f"row_{i:07d}",
                "reaction": reaction,
                "input": obj,
            }


def _iter_text(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as fh:
        n = 0
        for i, line in enumerate(fh):
            reaction = line.strip()
            if not reaction or reaction.startswith("#"):
                continue
            yield {
                "row_index": i,
                "case_id": f"row_{n:07d}",
                "reaction": reaction,
                "input": {"reaction": reaction},
            }
            n += 1


def load_rows(
    path: Path,
    reaction_column: str,
    case_id_column: str | None,
) -> list[dict[str, Any]]:
    ext = path.suffix.lower()
    if ext in {".jsonl", ".ndjson"}:
        rows = list(_iter_jsonl(path, reaction_column, case_id_column))
    elif ext in {".txt", ".smi", ".smiles"}:
        rows = list(_iter_text(path))
    else:
        rows = list(_iter_table(path, reaction_column, case_id_column))
    if not rows:
        raise SystemExit(f"No non-empty reactions found in {path}")
    return rows


# ---------------------------------------------------------------------------
# CLI execution
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ToolSpec:
    name: str
    script: Path
    extra_args: tuple[str, ...]
    reaction_flag: str


def _try_json(text: str) -> Any | None:
    s = text.strip()
    if not s:
        return None

    # Exact JSON output.
    try:
        return json.loads(s)
    except Exception:
        pass

    # Some CLIs print one or more informational lines before JSON.
    for pos, ch in enumerate(s):
        if ch not in "[{":
            continue
        try:
            return json.loads(s[pos:])
        except Exception:
            continue
    return None


def _records_from(value: Any, key: str) -> list[dict[str, Any]]:
    if not isinstance(value, dict):
        return []
    obj = value.get(key)
    if isinstance(obj, list):
        return [x for x in obj if isinstance(x, dict)]
    return []


def _semantics_root(parsed: Any) -> dict[str, Any]:
    if not isinstance(parsed, dict):
        return {}

    # Common wrappers used by the project over time.
    for key in ("semantics", "canonical_semantics", "result"):
        x = parsed.get(key)
        if isinstance(x, dict):
            if key == "result":
                for subkey in ("semantics", "canonical_semantics"):
                    y = x.get(subkey)
                    if isinstance(y, dict):
                        return y
            else:
                return x

    return parsed


def _event_names(records: list[dict[str, Any]]) -> list[str]:
    ans: list[str] = []
    for rec in records:
        name = rec.get("event")
        if name is not None and str(name).strip():
            ans.append(str(name))
    return ans


def _profile_labels(sem: dict[str, Any]) -> list[str]:
    p = sem.get("semantic_profile")
    if isinstance(p, dict):
        return sorted(str(k) for k in p.keys())
    if isinstance(p, list):
        return [str(x) for x in p]
    return []


def compact_fields(parsed: Any) -> dict[str, Any]:
    """Best-effort extraction only; full parsed JSON remains preserved."""
    sem = _semantics_root(parsed)

    return {
        "primary_events": _event_names(_records_from(sem, "primary_events")),
        "specific_events": _event_names(_records_from(sem, "specific_events")),
        "local_events": _event_names(_records_from(sem, "local_events")),
        "directional_conversions": _event_names(
            _records_from(sem, "directional_conversions")
        ),
        "redox_axis": _event_names(_records_from(sem, "redox_axis")),
        "formed_bond_classes": (
            sem.get("formed_bond_classes")
            if isinstance(sem.get("formed_bond_classes"), list)
            else []
        ),
        "semantic_profile_labels": _profile_labels(sem),
    }


def run_one(
    tool: ToolSpec,
    row: dict[str, Any],
    python_exe: str,
    timeout_s: float | None,
    env: dict[str, str],
) -> dict[str, Any]:
    cmd = [
        python_exe,
        str(tool.script),
        *tool.extra_args,
        tool.reaction_flag,
        row["reaction"],
    ]

    t0 = time.perf_counter()
    try:
        cp = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            env=env,
            shell=False,
        )
        elapsed = time.perf_counter() - t0
        stdout = cp.stdout
        stderr = cp.stderr
        parsed = _try_json(stdout)
        compact = compact_fields(parsed) if parsed is not None else {}
        return {
            "record_type": "reaction_cli_batch_result",
            "tool": tool.name,
            "case_id": row["case_id"],
            "row_index": row["row_index"],
            "reaction": row["reaction"],
            "input": row["input"],
            "ok": cp.returncode == 0,
            "exit_code": cp.returncode,
            "timed_out": False,
            "elapsed_s": round(elapsed, 6),
            "command_argv": cmd,
            "stdout": stdout,
            "stderr": stderr,
            "parsed_json": parsed,
            "compact": compact,
        }

    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - t0
        return {
            "record_type": "reaction_cli_batch_result",
            "tool": tool.name,
            "case_id": row["case_id"],
            "row_index": row["row_index"],
            "reaction": row["reaction"],
            "input": row["input"],
            "ok": False,
            "exit_code": None,
            "timed_out": True,
            "elapsed_s": round(elapsed, 6),
            "command_argv": cmd,
            "stdout": exc.stdout or "",
            "stderr": exc.stderr or "",
            "parsed_json": None,
            "compact": {},
            "error": f"timeout after {timeout_s}s",
        }

    except Exception as exc:
        elapsed = time.perf_counter() - t0
        return {
            "record_type": "reaction_cli_batch_result",
            "tool": tool.name,
            "case_id": row["case_id"],
            "row_index": row["row_index"],
            "reaction": row["reaction"],
            "input": row["input"],
            "ok": False,
            "exit_code": None,
            "timed_out": False,
            "elapsed_s": round(elapsed, 6),
            "command_argv": cmd,
            "stdout": "",
            "stderr": "",
            "parsed_json": None,
            "compact": {},
            "error": f"{type(exc).__name__}: {exc}",
        }


# ---------------------------------------------------------------------------
# Output / resume
# ---------------------------------------------------------------------------


def _resume_keys(path: Path) -> set[tuple[str, str]]:
    done: set[tuple[str, str]] = set()
    if not path.exists():
        return done
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue
            tool = str(obj.get("tool") or "")
            case_id = str(obj.get("case_id") or "")
            if tool and case_id:
                done.add((tool, case_id))
    return done


def _join(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, list):
        return ";".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v, ensure_ascii=False, sort_keys=True)
    return str(v)


def append_summary(writer: csv.DictWriter, result: dict[str, Any]) -> None:
    c = result.get("compact") or {}
    writer.writerow(
        {
            "tool": result.get("tool", ""),
            "case_id": result.get("case_id", ""),
            "row_index": result.get("row_index", ""),
            "ok": result.get("ok", False),
            "exit_code": result.get("exit_code", ""),
            "timed_out": result.get("timed_out", False),
            "elapsed_s": result.get("elapsed_s", ""),
            "primary_events": _join(c.get("primary_events")),
            "specific_events": _join(c.get("specific_events")),
            "local_events": _join(c.get("local_events")),
            "directional_conversions": _join(c.get("directional_conversions")),
            "redox_axis": _join(c.get("redox_axis")),
            "formed_bond_classes": _join(c.get("formed_bond_classes")),
            "semantic_profile_labels": _join(c.get("semantic_profile_labels")),
            "reaction": result.get("reaction", ""),
            "error": result.get("error", ""),
            "stderr_excerpt": str(result.get("stderr") or "").replace("\n", " ")[:500],
        }
    )


SUMMARY_FIELDS = [
    "tool",
    "case_id",
    "row_index",
    "ok",
    "exit_code",
    "timed_out",
    "elapsed_s",
    "primary_events",
    "specific_events",
    "local_events",
    "directional_conversions",
    "redox_axis",
    "formed_bond_classes",
    "semantic_profile_labels",
    "reaction",
    "error",
    "stderr_excerpt",
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Batch-run explicitly selected RXNI analysis/review scripts over a reaction list "
            "without changing their single-reaction chemistry."
        )
    )

    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--reaction-column", default="reaction")
    ap.add_argument("--case-id-column", default="case_id")

    ap.add_argument("--cli-script", type=Path)
    ap.add_argument(
        "--cli-args",
        default="",
        help='Arguments placed before --reaction, e.g. "analyze --format json"',
    )

    ap.add_argument("--review-script", type=Path)
    ap.add_argument(
        "--review-args",
        default="",
        help='Arguments placed before --reaction, e.g. "--all"',
    )

    ap.add_argument(
        "--reaction-flag",
        default="--reaction",
        help="Reaction option used by the target CLIs (default: --reaction)",
    )
    ap.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used for target scripts (default: current interpreter)",
    )
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=0.0)
    ap.add_argument("--start-row", type=int, default=0)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--fail-fast", action="store_true")
    ap.add_argument(
        "--env",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Additional environment variable; repeatable",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()

    if not args.cli_script and not args.review_script:
        raise SystemExit(
            "Configure at least one target: --cli-script and/or --review-script"
        )
    if args.workers < 1:
        raise SystemExit("--workers must be >= 1")

    tools: list[ToolSpec] = []
    if args.cli_script:
        if not args.cli_script.is_file():
            raise SystemExit(f"CLI script not found: {args.cli_script}")
        tools.append(
            ToolSpec(
                name="analysis",
                script=args.cli_script.resolve(),
                extra_args=tuple(shlex.split(args.cli_args)),
                reaction_flag=args.reaction_flag,
            )
        )
    if args.review_script:
        if not args.review_script.is_file():
            raise SystemExit(f"Review CLI script not found: {args.review_script}")
        tools.append(
            ToolSpec(
                name="review",
                script=args.review_script.resolve(),
                extra_args=tuple(shlex.split(args.review_args)),
                reaction_flag=args.reaction_flag,
            )
        )

    rows = load_rows(
        args.input,
        args.reaction_column,
        args.case_id_column or None,
    )
    if args.start_row:
        rows = rows[args.start_row :]
    if args.max_rows:
        rows = rows[: args.max_rows]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = args.output_dir / "batch_results.jsonl"
    tsv_path = args.output_dir / "batch_summary.tsv"
    run_path = args.output_dir / "run_manifest.json"

    done = _resume_keys(jsonl_path) if args.resume else set()

    jobs: list[tuple[ToolSpec, dict[str, Any]]] = []
    for row in rows:
        for tool in tools:
            if (tool.name, row["case_id"]) not in done:
                jobs.append((tool, row))

    env = os.environ.copy()
    for item in args.env:
        if "=" not in item:
            raise SystemExit(f"Bad --env value {item!r}; expected KEY=VALUE")
        k, v = item.split("=", 1)
        env[k] = v

    manifest = {
        "input": str(args.input.resolve()),
        "reaction_column": args.reaction_column,
        "case_id_column": args.case_id_column,
        "rows_selected": len(rows),
        "tools": [
            {
                "name": t.name,
                "script": str(t.script),
                "extra_args": list(t.extra_args),
                "reaction_flag": t.reaction_flag,
            }
            for t in tools
        ],
        "python": args.python,
        "workers": args.workers,
        "timeout_s": args.timeout or None,
        "resume": bool(args.resume),
        "jobs_to_run": len(jobs),
    }
    run_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    if not jobs:
        print("Nothing to run.")
        print(f"JSONL: {jsonl_path}")
        return 0

    tsv_exists = tsv_path.exists() and args.resume
    jsonl_mode = "a" if args.resume else "w"
    tsv_mode = "a" if tsv_exists else "w"

    ok_count = 0
    fail_count = 0
    t_start = time.perf_counter()

    with (
        jsonl_path.open(jsonl_mode, encoding="utf-8") as jout,
        tsv_path.open(tsv_mode, encoding="utf-8", newline="") as tout,
    ):

        writer = csv.DictWriter(tout, delimiter="\t", fieldnames=SUMMARY_FIELDS)
        if not tsv_exists:
            writer.writeheader()

        def consume(result: dict[str, Any]) -> None:
            nonlocal ok_count, fail_count
            jout.write(json.dumps(result, ensure_ascii=False) + "\n")
            jout.flush()
            append_summary(writer, result)
            tout.flush()

            if result.get("ok"):
                ok_count += 1
            else:
                fail_count += 1

            n = ok_count + fail_count
            print(
                f"[{n:>6}/{len(jobs)}] "
                f"tool={result.get('tool')} "
                f"case={result.get('case_id')} "
                f"ok={result.get('ok')} "
                f"t={result.get('elapsed_s')}s"
            )

        if args.workers == 1:
            for tool, row in jobs:
                result = run_one(
                    tool,
                    row,
                    args.python,
                    args.timeout or None,
                    env,
                )
                consume(result)
                if args.fail_fast and not result.get("ok"):
                    break
        else:
            # Thread workers orchestrate independent subprocesses.
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
                futs = [
                    ex.submit(
                        run_one,
                        tool,
                        row,
                        args.python,
                        args.timeout or None,
                        env,
                    )
                    for tool, row in jobs
                ]
                for fut in as_completed(futs):
                    result = fut.result()
                    consume(result)
                    if args.fail_fast and not result.get("ok"):
                        for f in futs:
                            f.cancel()
                        break

    elapsed = time.perf_counter() - t_start
    final = {
        **manifest,
        "completed_jobs": ok_count + fail_count,
        "ok": ok_count,
        "failed": fail_count,
        "elapsed_s": round(elapsed, 6),
        "results_jsonl": str(jsonl_path),
        "summary_tsv": str(tsv_path),
    }
    run_path.write_text(json.dumps(final, indent=2) + "\n", encoding="utf-8")

    print()
    print(json.dumps(final, indent=2))
    return 0 if fail_count == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
