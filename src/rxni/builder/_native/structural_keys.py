"""Stage 2: derive authoritative Runtime-v6.5 structural keys.

All hashes are recomputed from forward_reaction_smiles. No old derived hash is
read from the source membership file.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sqlite3
import time
from pathlib import Path
from typing import Any, Iterator

import pyarrow as pa
import pyarrow.parquet as pq

from rxni.retrieval._core.algorithms.query_structural_keys import derive_query_structural_keys

from .common import ensure_dir, require_u64, sql_to_u64, u64_to_sql, write_json
from .schema import (
    BUILDER_VERSION,
    LIBRARY_SCHEMA_VERSION,
    REACTION_CORE_BASE_FILENAME,
    REACTION_CORE_FILENAME,
    REACTION_KEY_VERSION,
    STAGING_DB,
    STRUCTURAL_HASH_FIELDS,
    STRUCTURAL_KEY_VERSION,
    STRUCTURAL_KEYS_FILENAME,
    STRUCTURAL_STATS,
)


def _build_one(task: tuple[str, str]) -> dict[str, Any]:
    key, reaction = task
    try:
        d = derive_query_structural_keys(reaction, store_readable=False)
        out: dict[str, Any] = {
            "ok": True,
            "reaction_key": key,
            "reaction": reaction,
        }
        for field in STRUCTURAL_HASH_FIELDS:
            out[field] = require_u64(d.get(field), field)
        out["forward_reactant_component_hashes"] = d.get("forward_reactant_component_hashes") or []
        out["forward_product_component_hashes"] = d.get("forward_product_component_hashes") or []
        for field in ("i_forward_raw_hash", "i_forward_edit_hash", "i_forward_component_hash"):
            out[field] = require_u64(d.get(field), field) if d.get(field) is not None else None
        return out
    except Exception as exc:
        return {
            "ok": False,
            "reaction_key": key,
            "reaction": reaction,
            "error_type": type(exc).__name__,
            "error": repr(exc),
        }


def _iter_tasks(core_base: Path, batch_size: int, counters: dict[str, int]) -> Iterator[tuple[str, str]]:
    pf = pq.ParquetFile(core_base)
    for batch in pf.iter_batches(columns=["reaction_key", "forward_reaction_smiles"], batch_size=batch_size):
        d = batch.to_pydict()
        n = len(d["reaction_key"])
        counters["core_rows"] += n
        for i in range(n):
            reaction = d["forward_reaction_smiles"][i]
            if not reaction:
                counters["missing_reaction"] += 1
                continue
            yield str(d["reaction_key"][i]), str(reaction)


def _structural_schema() -> pa.Schema:
    return pa.schema([
        ("reaction_key", pa.string()),
        ("forward_edit_signature_hash", pa.uint64()),
        ("forward_raw_center_hash_r0", pa.uint64()),
        ("forward_raw_center_hash_r1", pa.uint64()),
        ("forward_raw_center_hash_r2", pa.uint64()),
        ("forward_component_topology_hash", pa.uint64()),
        ("forward_component_evidence_hash", pa.uint64()),
        ("forward_reactant_component_hashes_json", pa.string()),
        ("forward_product_component_hashes_json", pa.string()),
        ("i_forward_raw_hash", pa.uint64()),
        ("i_forward_edit_hash", pa.uint64()),
        ("i_forward_component_hash", pa.uint64()),
        ("structural_key_version", pa.string()),
        ("library_schema_version", pa.string()),
    ])


def _core_schema() -> pa.Schema:
    return pa.schema([
        ("reaction_key", pa.string()),
        ("forward_reaction_hash", pa.string()),
        ("forward_reaction_smiles", pa.string()),
        ("forward_reaction_valid", pa.bool_()),
        ("feature_error", pa.bool_()),
        ("structural_error", pa.string()),
        ("forward_edit_signature_hash", pa.uint64()),
        ("forward_raw_center_hash_r0", pa.uint64()),
        ("forward_raw_center_hash_r1", pa.uint64()),
        ("forward_raw_center_hash_r2", pa.uint64()),
        ("forward_component_topology_hash", pa.uint64()),
        ("forward_component_evidence_hash", pa.uint64()),
        ("source_first_reaction_id", pa.string()),
        ("source_occurrence_count", pa.int64()),
        ("library_schema_version", pa.string()),
        ("reaction_key_version", pa.string()),
        ("structural_key_version", pa.string()),
        ("builder_version", pa.string()),
    ])


def _export_final_core(conn: sqlite3.Connection, path: Path, batch_size: int) -> int:
    schema = _core_schema()
    cur = conn.execute(
        """
        SELECT reaction_key, forward_reaction_hash, forward_reaction_smiles,
               forward_reaction_valid, feature_error, structural_error,
               forward_edit_signature_hash, forward_raw_center_hash_r0,
               forward_raw_center_hash_r1, forward_raw_center_hash_r2,
               forward_component_topology_hash, forward_component_evidence_hash,
               source_first_reaction_id, source_occurrence_count
        FROM reaction_core ORDER BY reaction_key
        """
    )
    writer = pq.ParquetWriter(path, schema, compression="zstd")
    total = 0
    try:
        while True:
            rows = cur.fetchmany(batch_size)
            if not rows:
                break
            payload = []
            for r in rows:
                payload.append({
                    "reaction_key": r[0],
                    "forward_reaction_hash": r[1],
                    "forward_reaction_smiles": r[2],
                    "forward_reaction_valid": bool(r[3]),
                    "feature_error": bool(r[4]),
                    "structural_error": r[5],
                    "forward_edit_signature_hash": sql_to_u64(r[6], "forward_edit_signature_hash"),
                    "forward_raw_center_hash_r0": sql_to_u64(r[7], "forward_raw_center_hash_r0"),
                    "forward_raw_center_hash_r1": sql_to_u64(r[8], "forward_raw_center_hash_r1"),
                    "forward_raw_center_hash_r2": sql_to_u64(r[9], "forward_raw_center_hash_r2"),
                    "forward_component_topology_hash": sql_to_u64(r[10], "forward_component_topology_hash"),
                    "forward_component_evidence_hash": sql_to_u64(r[11], "forward_component_evidence_hash"),
                    "source_first_reaction_id": r[12],
                    "source_occurrence_count": int(r[13]),
                    "library_schema_version": LIBRARY_SCHEMA_VERSION,
                    "reaction_key_version": REACTION_KEY_VERSION,
                    "structural_key_version": STRUCTURAL_KEY_VERSION,
                    "builder_version": BUILDER_VERSION,
                })
            writer.write_table(pa.Table.from_pylist(payload, schema=schema), row_group_size=len(payload))
            total += len(payload)
    finally:
        writer.close()
    return total


def build_structural_keys(
    output_dir: str | Path,
    *,
    core_base_path: str | Path | None = None,
    workers: int | None = None,
    batch_size: int = 8192,
    row_group_size: int = 4096,
    chunksize: int = 32,
    start_method: str = "fork",
    overwrite: bool = False,
    progress_every: int = 5000,
) -> dict[str, Any]:
    outdir = ensure_dir(output_dir)
    core_base = Path(core_base_path) if core_base_path else outdir / REACTION_CORE_BASE_FILENAME
    staging = outdir / STAGING_DB
    if not core_base.exists():
        raise FileNotFoundError(core_base)
    if not staging.exists():
        raise FileNotFoundError(staging)

    struct_path = outdir / STRUCTURAL_KEYS_FILENAME
    core_path = outdir / REACTION_CORE_FILENAME
    for p in (struct_path, core_path):
        if p.exists():
            if not overwrite:
                raise FileExistsError(p)
            p.unlink()

    workers = workers or max(1, (mp.cpu_count() or 2) - 1)
    workers = max(1, int(workers))
    conn = sqlite3.connect(staging)

    writer = pq.ParquetWriter(struct_path, _structural_schema(), compression="zstd")
    counters = {
        "core_rows": 0,
        "missing_reaction": 0,
        "processed": 0,
        "written": 0,
        "errors": 0,
    }
    t0 = time.perf_counter()
    buffer: list[dict[str, Any]] = []
    update_rows = []

    update_sql = """
        UPDATE reaction_core SET
            forward_reaction_valid=?, feature_error=?, structural_error=?,
            forward_edit_signature_hash=?, forward_raw_center_hash_r0=?,
            forward_raw_center_hash_r1=?, forward_raw_center_hash_r2=?,
            forward_component_topology_hash=?, forward_component_evidence_hash=?
        WHERE reaction_key=?
    """

    def flush() -> None:
        nonlocal buffer, update_rows
        if buffer:
            writer.write_table(pa.Table.from_pylist(buffer, schema=_structural_schema()), row_group_size=len(buffer))
            counters["written"] += len(buffer)
            buffer = []
        if update_rows:
            conn.executemany(update_sql, update_rows)
            conn.commit()
            update_rows = []

    tasks = _iter_tasks(core_base, batch_size, counters)
    pool = None
    try:
        if workers == 1:
            results = map(_build_one, tasks)
        else:
            ctx = mp.get_context(start_method)
            pool = ctx.Pool(processes=workers)
            results = pool.imap(_build_one, tasks, chunksize=chunksize)

        for item in results:
            counters["processed"] += 1
            if item["ok"]:
                payload: dict[str, Any] = {
                    "reaction_key": item["reaction_key"],
                    "forward_reactant_component_hashes_json": json.dumps(item["forward_reactant_component_hashes"], separators=(",", ":")),
                    "forward_product_component_hashes_json": json.dumps(item["forward_product_component_hashes"], separators=(",", ":")),
                    "i_forward_raw_hash": item["i_forward_raw_hash"],
                    "i_forward_edit_hash": item["i_forward_edit_hash"],
                    "i_forward_component_hash": item["i_forward_component_hash"],
                    "structural_key_version": STRUCTURAL_KEY_VERSION,
                    "library_schema_version": LIBRARY_SCHEMA_VERSION,
                }
                for f in STRUCTURAL_HASH_FIELDS:
                    payload[f] = item[f]
                buffer.append(payload)
                update_rows.append((
                    1, 0, None,
                    u64_to_sql(item["forward_edit_signature_hash"], "forward_edit_signature_hash"),
                    u64_to_sql(item["forward_raw_center_hash_r0"], "forward_raw_center_hash_r0"),
                    u64_to_sql(item["forward_raw_center_hash_r1"], "forward_raw_center_hash_r1"),
                    u64_to_sql(item["forward_raw_center_hash_r2"], "forward_raw_center_hash_r2"),
                    u64_to_sql(item["forward_component_topology_hash"], "forward_component_topology_hash"),
                    u64_to_sql(item["forward_component_evidence_hash"], "forward_component_evidence_hash"),
                    item["reaction_key"],
                ))
            else:
                counters["errors"] += 1
                update_rows.append((0, 1, f"{item['error_type']}: {item['error']}", None, None, None, None, None, None, item["reaction_key"]))

            if len(buffer) >= row_group_size or len(update_rows) >= row_group_size:
                flush()
            if progress_every > 0 and counters["processed"] % progress_every == 0:
                print(
                    f"[structural] processed={counters['processed']} written={counters['written']} "
                    f"errors={counters['errors']} elapsed={time.perf_counter()-t0:.1f}s",
                    flush=True,
                )
        flush()
    finally:
        if pool is not None:
            pool.close()
            pool.join()
        writer.close()

    core_rows = _export_final_core(conn, core_path, batch_size)
    conn.close()
    stats = {
        "stage": "build_structural_keys",
        **counters,
        "reaction_core_rows": core_rows,
        "workers": workers,
        "structural_key_version": STRUCTURAL_KEY_VERSION,
        "elapsed_seconds": time.perf_counter() - t0,
    }
    write_json(outdir / STRUCTURAL_STATS, stats)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--core-base")
    ap.add_argument("--workers", type=int)
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--row-group-size", type=int, default=4096)
    ap.add_argument("--chunksize", type=int, default=32)
    ap.add_argument("--start-method", default="fork")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    print(json.dumps(build_structural_keys(
        args.output_dir,
        core_base_path=args.core_base,
        workers=args.workers,
        batch_size=args.batch_size,
        row_group_size=args.row_group_size,
        chunksize=args.chunksize,
        start_method=args.start_method,
        overwrite=args.overwrite,
    ), indent=2))


if __name__ == "__main__":
    main()
