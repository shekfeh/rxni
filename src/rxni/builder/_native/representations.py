"""Stage 3: precompute canonical Runtime-v6.5 System-I retrieval records."""

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

from rxni.retrieval.record import build_retrieval_record
from rxni.retrieval._core.algorithms.stored_record_codec import retrieval_record_to_storage_row

from .common import ensure_dir, u64_to_sql, write_json
from .schema import (
    LIBRARY_SCHEMA_VERSION,
    LOCATOR_VERSION,
    REACTION_CORE_FILENAME,
    REACTION_KEY_VERSION,
    REPRESENTATION_ERRORS_FILENAME,
    REPRESENTATION_LOCATOR_DB,
    REPRESENTATION_LOCATOR_PARQUET,
    REPRESENTATION_STATS,
    REPRESENTATIONS_FILENAME,
    SYSTEM_I_REPRESENTATION_VERSION,
)


def _build_one(task: tuple[str, str, int | None, int, int]) -> dict[str, Any]:
    key, reaction, edit_hash, radius, fp_size = task
    try:
        record = build_retrieval_record(reaction, radius=radius, fp_size=fp_size, record_id=key)
        stored = retrieval_record_to_storage_row(record)
        return {
            "ok": True,
            "reaction_key": key,
            "forward_edit_signature_hash": edit_hash,
            **stored,
        }
    except Exception as exc:
        return {
            "ok": False,
            "reaction_key": key,
            "reaction": reaction,
            "forward_edit_signature_hash": edit_hash,
            "error_type": type(exc).__name__,
            "error": repr(exc),
        }


def _rep_schema() -> pa.Schema:
    return pa.schema([
        ("reaction_key", pa.string()),
        ("forward_edit_signature_hash", pa.uint64()),
        ("record_id", pa.string()),
        ("reaction", pa.string()),
        ("semantic_profile_json", pa.string()),
        ("semantics_json", pa.string()),
        ("reaction_centers_json", pa.string()),
        ("localized_fingerprints_json", pa.string()),
        ("multichannel_fingerprint_json", pa.string()),
        ("library_schema_version", pa.string()),
        ("reaction_key_version", pa.string()),
        ("system_i_representation_version", pa.string()),
    ])


def _locator_schema() -> pa.Schema:
    return pa.schema([
        ("reaction_key", pa.string()),
        ("row_group_id", pa.int32()),
        ("row_in_group", pa.int32()),
        ("forward_edit_signature_hash", pa.uint64()),
        ("locator_version", pa.string()),
    ])


def _iter_tasks(core: Path, batch_size: int, radius: int, fp_size: int, counters: dict[str, int]) -> Iterator[tuple[str, str, int | None, int, int]]:
    cols = ["reaction_key", "forward_reaction_smiles", "forward_reaction_valid", "feature_error", "forward_edit_signature_hash"]
    for batch in pq.ParquetFile(core).iter_batches(columns=cols, batch_size=batch_size):
        d = batch.to_pydict()
        n = len(d["reaction_key"])
        counters["core_rows"] += n
        for i in range(n):
            if not d["forward_reaction_valid"][i] or d["feature_error"][i]:
                counters["skipped_invalid_or_error"] += 1
                continue
            reaction = d["forward_reaction_smiles"][i]
            if not reaction:
                counters["skipped_missing_reaction"] += 1
                continue
            counters["eligible"] += 1
            yield (
                str(d["reaction_key"][i]),
                str(reaction),
                int(d["forward_edit_signature_hash"][i]) if d["forward_edit_signature_hash"][i] is not None else None,
                radius,
                fp_size,
            )


def build_system_i_representations(
    output_dir: str | Path,
    *,
    core_path: str | Path | None = None,
    workers: int | None = None,
    batch_size: int = 8192,
    row_group_size: int = 4096,
    chunksize: int = 32,
    radius: int = 2,
    fp_size: int = 2048,
    start_method: str = "fork",
    overwrite: bool = False,
    progress_every: int = 5000,
) -> dict[str, Any]:
    outdir = ensure_dir(output_dir)
    core = Path(core_path) if core_path else outdir / REACTION_CORE_FILENAME
    if not core.exists():
        raise FileNotFoundError(core)

    rep_path = outdir / REPRESENTATIONS_FILENAME
    locator_db = outdir / REPRESENTATION_LOCATOR_DB
    locator_pq = outdir / REPRESENTATION_LOCATOR_PARQUET
    errors = outdir / REPRESENTATION_ERRORS_FILENAME
    for p in (rep_path, locator_db, locator_pq, errors):
        if p.exists():
            if not overwrite:
                raise FileExistsError(p)
            p.unlink()

    workers = workers or max(1, (mp.cpu_count() or 2) - 1)
    workers = max(1, int(workers))

    conn = sqlite3.connect(locator_db)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.executescript(
        """
        CREATE TABLE locator (
            reaction_key TEXT PRIMARY KEY,
            row_group_id INTEGER NOT NULL,
            row_in_group INTEGER NOT NULL,
            forward_edit_signature_hash TEXT
        );
        CREATE INDEX idx_locator_row_group ON locator(row_group_id);
        CREATE INDEX idx_locator_edit_hash ON locator(forward_edit_signature_hash);
        """
    )

    rep_writer = pq.ParquetWriter(rep_path, _rep_schema(), compression="zstd")
    loc_writer = pq.ParquetWriter(locator_pq, _locator_schema(), compression="zstd")
    err_fh = errors.open("w", encoding="utf-8")
    counters = {"core_rows": 0, "eligible": 0, "skipped_invalid_or_error": 0, "skipped_missing_reaction": 0, "representations_written": 0, "errors": 0}
    row_group_id = 0
    buffer: list[dict[str, Any]] = []
    t0 = time.perf_counter()

    def flush() -> None:
        nonlocal row_group_id, buffer
        if not buffer:
            return
        payload = []
        loc_payload = []
        loc_sql = []
        for row_in_group, item in enumerate(buffer):
            payload.append({
                "reaction_key": item["reaction_key"],
                "forward_edit_signature_hash": item["forward_edit_signature_hash"],
                "record_id": item["record_id"],
                "reaction": item["reaction"],
                "semantic_profile_json": item["semantic_profile_json"],
                "semantics_json": item["semantics_json"],
                "reaction_centers_json": item["reaction_centers_json"],
                "localized_fingerprints_json": item["localized_fingerprints_json"],
                "multichannel_fingerprint_json": item["multichannel_fingerprint_json"],
                "library_schema_version": LIBRARY_SCHEMA_VERSION,
                "reaction_key_version": REACTION_KEY_VERSION,
                "system_i_representation_version": SYSTEM_I_REPRESENTATION_VERSION,
            })
            loc_payload.append({
                "reaction_key": item["reaction_key"],
                "row_group_id": row_group_id,
                "row_in_group": row_in_group,
                "forward_edit_signature_hash": item["forward_edit_signature_hash"],
                "locator_version": LOCATOR_VERSION,
            })
            loc_sql.append((item["reaction_key"], row_group_id, row_in_group, u64_to_sql(item["forward_edit_signature_hash"], "forward_edit_signature_hash")))
        rep_writer.write_table(pa.Table.from_pylist(payload, schema=_rep_schema()), row_group_size=len(payload))
        loc_writer.write_table(pa.Table.from_pylist(loc_payload, schema=_locator_schema()), row_group_size=len(loc_payload))
        conn.executemany("INSERT INTO locator VALUES (?, ?, ?, ?)", loc_sql)
        conn.commit()
        counters["representations_written"] += len(payload)
        row_group_id += 1
        buffer = []

    tasks = _iter_tasks(core, batch_size, radius, fp_size, counters)
    pool = None
    try:
        if workers == 1:
            results = map(_build_one, tasks)
        else:
            ctx = mp.get_context(start_method)
            pool = ctx.Pool(processes=workers)
            results = pool.imap(_build_one, tasks, chunksize=chunksize)
        processed = 0
        for item in results:
            processed += 1
            if item["ok"]:
                buffer.append(item)
                if len(buffer) >= row_group_size:
                    flush()
            else:
                counters["errors"] += 1
                err_fh.write(json.dumps(item, sort_keys=True) + "\n")
            if progress_every > 0 and processed % progress_every == 0:
                print(f"[representations] processed={processed} written={counters['representations_written']} errors={counters['errors']} elapsed={time.perf_counter()-t0:.1f}s", flush=True)
        flush()
    finally:
        if pool is not None:
            pool.close()
            pool.join()
        rep_writer.close()
        loc_writer.close()
        err_fh.close()
        conn.close()

    stats = {
        "stage": "build_system_i_representations",
        **counters,
        "workers": workers,
        "radius": radius,
        "fp_size": fp_size,
        "system_i_representation_version": SYSTEM_I_REPRESENTATION_VERSION,
        "elapsed_seconds": time.perf_counter() - t0,
    }
    write_json(outdir / REPRESENTATION_STATS, stats)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--core")
    ap.add_argument("--workers", type=int)
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--row-group-size", type=int, default=4096)
    ap.add_argument("--chunksize", type=int, default=32)
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--fp-size", type=int, default=2048)
    ap.add_argument("--start-method", default="fork")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    print(json.dumps(build_system_i_representations(
        args.output_dir, core_path=args.core, workers=args.workers, batch_size=args.batch_size,
        row_group_size=args.row_group_size, chunksize=args.chunksize, radius=args.radius,
        fp_size=args.fp_size, start_method=args.start_method, overwrite=args.overwrite,
    ), indent=2))


if __name__ == "__main__":
    main()
