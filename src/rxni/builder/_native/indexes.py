"""Stage 6: build exact-edit and structural lookup indexes from native V5 core."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .common import ensure_dir, u64_to_sql, write_json
from .schema import (
    EDIT_INDEX_DB,
    EDIT_INDEX_FILENAME,
    EDIT_INDEX_STATS,
    EDIT_INDEX_VERSION,
    REACTION_CORE_FILENAME,
    STRUCTURAL_INDEX_DB,
    STRUCTURAL_INDEX_VERSION,
)


def build_indexes(
    output_dir: str | Path,
    *,
    core_path: str | Path | None = None,
    batch_size: int = 8192,
    overwrite: bool = False,
) -> dict[str, Any]:
    outdir = ensure_dir(output_dir)
    core = Path(core_path) if core_path else outdir / REACTION_CORE_FILENAME
    edit_pq = outdir / EDIT_INDEX_FILENAME
    edit_db = outdir / EDIT_INDEX_DB
    structural_db = outdir / STRUCTURAL_INDEX_DB
    if not core.exists():
        raise FileNotFoundError(core)
    for p in (edit_pq, edit_db, structural_db):
        if p.exists():
            if not overwrite:
                raise FileExistsError(p)
            p.unlink()

    schema = pa.schema([
        ("forward_edit_signature_hash", pa.uint64()),
        ("reaction_key", pa.string()),
        ("forward_raw_center_hash_r0", pa.uint64()),
        ("forward_raw_center_hash_r1", pa.uint64()),
        ("forward_raw_center_hash_r2", pa.uint64()),
        ("edit_index_version", pa.string()),
    ])
    writer = pq.ParquetWriter(edit_pq, schema, compression="zstd")

    econn = sqlite3.connect(edit_db)
    econn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE edit_index (
            forward_edit_signature_hash TEXT NOT NULL,
            reaction_key TEXT NOT NULL,
            PRIMARY KEY (forward_edit_signature_hash, reaction_key)
        );
        CREATE INDEX idx_edit_hash ON edit_index(forward_edit_signature_hash);
        """
    )
    sconn = sqlite3.connect(structural_db)
    sconn.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE structural_index (
            reaction_key TEXT PRIMARY KEY,
            forward_edit_signature_hash TEXT,
            forward_raw_center_hash_r0 TEXT,
            forward_raw_center_hash_r1 TEXT,
            forward_raw_center_hash_r2 TEXT,
            forward_component_topology_hash TEXT,
            forward_component_evidence_hash TEXT
        );
        CREATE INDEX idx_struct_edit ON structural_index(forward_edit_signature_hash);
        CREATE INDEX idx_struct_r0 ON structural_index(forward_raw_center_hash_r0);
        CREATE INDEX idx_struct_r1 ON structural_index(forward_raw_center_hash_r1);
        CREATE INDEX idx_struct_r2 ON structural_index(forward_raw_center_hash_r2);
        CREATE INDEX idx_struct_ct ON structural_index(forward_component_topology_hash);
        CREATE INDEX idx_struct_ce ON structural_index(forward_component_evidence_hash);
        """
    )

    cols = [
        "reaction_key", "forward_reaction_valid", "feature_error",
        "forward_edit_signature_hash", "forward_raw_center_hash_r0",
        "forward_raw_center_hash_r1", "forward_raw_center_hash_r2",
        "forward_component_topology_hash", "forward_component_evidence_hash",
    ]
    total = 0
    indexed = 0
    missing_edit = 0
    try:
        for batch in pq.ParquetFile(core).iter_batches(columns=cols, batch_size=batch_size):
            d = batch.to_pydict()
            payload = []
            edit_rows = []
            struct_rows = []
            for i in range(batch.num_rows):
                total += 1
                if not d["forward_reaction_valid"][i] or d["feature_error"][i]:
                    continue
                key = str(d["reaction_key"][i])
                edit = d["forward_edit_signature_hash"][i]
                struct_rows.append((
                    key,
                    u64_to_sql(edit, "forward_edit_signature_hash"),
                    u64_to_sql(d["forward_raw_center_hash_r0"][i], "forward_raw_center_hash_r0"),
                    u64_to_sql(d["forward_raw_center_hash_r1"][i], "forward_raw_center_hash_r1"),
                    u64_to_sql(d["forward_raw_center_hash_r2"][i], "forward_raw_center_hash_r2"),
                    u64_to_sql(d["forward_component_topology_hash"][i], "forward_component_topology_hash"),
                    u64_to_sql(d["forward_component_evidence_hash"][i], "forward_component_evidence_hash"),
                ))
                if edit is None:
                    missing_edit += 1
                    continue
                payload.append({
                    "forward_edit_signature_hash": int(edit),
                    "reaction_key": key,
                    "forward_raw_center_hash_r0": d["forward_raw_center_hash_r0"][i],
                    "forward_raw_center_hash_r1": d["forward_raw_center_hash_r1"][i],
                    "forward_raw_center_hash_r2": d["forward_raw_center_hash_r2"][i],
                    "edit_index_version": EDIT_INDEX_VERSION,
                })
                edit_rows.append((u64_to_sql(edit, "forward_edit_signature_hash"), key))
                indexed += 1
            if payload:
                writer.write_table(pa.Table.from_pylist(payload, schema=schema), row_group_size=len(payload))
            if edit_rows:
                econn.executemany("INSERT INTO edit_index VALUES (?, ?)", edit_rows)
                econn.commit()
            if struct_rows:
                sconn.executemany("INSERT INTO structural_index VALUES (?, ?, ?, ?, ?, ?, ?)", struct_rows)
                sconn.commit()
    finally:
        writer.close()
        econn.close()
        sconn.close()

    stats = {
        "stage": "build_indexes",
        "core_rows": total,
        "indexed_reactions": indexed,
        "missing_edit_hash": missing_edit,
        "edit_index_version": EDIT_INDEX_VERSION,
        "structural_index_version": STRUCTURAL_INDEX_VERSION,
        "edit_index": str(edit_pq),
        "edit_index_db": str(edit_db),
        "structural_index_db": str(structural_db),
    }
    write_json(outdir / EDIT_INDEX_STATS, stats)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--core")
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    print(json.dumps(build_indexes(args.output_dir, core_path=args.core, batch_size=args.batch_size, overwrite=args.overwrite), indent=2))


if __name__ == "__main__":
    main()
