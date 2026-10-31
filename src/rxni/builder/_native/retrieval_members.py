"""Stage 4: materialize runtime-compatible retrieval member rows from native V5 data."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .common import ensure_dir, sql_to_u64, write_json
from .schema import LIBRARY_SCHEMA_VERSION, REACTION_KEY_VERSION, RETRIEVAL_MEMBERS_FILENAME, STAGING_DB


def build_retrieval_members(output_dir: str | Path, *, batch_size: int = 8192, overwrite: bool = False) -> dict:
    outdir = ensure_dir(output_dir)
    staging = outdir / STAGING_DB
    out = outdir / RETRIEVAL_MEMBERS_FILENAME
    if not staging.exists():
        raise FileNotFoundError(staging)
    if out.exists():
        if not overwrite:
            raise FileExistsError(out)
        out.unlink()

    schema = pa.schema([
        ("retro_template_hash", pa.string()),
        ("retro_template_number", pa.string()),
        ("reaction_id", pa.string()),
        ("forward_reaction_hash", pa.string()),
        ("forward_reaction_smiles", pa.string()),
        ("forward_reaction_valid", pa.bool_()),
        ("feature_error", pa.bool_()),
        ("forward_edit_signature_hash", pa.uint64()),
        ("forward_raw_center_hash_r0", pa.uint64()),
        ("forward_raw_center_hash_r1", pa.uint64()),
        ("forward_raw_center_hash_r2", pa.uint64()),
        ("reaction_key", pa.string()),
        ("membership_occurrence_count", pa.int64()),
        ("library_schema_version", pa.string()),
        ("reaction_key_version", pa.string()),
    ])
    conn = sqlite3.connect(staging)
    cur = conn.execute(
        """
        SELECT m.retro_template_hash, m.retro_template_number,
               m.source_sample_reaction_id, c.forward_reaction_hash,
               c.forward_reaction_smiles, c.forward_reaction_valid,
               c.feature_error, c.forward_edit_signature_hash,
               c.forward_raw_center_hash_r0, c.forward_raw_center_hash_r1,
               c.forward_raw_center_hash_r2, c.reaction_key,
               m.membership_occurrence_count
        FROM membership AS m
        JOIN reaction_core AS c USING (reaction_key)
        ORDER BY m.retro_template_hash, c.reaction_key
        """
    )
    writer = pq.ParquetWriter(out, schema, compression="zstd")
    total = 0
    try:
        while True:
            rows = cur.fetchmany(batch_size)
            if not rows:
                break
            payload = []
            for r in rows:
                payload.append({
                    "retro_template_hash": r[0],
                    "retro_template_number": r[1],
                    "reaction_id": r[2],
                    "forward_reaction_hash": r[3],
                    "forward_reaction_smiles": r[4],
                    "forward_reaction_valid": bool(r[5]),
                    "feature_error": bool(r[6]),
                    "forward_edit_signature_hash": sql_to_u64(r[7], "forward_edit_signature_hash"),
                    "forward_raw_center_hash_r0": sql_to_u64(r[8], "forward_raw_center_hash_r0"),
                    "forward_raw_center_hash_r1": sql_to_u64(r[9], "forward_raw_center_hash_r1"),
                    "forward_raw_center_hash_r2": sql_to_u64(r[10], "forward_raw_center_hash_r2"),
                    "reaction_key": r[11],
                    "membership_occurrence_count": int(r[12]),
                    "library_schema_version": LIBRARY_SCHEMA_VERSION,
                    "reaction_key_version": REACTION_KEY_VERSION,
                })
            writer.write_table(pa.Table.from_pylist(payload, schema=schema), row_group_size=len(payload))
            total += len(payload)
    finally:
        writer.close()
        conn.close()

    stats = {"stage": "build_retrieval_members", "rows": total, "output": str(out)}
    write_json(outdir / "retrieval_members_stats.json", stats)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    print(json.dumps(build_retrieval_members(args.output_dir, batch_size=args.batch_size, overwrite=args.overwrite), indent=2))


if __name__ == "__main__":
    main()
