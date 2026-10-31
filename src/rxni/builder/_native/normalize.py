"""Stage 1: normalize old source identity/membership/provenance only.

No structural hash, semantic label, fingerprint, class/family, or System-I
feature is imported from the source file.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .common import ensure_dir, reaction_key, sha256_text, write_json
from .schema import (
    BUILDER_VERSION,
    LIBRARY_SCHEMA_VERSION,
    REACTION_KEY_VERSION,
    MEMBERSHIP_REL_FILENAME,
    PROVENANCE_FILENAME,
    REACTION_CORE_BASE_FILENAME,
    SOURCE_NORMALIZE_STATS,
    STAGING_DB,
)
from .source_adapter import derived_columns, iter_batches, resolve_membership_mapping, source_schema


def _safe_str(value: Any) -> str | None:
    if value is None:
        return None
    s = str(value)
    return s if s else None


def _open_staging(path: Path, *, overwrite: bool) -> sqlite3.Connection:
    if path.exists():
        if not overwrite:
            raise FileExistsError(f"Staging database already exists: {path}")
        path.unlink()
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.executescript(
        """
        CREATE TABLE reaction_core (
            reaction_key TEXT PRIMARY KEY,
            forward_reaction_hash TEXT NOT NULL,
            forward_reaction_smiles TEXT NOT NULL,
            source_first_reaction_id TEXT,
            source_occurrence_count INTEGER NOT NULL DEFAULT 1,
            forward_reaction_valid INTEGER,
            feature_error INTEGER,
            structural_error TEXT,
            forward_edit_signature_hash TEXT,
            forward_raw_center_hash_r0 TEXT,
            forward_raw_center_hash_r1 TEXT,
            forward_raw_center_hash_r2 TEXT,
            forward_component_topology_hash TEXT,
            forward_component_evidence_hash TEXT
        );

        CREATE TABLE membership (
            reaction_key TEXT NOT NULL,
            retro_template_hash TEXT NOT NULL,
            retro_template_number TEXT,
            source_sample_reaction_id TEXT,
            membership_occurrence_count INTEGER NOT NULL DEFAULT 1,
            PRIMARY KEY (reaction_key, retro_template_hash)
        );

        CREATE INDEX idx_membership_template
            ON membership(retro_template_hash, reaction_key);
        """
    )
    return conn


def _export_core_base(conn: sqlite3.Connection, out: Path, batch_size: int) -> int:
    schema = pa.schema([
        ("reaction_key", pa.string()),
        ("forward_reaction_hash", pa.string()),
        ("forward_reaction_smiles", pa.string()),
        ("source_first_reaction_id", pa.string()),
        ("source_occurrence_count", pa.int64()),
        ("library_schema_version", pa.string()),
        ("reaction_key_version", pa.string()),
        ("builder_version", pa.string()),
    ])
    cur = conn.execute(
        """
        SELECT reaction_key, forward_reaction_hash, forward_reaction_smiles,
               source_first_reaction_id, source_occurrence_count
        FROM reaction_core ORDER BY reaction_key
        """
    )
    writer = pq.ParquetWriter(out, schema, compression="zstd")
    total = 0
    try:
        while True:
            rows = cur.fetchmany(batch_size)
            if not rows:
                break
            payload = [{
                "reaction_key": r[0],
                "forward_reaction_hash": r[1],
                "forward_reaction_smiles": r[2],
                "source_first_reaction_id": r[3],
                "source_occurrence_count": int(r[4]),
                "library_schema_version": LIBRARY_SCHEMA_VERSION,
                "reaction_key_version": REACTION_KEY_VERSION,
                "builder_version": BUILDER_VERSION,
            } for r in rows]
            writer.write_table(pa.Table.from_pylist(payload, schema=schema), row_group_size=len(payload))
            total += len(payload)
    finally:
        writer.close()
    return total


def _export_membership(conn: sqlite3.Connection, out: Path, batch_size: int) -> int:
    schema = pa.schema([
        ("reaction_key", pa.string()),
        ("retro_template_hash", pa.string()),
        ("retro_template_number", pa.string()),
        ("membership_occurrence_count", pa.int64()),
        ("source_sample_reaction_id", pa.string()),
        ("library_schema_version", pa.string()),
        ("reaction_key_version", pa.string()),
    ])
    cur = conn.execute(
        """
        SELECT reaction_key, retro_template_hash, retro_template_number,
               membership_occurrence_count, source_sample_reaction_id
        FROM membership ORDER BY retro_template_hash, reaction_key
        """
    )
    writer = pq.ParquetWriter(out, schema, compression="zstd")
    total = 0
    try:
        while True:
            rows = cur.fetchmany(batch_size)
            if not rows:
                break
            payload = [{
                "reaction_key": r[0],
                "retro_template_hash": r[1],
                "retro_template_number": r[2],
                "membership_occurrence_count": int(r[3]),
                "source_sample_reaction_id": r[4],
                "library_schema_version": LIBRARY_SCHEMA_VERSION,
                "reaction_key_version": REACTION_KEY_VERSION,
            } for r in rows]
            writer.write_table(pa.Table.from_pylist(payload, schema=schema), row_group_size=len(payload))
            total += len(payload)
    finally:
        writer.close()
    return total


def normalize_source_membership(
    membership_source: str | Path,
    output_dir: str | Path,
    *,
    batch_size: int = 8192,
    overwrite: bool = False,
    progress_every_batches: int = 25,
) -> dict[str, Any]:
    source = Path(membership_source)
    outdir = ensure_dir(output_dir)
    if not source.exists():
        raise FileNotFoundError(source)

    mapping = resolve_membership_mapping(source)
    schema = source_schema(source)
    ignored = derived_columns(set(schema.names))

    staging = outdir / STAGING_DB
    conn = _open_staging(staging, overwrite=overwrite)

    provenance = outdir / PROVENANCE_FILENAME
    if provenance.exists():
        if not overwrite:
            raise FileExistsError(provenance)
        provenance.unlink()

    # Provenance intentionally preserves the original row plus the native
    # reaction_key. Derived source fields may appear here only as history.
    prov_writer: pq.ParquetWriter | None = None
    source_rows = 0
    missing_reaction_rows = 0
    missing_template_rows = 0
    batch_no = 0

    core_sql = """
        INSERT INTO reaction_core (
            reaction_key, forward_reaction_hash, forward_reaction_smiles,
            source_first_reaction_id, source_occurrence_count
        ) VALUES (?, ?, ?, ?, 1)
        ON CONFLICT(reaction_key) DO UPDATE SET
            source_occurrence_count = reaction_core.source_occurrence_count + 1,
            source_first_reaction_id = COALESCE(reaction_core.source_first_reaction_id, excluded.source_first_reaction_id)
    """
    member_sql = """
        INSERT INTO membership (
            reaction_key, retro_template_hash, retro_template_number,
            source_sample_reaction_id, membership_occurrence_count
        ) VALUES (?, ?, ?, ?, 1)
        ON CONFLICT(reaction_key, retro_template_hash) DO UPDATE SET
            retro_template_number = COALESCE(membership.retro_template_number, excluded.retro_template_number),
            source_sample_reaction_id = COALESCE(membership.source_sample_reaction_id, excluded.source_sample_reaction_id),
            membership_occurrence_count = membership.membership_occurrence_count + 1
    """

    required_columns = [mapping.reaction, mapping.template_hash]
    for c in (mapping.template_number, mapping.reaction_id, mapping.reaction_hash):
        if c:
            required_columns.append(c)
    required_columns = list(dict.fromkeys(c for c in required_columns if c))

    try:
        for batch in iter_batches(source, batch_size=batch_size):
            batch_no += 1
            table = pa.Table.from_batches([batch])
            data = table.to_pydict()
            n = table.num_rows
            source_rows += n
            keys: list[str | None] = []
            core_rows = []
            member_rows = []

            for i in range(n):
                rv = data[mapping.reaction][i]  # type: ignore[index]
                if rv is None or not str(rv):
                    keys.append(None)
                    missing_reaction_rows += 1
                    continue
                rxn = str(rv)
                rkey = reaction_key(rxn)
                keys.append(rkey)
                rid = _safe_str(data[mapping.reaction_id][i]) if mapping.reaction_id else None
                supplied_hash = _safe_str(data[mapping.reaction_hash][i]) if mapping.reaction_hash else None
                # Never trust source reaction hash for identity; keep the native
                # hash deterministic from the exact normalized source string.
                rhash = sha256_text(rxn)
                core_rows.append((rkey, rhash, rxn, rid))

                tv = data[mapping.template_hash][i]  # type: ignore[index]
                if tv is None or not str(tv):
                    missing_template_rows += 1
                    continue
                tnum = _safe_str(data[mapping.template_number][i]) if mapping.template_number else None
                member_rows.append((rkey, str(tv), tnum, rid))

            if core_rows:
                conn.executemany(core_sql, core_rows)
            if member_rows:
                conn.executemany(member_sql, member_rows)
            conn.commit()

            if "builder_v5_reaction_key" in table.column_names:
                raise ValueError("Source already contains builder_v5_reaction_key")
            prov_table = table.append_column("builder_v5_reaction_key", pa.array(keys, type=pa.string()))
            if prov_writer is None:
                prov_writer = pq.ParquetWriter(provenance, prov_table.schema, compression="zstd")
            prov_writer.write_table(prov_table, row_group_size=max(1, prov_table.num_rows))

            if progress_every_batches > 0 and batch_no % progress_every_batches == 0:
                ur = conn.execute("SELECT COUNT(*) FROM reaction_core").fetchone()[0]
                up = conn.execute("SELECT COUNT(*) FROM membership").fetchone()[0]
                print(f"[normalize] batches={batch_no} rows={source_rows} unique_reactions={ur} unique_pairs={up}", flush=True)
    finally:
        if prov_writer is not None:
            prov_writer.close()

    unique_reactions = int(conn.execute("SELECT COUNT(*) FROM reaction_core").fetchone()[0])
    unique_pairs = int(conn.execute("SELECT COUNT(*) FROM membership").fetchone()[0])
    occurrence_sum = int(conn.execute("SELECT COALESCE(SUM(membership_occurrence_count),0) FROM membership").fetchone()[0])

    base_rows = _export_core_base(conn, outdir / REACTION_CORE_BASE_FILENAME, batch_size)
    membership_rows = _export_membership(conn, outdir / MEMBERSHIP_REL_FILENAME, batch_size)
    conn.close()

    stats = {
        "stage": "normalize_source_membership",
        "source": str(source),
        "source_rows": source_rows,
        "unique_reactions": unique_reactions,
        "unique_reaction_template_pairs": unique_pairs,
        "membership_occurrence_sum": occurrence_sum,
        "missing_reaction_rows": missing_reaction_rows,
        "missing_template_rows": missing_template_rows,
        "reaction_core_base_rows": base_rows,
        "membership_rows": membership_rows,
        "mapping": mapping.as_dict(),
        "ignored_derived_source_columns": ignored,
        "source_reaction_hash_is_authoritative": False,
    }
    write_json(outdir / SOURCE_NORMALIZE_STATS, stats)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership-source", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    print(json.dumps(normalize_source_membership(
        args.membership_source,
        args.output_dir,
        batch_size=args.batch_size,
        overwrite=args.overwrite,
    ), indent=2))


if __name__ == "__main__":
    main()
