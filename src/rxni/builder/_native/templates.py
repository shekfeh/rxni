"""Stage 5: aggregate native V5 reaction evidence to one row per template.

Optional old template-source files contribute identity/SMARTS metadata only.
No old structural, semantic, fingerprint, class, or family field is used.
"""

from __future__ import annotations
from .smarts_integrity import require_source_template_smarts, require_joined_template_smarts

import argparse
import json
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .common import ensure_dir, sql_to_u64, write_json
from .schema import (
    LIBRARY_SCHEMA_VERSION,
    STRUCTURAL_KEY_VERSION,
    TEMPLATE_LIBRARY_FILENAME,
    TEMPLATE_STATS,
    STAGING_DB,
)
from .source_adapter import iter_batches, resolve_template_mapping


HASH_FIELDS = (
    "forward_edit_signature_hash",
    "forward_raw_center_hash_r0",
    "forward_raw_center_hash_r1",
    "forward_raw_center_hash_r2",
    "forward_component_topology_hash",
    "forward_component_evidence_hash",
)


def _mode(counter: Counter[int]) -> tuple[int | None, int, int]:
    if not counter:
        return None, 0, 0
    # deterministic tie-break: highest count, then lowest exact uint64
    value, count = min(counter.items(), key=lambda kv: (-kv[1], kv[0]))
    return int(value), int(count), len(counter)


def _load_template_metadata(path: str | Path | None, batch_size: int) -> tuple[dict[str, dict[str, str | None]], dict[str, Any]]:
    require_source_template_smarts(path)
    if path is None:
        return {}, {"template_source": None, "rows": 0}
    p = Path(path)
    mapping = resolve_template_mapping(p)
    wanted = [mapping.template_hash]
    for c in (mapping.template_number, mapping.template_smarts):
        if c:
            wanted.append(c)
    wanted = [c for c in dict.fromkeys(wanted) if c]
    meta: dict[str, dict[str, str | None]] = {}
    rows = 0
    duplicate_hash_rows = 0
    for batch in iter_batches(p, batch_size=batch_size, columns=wanted):
        d = batch.to_pydict()
        n = batch.num_rows
        for i in range(n):
            rows += 1
            hv = d[mapping.template_hash][i]  # type: ignore[index]
            if hv is None or not str(hv):
                continue
            h = str(hv)
            rec = {
                "retro_template_number": str(d[mapping.template_number][i]) if mapping.template_number and d[mapping.template_number][i] is not None else None,
                "retro_template_smarts": str(d[mapping.template_smarts][i]) if mapping.template_smarts and d[mapping.template_smarts][i] is not None else None,
            }
            if h in meta:
                duplicate_hash_rows += 1
                old = meta[h]
                if old.get("retro_template_smarts") and rec.get("retro_template_smarts") and old["retro_template_smarts"] != rec["retro_template_smarts"]:
                    raise ValueError(f"Conflicting SMARTS for template hash {h}")
                if old.get("retro_template_number") is None:
                    old["retro_template_number"] = rec.get("retro_template_number")
                if old.get("retro_template_smarts") is None:
                    old["retro_template_smarts"] = rec.get("retro_template_smarts")
            else:
                meta[h] = rec
    return meta, {
        "template_source": str(p),
        "rows": rows,
        "unique_template_hashes": len(meta),
        "duplicate_hash_rows": duplicate_hash_rows,
        "mapping": mapping.as_dict(),
    }


def _schema() -> pa.Schema:
    fields: list[tuple[str, pa.DataType]] = [
        ("retro_template_hash", pa.string()),
        ("retro_template_number", pa.string()),
        ("retro_template_smarts", pa.string()),
        ("n_members", pa.int64()),
        ("membership_occurrence_count", pa.int64()),
        ("representative_reaction_key", pa.string()),
        ("representative_reaction_smiles", pa.string()),
    ]
    for f in HASH_FIELDS:
        fields.extend([
            (f, pa.uint64()),
            (f + "_consensus_fraction", pa.float64()),
            (f + "_n_distinct", pa.int64()),
        ])
    fields.extend([
        ("library_schema_version", pa.string()),
        ("structural_key_version", pa.string()),
    ])
    return pa.schema(fields)


def build_template_library(
    output_dir: str | Path,
    *,
    template_source: str | Path | None = None,
    batch_size: int = 8192,
    row_group_size: int = 4096,
    overwrite: bool = False,
) -> dict[str, Any]:
    outdir = ensure_dir(output_dir)
    staging = outdir / STAGING_DB
    out = outdir / TEMPLATE_LIBRARY_FILENAME
    if not staging.exists():
        raise FileNotFoundError(staging)
    if out.exists():
        if not overwrite:
            raise FileExistsError(out)
        out.unlink()

    metadata, meta_stats = _load_template_metadata(template_source, batch_size)
    conn = sqlite3.connect(staging)
    cur = conn.execute(
        """
        SELECT m.retro_template_hash, m.retro_template_number,
               m.membership_occurrence_count,
               c.reaction_key, c.forward_reaction_smiles,
               c.forward_edit_signature_hash,
               c.forward_raw_center_hash_r0,
               c.forward_raw_center_hash_r1,
               c.forward_raw_center_hash_r2,
               c.forward_component_topology_hash,
               c.forward_component_evidence_hash
        FROM membership AS m
        JOIN reaction_core AS c USING (reaction_key)
        ORDER BY m.retro_template_hash, c.reaction_key
        """
    )

    writer = pq.ParquetWriter(out, _schema(), compression="zstd")
    buffer: list[dict[str, Any]] = []
    template_count = 0
    member_count = 0
    occurrence_sum = 0

    current_hash: str | None = None
    current_number: str | None = None
    current_members = 0
    current_occ = 0
    representative_key: str | None = None
    representative_rxn: str | None = None
    counters = {f: Counter() for f in HASH_FIELDS}
    nonnull_counts = {f: 0 for f in HASH_FIELDS}

    def emit_current() -> None:
        nonlocal current_hash, current_number, current_members, current_occ
        nonlocal representative_key, representative_rxn, counters, nonnull_counts, template_count, buffer
        if current_hash is None:
            return
        meta = metadata.get(current_hash, {})
        row: dict[str, Any] = {
            "retro_template_hash": current_hash,
            "retro_template_number": current_number or meta.get("retro_template_number"),
            "retro_template_smarts": require_joined_template_smarts(metadata, current_hash),
            "n_members": current_members,
            "membership_occurrence_count": current_occ,
            "representative_reaction_key": representative_key,
            "representative_reaction_smiles": representative_rxn,
            "library_schema_version": LIBRARY_SCHEMA_VERSION,
            "structural_key_version": STRUCTURAL_KEY_VERSION,
        }
        for f in HASH_FIELDS:
            value, count, ndist = _mode(counters[f])
            row[f] = value
            row[f + "_consensus_fraction"] = (count / nonnull_counts[f]) if nonnull_counts[f] else 0.0
            row[f + "_n_distinct"] = ndist
        buffer.append(row)
        template_count += 1
        if len(buffer) >= row_group_size:
            writer.write_table(pa.Table.from_pylist(buffer, schema=_schema()), row_group_size=len(buffer))
            buffer.clear()

        current_hash = None
        current_number = None
        current_members = 0
        current_occ = 0
        representative_key = None
        representative_rxn = None
        counters = {f: Counter() for f in HASH_FIELDS}
        nonnull_counts = {f: 0 for f in HASH_FIELDS}

    try:
        while True:
            rows = cur.fetchmany(batch_size)
            if not rows:
                break
            for r in rows:
                h = str(r[0])
                if current_hash is not None and h != current_hash:
                    emit_current()
                if current_hash is None:
                    current_hash = h
                    current_number = r[1]
                    representative_key = r[3]
                    representative_rxn = r[4]
                elif current_number is None and r[1] is not None:
                    current_number = r[1]
                current_members += 1
                current_occ += int(r[2])
                member_count += 1
                occurrence_sum += int(r[2])
                for idx, f in enumerate(HASH_FIELDS, start=5):
                    v = sql_to_u64(r[idx], f)
                    if v is not None:
                        counters[f][v] += 1
                        nonnull_counts[f] += 1
        emit_current()
        if buffer:
            writer.write_table(pa.Table.from_pylist(buffer, schema=_schema()), row_group_size=len(buffer))
    finally:
        writer.close()
        conn.close()

    stats = {
        "stage": "build_template_library",
        "templates": template_count,
        "member_rows": member_count,
        "membership_occurrence_sum": occurrence_sum,
        "template_metadata": meta_stats,
        "output": str(out),
    }
    write_json(outdir / TEMPLATE_STATS, stats)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--template-source")
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--row-group-size", type=int, default=4096)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    print(json.dumps(build_template_library(
        args.output_dir, template_source=args.template_source, batch_size=args.batch_size,
        row_group_size=args.row_group_size, overwrite=args.overwrite,
    ), indent=2))


if __name__ == "__main__":
    main()
