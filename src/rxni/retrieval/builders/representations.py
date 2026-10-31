"""Restartable builder for Runtime-v6.6 stored retrieval representation v2.

Only derived retrieval representations are rebuilt. Library-v5 reaction core,
membership, provenance, and exact structural indices are read-only inputs.

Build lifecycle:
  reaction_core_v5.parquet
    -> deterministic chunks
    -> parts/part-XXXXXX.parquet
    -> finalize
    -> retrieval_representations_v2.parquet
    -> representation_locator_v2.sqlite
    -> representation_locator_v2.parquet
    -> retrieval_v2_manifest.json
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from ..storage.codec import (
    retrieval_storage_schema,
    retrieval_record_to_storage_row_v2,
)
from ..record import build_retrieval_record
from ..storage.schema_profiles import CANONICAL_SCHEMA_PROFILE
from .source_schema import (
    load_profile, source_columns, iter_source, sha256 as source_sha256,
)
from ..version import (
    RETRIEVAL_API_VERSION,
    STORED_REPRESENTATION_GENERATION,
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def _worker(item):
    ordinal, reaction_key, reaction, build_id = item
    try:
        rec = build_retrieval_record(
            reaction,
            record_id=reaction_key,
        )
        row = retrieval_record_to_storage_row_v2(
            rec,
            build_id=build_id,
        )
        # The V5 codec normally persists record_id.  The v2 physical contract
        # requires reaction_key explicitly for locator validation.
        row["reaction_key"] = reaction_key
        return {
            "ok": True,
            "ordinal": ordinal,
            "reaction_key": reaction_key,
            "row": row,
        }
    except Exception as exc:
        return {
            "ok": False,
            "ordinal": ordinal,
            "reaction_key": reaction_key,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _core_columns(core: Path):
    schema = pq.read_schema(core)
    names = schema.names
    key = next(
        c
        for c in (
            "reaction_key",
            "forward_reaction_hash",
            "reaction_hash",
        )
        if c in names
    )
    rxn = next(
        c
        for c in (
            "forward_reaction_smiles",
            "reaction",
            "reaction_smiles",
            "mapped_reaction",
        )
        if c in names
    )
    valid = "forward_reaction_valid" if "forward_reaction_valid" in names else None
    ferr = "feature_error" if "feature_error" in names else None
    cols = [key, rxn]
    if valid:
        cols.append(valid)
    if ferr:
        cols.append(ferr)
    return key, rxn, valid, ferr, cols


def _iter_core(core: Path, limit: int | None):
    key, rxn, valid, ferr, cols = _core_columns(core)
    pf = pq.ParquetFile(core)
    ordinal = 0
    emitted = 0
    for batch in pf.iter_batches(columns=cols, batch_size=8192):
        for row in batch.to_pylist():
            if valid and not bool(row.get(valid)):
                continue
            if ferr and str(row.get(ferr) or "").strip():
                continue
            if not row.get(key) or not row.get(rxn):
                continue
            yield ordinal, str(row[key]), str(row[rxn])
            ordinal += 1
            emitted += 1
            if limit is not None and emitted >= limit:
                return


def _storage_schema(library: Path) -> pa.Schema:
    """Build the V2 physical schema directly from its own storage codec.

    The ``library`` argument remains for internal call-site compatibility.
    A zero-row V1 representation is no longer required or consulted.
    """
    return retrieval_storage_schema()


def _part_paths(parts: Path, part_id: int):
    stem = f"part-{part_id:06d}"
    return (
        parts / f"{stem}.parquet",
        parts / f"{stem}.json",
        parts / f"{stem}.errors.jsonl",
    )


def build_parts(
    *,
    library: Path,
    output_dir: Path,
    build_id: str,
    workers: int,
    chunk_size: int,
    limit: int | None,
    resume: bool,
    source_core: Path | None = None,
    source_profile=None,
):
    parts = output_dir / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    schema = _storage_schema(library)
    core = source_core or (library / "reaction_core_v5.parquet")
    source_rows = (iter_source(core, source_profile, limit) if source_profile is not None
                   else _iter_core(core, limit))

    total_ok = total_err = 0
    part_id = 0
    chunk = []

    executor = cf.ProcessPoolExecutor(max_workers=workers)
    try:
        for ordinal, key, rxn in source_rows:
            chunk.append((ordinal, key, rxn, build_id))
            if len(chunk) >= chunk_size:
                ok, err = _build_one_part(
                    executor, chunk, schema, parts, part_id, resume
                )
                total_ok += ok
                total_err += err
                part_id += 1
                chunk = []
        if chunk:
            ok, err = _build_one_part(executor, chunk, schema, parts, part_id, resume)
            total_ok += ok
            total_err += err
            part_id += 1
    finally:
        executor.shutdown(wait=True)

    return {
        "parts": part_id,
        "represented": total_ok,
        "errors": total_err,
    }


def _build_one_part(executor, chunk, schema, parts, part_id, resume):
    parquet_path, meta_path, error_path = _part_paths(parts, part_id)
    expected_keys = [x[1] for x in chunk]

    if resume and parquet_path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if (
            meta.get("input_reaction_keys_sha256")
            == hashlib.sha256("\n".join(expected_keys).encode()).hexdigest()
            and meta.get("status") == "complete"
        ):
            print(
                f"[part {part_id:06d}] resume-skip "
                f"ok={meta['represented']} errors={meta['errors']}"
            )
            return int(meta["represented"]), int(meta["errors"])

    results = list(
        executor.map(
            _worker,
            chunk,
            chunksize=max(1, min(32, len(chunk) // max(1, os.cpu_count() or 1))),
        )
    )
    results.sort(key=lambda x: int(x["ordinal"]))
    good = [x for x in results if x["ok"]]
    bad = [x for x in results if not x["ok"]]

    rows = [x["row"] for x in good]
    table = pa.Table.from_pylist(rows, schema=schema)
    pq.write_table(
        table,
        parquet_path,
        compression="zstd",
        row_group_size=max(1, len(rows)),
    )

    with error_path.open("w", encoding="utf-8") as f:
        for x in bad:
            f.write(json.dumps(x, sort_keys=True) + "\n")

    meta = {
        "status": "complete",
        "part_id": part_id,
        "input_rows": len(chunk),
        "represented": len(good),
        "errors": len(bad),
        "first_ordinal": chunk[0][0],
        "last_ordinal": chunk[-1][0],
        "input_reaction_keys_sha256": hashlib.sha256(
            "\n".join(expected_keys).encode()
        ).hexdigest(),
        "parquet_sha256": _sha256(parquet_path),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    print(
        f"[part {part_id:06d}] " f"input={len(chunk)} ok={len(good)} errors={len(bad)}"
    )
    return len(good), len(bad)


def finalize(
    *,
    library: Path,
    output_dir: Path,
    build_id: str,
    source_core: Path | None = None,
    source_binding: dict | None = None,
):
    parts = output_dir / "parts"
    part_files = sorted(parts.glob("part-*.parquet"))
    if not part_files:
        raise SystemExit("No part files to finalize")

    out_parquet = output_dir / "retrieval_representations_v2.parquet"
    out_locator = output_dir / "representation_locator_v2.sqlite"
    out_locator_pq = output_dir / "representation_locator_v2.parquet"
    manifest_path = output_dir / "retrieval_v2_manifest.json"

    if out_parquet.exists():
        out_parquet.unlink()
    if out_locator.exists():
        out_locator.unlink()
    if out_locator_pq.exists():
        out_locator_pq.unlink()

    writer = None
    locator_rows = []
    row_group_id = 0
    total = 0

    try:
        for part in part_files:
            table = pq.read_table(part)
            if writer is None:
                writer = pq.ParquetWriter(
                    out_parquet,
                    table.schema,
                    compression="zstd",
                )
            # Force exactly one final row group per deterministic part.
            writer.write_table(
                table,
                row_group_size=max(1, table.num_rows),
            )

            keys = table["reaction_key"].to_pylist()
            for i, key in enumerate(keys):
                locator_rows.append(
                    {
                        "reaction_key": str(key),
                        "row_group_id": row_group_id,
                        "row_in_group": i,
                        "build_id": build_id,
                        "representation_generation": STORED_REPRESENTATION_GENERATION,
                    }
                )
            total += table.num_rows
            row_group_id += 1
    finally:
        if writer is not None:
            writer.close()

    conn = sqlite3.connect(out_locator)
    try:
        conn.execute(
            """
            CREATE TABLE representation_locator_v2 (
                reaction_key TEXT PRIMARY KEY,
                row_group_id INTEGER NOT NULL,
                row_in_group INTEGER NOT NULL,
                build_id TEXT NOT NULL,
                representation_generation TEXT NOT NULL
            )
            """
        )
        conn.executemany(
            """
            INSERT INTO representation_locator_v2
            (reaction_key,row_group_id,row_in_group,build_id,representation_generation)
            VALUES (:reaction_key,:row_group_id,:row_in_group,:build_id,:representation_generation)
            """,
            locator_rows,
        )
        conn.execute(
            """
            CREATE INDEX locator_v2_row_group
            ON representation_locator_v2(row_group_id)
            """
        )
        conn.commit()
    finally:
        conn.close()

    pq.write_table(
        pa.Table.from_pylist(locator_rows),
        out_locator_pq,
        compression="zstd",
    )

    error_count = 0
    for ep in sorted(parts.glob("part-*.errors.jsonl")):
        with ep.open() as f:
            error_count += sum(1 for line in f if line.strip())

    manifest = {
        "schema_version": "v660_retrieval_v2_manifest_v1",
        "physical_schema_profile": CANONICAL_SCHEMA_PROFILE,
        "build_id": build_id,
        "representation_generation": STORED_REPRESENTATION_GENERATION,
        "retrieval_api_version": RETRIEVAL_API_VERSION,
        "source_library": str(library),
        "source_reaction_core": str(source_core or (library / "reaction_core_v5.parquet")),
        "source_reaction_core_sha256": _sha256(source_core or (library / "reaction_core_v5.parquet")),
        "rows": total,
        "row_groups": row_group_id,
        "errors": error_count,
        "representations_sha256": _sha256(out_parquet),
        "locator_sqlite_sha256": _sha256(out_locator),
        "locator_parquet_sha256": _sha256(out_locator_pq),
        "scientific_contract": {
            "semantic_authority": "runtime_v6_7.semantic validated canonical backend",
            "structural_builder": "v6.5 structural compatibility bridge under v6.6 retrieval API",
            "scoring": "frozen v6.5 scoring equations exposed through runtime_v6_6.retrieval",
            "candidate_core_membership_provenance": "unchanged Library-v5 foundation",
        },
    }
    if source_binding is not None:
        # Optional provenance extension; all frozen Retrieval-V2 ABI keys and
        # physical record schemas remain unchanged.
        manifest["source_schema_adapter"] = {
            "schema_version": 1,
            "profile_id": source_binding["profile"]["profile_id"],
            "profile_sha256": source_binding["profile_sha256"],
            "source_sha256": source_binding["source_sha256"],
        }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    print(json.dumps(manifest, indent=2))
    return manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--input-core", default=None,
                    help="Alternate Parquet/CSV reaction source (requires --input-schema)")
    ap.add_argument("--input-schema", default=None,
                    help="Versioned JSON source-schema profile; never changes V2 storage ABI")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--chunk-size", type=int, default=2048)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--finalize-only", action="store_true")
    ap.add_argument("--build-id", default=None)
    args = ap.parse_args()

    library = Path(args.library)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    source_profile = None
    source_core = None
    binding = None
    binding_file = output_dir / "source_binding_v1.json"
    if args.input_core and not args.input_schema:
        raise SystemExit("--input-core requires an explicit --input-schema profile")
    if args.input_schema:
        profile_path = Path(args.input_schema).expanduser().resolve()
        source_profile = load_profile(profile_path)
        source_core = Path(args.input_core or (library / "reaction_core_v5.parquet")).expanduser().resolve()
        source_columns(source_core, source_profile)
        binding = {
            "schema_version": 1,
            "source": str(source_core),
            "source_sha256": source_sha256(source_core),
            "profile_sha256": source_sha256(profile_path),
            "profile": source_profile.as_dict(),
        }
        if binding_file.is_file():
            existing = json.loads(binding_file.read_text(encoding="utf-8"))
            if existing != binding:
                raise SystemExit("Source profile or file changed since first build; use a new output-dir")
        else:
            if (output_dir / "parts").exists() or (output_dir / "BUILD_ID").exists():
                raise SystemExit("Existing output has no source binding; use a clean output-dir")
            binding_file.write_text(json.dumps(binding, indent=2, sort_keys=True) + "\n")
    elif binding_file.exists():
        raise SystemExit("Existing output uses source profile; supply matching --input-schema")

    build_id_file = output_dir / "BUILD_ID"
    if args.build_id:
        build_id = args.build_id
        if build_id_file.exists() and build_id_file.read_text().strip() != build_id:
            raise SystemExit("Existing BUILD_ID differs from --build-id")
        build_id_file.write_text(build_id + "\n")
    elif build_id_file.exists():
        build_id = build_id_file.read_text().strip()
    else:
        build_id = datetime.now(timezone.utc).strftime(
            "v66_retrieval_v2_%Y%m%dT%H%M%SZ"
        )
        build_id_file.write_text(build_id + "\n")

    print(f"build_id={build_id}")

    if not args.finalize_only:
        stats = build_parts(
            library=library,
            output_dir=output_dir,
            build_id=build_id,
            workers=args.workers,
            chunk_size=args.chunk_size,
            limit=args.limit,
            resume=args.resume,
            source_core=source_core,
            source_profile=source_profile,
        )
        print("PART BUILD")
        print(json.dumps(stats, indent=2))

    finalize(
        library=library,
        output_dir=output_dir,
        build_id=build_id,
        source_core=source_core,
        source_binding=binding,
    )


if __name__ == "__main__":
    main()
