"""Stage 7: validation for a native Runtime-v6.5 V5 build."""

from __future__ import annotations

import argparse
import json
import random
import sqlite3
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from rxni.retrieval._core.algorithms.library_record_provider import LibraryRecordProvider
from rxni.retrieval._core.algorithms.query_structural_keys import derive_query_structural_keys
from rxni.retrieval._core.algorithms.retrieval_record import build_retrieval_record
from rxni.retrieval._core.algorithms.scoring import member_scores
from rxni.retrieval._core.algorithms.stored_record_codec import (
    retrieval_record_from_storage_row,
    retrieval_record_to_storage_row,
)

from .common import ensure_dir, u64_to_sql, write_json
from .schema import (
    EDIT_INDEX_DB,
    MEMBERSHIP_REL_FILENAME,
    PROVENANCE_FILENAME,
    REACTION_CORE_FILENAME,
    REPRESENTATION_ERRORS_FILENAME,
    REPRESENTATION_LOCATOR_DB,
    REPRESENTATIONS_FILENAME,
    RETRIEVAL_MEMBERS_FILENAME,
    STRUCTURAL_HASH_FIELDS,
    STRUCTURAL_INDEX_DB,
    STRUCTURAL_KEYS_FILENAME,
    TEMPLATE_LIBRARY_FILENAME,
    VALIDATION_FILENAME,
)


def _parquet_rows(path: Path) -> int:
    return int(pq.ParquetFile(path).metadata.num_rows)


def _count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8") as fh:
        return sum(1 for line in fh if line.strip())


def _membership_occurrence_sum(path: Path, batch_size: int = 65536) -> int:
    total = 0
    for batch in pq.ParquetFile(path).iter_batches(columns=["membership_occurrence_count"], batch_size=batch_size):
        total += sum(int(v) for v in batch.column(0).to_pylist() if v is not None)
    return total


def _eligible_count(path: Path, batch_size: int = 65536) -> int:
    total = 0
    cols = ["forward_reaction_valid", "feature_error", "forward_reaction_smiles"]
    for batch in pq.ParquetFile(path).iter_batches(columns=cols, batch_size=batch_size):
        d = batch.to_pydict()
        for i in range(batch.num_rows):
            if d["forward_reaction_valid"][i] and not d["feature_error"][i] and d["forward_reaction_smiles"][i]:
                total += 1
    return total


def _reservoir_sample(path: Path, n: int, seed: int = 6505, batch_size: int = 8192) -> list[dict[str, Any]]:
    rng = random.Random(seed)
    cols = [
        "reaction_key", "forward_reaction_smiles", "forward_reaction_valid", "feature_error",
        *STRUCTURAL_HASH_FIELDS,
    ]
    sample: list[dict[str, Any]] = []
    seen = 0
    for batch in pq.ParquetFile(path).iter_batches(columns=cols, batch_size=batch_size):
        for row in batch.to_pylist():
            if not row["forward_reaction_valid"] or row["feature_error"] or not row["forward_reaction_smiles"]:
                continue
            seen += 1
            if len(sample) < n:
                sample.append(row)
            else:
                j = rng.randrange(seen)
                if j < n:
                    sample[j] = row
    return sample


def validate_library_v5(output_dir: str | Path, *, sample_size: int = 10000, seed: int = 6505) -> dict[str, Any]:
    outdir = ensure_dir(output_dir)
    paths = {
        "core": outdir / REACTION_CORE_FILENAME,
        "membership": outdir / MEMBERSHIP_REL_FILENAME,
        "retrieval": outdir / RETRIEVAL_MEMBERS_FILENAME,
        "provenance": outdir / PROVENANCE_FILENAME,
        "representations": outdir / REPRESENTATIONS_FILENAME,
        "locator_db": outdir / REPRESENTATION_LOCATOR_DB,
        "errors": outdir / REPRESENTATION_ERRORS_FILENAME,
        "template": outdir / TEMPLATE_LIBRARY_FILENAME,
        "edit_db": outdir / EDIT_INDEX_DB,
        "structural_db": outdir / STRUCTURAL_INDEX_DB,
        "structural_keys": outdir / STRUCTURAL_KEYS_FILENAME,
    }
    missing = [str(p) for k, p in paths.items() if k != "errors" and not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing native V5 artifacts: {missing}")

    failures: list[str] = []
    checks: dict[str, Any] = {}

    core_rows = _parquet_rows(paths["core"])
    membership_rows = _parquet_rows(paths["membership"])
    retrieval_rows = _parquet_rows(paths["retrieval"])
    provenance_rows = _parquet_rows(paths["provenance"])
    representation_rows = _parquet_rows(paths["representations"])
    template_rows = _parquet_rows(paths["template"])
    structural_key_rows = _parquet_rows(paths["structural_keys"])
    representation_errors = _count_lines(paths["errors"])
    eligible = _eligible_count(paths["core"])
    occurrence_sum = _membership_occurrence_sum(paths["membership"])

    lconn = sqlite3.connect(paths["locator_db"])
    locator_rows = int(lconn.execute("SELECT COUNT(*) FROM locator").fetchone()[0])
    lconn.close()
    sconn = sqlite3.connect(paths["structural_db"])
    structural_index_rows = int(sconn.execute("SELECT COUNT(*) FROM structural_index").fetchone()[0])
    sconn.close()
    econn = sqlite3.connect(paths["edit_db"])
    edit_index_rows = int(econn.execute("SELECT COUNT(*) FROM edit_index").fetchone()[0])
    econn.close()

    checks["counts"] = {
        "reaction_core_rows": core_rows,
        "eligible_reaction_rows": eligible,
        "membership_rows": membership_rows,
        "retrieval_member_rows": retrieval_rows,
        "provenance_rows": provenance_rows,
        "membership_occurrence_sum": occurrence_sum,
        "representation_rows": representation_rows,
        "representation_error_rows": representation_errors,
        "locator_rows": locator_rows,
        "structural_key_rows": structural_key_rows,
        "structural_index_rows": structural_index_rows,
        "edit_index_rows": edit_index_rows,
        "template_rows": template_rows,
    }

    if membership_rows != retrieval_rows:
        failures.append("membership row count != retrieval member row count")
    if occurrence_sum > provenance_rows:
        failures.append("membership occurrence sum exceeds provenance row count")
    if structural_key_rows != eligible:
        failures.append("structural-key row count != eligible core row count")
    if representation_rows != eligible:
        failures.append("representation row count != eligible core row count")
    if locator_rows != representation_rows:
        failures.append("locator row count != representation row count")
    if structural_index_rows != eligible:
        failures.append("structural index row count != eligible core row count")
    if edit_index_rows > eligible:
        failures.append("edit index row count exceeds eligible core row count")
    if representation_errors:
        failures.append(f"representation build errors: {representation_errors}")

    sample = _reservoir_sample(paths["core"], sample_size, seed=seed)
    provider = LibraryRecordProvider(paths["representations"], paths["locator_db"], fallback_provider=None)

    structural_mismatches = {f: 0 for f in STRUCTURAL_HASH_FIELDS}
    storage_mismatches = 0
    scoring_mismatches = 0
    missing_stored = 0
    codec_normalizations = 0
    locator_edit_mismatches = 0
    edit_index_missing = 0
    structural_index_mismatches = 0
    live_records: list[dict[str, Any]] = []
    stored_records: list[dict[str, Any]] = []
    sample_by_key = {str(r["reaction_key"]): r for r in sample}

    lconn = sqlite3.connect(paths["locator_db"])
    econn = sqlite3.connect(paths["edit_db"])
    sconn = sqlite3.connect(paths["structural_db"])

    for row in sample:
        key = str(row["reaction_key"])
        reaction = str(row["forward_reaction_smiles"])
        live_keys = derive_query_structural_keys(reaction, store_readable=False)
        for field in STRUCTURAL_HASH_FIELDS:
            if live_keys.get(field) != row.get(field):
                structural_mismatches[field] += 1

        loc = lconn.execute("SELECT forward_edit_signature_hash FROM locator WHERE reaction_key=?", (key,)).fetchone()
        if loc is None or (None if loc[0] is None else int(loc[0])) != row.get("forward_edit_signature_hash"):
            locator_edit_mismatches += 1

        edit = row.get("forward_edit_signature_hash")
        if edit is not None:
            found = econn.execute(
                "SELECT 1 FROM edit_index WHERE forward_edit_signature_hash=? AND reaction_key=?",
                (u64_to_sql(edit, "forward_edit_signature_hash"), key),
            ).fetchone()
            if found is None:
                edit_index_missing += 1

        sr = sconn.execute(
            """
            SELECT forward_edit_signature_hash, forward_raw_center_hash_r0,
                   forward_raw_center_hash_r1, forward_raw_center_hash_r2,
                   forward_component_topology_hash, forward_component_evidence_hash
            FROM structural_index WHERE reaction_key=?
            """,
            (key,),
        ).fetchone()
        expected = [row.get(f) for f in STRUCTURAL_HASH_FIELDS]
        got = None if sr is None else [None if v is None else int(v) for v in sr]
        if got != expected:
            structural_index_mismatches += 1

        live = build_retrieval_record(reaction, record_id=key)
        try:
            stored = provider.get(reaction, record_id=key)
        except KeyError:
            missing_stored += 1
            continue
        canonical_live = retrieval_record_from_storage_row(retrieval_record_to_storage_row(live))
        if canonical_live != live:
            codec_normalizations += 1
        if canonical_live != stored:
            storage_mismatches += 1
        live_records.append(live)
        stored_records.append(stored)

    lconn.close()
    econn.close()
    sconn.close()

    # Validate retrieval member copies for sampled keys in one sequential scan.
    retrieval_mismatches = 0
    retrieval_seen: set[str] = set()
    rfields = ("forward_edit_signature_hash", "forward_raw_center_hash_r0", "forward_raw_center_hash_r1", "forward_raw_center_hash_r2")
    for batch in pq.ParquetFile(paths["retrieval"]).iter_batches(columns=["reaction_key", *rfields], batch_size=65536):
        for rr in batch.to_pylist():
            key = str(rr["reaction_key"])
            core = sample_by_key.get(key)
            if core is None:
                continue
            retrieval_seen.add(key)
            if any(rr.get(f) != core.get(f) for f in rfields):
                retrieval_mismatches += 1
    retrieval_sample_missing = len(set(sample_by_key) - retrieval_seen)

    for i in range(max(0, min(len(live_records), len(stored_records)) - 1)):
        if member_scores(live_records[i], live_records[i + 1]) != member_scores(stored_records[i], stored_records[i + 1]):
            scoring_mismatches += 1

    checks["sample_parity"] = {
        "requested": sample_size,
        "sampled": len(sample),
        "seed": seed,
        "structural_key_mismatches": structural_mismatches,
        "stored_record_missing": missing_stored,
        "retrieval_record_mismatches": storage_mismatches,
        "scoring_mismatches": scoring_mismatches,
        "storage_codec_normalizations": codec_normalizations,
        "locator_edit_mismatches": locator_edit_mismatches,
        "edit_index_missing": edit_index_missing,
        "structural_index_mismatches": structural_index_mismatches,
        "retrieval_member_mismatches": retrieval_mismatches,
        "retrieval_sample_missing": retrieval_sample_missing,
        "provider_snapshot": provider.snapshot(),
    }

    total_structural = sum(structural_mismatches.values())
    if total_structural:
        failures.append(f"live/stored structural-key mismatches: {structural_mismatches}")
    if missing_stored:
        failures.append(f"sample stored representations missing: {missing_stored}")
    if storage_mismatches:
        failures.append(f"stored/live canonical retrieval-record mismatches: {storage_mismatches}")
    if scoring_mismatches:
        failures.append(f"stored/live scoring mismatches: {scoring_mismatches}")
    if locator_edit_mismatches:
        failures.append(f"locator edit-hash mismatches: {locator_edit_mismatches}")
    if edit_index_missing:
        failures.append(f"sample edit-index entries missing: {edit_index_missing}")
    if structural_index_mismatches:
        failures.append(f"structural-index mismatches: {structural_index_mismatches}")
    if retrieval_mismatches or retrieval_sample_missing:
        failures.append(f"retrieval-member parity failure: mismatches={retrieval_mismatches}, missing={retrieval_sample_missing}")

    result = {
        "status": "PASS" if not failures else "FAIL",
        "output_dir": str(outdir),
        "checks": checks,
        "failures": failures,
    }
    write_json(outdir / VALIDATION_FILENAME, result)
    if failures:
        raise RuntimeError("Native Library-v5 validation failed: " + "; ".join(failures))
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--sample-size", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=6505)
    args = ap.parse_args()
    print(json.dumps(validate_library_v5(args.output_dir, sample_size=args.sample_size, seed=args.seed), indent=2))


if __name__ == "__main__":
    main()
