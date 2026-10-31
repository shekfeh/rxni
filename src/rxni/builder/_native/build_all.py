"""Orchestrate a native Runtime-v6.5 Template Library v5 build.

No V4 library artifact is required. Old files are used only for source facts:
reaction/template identity, template SMARTS/number, membership, provenance.
Every structural/System-I field is recomputed by Runtime-v6.5.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .common import ensure_dir
from .indexes import build_indexes
from .manifest import write_manifest
from .normalize import normalize_source_membership
from .representations import build_system_i_representations
from .retrieval_members import build_retrieval_members
from .schema import (
    EDIT_INDEX_FILENAME,
    MEMBERSHIP_REL_FILENAME,
    REACTION_CORE_BASE_FILENAME,
    REACTION_CORE_FILENAME,
    REPRESENTATIONS_FILENAME,
    REPRESENTATION_LOCATOR_DB,
    REPRESENTATION_LOCATOR_PARQUET,
    REPRESENTATION_ERRORS_FILENAME,
    RETRIEVAL_MEMBERS_FILENAME,
    SOURCE_NORMALIZE_STATS,
    STRUCTURAL_KEYS_FILENAME,
    STRUCTURAL_STATS,
    TEMPLATE_LIBRARY_FILENAME,
    TEMPLATE_STATS,
    PROVENANCE_FILENAME,
    STAGING_DB,
    EDIT_INDEX_DB,
    STRUCTURAL_INDEX_DB,
)
from .structural_keys import build_structural_keys
from .templates import build_template_library
from .validate import validate_library_v5


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _stage_done(outdir: Path, stats_name: str, required: list[str]) -> bool:
    return (outdir / stats_name).exists() and all((outdir / name).exists() for name in required)


def build_all(
    *,
    membership_source: str | Path,
    output_dir: str | Path,
    template_source: str | Path | None = None,
    workers: int | None = None,
    batch_size: int = 8192,
    row_group_size: int = 4096,
    chunksize: int = 32,
    radius: int = 2,
    fp_size: int = 2048,
    start_method: str = "fork",
    validation_sample: int = 10000,
    validation_seed: int = 6505,
    overwrite: bool = False,
    resume: bool = False,
    hash_sources: bool = True,
) -> dict[str, Any]:
    if overwrite and resume:
        raise ValueError("--overwrite and --resume are mutually exclusive")
    source = Path(membership_source)
    if not source.exists():
        raise FileNotFoundError(source)
    tsource = Path(template_source) if template_source else None
    if tsource is not None and not tsource.exists():
        raise FileNotFoundError(tsource)
    outdir = ensure_dir(output_dir)

    results: dict[str, Any] = {}

    print("\n=== STAGE 1/7: INGEST + NORMALIZE SOURCE FACTS ===", flush=True)
    if resume and _stage_done(outdir, SOURCE_NORMALIZE_STATS, [REACTION_CORE_BASE_FILENAME, MEMBERSHIP_REL_FILENAME, PROVENANCE_FILENAME, STAGING_DB]):
        print("[resume] Stage 1 already complete", flush=True)
        results["normalize"] = _load_json(outdir / SOURCE_NORMALIZE_STATS)
    else:
        results["normalize"] = normalize_source_membership(
            source, outdir, batch_size=batch_size, overwrite=overwrite
        )

    print("\n=== STAGE 2/7: DERIVE RUNTIME-v6.5 STRUCTURAL KEYS ===", flush=True)
    if resume and _stage_done(outdir, STRUCTURAL_STATS, [REACTION_CORE_FILENAME, STRUCTURAL_KEYS_FILENAME]):
        print("[resume] Stage 2 already complete", flush=True)
        results["structural"] = _load_json(outdir / STRUCTURAL_STATS)
    else:
        results["structural"] = build_structural_keys(
            outdir,
            workers=workers,
            batch_size=batch_size,
            row_group_size=row_group_size,
            chunksize=chunksize,
            start_method=start_method,
            overwrite=overwrite,
        )

    print("\n=== STAGE 3/7: BUILD SYSTEM-I REPRESENTATIONS ===", flush=True)
    if resume and _stage_done(outdir, "representation_stats.json", [REPRESENTATIONS_FILENAME, REPRESENTATION_LOCATOR_DB, REPRESENTATION_LOCATOR_PARQUET, REPRESENTATION_ERRORS_FILENAME]):
        print("[resume] Stage 3 already complete", flush=True)
        results["representations"] = _load_json(outdir / "representation_stats.json")
    else:
        results["representations"] = build_system_i_representations(
            outdir,
            workers=workers,
            batch_size=batch_size,
            row_group_size=row_group_size,
            chunksize=chunksize,
            radius=radius,
            fp_size=fp_size,
            start_method=start_method,
            overwrite=overwrite,
        )

    print("\n=== STAGE 4/7: BUILD RETRIEVAL MEMBERS ===", flush=True)
    if resume and (outdir / RETRIEVAL_MEMBERS_FILENAME).exists() and (outdir / "retrieval_members_stats.json").exists():
        print("[resume] Stage 4 already complete", flush=True)
        results["retrieval_members"] = _load_json(outdir / "retrieval_members_stats.json")
    else:
        results["retrieval_members"] = build_retrieval_members(outdir, batch_size=batch_size, overwrite=overwrite)

    print("\n=== STAGE 5/7: BUILD TEMPLATE LIBRARY ===", flush=True)
    if resume and _stage_done(outdir, TEMPLATE_STATS, [TEMPLATE_LIBRARY_FILENAME]):
        print("[resume] Stage 5 already complete", flush=True)
        results["templates"] = _load_json(outdir / TEMPLATE_STATS)
    else:
        results["templates"] = build_template_library(
            outdir,
            template_source=tsource,
            batch_size=batch_size,
            row_group_size=row_group_size,
            overwrite=overwrite,
        )

    print("\n=== STAGE 6/7: BUILD STRUCTURAL INDEXES ===", flush=True)
    if resume and _stage_done(outdir, "edit_index_stats.json", [EDIT_INDEX_FILENAME, EDIT_INDEX_DB, STRUCTURAL_INDEX_DB]):
        print("[resume] Stage 6 already complete", flush=True)
        results["indexes"] = _load_json(outdir / "edit_index_stats.json")
    else:
        results["indexes"] = build_indexes(outdir, batch_size=batch_size, overwrite=overwrite)

    print("\n=== STAGE 7/7: VALIDATE ===", flush=True)
    validation = validate_library_v5(outdir, sample_size=validation_sample, seed=validation_seed)
    results["validation"] = validation

    config = {
        "membership_source": str(source),
        "template_source": None if tsource is None else str(tsource),
        "workers": workers,
        "batch_size": batch_size,
        "row_group_size": row_group_size,
        "chunksize": chunksize,
        "radius": radius,
        "fp_size": fp_size,
        "start_method": start_method,
        "validation_sample": validation_sample,
        "validation_seed": validation_seed,
        "resume": resume,
    }
    write_manifest(
        outdir,
        membership_source=source,
        template_source=tsource,
        build_config=config,
        stage_results={k: v for k, v in results.items() if k != "validation"},
        validation=validation,
        hash_sources=hash_sources,
    )

    return {"status": "PASS", "output_dir": str(outdir), "stages": results}


def main() -> None:
    ap = argparse.ArgumentParser(description="Build native Runtime-v6.5 Template Library v5")
    ap.add_argument("--membership-source", required=True, help="Old source containing reaction + template membership facts")
    ap.add_argument("--template-source", help="Optional old template registry; only hash/number/SMARTS are used")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--workers", type=int)
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--row-group-size", type=int, default=4096)
    ap.add_argument("--chunksize", type=int, default=32)
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--fp-size", type=int, default=2048)
    ap.add_argument("--start-method", default="fork")
    ap.add_argument("--validation-sample", type=int, default=10000)
    ap.add_argument("--validation-seed", type=int, default=6505)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--no-source-sha256", action="store_true")
    args = ap.parse_args()
    result = build_all(
        membership_source=args.membership_source,
        template_source=args.template_source,
        output_dir=args.output_dir,
        workers=args.workers,
        batch_size=args.batch_size,
        row_group_size=args.row_group_size,
        chunksize=args.chunksize,
        radius=args.radius,
        fp_size=args.fp_size,
        start_method=args.start_method,
        validation_sample=args.validation_sample,
        validation_seed=args.validation_seed,
        overwrite=args.overwrite,
        resume=args.resume,
        hash_sources=not args.no_source_sha256,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
