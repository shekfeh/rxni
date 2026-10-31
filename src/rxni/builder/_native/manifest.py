"""Build manifest/provenance for native V5."""

from __future__ import annotations

import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import pyarrow

from .common import file_sha256, write_json
from .schema import (
    BUILD_MANIFEST_FILENAME,
    BUILDER_VERSION,
    LIBRARY_SCHEMA_VERSION,
    REACTION_KEY_VERSION,
    STRUCTURAL_KEY_VERSION,
    SYSTEM_I_REPRESENTATION_VERSION,
)
from .source_adapter import inspect_source


def _git_commit(project_root: Path) -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=project_root, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def _rdkit_version() -> str | None:
    try:
        import rdkit
        return getattr(rdkit, "__version__", None)
    except Exception:
        return None


def write_manifest(
    output_dir: str | Path,
    *,
    membership_source: str | Path,
    template_source: str | Path | None,
    build_config: dict[str, Any],
    stage_results: dict[str, Any],
    validation: dict[str, Any] | None,
    project_root: str | Path | None = None,
    hash_sources: bool = True,
) -> dict[str, Any]:
    outdir = Path(output_dir)
    membership = Path(membership_source)
    template = Path(template_source) if template_source else None
    root = Path(project_root) if project_root else Path.cwd()
    sources: list[dict[str, Any]] = []
    for kind, p in (("membership", membership), ("template", template)):
        if p is None:
            continue
        rec = {
            "kind": kind,
            "path": str(p.resolve()),
            "size_bytes": p.stat().st_size,
            "inspection": inspect_source(p, kind=kind),
        }
        if hash_sources:
            rec["sha256"] = file_sha256(p)
        sources.append(rec)

    manifest = {
        "builder_version": BUILDER_VERSION,
        "library_schema_version": LIBRARY_SCHEMA_VERSION,
        "reaction_key_version": REACTION_KEY_VERSION,
        "structural_key_version": STRUCTURAL_KEY_VERSION,
        "system_i_representation_version": SYSTEM_I_REPRESENTATION_VERSION,
        "sources": sources,
        "build_config": build_config,
        "stage_results": stage_results,
        "validation_status": None if validation is None else validation.get("status"),
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "pyarrow": pyarrow.__version__,
            "rdkit": _rdkit_version(),
            "git_commit": _git_commit(root),
        },
    }
    write_json(outdir / BUILD_MANIFEST_FILENAME, manifest)
    return manifest
