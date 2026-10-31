"""Read-only RXNI library discovery and versioned artifact layout resolution.

No inference about chemical correctness or semantic profile equivalence is made
from filenames.  The read path currently supports the stored V2 codec only;
older V1-only sources require an explicit materialization/conversion step.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping

DESCRIPTOR_FILENAME = "rxni-library.json"
DESCRIPTOR_VERSION = 1
STORED_GENERATION = "retrieval_representations_v2"

ARTIFACT_KEYS = (
    "structural_index", "metadata_lookup", "representations", "locator", "retrieval_manifest"
)
LEGACY_V2_NAMES = {
    "structural_index": "structural_index_v1.sqlite",
    "metadata_lookup": "runtime_lookup_v1.sqlite",
    "representations": "retrieval_representations_v2.parquet",
    "locator": "representation_locator_v2.sqlite",
    "retrieval_manifest": "retrieval_v2_manifest.json",
}


class LibraryFormatError(ValueError):
    """The supplied library cannot be consumed by the V2 retrieval reader."""


@dataclass(frozen=True)
class LibraryArtifacts:
    root: Path
    structural_index: Path
    metadata_lookup: Path
    representations: Path
    locator: Path
    retrieval_manifest: Path
    layout: str
    foundation_profile: Path | None = None
    library_id: str | None = None
    profile_policy_claim: str | None = None
    semantic_engine_claim: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "root": str(self.root),
            "layout": self.layout,
            "library_id": self.library_id,
            "profile_policy_claim": self.profile_policy_claim,
            "semantic_engine_claim": self.semantic_engine_claim,
            "artifacts": {name: str(getattr(self, name)) for name in ARTIFACT_KEYS},
            **({"foundation_schema_profile": str(self.foundation_profile)} if self.foundation_profile else {}),
        }


def _read_json(path: Path) -> Mapping[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise LibraryFormatError(f"Invalid library JSON {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise LibraryFormatError(f"Expected a JSON object in {path}")
    return data


def _artifact_path(root: Path, name: str) -> Path:
    if not isinstance(name, str) or not name.strip():
        raise LibraryFormatError("Artifact paths must be nonempty relative strings")
    rel = Path(name)
    if rel.is_absolute() or any(part in ("..", ".") for part in rel.parts):
        raise LibraryFormatError(f"Artifact path must stay within the library: {name!r}")
    candidate = (root / rel).resolve()
    if not candidate.is_relative_to(root):
        raise LibraryFormatError(f"Artifact path escapes library root: {name!r}")
    return candidate


def _resolve(source: str | Path) -> LibraryArtifacts:
    source_path = Path(source).expanduser().resolve()
    if source_path.is_file():
        descriptor_path = source_path
        root = source_path.parent
    else:
        root = source_path
        descriptor_path = root / DESCRIPTOR_FILENAME
    if not root.is_dir():
        raise FileNotFoundError(f"Library directory does not exist: {root}")

    if descriptor_path.is_file():
        raw = _read_json(descriptor_path)
        if raw.get("schema_version") != DESCRIPTOR_VERSION:
            raise LibraryFormatError(
                f"Unsupported library descriptor schema_version={raw.get('schema_version')!r}; "
                f"expected {DESCRIPTOR_VERSION}"
            )
        if raw.get("representation_format") != STORED_GENERATION:
            raise LibraryFormatError(
                "Library does not declare stored V2 representations; "
                "convert or materialize its profiles before retrieval"
            )
        artifacts = raw.get("artifacts")
        if not isinstance(artifacts, dict) or set(artifacts) != set(ARTIFACT_KEYS):
            raise LibraryFormatError(f"Descriptor artifacts must contain exactly: {ARTIFACT_KEYS!r}")
        paths = {name: _artifact_path(root, artifacts[name]) for name in ARTIFACT_KEYS}
        def optional_text(key):
            value = raw.get(key)
            if value is None:
                return None
            if not isinstance(value, str):
                raise LibraryFormatError(f"Descriptor {key} must be a string when supplied")
            return value
        profile_rel = raw.get('foundation_schema_profile')
        profile_path = None
        if profile_rel is not None:
            profile_path = _artifact_path(root, profile_rel)
            if not profile_path.is_file():
                raise LibraryFormatError(f'Foundation schema profile missing: {profile_path}')
            from .foundation_schema import load_foundation_profile, FoundationSchemaError
            try:
                load_foundation_profile(profile_path)
            except FoundationSchemaError as exc:
                raise LibraryFormatError(f'Invalid foundation schema profile: {exc}') from exc
        out = LibraryArtifacts(root=root, **paths, layout="descriptor-v1",
                               foundation_profile=profile_path,
                               library_id=optional_text("library_id"),
                               profile_policy_claim=optional_text("profile_policy"),
                               semantic_engine_claim=optional_text("semantic_engine_version"))
    else:
        if source_path.is_file():
            raise LibraryFormatError(f"Descriptor file not found: {source_path}")
        paths = {key: root / basename for key, basename in LEGACY_V2_NAMES.items()}
        out = LibraryArtifacts(root=root, **paths, layout="legacy-v2-directory")

    missing = [key for key in ARTIFACT_KEYS if not getattr(out, key).is_file()]
    if missing:
        if out.layout == "legacy-v2-directory" and (
            (root / "system_i_representations_v1.parquet").is_file()
            or (root / "representation_locator_v1.sqlite").is_file()
        ):
            raise LibraryFormatError(
                "V1-only stored profiles are not V2 retrieval records. "
                "A separate explicit profile-conversion/build step is required. "
                f"Missing V2 artifacts: {missing!r}"
            )
        if out.layout == "legacy-v2-directory" and (
            (root / "reactions.parquet").is_file()
            or (root / "reaction_core_v5.parquet").is_file()
        ):
            raise LibraryFormatError(
                "Reaction Core detected without a complete retrieval auxiliary. "
                f"Build the index and stored representations; missing: {missing!r}"
            )
        raise LibraryFormatError(f"Missing required retrieval artifacts: {missing!r} in {root}")

    if out.foundation_profile is not None:
        from .foundation_schema import (load_foundation_profile,
                                        validate_foundation_schema, FoundationSchemaError)
        try:
            validate_foundation_schema(out.structural_index, out.metadata_lookup,
                                       load_foundation_profile(out.foundation_profile))
        except FoundationSchemaError as exc:
            raise LibraryFormatError(f'Foundation SQLite schema mismatch: {exc}') from exc
    manifest = _read_json(out.retrieval_manifest)
    if manifest.get("representation_generation") != STORED_GENERATION:
        raise LibraryFormatError(
            f"Stored representation generation {manifest.get('representation_generation')!r} "
            f"does not match {STORED_GENERATION!r}"
        )
    if not manifest.get("build_id"):
        raise LibraryFormatError("Retrieval manifest missing build_id")
    return out


def resolve_library(source: str | Path) -> LibraryArtifacts:
    """Locate V2 artifacts; preserve existing directory behavior by default.

    An optional ``rxni-library.json`` descriptor supports arbitrary relative
    artifact names without changing the validated stored V2 codec or ranking.
    """
    return _resolve(source)


def inspect_library(source: str | Path) -> dict[str, Any]:
    """Inspect declared V2 compatibility; does not validate chemical accuracy."""
    try:
        out = _resolve(source)
    except (LibraryFormatError, FileNotFoundError) as exc:
        return {"status": "NOT_RETRIEVAL_READY", "reason": str(exc)}
    manifest = _read_json(out.retrieval_manifest)
    return {
        "status": "V2_READABLE",
        **out.as_dict(),
        "representation_generation": manifest["representation_generation"],
        "build_id": str(manifest["build_id"]),
        "retrieval_api_claim": manifest.get("retrieval_api_version"),
        "stored_semantic_claim": (manifest.get("scientific_contract") or {}).get("semantic_authority")
            if isinstance(manifest.get("scientific_contract"), dict) else None,
        "scientific_validation": "NOT_PERFORMED",
    }


__all__ = ["LibraryArtifacts", "LibraryFormatError", "inspect_library", "resolve_library"]
