#!/usr/bin/env python3
"""Document, verify, and ZIP a FINALIZED RXNI USPTO P6/Ring031 library.

Does NOT modify the scientific payload or retrieval manifest. Does NOT approve publication.
Run AFTER finalize_uspto_historical.py using the desired output folder name.
Requires the RXNI v0.1.0 Python environment (pyarrow and rxni installed).
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import sys
import zipfile

import pyarrow.parquet as pq

PROJECT = "https://github.com/shekfeh/rxni"
DATA = {
    "structural_index.sqlite": "Indexed structural features for reaction candidate generation; reaction_key joins other runtime files.",
    "reaction_lookup.sqlite": "Per-reaction metadata, source reaction SMILES, template assignment, observation metadata.",
    "reaction_representations.parquet": "Stored P6/Ring031 semantic, reaction-centre and fingerprint representations for 812,916 reactions.",
    "representation_locator.sqlite": "Lookup of reaction_key to Parquet row group/row offset, guarded by build identity.",
    "retrieval_manifest.json": "Storage generation, build identity, cryptographic provenance, frozen retrieval scientific contract.",
    "template_catalog.parquet": "42,128 historical retrosynthetic SMARTS and per-template structural consensus/metadata.",
    "rxni-library.json": "RXNI library descriptor mapping the public filenames to the five mandatory runtime artifacts.",
    "SCIENTIFIC_PROVENANCE.json": "Derived-data scientific lineage, source hashes, scoped acceptance and outstanding scientific review queue.",
    "RELEASE_AUDIT.json": "Construction and technical validation results plus explicitly pending release gates.",
    "README.md": "Human-readable dataset description, complete file/column dictionary and usage examples.",
    "NESTED_JSON_FIELDS.md": "Observed nested JSON member paths in the stored representation payloads (numeric keys normalized).",
    "SCHEMA_INVENTORY.json": "Machine-readable field inventory extracted from the actual Parquet, SQLite and JSON files.",
    "SHA256SUMS": "SHA256 digest of every dataset file except SHA256SUMS itself, to verify after download.",
}
DBS = {
    "structural_index.sqlite": "structural_index",
    "reaction_lookup.sqlite": "member_lookup",
    "representation_locator.sqlite": "representation_locator_v2",
}
PARQUETS = ("reaction_representations.parquet", "template_catalog.parquet")
REQUIRED = (
    "structural_index.sqlite", "reaction_lookup.sqlite", "reaction_representations.parquet",
    "representation_locator.sqlite", "retrieval_manifest.json", "template_catalog.parquet",
    "rxni-library.json", "SCIENTIFIC_PROVENANCE.json", "RELEASE_AUDIT.json",
)
HASH_COLUMNS = (
    "forward_edit_signature_hash", "forward_raw_center_hash_r0", "forward_raw_center_hash_r1",
    "forward_raw_center_hash_r2", "forward_component_topology_hash", "forward_component_evidence_hash",
)
FIELD_TEXT = {
    "reaction_key": "Opaque stable key identifying the source reaction; exact join key across the index, lookup, locator and representation store.",
    "forward_reaction_smiles": "Mapped forward reaction SMILES (reactants>agents>products) preserved for retrieval metadata.",
    "reaction": "Atom-mapped reaction SMILES represented by this stored semantic record.",
    "reaction_id": "Original source-reaction identifier where available; source-specific and not necessarily globally unique.",
    "record_id": "Stored record identifier; for this dataset typically identical to reaction_key.",
    "retro_template_hash": "Original AiZynthTrain fingerprint-bit-difference template identifier, NOT SHA256(raw SMARTS).",
    "retro_template_number": "Optional template numeric/index label supplied by the source extraction pipeline.",
    "retro_template_smarts": "Original retrosynthetic reaction SMARTS, suitable for rxnutils/RDKit parsing; chemical applicability is not guaranteed.",
    "forward_edit_signature_hash": "Opaque hashed graph-edit signature for exact edit-channel candidate generation.",
    "forward_raw_center_hash_r0": "Hashed local reaction-centre structure at graph radius zero.",
    "forward_raw_center_hash_r1": "Hashed local reaction-centre structure at graph radius one.",
    "forward_raw_center_hash_r2": "Hashed local reaction-centre structure at graph radius two.",
    "forward_component_topology_hash": "Hashed graph-level topology of reaction components, used as a retrieval channel.",
    "forward_component_evidence_hash": "Hashed component-evidence representation, used as a retrieval channel.",
    "membership_occurrence_count": "Number of source membership occurrences represented by this lookup or template membership.",
    "n_members": "Number of distinct reaction/template members supporting this template row.",
    "representative_reaction_key": "Reaction key selected as the representative member of the template.",
    "representative_reaction_smiles": "Forward reaction SMILES selected as representative of the template.",
    "reaction_key_version": "Historical version label for reaction-key generation; preserve as stored.",
    "library_schema_version": "Version string for the native V5 library layout/representation bridge.",
    "structural_key_version": "Provenance/version string for the frozen structural key implementation.",
    "system_i_representation_version": "Legacy optional field from an earlier representation stage; all 812,916 values were observed null in the archived historical audit. Retained for physical ABI.",
    "semantic_profile_json": "Serialized JSON of canonical retrieval semantic axes, e.g. accepted primary events and other semantic tokens.",
    "semantics_json": "Serialized comprehensive canonical semantic analysis: events, graph evidence, detector decisions, diagnostics and validations.",
    "reaction_centers_json": "Serialized computed reaction centres and graph-edit neighbourhood evidence.",
    "localized_fingerprints_json": "Serialized local count fingerprints associated with graph/reaction centres.",
    "multichannel_fingerprint_json": "Serialized fingerprints and semantic channels used by the frozen retrieval ranker.",
    "representation_generation": "Physical storage ABI identifier; expected retrieval_representations_v2 despite clean external filenames.",
    "representation_build_id": "Opaque record build identity that must match the retrieval manifest and locator.",
    "build_id": "Opaque retrieval representation build identity; v66_ prefix is historical and does not determine the semantic authority.",
    "row_group_id": "Zero-based Parquet row-group number used by the locator.",
    "row_in_group": "Zero-based ordinal within the specified Parquet row group.",
    "index": "SQLite index or source template code label, if present.",
    "source_core_sha256": "SHA256 digest of original source reaction Core; source Core is not shipped in the minimal distribution.",
    "source_reaction_core_sha256": "SHA256 digest of original source reaction Core; source Core is not shipped.",
    "representations_sha256": "SHA256 digest of the stored reaction representations file.",
    "locator_sqlite_sha256": "SHA256 digest of the safely snapshotted locator SQLite file.",
    "rows": "Number of stored reaction representation records.",
    "row_groups": "Number of physical Parquet row groups.",
    "errors": "Number of representation build errors reported in the manifest.",
    "schema_version": "Version of this file's JSON schema or build manifest, not the public dataset release title.",
    "semantic_authority": "Scientific runtime/rule lineage that actually generated stored semantic profiles.",
    "retrieval_api_version": "Version/compatibility claim for the retrieval API used during materialization.",
    "scientific_contract": "Nested declaration of semantic authority, ranking lineage and structural compatibility constraints.",
    "distribution_metadata": "Additional metadata attached during portable packaging; does not change the represented chemistry.",
    "representation_format": "Physical stored-representation ABI. Must remain retrieval_representations_v2 for RXNI v0.1.0.",
    "library_id": "Human-readable library identifier used by the RXNI descriptor.",
    "semantic_engine_version": "Declared stored semantic-scientific lineage; not a guarantee of all-record classification correctness.",
    "artifacts": "Mapping from five logical retrieval asset roles to relative paths within the dataset directory.",
    "status": "Status code of the corresponding build, provenance or validation stage.",
    "source_pin": "Map of build input source files to SHA256 code digests pinned during construction.",
    "known_scientific_review_queue": "Known reaction keys requiring scientific adjudication; publication gate remains open if unresolved.",
    "publication_gates": "Declared unclosed release requirements: rights, chemical review and clean-install testing.",
    "counts": "Record counts measured by the corresponding audit.",
    "issues": "Population/key discrepancy counts from the corresponding audit.",
    "retrieval_smoke": "Bounded source-versus-package retrieval parity check; not independent chemistry validation.",
    "smarts_catalogue_audit": "RDKit SMARTS parse/Validate counts and warnings; not experimental applicability evidence.",
    "full_index_audit": "Exact key-population and Parquet locator-position verification report.",
    "file_sha256": "Per-artifact SHA256 digests for content integrity.",
    "retro_template_number": "Optional source template numbering; not interchangeable with template fingerprint hash.",
    "structural_builder": "Frozen structural-key producer compatibility lineage; separate from semantic classification.",
    "scoring": "Frozen fingerprint and semantic retrieval-scoring implementation lineage.",
    "candidate_core_membership_provenance": "Scientific source of candidate reactions and their original template associations.",
    "p6_scoped_acceptance": "P6 scoped scientific rule acceptance counts; not a whole-library accuracy claim.",
    "p6_installation_gate_report_sha256": "Digest of the accepted P6 installation scientific gate report.",
    "ring_v031_acceptance": "Ring v0.3.1 scoped scientific acceptance counts; not top-k retrieval validation.",
    "ring_v031_gate_report_sha256": "Digest of the Ring v0.3.1 scientific acceptance gate report.",
    "source_provenance_note": "Explanation of original build-code lineage and the semantic execution authority.",
    "v67_accepted_detector_sha256": "Digest of the accepted v6.7 non-FGI detector source file.",
    "v67_accepted_resolver_sha256": "Digest of the accepted v6.7 non-FGI resolver source file.",
    "retrieval_build_source_sha256": "Digest of the source retrieval representation build program.",
    "retrieval_record_source_sha256": "Digest of the source live retrieval-record construction program.",
    "physical_build_id_note": "Clarifies historical naming in opaque builder IDs.",
    "original_manifest_sha256": "SHA256 of the unmodified scientific-source retrieval manifest.",
    "source_retrieval_manifest_sha256": "SHA256 of the source Ring031 retrieval manifest before portable repackaging.",
    "source_truthful_provenance_sha256": "SHA256 of the original source semantic provenance receipt.",
    "source_template_catalogue_sha256": "SHA256 of the untouched source retro SMARTS catalogue.",
    "catalogue_sha256": "SHA256 of the included retro SMARTS catalogue.",
    "source_profile": "Runtime scientific scope of the selected source profiles.",
    "packaging_release": "Release label for this portable packaging operation.",
    "dataset_title": "Readable name for the dataset package; not an ABI identifier.",
    "locator_resnapshotted": "Records whether SQLite locator bytes were backed up and rehashed for the distribution.",
    "scientific_representation_unchanged": "Attests that the Parquet scientific representation payload was not altered during packaging.",
    "stored_representation_contract": "Frozen physical V2 persisted-reaction representation format.",
    "source_core_sha256": "Hash of the original source Core data; Core is not included in the release.",
    "source_paths_excluded": "Attests that local absolute development paths were removed from packaged provenance.",
    "scientific_source": "Original scientific candidate/build identifier selected as source.",
    "scientific_source_status": "Source candidate promotion/status, which may remain unapproved.",
    "accepted_semantic_authority": "V6.7/P6/Ring031 semantic source declared by the selected scientific manifest.",
    "semantic_source_pins_sha256": "Source code file/digest mapping preserved from the accepted scientific run.",
    "known_scientific_review_queue": "Reaction keys specifically flagged for further scientific inspection.",
    "published_scientific_validation": "Whether the source claims independent scientific acceptance; NOT_CLAIMED means no such certification.",
    "rights_clearance": "Whether derivative dataset redistribution rights have been attested.",
    "generated_utc": "UTC time when the packaging audit report was generated.",
    "rxni_library_compatibility": "Result of public RXNI's loader/descriptor structural compatibility inspection.",
    "retrieval_smoke": "Number and outcome of selected source/package retrieval-parity queries.",
    "historical_template_hash_replay": "Evidence status for rxnutils fingerprint-bit template-hash reproduction.",
    "clean_environment_acceptance": "Status of an independent fresh-install and data-only library acceptance check.",
    "corpus_redistribution_rights": "Status of the source and derivative redistribution-rights review.",
    "template_catalogue_smarts": "Number of original retrosynthetic reaction templates in the catalogue.",
    "source_templates": "Number of templates supplied to the native builder.",
    "built_templates": "Number of templates produced by the native builder.",
    "missing_or_changed_smarts": "Number of templates with missing or mutated source SMARTS after transformation.",
    "hashes_unmatched": "Number of source template hashes without corresponding native template catalogue rows.",
    "exact_source_string_equality": "Indicates whether original SMARTS text was byte-for-byte preserved.",
    "metadata_only_records": "Count of metadata lookup reactions excluded from the structural candidate population.",
    "incorrect_locator_coordinates": "Number of locator entries pointing to a different Parquet row/row-group than their reaction key.",
    "incorrect_locator_build": "Number of locator entries violating the expected representation generation/build ID.",
    "parquet_rows_position_verified": "Number of Parquet records whose locator positions were checked exhaustively.",
    "distinct_lookup_templates": "Number of unique template hashes referenced in the metadata lookup table.",
    "catalogue_coverage": "Number of lookup template hashes present with corresponding SMARTS entries.",
}

NESTED_FIELD_TEXT = {
    "semantic_profile_json": "Canonical semantic feature-token profile consumed by retrieval.",
    "semantics_json": "Complete canonical semantic-analysis data object.",
    "reaction_centers_json": "Localized atom/bond edit and reaction-centre representation.",
    "localized_fingerprints_json": "Hashed local reaction-centre feature maps and counts.",
    "multichannel_fingerprint_json": "Multi-axis fingerprint representation used by the frozen scorer.",
    "primary_events": "Accepted primary reaction-event labels.",
    "specific_events": "Specific mechanistic or transformation-family event labels.",
    "semantic_events": "Resolved transformation-event inventory and supporting evidence.",
    "directional_conversions": "Direction-sensitive functional group transformation evidence.",
    "local_events": "Local graph-edit or reaction-centre event data.",
    "primitives": "Low-level graph-change primitives derived from the mapped reaction.",
    "rejected_hypotheses": "Candidate detector interpretations rejected during semantic resolution.",
    "raw_hypotheses": "Unresolved candidates proposed by graph-event detectors.",
    "validations": "Semantic and structural validation results collected by the analysis runtime.",
    "mapping_qc": "Atom-mapping validation and quality-control diagnostics.",
    "graph_summary": "Summary of atom/bond graph changes and component structure.",
    "event_topology": "Graph topology associated with transformation events.",
    "topology_evidence": "Supporting topology-level graph evidence for classification.",
    "component_provenance": "Origin/ownership metadata for reaction components and graph events.",
    "diagnostics": "Non-authoritative warnings and diagnostic outputs from the pipeline.",
    "engine_version": "Semantic-engine implementation version recorded in the individual record.",
    "rules_version": "Semantic rule-set version recorded in the individual record.",
    "redox_axis": "Reaction redox-related semantic evidence and axes.",
    "formed_bond_classes": "Categorized bond types formed in the mapped reaction.",
}


def nested_explanation(path: str) -> str:
    parts = path.replace("[]", "").split(".")
    for part in reversed(parts):
        if part in NESTED_FIELD_TEXT:
            return NESTED_FIELD_TEXT[part]
    return "Observed JSON member under the documented parent; exact scientific semantics are determined by the pinned RXNI runtime."

for base in HASH_COLUMNS:
    FIELD_TEXT[base + "_consensus_fraction"] = "Share of nonnull template members agreeing with the selected modal value of " + base + "."
    FIELD_TEXT[base + "_n_distinct"] = "Number of distinct per-reaction " + base + " values among members of this template."


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def md_cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ").replace("`", "'")


def type_label(field: object) -> str:
    return str(field.type) + (" (nullable)" if field.nullable else " (required)")


def ro(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)


def local_fields_for_parquet(path: Path) -> dict:
    pf = pq.ParquetFile(path)
    return {
        "kind": "parquet", "row_count": pf.metadata.num_rows, "row_groups": pf.num_row_groups,
        "fields": [{"name": f.name, "type": str(f.type), "nullable": bool(f.nullable),
                    "description": FIELD_TEXT.get(f.name, "Physical data field retained from the historical build; interpret using the producing source schema.")}
                   for f in pf.schema_arrow]
    }


def local_fields_for_sqlite(path: Path, table: str) -> dict:
    cx = ro(path)
    try:
        tables = [name for (name,) in cx.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name")]
        out = {}
        for name in tables:
            rows = cx.execute(f'PRAGMA table_info("{name}")').fetchall()
            entries = [{"name": c[1], "type": c[2], "not_null": bool(c[3]), "primary_key_position": c[5],
                        "description": FIELD_TEXT.get(c[1], "Historical SQLite field retained for source metadata or retrieval bookkeeping.")}
                       for c in rows]
            index_rows = cx.execute(f'PRAGMA index_list("{name}")').fetchall()
            out[name] = {"fields": entries,
                         "indexes": [{"name": x[1], "unique": bool(x[2])} for x in index_rows],
                         "rows": int(cx.execute(f'SELECT COUNT(*) FROM "{name}"').fetchone()[0])}
        return {"kind": "sqlite", "primary_runtime_table": table, "tables": out,
                "integrity": cx.execute("PRAGMA quick_check").fetchone()[0]}
    finally:
        cx.close()


def normalized_key(key: str) -> str:
    if (re.fullmatch(r"\d+", key) or re.fullmatch(r"[0-9a-f]{20,}", key, re.I)
            or len(key) > 100 or key.count("[") + key.count("(") > 2):
        return "{dynamic_key}"
    return key


def json_paths(value: object, prefix: str, collector: dict, max_depth: int = 24) -> None:
    # Accumulate *types* and occurrences; avoid recording user/source values.
    typ = "object" if isinstance(value, dict) else "array" if isinstance(value, list) else (
        "null" if value is None else "boolean" if isinstance(value, bool) else
        "number" if isinstance(value, (int, float)) else "string")
    collector[prefix][typ] += 1
    if prefix.count(".") + prefix.count("[]") >= max_depth:
        return
    if isinstance(value, dict):
        for k, v in value.items():
            if len(collector) > 25000:
                raise ValueError("Nested field catalogue exceeded 25,000 distinct paths; investigate dynamic identifiers")
            json_paths(v, prefix + "." + normalized_key(str(k)), collector, max_depth)
    elif isinstance(value, list) and value:
        # Inspect all elements of heterogeneous arrays for schema union.
        for item in value:
            json_paths(item, prefix + "[]", collector, max_depth)


def collect_nested(pf: pq.ParquetFile, mode: str, sample: int) -> tuple[dict, dict]:
    cols = [x for x in ("semantic_profile_json", "semantics_json", "reaction_centers_json",
                        "localized_fingerprints_json", "multichannel_fingerprint_json") if x in pf.schema_arrow.names]
    data = defaultdict(Counter)
    total = 0
    invalid = 0
    if mode == "none":
        return {}, {"scan_mode": "none", "records_scanned": 0, "invalid_json_values": 0, "columns": cols}
    # Stream through Parquet; never read all records into RAM.
    for batch in pf.iter_batches(batch_size=256, columns=cols):
        d = batch.to_pydict()
        for i in range(batch.num_rows):
            for col in cols:
                raw = d[col][i]
                try:
                    obj = json.loads(raw) if isinstance(raw, str) else raw
                except (ValueError, TypeError):
                    invalid += 1
                    continue
                json_paths(obj, col, data)
            total += 1
            if mode == "sample" and total >= sample:
                break
        if mode == "sample" and total >= sample:
            break
        if mode == "all" and total % 100000 == 0:
            print(f"Nested JSON schema: examined {total:,} records", flush=True)
    return {k: dict(v) for k, v in sorted(data.items())}, {
        "scan_mode": mode, "records_scanned": total, "invalid_json_values": invalid, "columns": cols,
        "dynamic_keys": "Numeric and long hexadecimal keys normalized to {dynamic_key}",
        "max_path_depth": 24,
    }


def inventory_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf8"))
    paths = defaultdict(Counter)
    json_paths(obj, path.name, paths)
    return {"kind": "json", "fields": [
        {"path": key, "types": dict(types), "description": FIELD_TEXT.get(key.split(".")[-1],
         "Historical provenance, validation or configuration field; inspect its parent JSON context.")}
        for key, types in sorted(paths.items())]}


def validate_existing(root: Path) -> tuple[dict, dict]:
    from rxni.retrieval.libraries import inspect_library, resolve_library
    missing = [p for p in REQUIRED if not (root / p).is_file()]
    if missing:
        raise ValueError(f"Not a finalized six-file RXNI distribution: missing {missing}; run finalize_uspto_historical.py first")
    lib = inspect_library(root)
    if lib.get("status") != "V2_READABLE" or lib.get("layout") != "descriptor-v1":
        raise ValueError(f"RXNI cannot resolve renamed files: {lib}")
    paths = resolve_library(root)
    mf = json.loads(paths.retrieval_manifest.read_text())
    if mf.get("representation_generation") != "retrieval_representations_v2":
        raise ValueError("Physical representation ABI changed unexpectedly")
    sem = str((mf.get("scientific_contract") or {}).get("semantic_authority") or "")
    if "runtime_v6_7.semantic" not in sem or "Ring" not in sem:
        raise ValueError(f"Not the selected P6/Ring031 scientific library: {sem!r}")
    if sha256(paths.representations) != mf.get("representations_sha256"):
        raise ValueError("Representation Parquet fails original scientific manifest hash")
    if sha256(paths.locator) != mf.get("locator_sqlite_sha256"):
        raise ValueError("SQLite locator checksum and manifest do not agree")
    pf = pq.ParquetFile(paths.representations)
    if pf.metadata.num_rows != mf.get("rows") or pf.metadata.num_row_groups != mf.get("row_groups"):
        raise ValueError("Parquet row/row-group population differs from scientific manifest")
    if pf.metadata.num_rows != 812916:
        raise ValueError(f"Expected 812,916 reactions; got {pf.metadata.num_rows}")
    cat = pq.ParquetFile(root / "template_catalog.parquet")
    if cat.metadata.num_rows != 42128:
        raise ValueError(f"Expected 42,128 template rows; got {cat.metadata.num_rows}")
    audit = json.loads((root / "RELEASE_AUDIT.json").read_text())
    if audit.get("status") != "TECHNICAL_CANDIDATE_PASS_PENDING_SCIENTIFIC_AND_RIGHTS_REVIEW":
        raise ValueError("Missing expected technical candidate audit status")
    for path in DBS:
        con = ro(root / path)
        try:
            if con.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                raise ValueError(f"SQLite integrity failed: {path}")
        finally:
            con.close()
    # Existing hashes must hold before changing documentation.
    old = root / "SHA256SUMS"
    if old.exists():
        for line in old.read_text().splitlines():
            if not line.strip():
                continue
            expected, sep, name = line.partition("  ")
            if sep != "  ":
                raise ValueError(f"Malformed existing SHA256SUMS line: {line!r}")
            if name != "README.md" and sha256(root / name) != expected:
                raise ValueError(f"Old SHA256SUMS mismatch for {name}")
    return mf, audit


def select_dataset_files(root: Path) -> tuple[list[Path], list[str]]:
    """Keep transient SQLite WAL/SHM files out of the published archive.

    Fail closed if a WAL has frames: the base .sqlite file might not contain
    committed data, so it would be unsafe to archive it without a checkpoint.
    Do not unlink sidecars: that could interfere with live SQLite connections.
    """
    excluded = []
    unexpected = []
    for p in root.iterdir():
        if p.is_dir() or p.name.startswith("."):
            unexpected.append(p.name)
            continue
        if p.name.endswith(("-wal", "-shm")):
            base, extension = p.name.rsplit("-", 1)
            if base not in DBS:
                unexpected.append(p.name)
                continue
            if extension == "wal" and p.stat().st_size > 32:
                raise ValueError(
                    f"SQLite WAL contains possible uncheckpointed frames: {p} "
                    f"({p.stat().st_size} bytes). Do not delete it. Stop other "
                    "readers/writers, checkpoint it safely, and re-run final "
                    "SQLite integrity and checksum validation."
                )
            excluded.append(p.name)
    if unexpected:
        raise ValueError(f"Refusing nested/temporary/hidden dataset contents: {unexpected}")
    included = sorted(
        (p for p in root.iterdir() if p.is_file() and p.name not in excluded),
        key=lambda p: p.name,
    )
    return included, excluded


def render_readme(root: Path, inventory: dict, manifest: dict, audit: dict) -> str:
    lines = [
        "# RXNI USPTO 1976–2016 Reaction Library",
        "",
        "**Dataset format:** RXNI five-asset retrieval ABI + historical template catalogue and audit metadata.",
        "**Source period:** USPTO patent-reaction corpus conventionally identified as 1976–2016; this release is a processed subset. Per-reaction publication dates are not included in the minimal runtime files, so the year span is not independently verified by this archive.",
        f"**Software:** [RXNI (`shekfeh/rxni`)]({PROJECT}); tested public interface: RXNI 0.1.0.",
        "**Stored semantic authority:** `runtime_v6_7.semantic`, after accepted P6 and Ring v0.3.1 scoped patches.",
        "**Storage ABI:** `retrieval_representations_v2` (an internal format identifier, NOT the public dataset title).",
        f"**Opaque build ID:** `{manifest['build_id']}` (the `v66_` prefix is not the semantic authority).",
        f"**Reaction representations:** {manifest['rows']:,}; **row groups:** {manifest['row_groups']:,}; **template SMARTS:** 42,128.",
        "",
        "> **Release status:** Technically prepared, NOT automatically approved for public publication. The source provenance records an unpromoted scientific candidate and a pending chemical review case. Independently confirm corpus/derivative redistribution permissions before making a Zenodo record public. RXNI's Apache-2.0 software license does not license these dataset contents.",
        "",
        "## 1. Contents and purpose of every file",
        "",
        "| File | Description | Size (bytes) |",
        "|---|---|---:|",
    ]
    for name in sorted(DATA):
        if name == "SHA256SUMS":
            size = "—"
        elif name == "README.md":
            size = "— (this document)"
        elif (root / name).is_file():
            size = f"{(root / name).stat().st_size:,}"
        else:
            size = "not yet generated"
        lines.append(f"| `{name}` | {md_cell(DATA.get(name, 'Dataset release support file.'))} | {size} |")
    lines += ["",
              "## 2. Retrieval data flow", "",
              "`reaction SMILES → RXNI v6.7 query semantics / graph-edit keys → structural_index.sqlite → reaction_lookup.sqlite → representation_locator.sqlite → reaction_representations.parquet → frozen multi-channel ranker → ranked reactions/templates`.", "",
              "The separate `template_catalog.parquet` joins by `retro_template_hash` for retro SMARTS inspection; it is not automatically loaded into every retrieval response.", "",
              "For details of structural key derivation, fingerprints and the ranking contract, see the [RXNI repository documentation](" + PROJECT + "/tree/main/docs).", "",
              "## 3. Exact physical schemas (inspected from this distribution)", ""]
    for name in (*DBS, *PARQUETS):
        obj = inventory["files"][name]
        lines += [f"### `{name}`", ""]
        if obj["kind"] == "sqlite":
            lines += [f"SQLite integrity: `{obj['integrity']}`. The runtime's primary table is `{obj['primary_runtime_table']}`.", ""]
            for table, detail in obj["tables"].items():
                lines += [f"**Table `{table}`:** {detail['rows']:,} rows.", "", "| Column | SQLite type | Nullability | Purpose |", "|---|---|---|---|"]
                for field in detail["fields"]:
                    nullability = "NOT NULL" if field["not_null"] or field["primary_key_position"] else "nullable"
                    lines.append(f"| `{field['name']}` | `{field['type']}` | {nullability} | {md_cell(field['description'])} |")
                lines.extend(["", "**Indexes:** " + ", ".join("`" + x["name"] + "`" for x in detail["indexes"]), ""])
        else:
            lines += [f"Parquet rows: **{obj['row_count']:,}**. Row groups: **{obj['row_groups']:,}**.", "",
                      "| Column | Arrow physical type | Purpose |", "|---|---|---|"]
            for field in obj["fields"]:
                nullable = " (nullable)" if field["nullable"] else " (required)"
                lines.append(f"| `{field['name']}` | `{field['type']}`{nullable} | {md_cell(field['description'])} |")
            lines.append("")
    lines += [
        "## 4. JSON metadata file dictionaries", "",
        "The tables enumerate every observed JSON object field path in each included metadata file. Arrays are denoted `[]`; dynamic numeric/long hex keys are shown as `{dynamic_key}`. `SCHEMA_INVENTORY.json` contains the same machine-readable inventory.", "",
    ]
    for name in ("rxni-library.json", "retrieval_manifest.json", "SCIENTIFIC_PROVENANCE.json", "RELEASE_AUDIT.json"):
        fields = inventory["files"][name]["fields"]
        lines += [f"### `{name}`", "", "| Path | JSON type(s) | Explanation |", "|---|---|---|"]
        for field in fields:
            lines.append(f"| `{field['path']}` | `{','.join(field['types'])}` | {md_cell(field['description'])} |")
        lines.append("")
    lines += [
        "## 5. Nested scientific JSON stored in Parquet", "",
        "Five Parquet columns contain JSON-serialized scientific objects. The separate [`NESTED_JSON_FIELDS.md`](NESTED_JSON_FIELDS.md) documents observed paths and value types. Important top-level families include:", "",
        "- `semantics_json`: canonical primary/specific events, primitives, mapping QC, graph evidence, hypothesis acceptance/rejection, engine/rules version and validations.",
        "- `semantic_profile_json`: query-index semantic axes and event tokens used for retrieval matching.",
        "- `reaction_centers_json`: changed bonds and localized graph features around the reaction centre.",
        "- `localized_fingerprints_json`: reaction-centre-local hashed feature/count dictionaries.",
        "- `multichannel_fingerprint_json`: multiple structural and semantic fingerprint channels used in ranking.", "",
        "Do not interpret fingerprint bits as unique chemical substructures without reference to the producing algorithm. `system_i_representation_version` is a retained legacy schema column, not the current semantic engine.", "",
        "## 6. How to install and query with RXNI", "",
        "```bash", "git clone https://github.com/shekfeh/rxni.git", "cd rxni", "python3 -m venv .venv", "./.venv/bin/python -m pip install -e .", "export REACTION_INTELLIGENCE_LIBRARY=/absolute/path/to/uspto_1976_2016_rxni_library", "./.venv/bin/rxni --version", "./.venv/bin/rxni retrieve --help", "```", "",
        "Choose an actual reaction from the library, to ensure a meaningful smoke test:", "",
        "```bash", "./.venv/bin/python - <<'PY'", "from pathlib import Path", "import pyarrow.parquet as pq", "from rxni.retrieval.libraries import inspect_library", "root = Path('/absolute/path/to/uspto_1976_2016_rxni_library')", "print(inspect_library(root))", "pf = pq.ParquetFile(root/'reaction_representations.parquet')", "print('Example reaction:', pf.read_row_group(0, columns=['reaction'])['reaction'][0].as_py())", "PY", "```", "",
        "Then pass the example to `rxni retrieve --reaction 'REACTANTS>AGENTS>PRODUCTS' --library $REACTION_INTELLIGENCE_LIBRARY --top-k 5 --candidate-budget 500 --format json`.", "",
        "The repository provides [Tutorial 05]("+PROJECT+"/blob/main/docs/tutorials/05-download-and-use-library.md), [Tutorial 06]("+PROJECT+"/blob/main/docs/tutorials/06-build-library-from-sources.md) and [Tutorial 07]("+PROJECT+"/blob/main/docs/tutorials/07-template-catalogue-and-provenance.md).", "",
        "**Inspect the optional template catalogue:**", "",
        "```python", "import pyarrow.parquet as pq", "from pathlib import Path", "library = Path('/absolute/path/to/uspto_1976_2016_rxni_library')", "tab = pq.read_table(library/'template_catalog.parquet', columns=['retro_template_hash','retro_template_smarts'])", "print(tab.slice(0, 3).to_pylist())", "```", "",
        "## 7. Validate after downloading", "",
        "```bash", "cd /absolute/path/to/uspto_1976_2016_rxni_library", "sha256sum -c SHA256SUMS", "python /path/to/rxni/scripts/release/verify_published_library.py --library . --smoke", "```", "",
        "The descriptor `rxni-library.json` resolves the five required runtime files by their public names. Do **not** rename internal SQLite tables, Parquet columns, the `representation_generation` string, or the build ID.", "",
        "## 8. Dataset provenance and limitations", "",
        "- Dataset is derived from a historical USPTO patent-reaction extraction and the RXNI native V5 reaction foundation. The exact source publication span and redistribution rights must be confirmed from upstream source manifests before making the Zenodo deposit public.",
        "- Source semantic records were reprofiled under `runtime_v6_7.semantic`, post-P6/Ring v0.3.1. The structural features and score equations retain earlier, frozen compatibility contracts.",
        "- This is a **reaction retrieval** library, not a trained retrosynthesis predictor, mechanistic step generator or experimentally validated condition recommender.",
        "- The historical metadata lookup contains one extra metadata-only reaction key outside the 812,916 indexed reactions. It is unreachable through structural-index candidate generation. It was preserved for provenance.",
        "- All 42,128 retrospective SMARTS parsed with RDKit; the historical source replayed 42,128/42,128 fingerprint-derived template hashes with rxnutils and RDKit 2023.09.6. This does not establish that all templates apply successfully to their corresponding reaction products.",
        "- The Ring031 source is a scientific candidate with a known unresolved review case and no independently attested whole-library retrieval accuracy. Technical parity is not a claim of chemical accuracy.",
        "- Source reaction IDs and reaction chemistry may be covered by independent data licensing or database rights. No patent-source ownership is transferred by the RXNI Apache-2.0 code license.", "",
        "**Links:** [RXNI source code]({}) · [RXNI distribution guide]({}/blob/main/docs/library-distribution.md) · [Zenodo](https://zenodo.org/)".format(PROJECT, PROJECT), "",
        "**Citation:** Insert the final Zenodo dataset DOI and creator metadata after an authorized publication. No DOI has been allocated in this archive.", "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", type=Path, required=True, help="Already finalized, six-data-file RXNI directory")
    ap.add_argument("--archive", type=Path, required=True, help="Destination .zip, outside the dataset folder")
    ap.add_argument("--nested-json", choices=("all", "sample", "none"), default="all",
                    help="all: enumerate observed nested JSON keys across every stored reaction (slower but comprehensive)")
    ap.add_argument("--sample-rows", type=int, default=1024)
    ap.add_argument("--reuse-docs", action="store_true",
                    help="Reuse complete, validated generated schema documentation from a previous interrupted packaging attempt")
    ap.add_argument("--no-zip-crc", action="store_true", help="Skip post-write complete ZIP CRC verification")
    args = ap.parse_args()
    root = args.library.expanduser().resolve()
    archive = args.archive.expanduser().resolve()
    if root.name != "uspto_1976_2016_rxni_library":
        ap.error("Expected exact final library directory name uspto_1976_2016_rxni_library")
    if not root.is_dir():
        ap.error("Library directory not found")
    if archive.suffix.lower() != ".zip" or archive.exists():
        ap.error("Archive path must end in .zip and must not already exist")
    if archive.is_relative_to(root):
        ap.error("Archive must be created OUTSIDE the data directory")
    if args.sample_rows <= 0:
        ap.error("--sample-rows must be >0")
    mf, audit = validate_existing(root)
    print("PASS: existing RXNI library, SHA256 and semantic authority", flush=True)
    if args.reuse_docs:
        required_docs = ("README.md", "NESTED_JSON_FIELDS.md", "SCHEMA_INVENTORY.json")
        missing_docs = [name for name in required_docs if not (root / name).is_file()]
        if missing_docs:
            raise ValueError(f"Cannot reuse missing documentation: {missing_docs}")
        inventory = json.loads((root / "SCHEMA_INVENTORY.json").read_text(encoding="utf8"))
        nested_status = inventory.get("nested_json", {})
        pf = pq.ParquetFile(root / "reaction_representations.parquet")
        if nested_status.get("scan_mode") != args.nested_json:
            raise ValueError("Existing JSON documentation scan mode does not match --nested-json")
        if args.nested_json == "all" and nested_status.get("records_scanned") != pf.metadata.num_rows:
            raise ValueError("Existing nested JSON inventory does not cover all Parquet rows")
        if nested_status.get("invalid_json_values") != 0:
            raise ValueError("Existing JSON inventory reports invalid JSON values")
        if inventory.get("dataset") != root.name or not inventory.get("files"):
            raise ValueError("Existing schema inventory has wrong dataset name or no field schemas")
        print(f"PASS: reused {nested_status['records_scanned']:,}-record JSON schema and generated documentation", flush=True)
    else:
        inventory = {"dataset": root.name, "repo": PROJECT, "files": {},
                     "generated_utc": datetime.now(timezone.utc).isoformat()}
        for filename, table in DBS.items():
            item = local_fields_for_sqlite(root / filename, table)
            if item["integrity"] != "ok":
                raise ValueError(f"SQLite quick_check failed for {filename}")
            inventory["files"][filename] = item
        for filename in PARQUETS:
            inventory["files"][filename] = local_fields_for_parquet(root / filename)
        for filename in ("retrieval_manifest.json", "rxni-library.json", "SCIENTIFIC_PROVENANCE.json", "RELEASE_AUDIT.json"):
            inventory["files"][filename] = inventory_json(root / filename)
        pf = pq.ParquetFile(root / "reaction_representations.parquet")
        nested, mode = collect_nested(pf, args.nested_json, args.sample_rows)
        inventory["nested_json"] = mode
        if mode["invalid_json_values"]:
            raise ValueError(f"Found {mode['invalid_json_values']} invalid stored JSON strings")
        details = ["# RXNI stored scientific JSON — observed nested field inventory", "",
                   "This inventory is derived by JSON parsing of the actual `reaction_representations.parquet` payloads.",
                   f"Scope: **{mode['scan_mode']}**, examined **{mode['records_scanned']:,}** records; **{len(nested):,}** distinct normalized paths.",
                   "Numeric and long hexadecimal JSON keys are collapsed to `{dynamic_key}`; this is an inspection convention only and does not alter the stored data.",
                   "Array elements use `[]`; heterogeneous list members are unioned. Deep paths are bounded at 24 structural levels.", "",
                   "| JSON member path | Types observed | Occurrences | Interpretation |", "|---|---|---:|---|"]
        for path, types in nested.items():
            details.append("| `" + md_cell(path) + "` | " + ", ".join(types.keys()) + " | " + str(sum(types.values())) + " | " + md_cell(nested_explanation(path)) + " |")
        (root / "NESTED_JSON_FIELDS.md").write_text("\n".join(details) + "\n", encoding="utf8")
        (root / "SCHEMA_INVENTORY.json").write_text(json.dumps(inventory, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf8")
        (root / "README.md").write_text(render_readme(root, inventory, mf, audit), encoding="utf8")
    # No local paths or private provenance receipts should leak through JSON/Markdown.
    for name in ("README.md", "SCHEMA_INVENTORY.json", "NESTED_JSON_FIELDS.md", "retrieval_manifest.json", "rxni-library.json", "RELEASE_AUDIT.json", "SCIENTIFIC_PROVENANCE.json"):
        text = (root / name).read_text(encoding="utf8")
        for marker in ("/home/suhaib/", "/media/suhaib/"):
            if marker in text:
                raise ValueError(f"Local user path found in publication candidate {name}; refuse archive")
    selected, sidecars = select_dataset_files(root)
    files = [p for p in selected if p.name != "SHA256SUMS"]
    if sidecars:
        print(f"INFO: excluding transient SQLite sidecars from the ZIP: {sidecars}", flush=True)
    (root / "SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in files), encoding="utf8")
    print("PASS: detailed README, Parquet/SQLite schema and nested JSON paths generated", flush=True)
    archive.parent.mkdir(parents=True, exist_ok=True)
    try:
        included, sidecars_now = select_dataset_files(root)
        if set(sidecars_now) != set(sidecars):
            raise ValueError("SQLite sidecars changed during packaging; retry with no other database users")
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1, allowZip64=True,
                             strict_timestamps=False) as zf:
            for file in included:
                arcname = root.name + "/" + file.name
                zf.write(file, arcname)
                print(f"ZIP: {arcname}", flush=True)
        with zipfile.ZipFile(archive, "r") as zf:
            names = zf.namelist()
            expected = [root.name + "/" + p.name for p in included]
            if names != expected:
                raise ValueError("ZIP members do not exactly match the selected release files")
            if not args.no_zip_crc:
                invalid = zf.testzip()
                if invalid:
                    raise ValueError(f"ZIP CRC failure in {invalid}")
        archive_hash = sha256(archive)
        checksum_file = archive.with_name(archive.name + ".sha256")
        checksum_file.write_text(f"{archive_hash}  {archive.name}\n", encoding="utf8")
        print("\nARCHIVE CREATED:", archive)
        print("BYTES:", archive.stat().st_size)
        print("SHA256:", archive_hash)
        print("Archive SHA256 file:", checksum_file)
        print("STATUS: TECHNICAL CANDIDATE ONLY. DO NOT PUBLISH UNTIL SCIENTIFIC/RIGHTS REVIEW CLOSED.")
    except BaseException:
        archive.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    main()
