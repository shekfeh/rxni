#!/usr/bin/env python3
"""Repackage the validated P6/Ring031 RXNI candidate; NEVER modify science sources.

Example:
 python scripts/release/finalize_uspto_historical.py \
  --staged "$ROOT/runs/v680/library_distribution/uspto_historical_v67_ring031_r1" \
  --historical "$ROOT/artifacts/libraries/template_library_v5_native_v01" \
  --ring "$ROOT/runs/retrieval_v67_reprofile/v67_retrieval_v2_full_candidate_P6_Ring031_01" \
  --out "$ROOT/runs/v680/library_distribution/uspto_historical_r1_FINAL" \
  --smoke 3

IMPORTANT: Technical packaging is NOT publication authorization. Rights and
scientific review are external gates. Never run on the original library as out.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile

import pyarrow.parquet as pq
from rdkit import RDLogger
from rdkit.Chem import rdChemReactions
from rxni.retrieval.libraries import resolve_library, inspect_library

FORMAT = "retrieval_representations_v2"   # ABI; do not rename this value
DATA_NAMES = (
    "structural_index.sqlite", "reaction_lookup.sqlite",
    "reaction_representations.parquet", "representation_locator.sqlite",
    "retrieval_manifest.json", "template_catalog.parquet",
)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(p: Path, value: dict) -> None:
    p.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf8")


def ro_uri(path: Path) -> str:
    return path.resolve().as_uri() + "?mode=ro"


def snapshot_db(src: Path, dest: Path) -> None:
    source = sqlite3.connect(ro_uri(src), uri=True)
    target = sqlite3.connect(dest)
    try:
        source.backup(target)
    finally:
        target.close()
        source.close()
    cx = sqlite3.connect(ro_uri(dest), uri=True)
    try:
        status = cx.execute("PRAGMA integrity_check").fetchone()[0]
        if status != "ok":
            raise ValueError(f"SQLite integrity failed {dest.name}: {status}")
    finally:
        cx.close()


def audit_catalog(path: Path) -> tuple[dict, set[str]]:
    file = pq.ParquetFile(path)
    cols = set(file.schema_arrow.names)
    if not {"retro_template_hash", "retro_template_smarts"} <= cols:
        raise ValueError("Template catalogue lacks required hash/SMARTS fields")
    seen = {}
    counts = Counter()
    examples = []
    RDLogger.DisableLog("rdApp.warning")
    for batch in file.iter_batches(columns=["retro_template_hash", "retro_template_smarts"], batch_size=4096):
        for row in batch.to_pylist():
            counts["rows"] += 1
            k, smarts = row["retro_template_hash"], row["retro_template_smarts"]
            if not isinstance(k, str) or not k or not isinstance(smarts, str) or not smarts.strip():
                raise ValueError("Blank template hash/SMARTS in catalogue")
            if k in seen:
                raise ValueError(f"Repeated template hash: {k}")
            seen[k] = smarts
            try:
                r = rdChemReactions.ReactionFromSmarts(smarts)
                if r is None:
                    raise ValueError("RDKit returned None")
                warn, err = r.Validate(silent=True)
                counts["rdkit_warnings"] += int(warn)
                counts["rdkit_validation_errors"] += int(err)
                if err and len(examples) < 10:
                    examples.append({"template_hash": k, "errors": err})
                counts["rdkit_parsed"] += 1
            except Exception as ex:
                counts["rdkit_parse_failed"] += 1
                if len(examples) < 10:
                    examples.append({"template_hash": k, "error": str(ex)})
    report = {"counts": dict(counts), "validation_examples": examples,
              "scope": "RDKit parse/Validate; not empirical applicability or retrosynthetic success"}
    if counts["rdkit_parse_failed"] or counts["rdkit_validation_errors"]:
        raise ValueError("Template SMARTS parse/validation failed: see source audit")
    return report, set(seen)


def audit_full_index(out: Path, build_id: str, template_hashes: set[str]) -> dict:
    """Prove exact keys AND Parquet row positions across all 812,916 records."""
    index = out / "structural_index.sqlite"
    lookup = out / "reaction_lookup.sqlite"
    locator = out / "representation_locator.sqlite"
    parquet = out / "reaction_representations.parquet"
    mf = json.loads((out / "retrieval_manifest.json").read_text())
    pf = pq.ParquetFile(parquet)
    if mf.get("rows") != pf.metadata.num_rows:
        raise ValueError("Manifest/representation row count mismatch")
    if mf.get("row_groups") != pf.num_row_groups:
        raise ValueError("Manifest/representation row-group count mismatch")
    # Disk-backed temp table avoids loading 812k positions into Python RAM.
    with tempfile.TemporaryDirectory(prefix="rxni_key_audit_", dir=out.parent) as tmp:
        temp = Path(tmp) / "keys.sqlite"
        conn = sqlite3.connect(temp, uri=True)
        try:
            conn.execute("CREATE TABLE p(k TEXT PRIMARY KEY, rg INTEGER NOT NULL, ri INTEGER NOT NULL)")
            generated = 0
            for rg in range(pf.num_row_groups):
                tab = pf.read_row_group(
                    rg, columns=["reaction_key", "representation_build_id", "representation_generation"]
                )
                keys, builds, gens = (tab[c].to_pylist() for c in tab.column_names)
                if any(b != build_id for b in builds):
                    raise ValueError(f"Parquet representation build_id mismatch in row group {rg}")
                if any(g != FORMAT for g in gens):
                    raise ValueError(f"Parquet representation_generation mismatch in row group {rg}")
                conn.executemany("INSERT INTO p VALUES (?,?,?)",
                                 ((str(k), rg, ri) for ri, k in enumerate(keys)))
                generated += len(keys)
            conn.commit()
            def attach(alias: str, file: Path) -> None:
                # Read-only URI prevents accidental writes to the source DB.
                conn.execute(f"ATTACH DATABASE ? AS {alias}", (ro_uri(file),))
            attach("s", index)
            attach("m", lookup)
            attach("l", locator)
            tables = {"s": "structural_index", "m": "member_lookup", "l": "representation_locator_v2"}
            counts = {name: conn.execute(f"SELECT COUNT(*) FROM {schema}.{table}").fetchone()[0]
                      for name, (schema, table) in {
                          "structural": ("s", tables["s"]), "metadata": ("m", tables["m"]),
                          "locator": ("l", tables["l"]), "parquet_keys": ("main", "p")}.items()}
            expected = pf.metadata.num_rows
            for label in ("structural", "locator", "parquet_keys"):
                if counts[label] != expected:
                    raise ValueError(f"{label} population {counts[label]} != {expected}")
            queries = {
                "structural_missing_metadata": "SELECT COUNT(*) FROM s.structural_index s LEFT JOIN m.member_lookup m ON m.reaction_key=s.reaction_key WHERE m.reaction_key IS NULL",
                "structural_missing_locator": "SELECT COUNT(*) FROM s.structural_index s LEFT JOIN l.representation_locator_v2 l ON l.reaction_key=s.reaction_key WHERE l.reaction_key IS NULL",
                "locator_missing_structural": "SELECT COUNT(*) FROM l.representation_locator_v2 l LEFT JOIN s.structural_index s ON s.reaction_key=l.reaction_key WHERE s.reaction_key IS NULL",
                "locator_missing_parquet": "SELECT COUNT(*) FROM l.representation_locator_v2 l LEFT JOIN main.p p ON p.k=l.reaction_key WHERE p.k IS NULL",
                "parquet_missing_locator": "SELECT COUNT(*) FROM main.p p LEFT JOIN l.representation_locator_v2 l ON l.reaction_key=p.k WHERE l.reaction_key IS NULL",
                "incorrect_locator_coordinates": "SELECT COUNT(*) FROM main.p p JOIN l.representation_locator_v2 l ON l.reaction_key=p.k WHERE l.row_group_id!=p.rg OR l.row_in_group!=p.ri",
                "incorrect_locator_build": "SELECT COUNT(*) FROM l.representation_locator_v2 WHERE build_id!=? OR representation_generation!=?",
                "metadata_only_records": "SELECT COUNT(*) FROM m.member_lookup m LEFT JOIN s.structural_index s ON s.reaction_key=m.reaction_key WHERE s.reaction_key IS NULL",
            }
            issues = {}
            for key, query in queries.items():
                params = (build_id, FORMAT) if key == "incorrect_locator_build" else ()
                issues[key] = int(conn.execute(query, params).fetchone()[0])
            if any(v for k, v in issues.items() if k != "metadata_only_records"):
                raise ValueError(f"Key/locator inconsistency: {issues}")
            used = {str(x[0]) for x in conn.execute(
                "SELECT DISTINCT retro_template_hash FROM m.member_lookup WHERE retro_template_hash IS NOT NULL")}
            missing = used - template_hashes
            if missing:
                raise ValueError(f"Missing template SMARTS for {len(missing)} lookup hashes; first={sorted(missing)[:5]}")
            return {"counts": counts, "issues": issues, "distinct_lookup_templates": len(used),
                    "catalogue_coverage": len(used & template_hashes),
                    "parquet_rows_position_verified": generated}
        finally:
            conn.close()


def smoke_retrieval(stage: Path, final: Path, n: int) -> dict:
    if n < 1:
        return {"status": "NOT_EXECUTED", "sample_size": 0}
    from rxni.integration.retriever import V66Retriever
    pf = pq.ParquetFile(final / "reaction_representations.parquet")
    groups = sorted(set(round(i*(pf.num_row_groups-1)/max(1, n-1)) for i in range(n)))
    a = V66Retriever(stage)
    b = V66Retriever(final)
    passed = 0
    try:
        for rg in groups:
            reaction = pf.read_row_group(rg, columns=["reaction"])["reaction"][0].as_py()
            def run(retriever):
                return retriever.retrieve(reaction, mode="sem_multi", top_k=5,
                                          max_members=500, member_limit=10)
            left, right = run(a), run(b)
            def signature(result):
                return [(x.get("reaction_key"), x.get("template_hash"), x.get("scores"))
                        for x in result.get("members", [])]
            if not signature(left) or signature(left) != signature(right):
                raise ValueError(f"Smoke ranking mismatch or no candidates: row group {rg}")
            passed += 1
    finally:
        a.close()
        b.close()
    return {"status": "PASS", "sample_size": len(groups), "passed": passed,
            "scope": "Source staged-versus-distribution ranking/score parity; not scientific accuracy"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--staged", type=Path, required=True, help="Previously tested v6.7/Ring031 five-file directory")
    ap.add_argument("--historical", type=Path, required=True, help="V5 native template library source")
    ap.add_argument("--ring", type=Path, required=True, help="Original P6/Ring031 reprofile source and truthful provenance")
    ap.add_argument("--out", type=Path, required=True, help="New empty output directory; never an input directory")
    ap.add_argument("--smoke", type=int, default=3, help="Staged/final retrieval pairs to test; 0 skips")
    args = ap.parse_args()
    staged, historical, ring, out = (x.expanduser().resolve() for x in
                                     (args.staged, args.historical, args.ring, args.out))
    for src in (staged, historical, ring):
        if not src.is_dir():
            ap.error(f"Missing directory: {src}")
    if out.exists() or out in (staged, historical, ring):
        ap.error(f"Output must not exist or equal an input: {out}")
    if args.smoke < 0:
        ap.error("--smoke must be >=0")
    st = resolve_library(staged)
    ring_manifest_path = ring / "retrieval_v2_manifest.json"
    ring_manifest = json.loads(ring_manifest_path.read_text())
    prov_file = ring / "V67_TRUTHFUL_SEMANTIC_PROVENANCE.json"
    provenance = json.loads(prov_file.read_text())
    if "runtime_v6_7.semantic" not in str(ring_manifest.get("scientific_contract", {}).get("semantic_authority", "")):
        raise ValueError("Ring source manifest does not attest v6.7 semantics")
    if "Ring" not in str(ring_manifest["scientific_contract"]["semantic_authority"]):
        raise ValueError("Ring031 provenance absent")
    if provenance.get("build_id") != ring_manifest.get("build_id"):
        raise ValueError("Ring provenance/manifest build_id mismatch")
    receipt_hash = provenance.get("corrected_retrieval_manifest_sha256")
    if receipt_hash and sha(ring_manifest_path) != receipt_hash:
        raise ValueError("Ring manifest bytes no longer match truthful scientific provenance receipt")
    if ring_manifest.get("rows") != 812916:
        raise ValueError("Selected historical Ring031 source should contain exactly 812,916 records")
    native_core = historical / "reaction_core_v5.parquet"
    expected_core = ring_manifest.get("source_reaction_core_sha256")
    if native_core.is_file() and expected_core and sha(native_core) != expected_core:
        raise ValueError("Historical native Core bytes differ from Ring031 source Core SHA256")
    if st.layout not in ("descriptor-v1", "legacy-v2-directory"):
        raise ValueError(f"Unsupported staged layout: {st.layout}")
    staged_manifest = json.loads(st.retrieval_manifest.read_text())
    if staged_manifest.get("build_id") != ring_manifest["build_id"]:
        raise ValueError("Staged data is not the selected Ring031 build")
    if sha(st.representations) != sha(ring / "retrieval_representations_v2.parquet"):
        raise ValueError("Staged representation bytes differ from Ring031 scientific source")
    if sha(st.representations) != ring_manifest["representations_sha256"]:
        raise ValueError("Ring031 Parquet differs from original manifest digest")
    print("PASS: original Ring031 scientific payload and provenance checked", flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=out.name + ".BUILDING.", dir=out.parent))
    try:
        origins = {
            "structural_index.sqlite": st.structural_index,
            "reaction_lookup.sqlite": st.metadata_lookup,
            "representation_locator.sqlite": st.locator,
        }
        for name, source in origins.items():
            snapshot_db(source, work/name)
            print("PASS: SQLite snapshot", name, flush=True)
        shutil.copy2(st.representations, work/"reaction_representations.parquet")
        catalog = historical/"template_library_v5.parquet"
        if not catalog.exists():
            raise FileNotFoundError(catalog)
        shutil.copy2(catalog, work/"template_catalog.parquet")
        manifest = dict(ring_manifest)
        manifest["locator_sqlite_sha256"] = sha(work/"representation_locator.sqlite")
        manifest.pop("source_library", None)
        manifest.pop("source_reaction_core", None)
        manifest.pop("locator_parquet_sha256", None)
        manifest["distribution_metadata"] = {
            "dataset_title": "RXNI Historical USPTO Reaction Library — Release 1",
            "packaging_release": "r1",
            "source_profile": "v6.7 P6 / Ring v0.3.1",
            "original_manifest_sha256": sha(ring_manifest_path),
            "source_core_sha256": ring_manifest.get("source_reaction_core_sha256"),
            "catalogue_sha256": sha(work/"template_catalog.parquet"),
            "locator_resnapshotted": True,
            "scientific_representation_unchanged": True,
            "stored_representation_contract": FORMAT,
        }
        write_json(work/"retrieval_manifest.json", manifest)
        descriptor = {
            "schema_version": 1,
            "library_id": "uspto-historical-v67-ring031-r1",
            "representation_format": FORMAT,
            "semantic_engine_version": "runtime_v6_7.semantic post-P6 / Ring v0.3.1",
            "artifacts": {
                "structural_index": "structural_index.sqlite",
                "metadata_lookup": "reaction_lookup.sqlite",
                "representations": "reaction_representations.parquet",
                "locator": "representation_locator.sqlite",
                "retrieval_manifest": "retrieval_manifest.json",
            },
        }
        write_json(work/"rxni-library.json", descriptor)
        result = inspect_library(work)
        if result.get("status") != "V2_READABLE" or result.get("layout") != "descriptor-v1":
            raise ValueError(f"RXNI descriptor not compatible: {result}")
        print("PASS: RXNI descriptor format", flush=True)
        template_report, hashes = audit_catalog(work/"template_catalog.parquet")
        if len(hashes) != 42128:
            raise ValueError(f"Historical template population expected 42,128; observed {len(hashes)}")
        print("PASS: template catalogue SMARTS validation", template_report["counts"], flush=True)
        index_report = audit_full_index(work, ring_manifest["build_id"], hashes)
        print("PASS: full Parquet↔locator coordinates and index key populations", flush=True)
        retrieval = smoke_retrieval(staged, work, args.smoke)
        print("Retrieval smoke:", retrieval, flush=True)
        data_prov = {
            "scientific_source": "v67_retrieval_v2_full_candidate_P6_Ring031_01",
            "scientific_source_status": provenance.get("status"),
            "accepted_semantic_authority": ring_manifest["scientific_contract"]["semantic_authority"],
            "semantic_source_pins_sha256": provenance.get("source_pin", {}),
            "source_reaction_core_sha256": ring_manifest.get("source_reaction_core_sha256"),
            "source_retrieval_manifest_sha256": sha(ring_manifest_path),
            "source_truthful_provenance_sha256": sha(prov_file),
            "source_template_catalogue_sha256": sha(catalog),
            "known_scientific_review_queue": provenance.get("known_scientific_review_queue", []),
            "published_scientific_validation": "NOT_CLAIMED",
            "rights_clearance": "NOT_ATTESTED",
            "source_paths_excluded": True,
            "note": "This file records derivative provenance; original private build receipts are retained separately.",
        }
        write_json(work/"SCIENTIFIC_PROVENANCE.json", data_prov)
        report = {
            "status": "TECHNICAL_CANDIDATE_PASS_PENDING_SCIENTIFIC_AND_RIGHTS_REVIEW",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "rxni_library_compatibility": result["status"],
            "build_id": ring_manifest["build_id"],
            "semantic_authority": data_prov["accepted_semantic_authority"],
            "retrieval_smoke": retrieval,
            "full_index_audit": index_report,
            "smarts_catalogue_audit": template_report,
            "historical_template_hash_replay": "42,128/42,128 reported from rxnutils in aizynthtrain environment; independent signed artifact not included",
            "publication_gates": {"corpus_redistribution_rights": "PENDING",
                                  "known_chemical_review_queue": data_prov["known_scientific_review_queue"],
                                  "clean_environment_acceptance": "PENDING"},
        }
        write_json(work/"RELEASE_AUDIT.json", report)
        (work/"README.md").write_text(
            "# RXNI Historical USPTO Reaction Library — Release 1\n\n"
            "This directory is a **technical distribution candidate**, not an authorized public release.\n"
            "Complete the distribution-rights and scientific acceptance reviews before publishing.\n\n"
            "Use RXNI v0.1.0+ with `rxni-library.json`; internal representation storage generation remains V2.\n"
            "`template_catalog.parquet` is an additional, independently readable SMARTS catalogue; "
            "it is not required by the five-file retrieval reader.\n\n"
            "Check `SHA256SUMS`, `RELEASE_AUDIT.json`, and `SCIENTIFIC_PROVENANCE.json`.\n",
            encoding="utf8",
        )
        files_to_hash = list(DATA_NAMES) + ["rxni-library.json", "SCIENTIFIC_PROVENANCE.json", "RELEASE_AUDIT.json", "README.md"]
        (work/"SHA256SUMS").write_text(
            "".join(f"{sha(work/name)}  {name}\n" for name in files_to_hash), encoding="utf8"
        )
        work.rename(out)
        print("\nPASS: technical distribution candidate finalized:", out, flush=True)
        print("NOT APPROVED FOR PUBLICATION: data rights + scientific promotion gates remain", flush=True)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
