"""Installed, read-only checks for an RXNI historical-library distribution.

A V2_READABLE result, integrity checks and one self-query are technical checks;
none independently establishes chemical accuracy or source redistribution rights.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_checksums(root: Path) -> int:
    sums = root / "SHA256SUMS"
    if not sums.is_file():
        raise ValueError("Missing SHA256SUMS")
    count = 0
    seen: set[str] = set()
    for line in sums.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        m = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if m is None:
            raise ValueError(f"Malformed SHA256SUMS entry: {line[:80]!r}")
        expected, name = m.groups()
        rel = Path(name)
        if rel.is_absolute() or rel.name != name or name in seen or name == "SHA256SUMS":
            raise ValueError(f"Unsafe/duplicate SHA256SUMS entry: {name!r}")
        seen.add(name)
        target = root / rel
        if not target.is_file() or target.is_symlink() or not target.resolve().is_relative_to(root):
            raise ValueError(f"Missing or unsafe file: {name}")
        if sha256(target) != expected:
            raise ValueError(f"SHA256 mismatch: {name}")
        print("PASS checksum", name)
        count += 1
    if not count:
        raise ValueError("Empty SHA256SUMS")
    return count


def verify_library(root: Path, *, smoke: bool = False) -> dict:
    import pyarrow.parquet as pq
    from rxni.retrieval.libraries import inspect_library, resolve_library

    root = root.resolve(strict=True)
    if not root.is_dir():
        raise ValueError("Library must be a directory")
    num = verify_checksums(root)
    info = inspect_library(root)
    if info.get("status") != "V2_READABLE":
        raise ValueError(f"Not retrieval-ready: {info}")
    lib = resolve_library(root)
    manifest = json.loads(lib.retrieval_manifest.read_text(encoding="utf-8"))
    if manifest.get("representation_generation") != "retrieval_representations_v2":
        raise ValueError("Unsupported stored representation generation")
    if sha256(lib.representations) != manifest.get("representations_sha256"):
        raise ValueError("Representation SHA256 differs from manifest")
    if sha256(lib.locator) != manifest.get("locator_sqlite_sha256"):
        raise ValueError("Locator SHA256 differs from manifest")
    for db in (lib.structural_index, lib.metadata_lookup, lib.locator):
        conn = sqlite3.connect(db.resolve().as_uri() + "?mode=ro", uri=True)
        try:
            if conn.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                raise ValueError(f"SQLite integrity failed: {db.name}")
        finally:
            conn.close()
    pf = pq.ParquetFile(lib.representations)
    if pf.metadata.num_rows != manifest.get("rows") or pf.metadata.num_rows <= 0:
        raise ValueError("Representation count mismatch")
    catalog = root / "template_catalog.parquet"
    if not catalog.is_file():
        raise ValueError("Missing release template catalogue")
    ntemplates = pq.ParquetFile(catalog).metadata.num_rows
    if ntemplates <= 0:
        raise ValueError("Empty template catalogue")
    if smoke:
        from rxni.integration import ReactionRetriever
        rxn = pf.read_row_group(0, columns=["reaction"])["reaction"][0].as_py()
        rt = ReactionRetriever(root)
        try:
            result = rt.retrieve(rxn, mode="sem_multi", top_k=5, max_members=500, member_limit=10)
        finally:
            rt.close()
        if not result.get("members"):
            raise ValueError("Real retrieval produced no members")
        print("PASS retrieval members", len(result["members"]))
    result = {"status": "TECHNICAL_PASS", "checksum_files": num,
              "reactions": pf.metadata.num_rows, "templates": ntemplates,
              "real_retrieval_smoke": "PASS" if smoke else "NOT_RUN",
              "chemical_accuracy": "NOT_ASSESSED", "data_rights": "NOT_ASSESSED"}
    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify a downloaded RXNI library; no private build files required")
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--smoke", action="store_true", help="Retrieve one stored reaction")
    args = ap.parse_args()
    verify_library(args.library, smoke=args.smoke)


if __name__ == "__main__":
    main()
