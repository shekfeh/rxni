"""Post-build independent metadata integrity verification: zero silent SMARTS losses."""
from __future__ import annotations
import json
import sqlite3
import tempfile
from pathlib import Path
from .bridge import sha256_file


def verify(prepared, built, out_report=None, *, require_native_validation=True):
    import pyarrow.parquet as pq
    prep=Path(prepared).resolve()
    build=Path(built).resolve()
    bridge=json.loads((prep/"bridge_manifest.json").read_text())
    prep_templates=prep/"source_templates.parquet"
    native=build/"template_library_v5.parquet"
    if not native.exists():raise FileNotFoundError(native)
    if sha256_file(prep_templates)!=bridge["prepared_sha256"]["source_templates.parquet"]:
        raise RuntimeError("Source template file digest changed")
    vfile=build/"validation_v5.json"
    if require_native_validation:
        if not vfile.exists():raise FileNotFoundError(vfile)
        v=json.loads(vfile.read_text())
        if v.get("status")!="PASS":raise RuntimeError(f"Native validation not PASS: {v.get('status')}")
    with tempfile.TemporaryDirectory(prefix="rxni_lb_verify_") as td:
        conn=sqlite3.connect(str(Path(td)/"check.sqlite"))
        conn.execute("CREATE TABLE expected(hash TEXT PRIMARY KEY, smarts TEXT NOT NULL, checked INT NOT NULL DEFAULT 0)")
        count=0
        for batch in pq.ParquetFile(prep_templates).iter_batches(batch_size=8192,columns=["retro_template_hash","retro_template_smarts"]):
            for h,s in zip(batch.column(0).to_pylist(),batch.column(1).to_pylist()):
                if not h or not s:raise RuntimeError("Prepared source has missing SMARTS")
                conn.execute("INSERT INTO expected(hash,smarts) VALUES (?,?)",(h,s))
                count+=1
            conn.commit()
        observed=0
        column_names=pq.ParquetFile(native).schema_arrow.names
        for needed in ("retro_template_hash","retro_template_smarts"):
            if needed not in column_names:raise RuntimeError(f"Built template catalogue missing {needed}")
        for batch in pq.ParquetFile(native).iter_batches(batch_size=8192,columns=["retro_template_hash","retro_template_smarts"]):
            for h,s in zip(batch.column(0).to_pylist(),batch.column(1).to_pylist()):
                if not h or not isinstance(s,str) or not s.strip():
                    raise RuntimeError("Built template catalogue has null or empty SMARTS")
                ref=conn.execute("SELECT smarts,checked FROM expected WHERE hash=?",(h,)).fetchone()
                if ref is None:raise RuntimeError(f"Unexpected built template hash: {h}")
                if ref[0]!=s:raise RuntimeError(f"SMARTS changed for {h}")
                if ref[1]:raise RuntimeError(f"Duplicate built template hash: {h}")
                conn.execute("UPDATE expected SET checked=1 WHERE hash=?",(h,))
                observed+=1
            conn.commit()
        absent=conn.execute("SELECT COUNT(*) FROM expected WHERE checked=0").fetchone()[0]
        if absent or observed!=count:
            raise RuntimeError(f"Source/build incomplete: source={count}, built={observed}, missing={absent}")
        conn.close()
    report={"status":"PASS","source_templates":count,"built_templates":observed,
            "missing_or_changed_smarts":0,"hashes_unmatched":0,"exact_source_string_equality":True,
            "native_validation_status":"PASS" if require_native_validation else "NOT_CHECKED",
            "original_source_sha256":bridge["source_sha256"],
            "built_template_catalogue_sha256":sha256_file(native),
            "scope":"SMARTS identity and installed native validation; not v6.7 retrieval parity"}
    if out_report:
        dst=Path(out_report)
        dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():raise FileExistsError(dst)
        dst.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report
