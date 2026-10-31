"""Strict preparation of rxni_v5_style_1 exports for the historical native builder.

Only reaction/template membership, source identity, and authentic template SMARTS
are authoritative. Structural hashes/semantics are *never* imported.
"""
from __future__ import annotations
import hashlib
import json
import os
import re
import shutil
import sqlite3
import tempfile
from pathlib import Path

TEMPLATE_NAME = "template_library_v5.parquet"
MEMBERSHIP_NAME = "reaction_template_membership_v5.parquet"
CORE_NAME = "reaction_core_v5.parquet"
# Priority is explicit; never use original, unmapped reaction for structural derivation.
REACTION_CANDIDATES = (
    "native_cleaned_reaction", "reaction", "mapped_rxn", "mapped_reaction_smiles", "mapped_reaction",
    "forward_reaction_smiles", "mapped_reaction_string", "reaction", "mapping_input",
)
CORE_KEY_CANDIDATES = ("reaction_key", "source_reaction_key", "source_observation_key", "observation_key")
MEMBERSHIP_KEY_CANDIDATES = CORE_KEY_CANDIDATES
TEMPLATE_SMARTS_CANDIDATES = ("retro_smarts", "retro_template_smarts", "retro_template")
SOURCE_ID_CANDIDATES = ("source_sample_reaction_id", "source_reaction_id", "original_source_id", "reaction_id", "source_id", "observation_id")


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def pick(names, candidates, title, *, required=True):
    found = [a for a in candidates if a in names]
    if not found and required:
        raise ValueError(f"Missing {title}. Expected one of {candidates}; got {sorted(names)}")
    return found[0] if found else None


def clean(value):
    if value is None:
        return None
    s = str(value)
    return s if s.strip() else None


def get_schema(export_dir, *, reaction_column=None):
    import pyarrow.parquet as pq
    d = Path(export_dir)
    out = {}
    for key, filename in (("templates", TEMPLATE_NAME), ("memberships", MEMBERSHIP_NAME), ("core", CORE_NAME)):
        f = d / filename
        if not f.is_file():
            raise FileNotFoundError(f)
        pf = pq.ParquetFile(f)
        out[key] = {"path": str(f.resolve()), "columns": pf.schema_arrow.names, "rows": pf.metadata.num_rows}
    t, m, c = (out[x]["columns"] for x in ("templates", "memberships", "core"))
    out["resolved"] = {
        "template_hash": pick(t, ("retro_template_hash",), "template hash"),
        "template_smarts": pick(t, TEMPLATE_SMARTS_CANDIDATES, "template SMARTS"),
        "template_number": pick(t, ("retro_template_number",), "template number", required=False),
        "membership_hash": pick(m, ("retro_template_hash",), "membership template hash"),
        "membership_key": pick(m, MEMBERSHIP_KEY_CANDIDATES, "membership reaction key"),
        "membership_source_id": pick(m, SOURCE_ID_CANDIDATES, "membership source ID", required=False),
        "membership_number": pick(m, ("retro_template_number",), "membership template number", required=False),
        "core_key": pick(c, CORE_KEY_CANDIDATES, "core reaction key"),
        "core_mapped_rxn": (reaction_column if reaction_column is not None else pick(c, REACTION_CANDIDATES, "mapped forward reaction SMILES")),
    }
    if reaction_column is not None and reaction_column not in c:
        raise ValueError(f"Chosen mapped forward reaction column missing: {reaction_column}")
    return out


def _batch_rows(path, columns=None, batch_size=8192):
    import pyarrow.parquet as pq
    yield from pq.ParquetFile(path).iter_batches(batch_size=batch_size, columns=columns)


def _serialize(data):
    def convert(v):
        if v is None or isinstance(v, (bool, str, int, float)):
            return v
        if isinstance(v, dict):
            return {str(k): convert(w) for k,w in v.items()}
        if isinstance(v, (list, tuple)):
            return [convert(w) for w in v]
        return str(v)
    return json.dumps(convert(data), ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _split_forward_reaction_sections(rxn: str):
    """Split reaction separators, excluding RDKit dative-bond `->` tokens."""
    parts=[];start=0
    for i,char in enumerate(rxn):
        if char=='>' and (i==0 or rxn[i-1]!='-'):
            parts.append(rxn[start:i]);start=i+1
    parts.append(rxn[start:])
    if len(parts)!=3:
        raise ValueError(f"Expected 2 reaction separators outside dative bonds: {rxn[:160]}")
    return tuple(parts)


def _check_forward_smiles(value):
    rxn=clean(value)
    if rxn is None:raise ValueError('Missing mapped forward reaction')
    left,agents,product=_split_forward_reaction_sections(rxn)
    if not left or not product:raise ValueError(f'Empty reactant/product: {rxn[:120]}')
    if not re.search(r"\[[^\]]+:\d+\]",left) or not re.search(r"\[[^\]]+:\d+\]",product):
        raise ValueError(f'Mapped atoms required in reactants and products: {rxn[:120]}')
    return rxn


def _create_db(path):
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=DELETE")
    conn.execute("PRAGMA synchronous=FULL")
    conn.execute("CREATE TABLE core (key TEXT PRIMARY KEY, rxn TEXT NOT NULL, data TEXT NOT NULL)")
    conn.execute("CREATE TABLE templates (hash TEXT PRIMARY KEY, smarts TEXT NOT NULL, number TEXT)")
    conn.execute("CREATE TABLE observed (hash TEXT PRIMARY KEY)")
    return conn


def _write_table(writer, rows, schema):
    import pyarrow as pa
    if rows:
        writer.write_table(pa.Table.from_pylist(rows, schema=schema), row_group_size=len(rows))
        rows.clear()


def prepare(export_dir, out_dir, *, reaction_column=None, keep_core_json=True, batch_size=8192, max_bad_examples=8):
    """Fails closed if the entire source cannot be joined, or template SMARTS missing.

    Fresh exclusive output directory; on failure removes its own partial artifacts.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    src = Path(export_dir).resolve()
    out = Path(out_dir).resolve()
    if out.exists():
        raise FileExistsError(f"Refusing existing output directory: {out}")
    layout = get_schema(src, reaction_column=reaction_column)
    res = dict(layout["resolved"])
    if reaction_column is not None:
        if reaction_column not in layout["core"]["columns"]:
            raise ValueError(f"--reaction-column {reaction_column} absent from source core")
        res["core_mapped_rxn"] = reaction_column
    if reaction_column is None and res["core_mapped_rxn"] in ("reaction_smiles", "original_reaction"):
        raise ValueError("Refusing unqualified original reaction column")
    # Do not silently drop schema columns that might contain a higher-fidelity mapped rxn.
    out.mkdir(parents=True)
    dbfile = out / "bridge_index.sqlite"
    tplfile = out / "source_templates.parquet"
    memfile = out / "source_memberships.parquet"
    tplschema = pa.schema([
        ("retro_template_hash", pa.string()), ("retro_template_number", pa.string()),
        ("retro_template_smarts", pa.string()),
    ])
    memschema = pa.schema([
        ("forward_reaction_smiles", pa.string()), ("retro_template_hash", pa.string()),
        ("retro_template_number", pa.string()), ("reaction_id", pa.string()),
        ("source_observation_key", pa.string()), ("source_dataset", pa.string()),
        ("source_core_json", pa.string()), ("source_membership_json", pa.string()),
    ])
    stats = {"source_layout": layout, "resolved": res, "counts": {},
             "template_identity": "sha256 of exact retro SMARTS",
             "chemistry_policy": "all structural and System-I fields recomputed by installed native builder"}
    t_count = c_count = m_count = 0
    core_duplicates = 0
    tpl_duplicates = 0
    used_templates = set()  # Only unique hashes, <= ~260k; no chemistry fields.
    conn = None
    try:
        conn = _create_db(dbfile)
        with pq.ParquetWriter(tplfile, tplschema, compression="zstd") as writer:
            buf=[]
            for batch in _batch_rows(src/TEMPLATE_NAME,batch_size=batch_size):
                for row in batch.to_pylist():
                    h=clean(row[res["template_hash"]]); s=clean(row[res["template_smarts"]])
                    n=clean(row.get(res["template_number"])) if res["template_number"] else None
                    if not h or not s:
                        raise ValueError("Missing template hash/SMARTS in template catalogue")
                    if sha256_text(s)!=h:
                        raise ValueError(f"Template hash mismatch: {h} != SHA256(SMARTS)")
                    old=conn.execute("SELECT smarts,number FROM templates WHERE hash=?",(h,)).fetchone()
                    if old is not None:
                        if old[0]!=s:
                            raise ValueError(f"Duplicate template hash with conflicting SMARTS: {h}")
                        tpl_duplicates+=1
                        continue
                    conn.execute("INSERT INTO templates VALUES (?,?,?)",(h,s,n))
                    buf.append({"retro_template_hash":h,"retro_template_number":n,"retro_template_smarts":s})
                    t_count+=1
                    if len(buf)>=batch_size:_write_table(writer,buf,tplschema)
                conn.commit()
            _write_table(writer,buf,tplschema)
        if t_count==0:
            raise ValueError("Zero templates in catalogue")

        key_col=res["core_key"]
        rxn_col=res["core_mapped_rxn"]
        for batch in _batch_rows(src/CORE_NAME,batch_size=batch_size):
            for row in batch.to_pylist():
                key=clean(row.get(key_col))
                if not key:raise ValueError("Null core reaction key")
                rxn=_check_forward_smiles(row.get(rxn_col))
                if 'native_cleaned_reaction' in row:
                    ncr=_check_forward_smiles(row.get('native_cleaned_reaction'))
                    if rxn!=ncr:
                        raise ValueError('Core reaction is not identical to the native cleaned template-extraction reaction')
                old=conn.execute("SELECT rxn,data FROM core WHERE key=?",(key,)).fetchone()
                blob=_serialize(row) if keep_core_json else "{}"
                if old:
                    if old[0]!=rxn:raise ValueError(f"Same source-scoped reaction key maps to different reactions: {key}")
                    core_duplicates+=1
                else:
                    conn.execute("INSERT INTO core VALUES (?,?,?)",(key,rxn,blob))
                    c_count+=1
            conn.commit()
        if c_count==0:
            raise ValueError("Zero core reactions")
        with pq.ParquetWriter(memfile, memschema, compression="zstd") as writer:
            buf=[]
            for batch in _batch_rows(src/MEMBERSHIP_NAME,batch_size=batch_size):
                for row in batch.to_pylist():
                    key=clean(row.get(res["membership_key"]))
                    h=clean(row.get(res["membership_hash"]))
                    if not key or not h:raise ValueError("Membership with missing reaction key or template hash")
                    cr=conn.execute("SELECT rxn,data FROM core WHERE key=?",(key,)).fetchone()
                    if cr is None:raise ValueError(f"Orphan membership reaction: {key}")
                    tr=conn.execute("SELECT number FROM templates WHERE hash=?",(h,)).fetchone()
                    if tr is None:raise ValueError(f"Orphan membership template: {h}")
                    source_id=(clean(row.get(res["membership_source_id"])) if res["membership_source_id"] else None) or key
                    number=(clean(row.get(res["membership_number"])) if res["membership_number"] else None) or tr[0]
                    if tr[0] and number and tr[0]!=number:
                        raise ValueError(f"Inconsistent template number for hash {h}: {tr[0]} vs {number}")
                    buf.append({
                        "forward_reaction_smiles":cr[0],"retro_template_hash":h,
                        "retro_template_number":number,"reaction_id":source_id,
                        "source_observation_key":key,
                        "source_dataset":clean(row.get("source_dataset")) or clean(layout.get("source_dataset")) or src.parent.name,
                        "source_core_json":cr[1],
                        "source_membership_json":_serialize(row) if keep_core_json else "{}",
                    })
                    used_templates.add(h)
                    m_count+=1
                    if len(buf)>=batch_size:_write_table(writer,buf,memschema)
                if m_count % (batch_size*25)<batch_size:
                    print(f"[bridge] {m_count:,} memberships",flush=True)
            _write_table(writer,buf,memschema)
        conn.commit()
        if not m_count:raise ValueError("Zero memberships")
        if len(used_templates)!=t_count:
            raise ValueError(f"Template catalogue contains {t_count-len(used_templates)} templates with zero memberships; builder would drop them")
        stats["counts"]={"templates":t_count,"core_reactions":c_count,"memberships":m_count,
                         "core_duplicate_keys":core_duplicates,"duplicate_template_rows":tpl_duplicates,
                         "distinct_used_templates":len(used_templates)}
        from datetime import datetime,timezone
        stats["created_utc"]=datetime.now(timezone.utc).isoformat()
        stats["source_sha256"]={n:sha256_file(src/n) for n in (CORE_NAME,MEMBERSHIP_NAME,TEMPLATE_NAME)}
        stats["prepared_sha256"]={n:sha256_file(out/n) for n in ("source_templates.parquet","source_memberships.parquet")}
        (out/"bridge_manifest.json").write_text(json.dumps(stats,indent=2,sort_keys=True)+"\n")
        conn.close();conn=None
        # Retain the index to aid further diagnostics, but not authoritative chemistry.
        return stats
    except BaseException:
        if conn is not None:
            conn.close()
        shutil.rmtree(out,ignore_errors=True)
        raise


def sha256_file(path, chunksize=8<<20):
    h=hashlib.sha256()
    with open(path,"rb") as fh:
        while data:=fh.read(chunksize):h.update(data)
    return h.hexdigest()
