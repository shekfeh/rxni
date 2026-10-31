from __future__ import annotations
import sqlite3
from pathlib import Path
from typing import Iterable, Any
import pyarrow.parquet as pq

CACHE_NAME = "runtime_lookup_v1.sqlite"

def build_runtime_lookup(library: str | Path, *, overwrite: bool=False, batch_size: int=65536) -> dict[str, Any]:
    lib=Path(library); src=lib/"retrieval_members_v5.parquet"; db=lib/CACHE_NAME
    if not src.exists(): raise FileNotFoundError(src)
    if db.exists():
        if not overwrite: return {"status":"EXISTS","path":str(db)}
        db.unlink()
    conn=sqlite3.connect(db)
    conn.executescript("""
    PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL;
    CREATE TABLE member_lookup (
      reaction_key TEXT PRIMARY KEY,
      forward_reaction_smiles TEXT NOT NULL,
      retro_template_hash TEXT,
      retro_template_number TEXT,
      reaction_id TEXT,
      membership_occurrence_count INTEGER NOT NULL,
      forward_edit_signature_hash TEXT,
      forward_raw_center_hash_r0 TEXT,
      forward_raw_center_hash_r1 TEXT,
      forward_raw_center_hash_r2 TEXT
    );
    CREATE INDEX idx_member_template ON member_lookup(retro_template_hash);
    CREATE INDEX idx_member_edit ON member_lookup(forward_edit_signature_hash);
    """)
    cols=["reaction_key","forward_reaction_smiles","retro_template_hash","retro_template_number","reaction_id","membership_occurrence_count","forward_edit_signature_hash","forward_raw_center_hash_r0","forward_raw_center_hash_r1","forward_raw_center_hash_r2"]
    n=0
    for batch in pq.ParquetFile(src).iter_batches(columns=cols,batch_size=batch_size):
        rows=[]
        for r in batch.to_pylist():
            rows.append((str(r["reaction_key"]),str(r["forward_reaction_smiles"]),r.get("retro_template_hash"),r.get("retro_template_number"),r.get("reaction_id"),int(r.get("membership_occurrence_count") or 1),None if r.get("forward_edit_signature_hash") is None else str(int(r["forward_edit_signature_hash"])),None if r.get("forward_raw_center_hash_r0") is None else str(int(r["forward_raw_center_hash_r0"])),None if r.get("forward_raw_center_hash_r1") is None else str(int(r["forward_raw_center_hash_r1"])),None if r.get("forward_raw_center_hash_r2") is None else str(int(r["forward_raw_center_hash_r2"]))))
        conn.executemany("INSERT OR REPLACE INTO member_lookup VALUES (?,?,?,?,?,?,?,?,?,?)",rows); conn.commit(); n += len(rows)
    conn.close()
    return {"status":"BUILT","path":str(db),"rows":n}

def connect_lookup(library: str | Path, *, auto_build: bool=True) -> sqlite3.Connection:
    lib=Path(library); db=lib/CACHE_NAME
    if not db.exists():
        if not auto_build: raise FileNotFoundError(db)
        build_runtime_lookup(lib)
    return sqlite3.connect(db)

def fetch_members(conn: sqlite3.Connection, keys: Iterable[str]) -> dict[str, dict[str, Any]]:
    ks=list(dict.fromkeys(str(x) for x in keys))
    out={}
    # SQLite parameter limit safe batching.
    for i in range(0,len(ks),800):
        chunk=ks[i:i+800]; qs=','.join('?' for _ in chunk)
        cur=conn.execute(f"SELECT reaction_key,forward_reaction_smiles,retro_template_hash,retro_template_number,reaction_id,membership_occurrence_count,forward_edit_signature_hash,forward_raw_center_hash_r0,forward_raw_center_hash_r1,forward_raw_center_hash_r2 FROM member_lookup WHERE reaction_key IN ({qs})",chunk)
        for r in cur:
            out[r[0]]={"reaction_key":r[0],"forward_reaction_smiles":r[1],"retro_template_hash":r[2],"retro_template_number":r[3],"reaction_id":r[4],"membership_occurrence_count":int(r[5]),"forward_edit_signature_hash":None if r[6] is None else int(r[6]),"forward_raw_center_hash_r0":None if r[7] is None else int(r[7]),"forward_raw_center_hash_r1":None if r[8] is None else int(r[8]),"forward_raw_center_hash_r2":None if r[9] is None else int(r[9])}
    return out


def main():
    import argparse, json
    ap=argparse.ArgumentParser(description="Build the original RXNI runtime lookup SQLite file")
    ap.add_argument("--library", required=True)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--batch-size", type=int, default=65536)
    args=ap.parse_args()
    print(json.dumps(build_runtime_lookup(args.library, overwrite=args.overwrite, batch_size=args.batch_size),indent=2))

if __name__=="__main__":
    main()
