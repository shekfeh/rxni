#!/usr/bin/env python3
"""Read-only verification of version-free RXNI library archive after download.

  python scripts/release/verify_published_library.py --library /data/RXNI-USPTO --smoke

The required five-file reader is supported by installed RXNI 0.1.0.
Template catalogue is optional for the reader, but required by this release.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sqlite3

import pyarrow.parquet as pq
from rxni.retrieval.libraries import inspect_library, resolve_library


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--library',required=True,type=Path)
    ap.add_argument('--smoke',action='store_true',help='Execute one example reaction retrieval')
    args=ap.parse_args()
    root=args.library.resolve()
    lines=(root/'SHA256SUMS').read_text().splitlines()
    print('=== Checksums ===')
    for line in lines:
        if not line.strip(): continue
        expected,name=line.split('  ',1)
        name=name.strip()
        if Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError('Unsafe checksum filename')
        got=sha256(root/name)
        if got!=expected:
            raise ValueError(f'SHA256 mismatch: {name}')
        print('PASS',name)
    print('=== RXNI format ===')
    info=inspect_library(root)
    if info.get('status')!='V2_READABLE':
        raise ValueError(info)
    resolved=resolve_library(root)
    manifest=json.loads(resolved.retrieval_manifest.read_text())
    if manifest['representation_generation']!='retrieval_representations_v2':
        raise ValueError('Unexpected representation ABI')
    if sha256(resolved.representations)!=manifest['representations_sha256']:
        raise ValueError('Representation digest mismatch against manifest')
    if sha256(resolved.locator)!=manifest['locator_sqlite_sha256']:
        raise ValueError('Locator digest mismatch against manifest')
    parquet=pq.ParquetFile(resolved.representations)
    if parquet.metadata.num_rows!=manifest['rows']:
        raise ValueError('Manifest/Parquet rows mismatch')
    for db in (resolved.structural_index,resolved.metadata_lookup,resolved.locator):
        cx=sqlite3.connect(db.resolve().as_uri()+'?mode=ro',uri=True)
        try:
            if cx.execute('PRAGMA integrity_check').fetchone()[0]!='ok':
                raise ValueError(f'SQLite database integrity failed: {db}')
        finally:cx.close()
    catalog=root/'template_catalog.parquet'
    if not catalog.exists():
        raise ValueError('Missing historical SMARTS catalogue')
    cat=pq.ParquetFile(catalog)
    print('PASS',parquet.metadata.num_rows,'representations')
    print('PASS',cat.metadata.num_rows,'templates')
    print('PASS storage format and SQLite integrity')
    if args.smoke:
        from rxni.integration.retriever import V66Retriever
        reaction=parquet.read_row_group(0,columns=['reaction'])['reaction'][0].as_py()
        rt=V66Retriever(root)
        try:
            result=rt.retrieve(reaction,mode='sem_multi',top_k=5,max_members=500,member_limit=10)
        finally:rt.close()
        if not result.get('members'):
            raise ValueError('Retrieval returned no candidates')
        print('PASS real retrieval: candidates=',len(result['members']))
    print('PASS: downloaded directory is usable by RXNI')
    print('This does not certify chemical accuracy or data redistribution rights.')

if __name__=='__main__':
    main()
