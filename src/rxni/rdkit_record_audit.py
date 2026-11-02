"""Streaming record-level RDKit parse/sanitization diagnostic tool.

Produces a private report keyed by stored reaction_key; not an auto-cleaner.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rxni.retrieval.libraries import resolve_library
import pyarrow.parquet as pq


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--library',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--max-rows',type=int,default=5000)
    args=ap.parse_args(argv)
    if args.max_rows < 1: ap.error('--max-rows must be >0')
    library=resolve_library(args.library)
    pf=pq.ParquetFile(library.representations)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    count=0; failed=0
    with args.output.open('w',encoding='utf8') as out:
        for rg in range(pf.num_row_groups):
            if count>=args.max_rows: break
            table=pf.read_row_group(rg,columns=['reaction_key','reaction'])
            for offset,row in enumerate(table.to_pylist()):
                if count>=args.max_rows: break
                count+=1
                errs=[]
                reaction=row['reaction'] or ''
                parts=reaction.split('>')
                if len(parts)!=3: errs.append('reaction_syntax')
                else:
                    for side in (0,2):
                        for fragment in parts[side].split('.'):
                            with rdBase.BlockLogs():
                                try:
                                    mol=Chem.MolFromSmiles(fragment,sanitize=False)
                                    if mol is None: errs.append(f'side_{side}_parse_failure');continue
                                    status=Chem.SanitizeMol(mol,catchErrors=True)
                                    if status != Chem.SanitizeFlags.SANITIZE_NONE:
                                        errs.append(f'side_{side}_{status}')
                                    if any(a.GetAtomicNum()==15 for a in mol.GetAtoms()):
                                        errs.append(f'side_{side}_contains_phosphorus')
                                    if any(a.GetAtomicNum()==1 and a.GetDegree()==0 for a in mol.GetAtoms()):
                                        errs.append(f'side_{side}_isolated_hydrogen')
                                except Exception as exc:
                                    errs.append(f'side_{side}_{type(exc).__name__}')
                if errs:
                    failed+=1
                    out.write(json.dumps({'reaction_key':row['reaction_key'],
                          'row_group':rg,'row_in_group':offset,'observations':sorted(set(errs))})+'\n')
    print(json.dumps({'inspected':count,'flagged':failed,'report':str(args.output),
                      'scope':'parse/sanitization findings and structure flags, not an exclusion recommendation'}))

if __name__=='__main__':main()
