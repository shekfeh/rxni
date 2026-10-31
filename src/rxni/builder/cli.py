"""Stable CLI: version-independent Core and explicitly profiled Auxiliary."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from .knowledge import make_core,audit_core
from .auxiliary import runtime_inspection,build_auxiliary,assemble_candidate
from .bridge import get_schema
from .profiles import V67
from dataclasses import asdict

def main():
    p=argparse.ArgumentParser(description='RXNI Library Builder: Core + runtime-aware Retrieval Auxiliary')
    sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('profiles',help='Available runtime profiles (v6.7 only)')
    c=sub.add_parser('inspect',help='Read source Parquet schemas and field mappings')
    c.add_argument('--export-dir',required=True);c.add_argument('--reaction-column')
    c=sub.add_parser('build-core',help='Immutable Core Library only; NO retrieval generation')
    src=c.add_mutually_exclusive_group(required=True)
    src.add_argument('--export-dir');src.add_argument('--prepared')
    c.add_argument('--out-dir',required=True);c.add_argument('--dataset',required=True)
    c.add_argument('--reaction-column',default='reaction')
    c=sub.add_parser('validate-core',help='Full integrity check of all Core members and SMARTS')
    c.add_argument('--core',required=True);c.add_argument('--rdkit',action='store_true');c.add_argument('--report')
    c=sub.add_parser('inspect-runtime',help='Inspect installed native/v6.7 code producers')
    c.add_argument('--root',required=True);c.add_argument('--runtime',default='v6.7')
    c=sub.add_parser('build-auxiliary',help='Generate isolated prerequisite + v6.7 retrieval-v2 artifacts')
    c.add_argument('--core',required=True);c.add_argument('--work-dir',required=True)
    c.add_argument('--root',required=True);c.add_argument('--runtime',default='v6.7')
    c.add_argument('--workers',type=int,default=4)
    c.add_argument('--validation-sample',type=int,default=500)
    c.add_argument('--chunk-size',type=int,default=2048)
    c.add_argument('--through',choices=['native','v2'],default='v2')
    c.add_argument('--resume',action='store_true')
    c=sub.add_parser('assemble-candidate',help='Assemble five runtime files; NOT scientific release')
    c.add_argument('--core',required=True);c.add_argument('--work-dir',required=True)
    c.add_argument('--root',required=True);c.add_argument('--runtime',default='v6.7')
    c.add_argument('--lookup-db',required=True,help='Lookup DB built from SAME Core with installed producer')
    c.add_argument('--out-dir',required=True)
    a=p.parse_args()
    if a.command=='profiles':result={'v6.7':asdict(V67)}
    elif a.command=='inspect':result=get_schema(a.export_dir,reaction_column=a.reaction_column)
    elif a.command=='build-core':result=make_core(out_dir=a.out_dir,export_dir=a.export_dir,
                    prepared=a.prepared,reaction_column=a.reaction_column,dataset=a.dataset)
    elif a.command=='validate-core':result=audit_core(a.core,rdkit=a.rdkit,report=a.report)
    elif a.command=='inspect-runtime':result=runtime_inspection(a.root,a.runtime)
    elif a.command=='build-auxiliary':result=build_auxiliary(core=a.core,work_dir=a.work_dir,root=a.root,
                    runtime=a.runtime,workers=a.workers,validation_sample=a.validation_sample,
                    chunk_size=a.chunk_size,resume=a.resume,through=a.through)
    elif a.command=='assemble-candidate':result=assemble_candidate(core=a.core,work_dir=a.work_dir,
                    lookup_db=a.lookup_db,out_dir=a.out_dir,root=a.root,runtime=a.runtime)
    else:raise AssertionError(a.command)
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
