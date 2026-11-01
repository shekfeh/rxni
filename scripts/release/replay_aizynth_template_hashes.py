#!/usr/bin/env python3
"""Replay AiZynthTrain fingerprint-based template identities from a JSONL export.

Run with the pinned historical AiZynthTrain environment (rxnutils required):
 conda run --no-capture-output -n aizynthtrain python \
   scripts/release/replay_aizynth_template_hashes.py \
   --jsonl /path/to/template_hash_full_42128.jsonl \
   --report /path/to/template_hash_replay_report.json
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import rdkit
from rxnutils.chem.template import ReactionTemplate


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--jsonl',required=True,type=Path)
    ap.add_argument('--report',required=True,type=Path)
    args=ap.parse_args()
    args.report.parent.mkdir(parents=True,exist_ok=True)
    matched=failed=0
    errors=[]
    with args.jsonl.open(encoding='utf-8') as f:
        for i,line in enumerate(f,1):
            if not line.strip(): continue
            row=json.loads(line)
            try:
                actual=ReactionTemplate(row['retro_template_smarts'],direction='retro').hash_from_bits()
                if actual==row['retro_template_hash']:
                    matched+=1
                else:
                    failed+=1
                    if len(errors)<20: errors.append({'line':i,'hash':row['retro_template_hash'],'actual':actual})
            except Exception as e:
                failed+=1
                if len(errors)<20: errors.append({'line':i,'hash':row.get('retro_template_hash'),'error':str(e)})
    report={'hash_policy':'AiZynthTrain ReactionTemplate.hash_from_bits default radius=2 chirality=True',
            'rdkit_version':rdkit.__version__,'reaction_utils_version':importlib.metadata.version('reaction-utils'),
            'python_version':platform.python_version(),
            'input_jsonl_sha256':digest(args.jsonl),'matched':matched,'failed':failed,'examples':errors,
            'status':'PASS' if matched and failed==0 else 'FAIL',
            'scope':'template identity reproducibility, not chemical applicability'}
    args.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))
    if failed: raise SystemExit(1)

if __name__=='__main__':main()
