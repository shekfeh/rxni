from __future__ import annotations
import json
from pathlib import Path
from typing import Any, Iterable
import pyarrow.parquet as pq
from .schema import ALIASES, ConditionRecord

def _choose(names:set[str], aliases):
    return next((x for x in aliases if x in names),None)

def _float(v):
    if v is None or v=="": return None
    try: return float(v)
    except Exception: return None

def _record(row:dict[str,Any], mapping:dict[str,str|None]):
    def g(k):
        c=mapping.get(k); return row.get(c) if c else None
    return ConditionRecord(reaction_key=None if g("reaction_key") is None else str(g("reaction_key")),template_hash=None if g("template_hash") is None else str(g("template_hash")),reaction_id=None if g("reaction_id") is None else str(g("reaction_id")),catalyst=g("catalyst"),base=g("base"),solvent=g("solvent"),reagent=g("reagent"),temperature=_float(g("temperature")),time=_float(g("time")),condition_cluster=None if g("condition_cluster") is None else str(g("condition_cluster")),raw=row)

class ConditionEvidenceStore:
    def __init__(self, records:list[ConditionRecord]):
        self.records=records; self.by_reaction={}; self.by_template={}
        for r in records:
            if r.reaction_key: self.by_reaction.setdefault(r.reaction_key,[]).append(r)
            if r.template_hash: self.by_template.setdefault(r.template_hash,[]).append(r)
    @classmethod
    def load(cls,path:str|Path):
        p=Path(path); records=[]
        if p.suffix.lower()=='.parquet':
            pf=pq.ParquetFile(p); names=set(pf.schema_arrow.names); mapping={k:_choose(names,v) for k,v in ALIASES.items()}; cols=sorted({x for x in mapping.values() if x})
            for b in pf.iter_batches(columns=cols,batch_size=65536):
                for row in b.to_pylist(): records.append(_record(row,mapping))
        else:
            rows=[]
            with p.open(encoding='utf-8') as fh:
                for line in fh:
                    if line.strip(): rows.append(json.loads(line))
            names=set().union(*(r.keys() for r in rows)) if rows else set(); mapping={k:_choose(names,v) for k,v in ALIASES.items()}
            records=[_record(r,mapping) for r in rows]
        return cls(records)
    def for_reactions(self,keys:Iterable[str]):
        out=[]
        for k in keys: out.extend(self.by_reaction.get(str(k),[]))
        return out
    def for_template(self,template_hash:str|None):
        return [] if not template_hash else list(self.by_template.get(str(template_hash),[]))
