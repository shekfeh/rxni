from __future__ import annotations
from collections import Counter
from statistics import median
from typing import Any

def _tokens(v:Any):
    if v is None: return []
    if isinstance(v,(list,tuple,set)): return [str(x).strip() for x in v if str(x).strip()]
    s=str(v).strip()
    if not s: return []
    for sep in (';','|'):
        if sep in s: return [x.strip() for x in s.split(sep) if x.strip()]
    return [s]

def summarize_conditions(records, *, weights:dict[str,float]|None=None):
    weights=weights or {"catalyst":3.0,"base":2.0,"reagent":2.0,"solvent":1.0,"condition_cluster":2.0}
    counts={k:Counter() for k in weights}; temps=[]; times=[]
    for r in records:
        for field in weights:
            for x in _tokens(getattr(r,field,None)): counts[field][x]+=1
        if r.temperature is not None: temps.append(float(r.temperature))
        if r.time is not None: times.append(float(r.time))
    ranked={}
    n=max(1,len(records))
    for field,c in counts.items():
        ranked[field]=[{"value":v,"count":cnt,"fraction":cnt/n,"weighted_support":weights[field]*cnt} for v,cnt in c.most_common(10)]
    return {"n_records":len(records),"components":ranked,"temperature_median":median(temps) if temps else None,"temperature_n":len(temps),"time_median":median(times) if times else None,"time_n":len(times)}
