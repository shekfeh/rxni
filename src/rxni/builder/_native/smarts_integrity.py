"""Fail closed on incomplete template metadata. Can be copied into builder_v5."""
from pathlib import Path

HASHES=('retro_template_hash','template_hash','reaction_template_hash')
SMARTS=('retro_template_smarts','retro_template','template_smarts','reaction_smarts','smarts')


def require_source_template_smarts(path):
    import pyarrow.parquet as pq
    p=Path(path)
    fields=pq.read_schema(p).names
    hash_col=next((x for x in HASHES if x in fields),None)
    smarts_col=next((x for x in SMARTS if x in fields),None)
    if not hash_col or not smarts_col:
        raise ValueError(f'Template-complete build requires template hash and SMARTS; file={p}; available={fields}')
    tab=pq.read_table(p,columns=[hash_col,smarts_col])
    seen={}; bad=[]; conflicts=[]
    for i,(hash_,smarts) in enumerate(zip(tab.column(hash_col).to_pylist(),tab.column(smarts_col).to_pylist())):
        if not isinstance(smarts,str) or not smarts.strip():bad.append(i)
        if not hash_:raise ValueError(f'Empty template hash at source row {i}')
        if hash_ in seen and seen[hash_]!=smarts:conflicts.append(hash_)
        else:seen[hash_]=smarts
    if bad or conflicts:
        raise ValueError(f'Incomplete/contradictory source template SMARTS: empty={len(bad)} conflicting_hashes={len(set(conflicts))}; first_empty_rows={bad[:8]}')
    return {'rows':tab.num_rows,'distinct_hashes':len(seen),'nonempty_smarts':tab.num_rows}


def require_joined_template_smarts(metadata, hash_):
    meta=metadata.get(hash_)
    if meta is None:raise ValueError(f'MISSING_TEMPLATE_METADATA_JOIN: template hash {hash_!r}')
    smarts=meta.get('retro_template_smarts')
    if not isinstance(smarts,str) or not smarts.strip():
        raise ValueError(f'MISSING_TEMPLATE_SMARTS: template hash {hash_!r}')
    return smarts
