"""Minimal immutable Core Parquet writer, independent of retrieval generations.

Input: validated 3-file modern rxni_v5_style_1 export, or a previously
verified v0.1 bridge. Reuses disk SQLite only during normalization.
"""
from __future__ import annotations
import json, sqlite3, tempfile, shutil
from pathlib import Path
from datetime import datetime,timezone
from .bridge import prepare, sha256_file, sha256_text, _check_forward_smiles

FORMAT='rxni_core_1'
FILE_NAMES=('reactions.parquet','templates.parquet','memberships.parquet','evidence.parquet')

def _j(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,ensure_ascii=False)+'\n')

def _rows(file, columns=None):
    import pyarrow.parquet as pq
    for batch in pq.ParquetFile(file).iter_batches(batch_size=8192,columns=columns):
        yield from batch.to_pylist()

def _out(file,schema,rows):
    import pyarrow as pa
    import pyarrow.parquet as pq
    count=0;buffer=[]
    with pq.ParquetWriter(file,schema,compression='zstd') as w:
        for r in rows:
            buffer.append(r);count+=1
            if len(buffer)==8192:
                w.write_table(pa.Table.from_pylist(buffer,schema=schema));buffer.clear()
        if buffer:w.write_table(pa.Table.from_pylist(buffer,schema=schema))
    return count

def make_core(*,out_dir, export_dir=None, prepared=None, reaction_column='reaction',dataset='UNSPECIFIED'):
    """Writes no structural features; preserves each observation and exact template SMARTS."""
    if bool(export_dir)==bool(prepared):
        raise ValueError('Supply exactly one of --export-dir/--prepared')
    if prepared:
        raise ValueError('v0.2.1 strict Core requires original native-cleaned export for lineage attestation; --prepared legacy bypass disabled')
    import pyarrow as pa
    dest=Path(out_dir).resolve()
    if dest.exists():raise FileExistsError(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix=dest.name+'.BUILDING.',dir=dest.parent))
    try:
        bridge=Path(prepared).resolve() if prepared else tmp/'bridge'
        if export_dir:
            export_meta=Path(export_dir)/'manifest.json'
            if not export_meta.exists():raise ValueError('New Core requires a lineage-attested dataset export manifest')
            m=json.loads(export_meta.read_text())
            if m.get('reaction_policy')!='native_cleaned_template_extraction_input':
                raise ValueError('Refusing mapper-reaction export for new Core: rerun corrected dataset export')
            prepare(export_dir,bridge,reaction_column=reaction_column)
        bm=json.loads((bridge/'bridge_manifest.json').read_text())
        for f,h in bm['prepared_sha256'].items():
            if sha256_file(bridge/f)!=h:raise ValueError('Source bridge checksum mismatch: '+f)
        cn=sqlite3.connect(tmp/'normalize.sqlite')
        cn.executescript('''
        CREATE TABLE t(h TEXT PRIMARY KEY,smarts TEXT NOT NULL,num TEXT,native_hash TEXT,source_json TEXT);
        CREATE TABLE r(k TEXT PRIMARY KEY,mapped TEXT NOT NULL,n INTEGER NOT NULL);
        CREATE TABLE m(k TEXT,h TEXT,num TEXT,n INTEGER NOT NULL,PRIMARY KEY(k,h));
        ''')
        try:
            for t in _rows(bridge/'source_templates.parquet'):
                h=t.get('retro_template_hash');s=t.get('retro_template_smarts')
                if not h or not s or sha256_text(s)!=h:
                    raise ValueError('Original SMARTS lost or exact hash changed')
                cn.execute('INSERT INTO t(h,smarts,num) VALUES (?,?,?)',(h,s,t.get('retro_template_number')))
            if export_dir:
                for t in _rows(Path(export_dir)/'template_library_v5.parquet'):
                    h=t.get('retro_template_hash')
                    s=t.get('retro_smarts',t.get('retro_template_smarts'))
                    if not h or not s or sha256_text(s)!=h:
                        raise ValueError('Original export has missing/changed SMARTS')
                    cur=cn.execute('UPDATE t SET native_hash=?,source_json=? WHERE h=?',
                                   (str(t.get('native_template_hash')) if t.get('native_template_hash') is not None else None,
                                    json.dumps(t,sort_keys=True,default=str),h))
                    if cur.rowcount!=1:raise ValueError('Unexpected template hash in source export: '+h)
            cn.commit()
            es=pa.schema([('evidence_id',pa.string()),('reaction_key',pa.string()),
                          ('retro_template_hash',pa.string()),('source_observation_key',pa.string()),
                          ('source_reaction_id',pa.string()),('source_dataset',pa.string()),
                          ('source_core_json',pa.string()),('source_membership_json',pa.string())])
            def evidence_rows():
                for i,v in enumerate(_rows(bridge/'source_memberships.parquet')):
                    rxn=_check_forward_smiles(v.get('forward_reaction_smiles'))
                    h=v.get('retro_template_hash');k=sha256_text(rxn)
                    if not cn.execute('SELECT 1 FROM t WHERE h=?',(h,)).fetchone():
                        raise ValueError('Membership references missing template: '+str(h))
                    old=cn.execute('SELECT mapped FROM r WHERE k=?',(k,)).fetchone()
                    if old and old[0]!=rxn:raise ValueError('Reaction hash collision')
                    cn.execute('INSERT INTO r VALUES (?,?,1) ON CONFLICT(k) DO UPDATE SET n=n+1',(k,rxn))
                    num=v.get('retro_template_number')
                    cn.execute('INSERT INTO m VALUES (?,?,?,1) ON CONFLICT(k,h) DO UPDATE SET n=n+1',
                               (k,h,num))
                    sk=str(v.get('source_observation_key') or '')
                    src=str(v.get('source_dataset') or dataset)
                    evid=sha256_text('\0'.join((src,sk,h,str(i))))
                    yield dict(evidence_id=evid,reaction_key=k,retro_template_hash=h,
                               source_observation_key=sk,source_reaction_id=v.get('reaction_id'),
                               source_dataset=src,source_core_json=v.get('source_core_json'),
                               source_membership_json=v.get('source_membership_json'))
                    if (i+1)%50000==0:cn.commit();print(f'[core] source evidence {i+1:,}',flush=True)
            no=_out(tmp/'evidence.parquet',es,evidence_rows())
            cn.commit()
            rs=pa.schema([('reaction_key',pa.string()),('mapped_reaction_smiles',pa.string()),
                          ('evidence_occurrence_count',pa.int64())])
            ms=pa.schema([('reaction_key',pa.string()),('retro_template_hash',pa.string()),
                          ('retro_template_number',pa.string()),('evidence_occurrence_count',pa.int64())])
            ts=pa.schema([('retro_template_hash',pa.string()),('retro_template_smarts',pa.string()),
                          ('retro_template_number',pa.string()),('native_template_hash',pa.string()),
                          ('source_template_json',pa.string())])
            nr=_out(tmp/'reactions.parquet',rs,(
                dict(reaction_key=k,mapped_reaction_smiles=r,evidence_occurrence_count=n)
                for k,r,n in cn.execute('SELECT k,mapped,n FROM r ORDER BY k')))
            nm=_out(tmp/'memberships.parquet',ms,(
                dict(reaction_key=k,retro_template_hash=h,retro_template_number=num,
                     evidence_occurrence_count=n)
                for k,h,num,n in cn.execute('SELECT k,h,num,n FROM m ORDER BY h,k')))
            nt=_out(tmp/'templates.parquet',ts,(
                dict(retro_template_hash=h,retro_template_smarts=s,retro_template_number=num,
                     native_template_hash=native,source_template_json=original)
                for h,s,num,native,original in cn.execute('SELECT h,smarts,num,native_hash,source_json FROM t ORDER BY h')))
            if nt!=bm['counts']['templates'] or no!=bm['counts']['memberships']:
                raise ValueError('Core lost templates or source evidence')
            if cn.execute('SELECT COUNT(*) FROM t WHERE NOT EXISTS(SELECT 1 FROM m WHERE m.h=t.h)').fetchone()[0]:
                raise ValueError('Orphan template after normalization')
            if cn.execute('SELECT SUM(n) FROM m').fetchone()[0]!=no:
                raise ValueError('Aggregated occurrences do not equal source evidence')
        finally:cn.close()
        (tmp/'normalize.sqlite').unlink()
        if export_dir:shutil.rmtree(bridge)
        manifest={'schema_version':FORMAT,'status':'CORE_BUILT_PENDING_AUDIT','created_utc':datetime.now(timezone.utc).isoformat(),
                  'dataset':dataset,
                  'counts':{'reactions':nr,'templates':nt,'memberships':nm,'evidence':no},
                  'identity':{'reaction':'sha256 exact mapped forward reaction string, mapping-sensitive',
                              'template':'sha256 exact authentic retrosynthetic SMARTS',
                              'source':'source observation ID preserved separately'},
                  'source_sha256':bm.get('source_sha256',{}),
                  'scientific_policy':'No v6.5/v6.6/v6.7 structural or semantic representations in Core',
                  'file_sha256':{f:sha256_file(tmp/f) for f in FILE_NAMES}}
        _j(tmp/'manifest.json',manifest)
        tmp.rename(dest)
        return manifest
    except BaseException:
        shutil.rmtree(tmp,ignore_errors=True)
        raise

def audit_core(core,*,rdkit=False,report=None):
    """Independent full-file integrity with counted documentary membership equality."""
    core=Path(core).resolve()
    mf=json.loads((core/'manifest.json').read_text())
    if mf['schema_version']!=FORMAT:raise ValueError('Unsupported Core format')
    for f,d in mf['file_sha256'].items():
        if sha256_file(core/f)!=d:raise ValueError('Core file checksum mismatch: '+f)
    with tempfile.TemporaryDirectory(prefix='rxni_core_audit_') as d:
        cn=sqlite3.connect(str(Path(d)/'audit.sqlite'))
        cn.executescript('''
        CREATE TABLE r(k TEXT PRIMARY KEY);
        CREATE TABLE t(h TEXT PRIMARY KEY);
        CREATE TABLE m(k TEXT,h TEXT,n INTEGER NOT NULL,PRIMARY KEY(k,h));
        CREATE TABLE e(k TEXT,h TEXT,n INTEGER NOT NULL,PRIMARY KEY(k,h));
        CREATE TABLE ids(id TEXT PRIMARY KEY);
        ''')
        nr=nt=nm=ne=0
        try:
            for x in _rows(core/'reactions.parquet'):
                _check_forward_smiles(x['mapped_reaction_smiles'])
                if sha256_text(x['mapped_reaction_smiles'])!=x['reaction_key']:
                    raise ValueError('Wrong reaction identity')
                cn.execute('INSERT INTO r VALUES (?)',(x['reaction_key'],));nr+=1
            for x in _rows(core/'templates.parquet'):
                h=x['retro_template_hash'];s=x['retro_template_smarts']
                if not s or sha256_text(s)!=h:raise ValueError('SMARTS corruption')
                if rdkit:
                    from rdkit.Chem import rdChemReactions
                    if rdChemReactions.ReactionFromSmarts(s) is None:
                        raise ValueError('Invalid original reaction SMARTS: '+h)
                cn.execute('INSERT INTO t VALUES (?)',(h,));nt+=1
            for x in _rows(core/'memberships.parquet'):
                k=x['reaction_key'];h=x['retro_template_hash'];n=x['evidence_occurrence_count']
                if not isinstance(n,int) or n<1:raise ValueError('Invalid occurrence')
                cn.execute('INSERT INTO m VALUES (?,?,?)',(k,h,n));nm+=1
            for x in _rows(core/'evidence.parquet'):
                k=x['reaction_key'];h=x['retro_template_hash']
                cn.execute('INSERT INTO ids VALUES (?)',(x['evidence_id'],))
                cn.execute('INSERT INTO e VALUES (?,?,1) ON CONFLICT(k,h) DO UPDATE SET n=n+1',(k,h));ne+=1
            for sql in ('SELECT COUNT(*) FROM m WHERE k NOT IN (SELECT k FROM r)',
                        'SELECT COUNT(*) FROM m WHERE h NOT IN (SELECT h FROM t)',
                        'SELECT COUNT(*) FROM m LEFT JOIN e USING (k,h) WHERE m.n!=COALESCE(e.n,0)'):
                if cn.execute(sql).fetchone()[0]:raise ValueError('Orphan or inconsistent evidence/membership')
            cnt={'reactions':nr,'templates':nt,'memberships':nm,'evidence':ne}
            if cnt!=mf['counts']:raise ValueError('Count mismatch: '+str(cnt))
            result={'status':'PASS','counts':cnt,'all_original_SMARTS_sha256_verified':True,
                    'orphan_reaction_or_template_keys':0,'documentary_evidence_counts_verified':True,
                    'rdkit_parseability_tested':rdkit}
            if report:_j(report,result)
            return result
        finally:cn.close()
