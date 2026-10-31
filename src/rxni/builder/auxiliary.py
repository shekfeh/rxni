"""v6.7-aware Auxiliary generator: isolated scientific build, no legacy overwrite.

v6.7 shares v6.6 stored retrieval-v2 & ranker but has an updated semantic
query policy. Do not infer library-side semantic compatibility from filenames.
"""
from __future__ import annotations
import json, shutil, sqlite3, subprocess, sys, tempfile
from pathlib import Path
from dataclasses import asdict
from .bridge import sha256_file
from .knowledge import audit_core, _rows, _out, _j
from .profiles import select
from .native_runtime import run_build
from .native_verify import verify

def runtime_inspection(root,runtime='v6.7'):
    prof=select(runtime);r=Path(root).resolve()
    public_package=Path(__file__).resolve().parent
    actual={'native_builder':public_package/'_native/build_all.py',
            'retrieval_v2_builder':public_package.parent/'retrieval/builders/representations.py',
            'lookup_builder':public_package/'lookup_cache.py'}
    return {'profile':asdict(prof),'producer_paths':{k:str(v) for k,v in actual.items()},
            'missing':{k:str(v) for k,v in actual.items() if not v.is_file()},
            'runtime_lookup_invocation':'inspect installed producer help; no fabricated CLI',
            'semantic_compatibility':'requires library-side profile and executed retrieval fixture audit'}

def prepare_runtime_inputs(core_dir,out_dir):
    import pyarrow as pa
    core=Path(core_dir).resolve();out=Path(out_dir).resolve()
    if out.exists():raise FileExistsError(out)
    audit=audit_core(core)
    out.mkdir(parents=True)
    cx=sqlite3.connect(out/'reaction_join.sqlite')
    try:
        cx.execute('CREATE TABLE r(k TEXT PRIMARY KEY,rxn TEXT NOT NULL)')
        for x in _rows(core/'reactions.parquet'):
            cx.execute('INSERT INTO r VALUES (?,?)',(x['reaction_key'],x['mapped_reaction_smiles']))
        cx.commit()
        ts=pa.schema([('retro_template_hash',pa.string()),('retro_template_number',pa.string()),
                      ('retro_template_smarts',pa.string())])
        _out(out/'source_templates.parquet',ts,(
            dict(retro_template_hash=t['retro_template_hash'],retro_template_number=t['retro_template_number'],
                 retro_template_smarts=t['retro_template_smarts'])
            for t in _rows(core/'templates.parquet')))
        tm={x['retro_template_hash']:x['retro_template_number'] for x in _rows(core/'templates.parquet')}
        ms=pa.schema([('forward_reaction_smiles',pa.string()),('retro_template_hash',pa.string()),
                      ('retro_template_number',pa.string()),('reaction_id',pa.string()),
                      ('source_observation_key',pa.string()),('source_dataset',pa.string()),
                      ('source_core_json',pa.string()),('source_membership_json',pa.string())])
        def gen():
            for r in _rows(core/'evidence.parquet'):
                found=cx.execute('SELECT rxn FROM r WHERE k=?',(r['reaction_key'],)).fetchone()
                if found is None:raise ValueError('Missing reaction in Core evidence')
                h=r['retro_template_hash']
                if h not in tm:raise ValueError('Missing template in Core evidence')
                yield {'forward_reaction_smiles':found[0],'retro_template_hash':h,
                       'retro_template_number':tm[h],'reaction_id':r['source_reaction_id'],
                       'source_observation_key':r['source_observation_key'],
                       'source_dataset':r['source_dataset'],
                       'source_core_json':r['source_core_json'],
                       'source_membership_json':r['source_membership_json']}
        no=_out(out/'source_memberships.parquet',ms,gen())
        if no!=audit['counts']['evidence']:raise ValueError('Lost documentary evidence')
        bridge={'counts':{'templates':audit['counts']['templates'],'memberships':no},
                'prepared_sha256':{n:sha256_file(out/n) for n in
                    ('source_memberships.parquet','source_templates.parquet')},
                'source_sha256':json.loads((core/'manifest.json').read_text())['file_sha256'],
                'core_manifest_sha256':sha256_file(core/'manifest.json')}
        _j(out/'bridge_manifest.json',bridge)
        return bridge
    except BaseException:
        cx.close();shutil.rmtree(out,ignore_errors=True);raise
    finally:
        try:cx.close()
        except Exception:pass
        (out/'reaction_join.sqlite').unlink(missing_ok=True)

def build_auxiliary(*,core,work_dir,root,runtime='v6.7',workers=4,
                    validation_sample=500,chunk_size=2048,resume=False,through='v2'):
    if through not in ('native','v2'):raise ValueError('--through must be native or v2')
    prof=select(runtime)
    root=Path(root).resolve() if root else Path.cwd().resolve();core=Path(core).resolve();work=Path(work_dir).resolve()
    # Production is resolved from installed public modules, not the private checkout.
    inspected=runtime_inspection(root,runtime)
    if inspected['missing']:raise FileNotFoundError(inspected['missing'])
    mf=json.loads((core/'manifest.json').read_text())
    lock={'core_manifest_sha256':sha256_file(core/'manifest.json'),
          'core_data_sha256':mf['file_sha256'],'profile':asdict(prof),
          'workers':workers,'validation_sample':validation_sample,'chunk_size':chunk_size}
    if work.exists():
        if not resume or not (work/'build_request.json').exists() or json.loads((work/'build_request.json').read_text())!=lock:
            raise ValueError('Workspace exists or unsafe resume: '+str(work))
    else:
        if resume:raise FileNotFoundError(work)
        work.mkdir(parents=True)
        _j(work/'build_request.json',lock)
    inp=work/'inputs';nat=work/'native';v2=work/'retrieval_v2'
    if not inp.exists():prepare_runtime_inputs(core,inp)
    else:
        chk=json.loads((inp/'bridge_manifest.json').read_text())
        if chk['core_manifest_sha256']!=lock['core_manifest_sha256']:raise ValueError('Cached inputs use another Core')
        for f,d in chk['prepared_sha256'].items():
            if sha256_file(inp/f)!=d:raise ValueError('Cached source mutated: '+f)
    run_build(inp,nat,root=root,workers=workers,validation_sample=validation_sample,resume=resume and nat.exists())
    verified=verify(inp,nat)
    if verified['status']!='PASS':raise ValueError('Native template preservation failed')
    _j(work/'native_verification.json',verified)
    if through=='native':
        return {'status':'NATIVE_PREREQUISITE_PASS','semantic_profile':'NOT_ATTESTED',
                'source_templates':verified['source_templates'],'workspace':str(work)}
    # This is the real installed v6.7 representation materializer, not a new implementation.
    if not (resume and (v2/'retrieval_v2_manifest.json').is_file()):
        cmd=[sys.executable,'-m',prof.v2_module,'--library',str(nat),'--output-dir',str(v2),
             '--workers',str(workers),'--chunk-size',str(chunk_size)]
        if resume and v2.exists():cmd.append('--resume')
        subprocess.run(cmd,cwd=root,check=True)
    for f in ('retrieval_representations_v2.parquet','representation_locator_v2.sqlite','retrieval_v2_manifest.json'):
        if not (v2/f).is_file():raise FileNotFoundError(v2/f)
    # Reuse the unchanged V5 retrieval member-to-lookup producer on this exact native build.
    lookup_db=nat/'runtime_lookup_v1.sqlite'
    if not lookup_db.is_file():
        subprocess.run([sys.executable,'-m','rxni.builder.lookup_cache',
                        '--library',str(nat)],check=True)
    if not _check_sqlite(lookup_db):raise RuntimeError('Invalid newly built lookup SQLite')
    report={'status':'CANDIDATE_SOURCES_READY_PENDING_PARITY','runtime':runtime,
            'workspace':str(work),'core_manifest_sha256':lock['core_manifest_sha256'],
            'semantic_profile_compatibility':'NOT_VALIDATED',
            'runtime_lookup':'BUILT_FROM_CURRENT_NATIVE',
            'runtime_lookup_sha256':sha256_file(lookup_db),
            'functional_ranking_parity':'NOT_TESTED',
            'v2_manifest_sha256':sha256_file(v2/'retrieval_v2_manifest.json')}
    _j(work/'auxiliary_build_status.json',report)
    return report

def _check_sqlite(p):
    cx=sqlite3.connect(f'file:{Path(p).resolve()}?mode=ro',uri=True)
    try:
        return cx.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    finally:cx.close()

def _sqlite_snapshot(source, destination):
    # Safe snapshot even when upstream source uses WAL; do not copy bare .sqlite.
    src=sqlite3.connect(f'file:{Path(source).resolve()}?mode=ro',uri=True)
    dst=sqlite3.connect(destination)
    try:src.backup(dst)
    finally:dst.close();src.close()
    if not _check_sqlite(destination):raise ValueError('Snapshot SQLite integrity check failed')

def assemble_candidate(*,core,work_dir,lookup_db,out_dir,root,runtime='v6.7'):
    """Non-production 5-artifact candidate; does not assert semantic/ranking validation."""
    prof=select(runtime)
    core=Path(core).resolve();work=Path(work_dir).resolve();out=Path(out_dir).resolve()
    if out.exists():raise FileExistsError(out)
    if not (work/'auxiliary_build_status.json').is_file():raise FileNotFoundError('Run build-auxiliary --through v2')
    req=json.loads((work/'build_request.json').read_text())
    if req['core_manifest_sha256']!=sha256_file(core/'manifest.json'):
        raise ValueError('Core modified after auxiliary build')
    lookup=Path(lookup_db).resolve()
    aux_status=json.loads((work/'auxiliary_build_status.json').read_text())
    if aux_status.get('runtime_lookup_sha256')!=sha256_file(lookup):
        raise ValueError('Lookup is not the locked output from this auxiliary build')
    if not lookup.is_relative_to(work):
        raise ValueError('Runtime lookup must be produced inside this dataset workspace; refusing historical/external database')
    sources={'structural_index_v1.sqlite':work/'native/structural_index_v1.sqlite',
             'runtime_lookup_v1.sqlite':lookup,
             'retrieval_representations_v2.parquet':work/'retrieval_v2/retrieval_representations_v2.parquet',
             'representation_locator_v2.sqlite':work/'retrieval_v2/representation_locator_v2.sqlite',
             'retrieval_v2_manifest.json':work/'retrieval_v2/retrieval_v2_manifest.json'}
    for name,f in sources.items():
        if not f.is_file():raise FileNotFoundError(f)
        if f.suffix=='.sqlite' and not _check_sqlite(f):raise ValueError('SQLite integrity failed: '+name)
    out.parent.mkdir(parents=True,exist_ok=True)
    temp=Path(tempfile.mkdtemp(prefix=out.name+'.BUILDING.',dir=out.parent))
    try:
        for name,f in sources.items():
            if f.suffix=='.sqlite':_sqlite_snapshot(f,temp/name)
            else:shutil.copy2(f,temp/name)
        report={'status':'CANDIDATE_PENDING_SEMANTIC_AND_FUNCTIONAL_VALIDATION',
                'runtime':runtime,'profile':asdict(prof),
                'core_manifest_sha256':req['core_manifest_sha256'],
                'artifacts_sha256':{n:sha256_file(temp/n) for n in sources},
                'retrieval_functional_fixtures':'NOT_EXECUTED',
                'semantic_profile_parity':'NOT_ESTABLISHED',
                'lookup_origin':str(sources['runtime_lookup_v1.sqlite'])}
        _j(temp/'candidate_manifest.json',report)
        temp.rename(out)
        return report
    except BaseException:
        shutil.rmtree(temp,ignore_errors=True);raise
