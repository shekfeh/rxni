"""Public RXNI native producer adapter; no private runtime_v6_x imports."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
from .bridge import sha256_file

def guard_builder(prepared):
    from ._native.source_adapter import resolve_template_mapping
    from ._native.templates import _load_template_metadata
    from ._native.smarts_integrity import require_source_template_smarts, require_joined_template_smarts
    if not callable(require_source_template_smarts) or not callable(require_joined_template_smarts):
        raise RuntimeError("Missing original SMARTS integrity validators")
    p=Path(prepared)/"source_templates.parquet"
    if not p.exists(): raise FileNotFoundError(p)
    mapped=resolve_template_mapping(str(p))
    if mapped.template_smarts!="retro_template_smarts":
        raise RuntimeError(f"Unexpected source SMARTS field: {mapped.template_smarts}")
    meta, stats=_load_template_metadata(str(p),4096)
    bad=[h for h,x in meta.items() if not x.get("retro_template_smarts") or not x["retro_template_smarts"].strip()]
    if bad: raise RuntimeError(f"{len(bad)} missing source SMARTS; examples={bad[:3]}")
    expected=json.loads((Path(prepared)/"bridge_manifest.json").read_text())["counts"]["templates"]
    if len(meta)!=expected: raise RuntimeError(f"Source templates {len(meta)} != bridge {expected}")
    return {"loaded_templates":len(meta),"resolved_smarts_column":mapped.template_smarts,"guard":"public_native_builder"}

def run_build(prepared,out,*,root=None,workers=4,row_group_size=4096,validation_sample=10000,resume=False):
    from ._native import build_all
    prep=Path(prepared).resolve(); out=Path(out).resolve()
    mf=prep/"bridge_manifest.json"
    if not mf.is_file(): raise FileNotFoundError(mf)
    m=json.loads(mf.read_text())
    for name,sha in m["prepared_sha256"].items():
        if sha256_file(prep/name)!=sha: raise RuntimeError("Prepared source changed: "+name)
    src=Path(build_all.__file__).resolve()
    lock=out/"rxni_library_builder_request.json"
    cfg={"bridge_manifest_sha256":sha256_file(mf),
         "membership_sha256":m["prepared_sha256"]["source_memberships.parquet"],
         "template_sha256":m["prepared_sha256"]["source_templates.parquet"],
         "workers":workers,"row_group_size":row_group_size,"validation_sample":validation_sample,
         "builder_module":"rxni.builder._native.build_all", "builder_sha256":sha256_file(src)}
    if out.exists():
        if not resume: raise FileExistsError(out)
        if not lock.is_file() or json.loads(lock.read_text())!=cfg:
            raise RuntimeError("Unsafe resume: provenance or parameters changed")
    else:
        if resume: raise FileNotFoundError(out)
        out.mkdir(parents=True)
        lock.write_text(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
    guard=guard_builder(prep)
    cmd=[sys.executable,"-m","rxni.builder._native.build_all",
         "--membership-source",str(prep/"source_memberships.parquet"),
         "--template-source",str(prep/"source_templates.parquet"),
         "--output-dir",str(out),"--workers",str(workers),
         "--row-group-size",str(row_group_size),
         "--validation-sample",str(validation_sample)]
    if resume: cmd.append("--resume")
    subprocess.run(cmd, check=True)
    return {"status":"BUILDER_PROCESS_EXIT_0_VERIFY_REQUIRED","output":str(out),"guard":guard}
