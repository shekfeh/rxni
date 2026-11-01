from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from rxni.public import (
    AnalysisOptions,
    ChemKNScope,
    ConditionInput,
    ReactionIntelligenceService,
    model_dump_compat,
)

from .render import (
    render_analysis,
    render_chemkn,
    render_conditions,
    render_diagnostics,
    render_json,
    render_methods,
    render_retrieval,
    render_semantics,
    render_summary,
)

CLI_VERSION = "rxni-cli/0.1.1"

app = typer.Typer(
    name="ri",
    no_args_is_help=True,
    add_completion=False,
    help="RXNI reaction analysis, retrieval and reaction-method CLI.",
)
console = Console()


@dataclass
class AppContext:
    library: Path
    mapping_control: str
    retrieval_semantic_policy: str | None
    service: ReactionIntelligenceService | None = None


def _resolve_library(value: Optional[Path]) -> Path:
    if value is not None:
        return value.expanduser()
    env = os.environ.get("REACTION_INTELLIGENCE_LIBRARY") or os.environ.get("V5")
    if env:
        return Path(env).expanduser()
    return Path("artefacts/libraries/template_library_v5_native_v01")


@app.callback()
def root(
    ctx: typer.Context,
    library: Optional[Path] = typer.Option(None, "--library", help="Library-v5 directory. Defaults to $REACTION_INTELLIGENCE_LIBRARY, $V5, then the project-relative V5 path."),
    mapping_control: str = typer.Option("conversion_validated", "--mapping-control", help="Semantic mapping-control policy."),
    retrieval_semantic_policy: Optional[str] = typer.Option(None, "--retrieval-semantic-policy", help="Stored Retrieval-V2 semantic lineage, e.g. strict."),
    version: bool = typer.Option(False, "--version", is_eager=True, help="Show CLI version and exit."),
) -> None:
    if version:
        console.print(CLI_VERSION)
        raise typer.Exit()
    ctx.obj = AppContext(
        library=_resolve_library(library),
        mapping_control=mapping_control,
        retrieval_semantic_policy=(
            retrieval_semantic_policy
            or os.environ.get("RI_RETRIEVAL_REPRESENTATION_SEMANTIC_POLICY")
        ),
    )


def _ctx(ctx: typer.Context) -> AppContext:
    if not isinstance(ctx.obj, AppContext):
        raise RuntimeError("CLI context not initialized")
    return ctx.obj


def _service(ctx: typer.Context) -> ReactionIntelligenceService:
    cfg = _ctx(ctx)
    if cfg.service is None:
        cfg.service = ReactionIntelligenceService(
            library=cfg.library,
            mapping_control=cfg.mapping_control,
            retrieval_representation_semantic_policy=cfg.retrieval_semantic_policy,
        )
    return cfg.service


def _reaction(reaction: Optional[str], reaction_file: Optional[Path]) -> str:
    if reaction:
        return reaction.strip()
    if reaction_file:
        for line in reaction_file.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                return line
        raise typer.BadParameter(f"No reaction found in {reaction_file}")
    raise typer.BadParameter("Provide --reaction/-r or --reaction-file")


def _split(values: list[str]) -> list[str]:
    out: list[str] = []
    for value in values:
        out.extend(x.strip() for x in value.split(",") if x.strip())
    return out


def _conditions(
    reagents: list[str], solvents: list[str], catalysts: list[str], bases: list[str],
    additives: list[str], acids: list[str], temperature_c: Optional[float],
    reaction_time_h: Optional[float], atmosphere: Optional[str], pressure: Optional[str],
) -> ConditionInput:
    return ConditionInput(
        reagents=_split(reagents), solvents=_split(solvents), catalysts=_split(catalysts),
        bases=_split(bases), additives=_split(additives), acids=_split(acids),
        temperature_c=temperature_c, reaction_time_h=reaction_time_h,
        atmosphere=atmosphere, pressure=pressure,
    )


def _options(
    *, mode: str = "sem_multi", top_k: int = 20, candidate_budget: int = 3000,
    member_limit: Optional[int] = None, reagent_weight: float = 0.20,
    profile_methods: bool = False, diagnostics: bool = False,
    include_chemkn: bool = False, chemkn_scope: str = "all", chemkn_card_limit: int = 8,
) -> AnalysisOptions:
    return AnalysisOptions(
        mode=mode,
        top_k=top_k,
        candidate_budget=candidate_budget,
        member_limit=member_limit,
        reagent_weight=reagent_weight,
        profile_methods=profile_methods,
        include_diagnostics=diagnostics,
        include_chemkn=include_chemkn,
        chemkn_scope=ChemKNScope(chemkn_scope),
        chemkn_card_limit=chemkn_card_limit,
    )


def _write_json(path: Optional[Path], result) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(model_dump_compat(result), indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    console.print(f"[dim]Wrote JSON: {path}[/dim]")


@app.command()
def analyze(
    ctx: typer.Context,
    reaction: Optional[str] = typer.Option(None, "--reaction", "-r"),
    reaction_file: Optional[Path] = typer.Option(None, "--reaction-file"),
    mode: str = typer.Option("sem_multi", "--mode"),
    top_k: int = typer.Option(20, "--top-k", min=1),
    candidate_budget: int = typer.Option(3000, "--candidate-budget", min=1),
    member_limit: Optional[int] = typer.Option(None, "--member-limit", min=0),
    reagent: list[str] = typer.Option([], "--reagent"),
    solvent: list[str] = typer.Option([], "--solvent"),
    catalyst: list[str] = typer.Option([], "--catalyst"),
    base: list[str] = typer.Option([], "--base"),
    additive: list[str] = typer.Option([], "--additive"),
    acid: list[str] = typer.Option([], "--acid"),
    temperature_c: Optional[float] = typer.Option(None, "--temperature-c"),
    reaction_time_h: Optional[float] = typer.Option(None, "--reaction-time-h"),
    atmosphere: Optional[str] = typer.Option(None, "--atmosphere"),
    pressure: Optional[str] = typer.Option(None, "--pressure"),
    reagent_weight: float = typer.Option(0.20, "--reagent-weight", min=0.0, max=1.0),
    methods: bool = typer.Option(False, "--methods", help="Build and show method landscapes."),
    chemkn: bool = typer.Option(False, "--chemkn", help="Enrich method landscape with ChemKN cards."),
    chemkn_scope: str = typer.Option("all", "--chemkn-scope", help="all|exact|related|none"),
    chemkn_card_limit: int = typer.Option(8, "--chemkn-card-limit", min=0),
    full: bool = typer.Option(False, "--full", help="Show conditions, graph/semantic detail and diagnostics."),
    candidates: int = typer.Option(0, "--candidates", min=0, help="Number of candidate rows to display."),
    templates: int = typer.Option(10, "--templates", min=0, help="Number of template rows to display."),
    json_output: bool = typer.Option(False, "--json", help="Print typed JSON instead of Rich report."),
    json_out: Optional[Path] = typer.Option(None, "--json-out"),
) -> None:
    rxn = _reaction(reaction, reaction_file)
    cond = _conditions(reagent, solvent, catalyst, base, additive, acid, temperature_c, reaction_time_h, atmosphere, pressure)
    opts = _options(
        mode=mode, top_k=top_k, candidate_budget=candidate_budget, member_limit=member_limit,
        reagent_weight=reagent_weight, profile_methods=(methods or chemkn), diagnostics=full,
        include_chemkn=chemkn, chemkn_scope=chemkn_scope, chemkn_card_limit=chemkn_card_limit,
    )
    result = _service(ctx).analyze(rxn, conditions=cond, options=opts)
    _write_json(json_out, result)
    if json_output:
        render_json(console, result)
    else:
        render_analysis(console, result, full=full, template_limit=templates, candidate_limit=candidates, methods=methods or chemkn, chemkn=chemkn)


@app.command()
def semantics(
    ctx: typer.Context,
    reaction: Optional[str] = typer.Option(None, "--reaction", "-r"),
    reaction_file: Optional[Path] = typer.Option(None, "--reaction-file"),
    deep: bool = typer.Option(False, "--deep", help="Include graph/topology/validation fields."),
    raw: bool = typer.Option(False, "--raw", help="Also print complete canonical semantics JSON and detector outputs."),
) -> None:
    result = _service(ctx).analyze(
        _reaction(reaction, reaction_file),
        options=_options(top_k=20, candidate_budget=3000, diagnostics=(deep or raw)),
    )
    render_summary(console, result)
    render_semantics(console, result, deep=deep or raw)
    if raw:
        render_diagnostics(console, result, raw_semantics=True)


@app.command()
def retrieve(
    ctx: typer.Context,
    reaction: Optional[str] = typer.Option(None, "--reaction", "-r"),
    reaction_file: Optional[Path] = typer.Option(None, "--reaction-file"),
    mode: str = typer.Option("sem_multi", "--mode"),
    top_k: int = typer.Option(20, "--top-k", min=1),
    candidate_budget: int = typer.Option(3000, "--candidate-budget", min=1),
    templates: int = typer.Option(20, "--templates", min=0),
    candidates: int = typer.Option(20, "--candidates", min=0),
) -> None:
    result = _service(ctx).analyze(
        _reaction(reaction, reaction_file),
        options=_options(mode=mode, top_k=top_k, candidate_budget=candidate_budget),
    )
    render_summary(console, result)
    render_retrieval(console, result, template_limit=templates, candidate_limit=candidates)


def _method_result(
    ctx: typer.Context, reaction: Optional[str], reaction_file: Optional[Path],
    reagent: list[str], solvent: list[str], catalyst: list[str], base: list[str],
    additive: list[str], acid: list[str], temperature_c: Optional[float],
    reaction_time_h: Optional[float], atmosphere: Optional[str], pressure: Optional[str],
    *, chemkn: bool, chemkn_scope: str = "all", chemkn_card_limit: int = 8,
):
    cond = _conditions(reagent, solvent, catalyst, base, additive, acid, temperature_c, reaction_time_h, atmosphere, pressure)
    opts = _options(profile_methods=True, diagnostics=True, include_chemkn=chemkn, chemkn_scope=chemkn_scope, chemkn_card_limit=chemkn_card_limit)
    return _service(ctx).analyze(_reaction(reaction, reaction_file), conditions=cond, options=opts)


@app.command()
def methods(
    ctx: typer.Context,
    reaction: Optional[str] = typer.Option(None, "--reaction", "-r"),
    reaction_file: Optional[Path] = typer.Option(None, "--reaction-file"),
    reagent: list[str] = typer.Option([], "--reagent"), solvent: list[str] = typer.Option([], "--solvent"),
    catalyst: list[str] = typer.Option([], "--catalyst"), base: list[str] = typer.Option([], "--base"),
    additive: list[str] = typer.Option([], "--additive"), acid: list[str] = typer.Option([], "--acid"),
    temperature_c: Optional[float] = typer.Option(None, "--temperature-c"), reaction_time_h: Optional[float] = typer.Option(None, "--reaction-time-h"),
    atmosphere: Optional[str] = typer.Option(None, "--atmosphere"), pressure: Optional[str] = typer.Option(None, "--pressure"),
    limit: int = typer.Option(15, "--limit", min=1),
) -> None:
    result = _method_result(ctx, reaction, reaction_file, reagent, solvent, catalyst, base, additive, acid, temperature_c, reaction_time_h, atmosphere, pressure, chemkn=False)
    render_summary(console, result)
    render_conditions(console, result)
    render_methods(console, result, limit=limit)


@app.command()
def chemkn(
    ctx: typer.Context,
    reaction: Optional[str] = typer.Option(None, "--reaction", "-r"),
    reaction_file: Optional[Path] = typer.Option(None, "--reaction-file"),
    reagent: list[str] = typer.Option([], "--reagent"), solvent: list[str] = typer.Option([], "--solvent"),
    catalyst: list[str] = typer.Option([], "--catalyst"), base: list[str] = typer.Option([], "--base"),
    additive: list[str] = typer.Option([], "--additive"), acid: list[str] = typer.Option([], "--acid"),
    temperature_c: Optional[float] = typer.Option(None, "--temperature-c"), reaction_time_h: Optional[float] = typer.Option(None, "--reaction-time-h"),
    atmosphere: Optional[str] = typer.Option(None, "--atmosphere"), pressure: Optional[str] = typer.Option(None, "--pressure"),
    scope: str = typer.Option("all", "--scope", help="all|exact|related|none"),
    limit: int = typer.Option(8, "--limit", min=0),
) -> None:
    result = _method_result(ctx, reaction, reaction_file, reagent, solvent, catalyst, base, additive, acid, temperature_c, reaction_time_h, atmosphere, pressure, chemkn=True, chemkn_scope=scope, chemkn_card_limit=limit)
    render_summary(console, result)
    render_methods(console, result)
    render_chemkn(console, result, limit=limit)


@app.command()
def inspect(
    ctx: typer.Context,
    reaction: Optional[str] = typer.Option(None, "--reaction", "-r"),
    reaction_file: Optional[Path] = typer.Option(None, "--reaction-file"),
    methods: bool = typer.Option(True, "--methods/--no-methods"),
    chemkn: bool = typer.Option(False, "--chemkn"),
    raw: bool = typer.Option(False, "--raw", help="Include full canonical semantics and detector output."),
    candidates: int = typer.Option(20, "--candidates", min=0),
) -> None:
    result = _service(ctx).analyze(
        _reaction(reaction, reaction_file),
        options=_options(profile_methods=(methods or chemkn), diagnostics=True, include_chemkn=chemkn),
    )
    render_summary(console, result)
    render_semantics(console, result, deep=True)
    render_retrieval(console, result, template_limit=20, candidate_limit=candidates)
    render_conditions(console, result)
    if methods or chemkn:
        render_methods(console, result, limit=20)
    if chemkn:
        render_chemkn(console, result, limit=12)
    render_diagnostics(console, result, raw_semantics=raw)


@app.command()
def doctor(ctx: typer.Context) -> None:
    cfg = _ctx(ctx)
    rows: list[tuple[str, str, bool]] = []
    rows.append(("Library", str(cfg.library), cfg.library.exists()))
    for name in ("structural_index_v1.sqlite", "runtime_lookup_v1.sqlite", "retrieval_representations_v2.parquet", "retrieval_v2_manifest.json"):
        p = cfg.library / name
        rows.append((name, str(p), p.exists()))
    try:
        import fastapi  # noqa: F401
        rows.append(("FastAPI", "import", True))
    except Exception:
        rows.append(("FastAPI", "import", False))
    try:
        import typer as _typer  # noqa: F401
        import rich as _rich  # noqa: F401
        rows.append(("Typer + Rich", "import", True))
    except Exception:
        rows.append(("Typer + Rich", "import", False))
    web = Path.cwd() / "webapp" / "package.json"
    rows.append(("React webapp", str(web), web.exists()))

    table = Table(title="Reaction Intelligence doctor")
    table.add_column("Check")
    table.add_column("Location")
    table.add_column("Status")
    failed = False
    for name, loc, ok in rows:
        table.add_row(name, loc, "[green]PASS[/green]" if ok else "[red]MISSING[/red]")
        failed |= not ok
    console.print(table)
    if failed:
        raise typer.Exit(code=1)


@app.command()
def serve(
    ctx: typer.Context,
    host: str = typer.Option("127.0.0.1", "--host"),
    port: int = typer.Option(8000, "--port"),
) -> None:
    cfg = _ctx(ctx)
    os.environ["REACTION_INTELLIGENCE_LIBRARY"] = str(cfg.library)
    os.environ["RUNTIME_V67_MAPPING_CONTROL"] = cfg.mapping_control
    os.environ.setdefault("RI_SERIALIZE_ANALYSIS", "1")
    if cfg.retrieval_semantic_policy:
        os.environ["RI_RETRIEVAL_REPRESENTATION_SEMANTIC_POLICY"] = cfg.retrieval_semantic_policy
    try:
        import uvicorn
    except Exception as exc:
        console.print(Panel(str(exc), title="Uvicorn unavailable", border_style="red"))
        raise typer.Exit(1)
    uvicorn.run("rxni.api.app:app", host=host, port=port, workers=1)


@app.command()
def web(
    directory: Path = typer.Option(Path("webapp"), "--directory"),
    host: str = typer.Option("127.0.0.1", "--host"),
    port: int = typer.Option(5173, "--port"),
) -> None:
    if not (directory / "package.json").exists():
        raise typer.BadParameter(f"No webapp/package.json under {directory}")
    cmd = ["npm", "run", "dev", "--", "--host", host, "--port", str(port)]
    console.print("[dim]" + " ".join(cmd) + "[/dim]")
    raise typer.Exit(subprocess.call(cmd, cwd=directory))


if __name__ == "__main__":
    app()
