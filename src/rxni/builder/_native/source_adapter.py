"""Streaming adapters for old reaction/template source files.

Supported inputs: Parquet, CSV/CSV.GZ, TSV/TSV.GZ.
Only identity/membership/provenance columns are authoritative. Any old derived
chemistry columns are reported and ignored by the native V5 build.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from .schema import FORBIDDEN_DERIVED_TOKENS, MEMBERSHIP_ALIASES, TEMPLATE_SOURCE_ALIASES


@dataclass(frozen=True)
class SourceMapping:
    reaction: str | None = None
    template_hash: str | None = None
    template_number: str | None = None
    reaction_id: str | None = None
    reaction_hash: str | None = None
    template_smarts: str | None = None

    def as_dict(self) -> dict[str, str | None]:
        return {
            "reaction": self.reaction,
            "template_hash": self.template_hash,
            "template_number": self.template_number,
            "reaction_id": self.reaction_id,
            "reaction_hash": self.reaction_hash,
            "template_smarts": self.template_smarts,
        }


def _choose(names: set[str], aliases: tuple[str, ...], *, required: bool) -> str | None:
    for name in aliases:
        if name in names:
            return name
    if required:
        raise KeyError(f"None of required aliases found: {aliases}; columns={sorted(names)}")
    return None


def derived_columns(names: set[str]) -> list[str]:
    out = []
    for name in sorted(names):
        low = name.lower()
        if any(token in low for token in FORBIDDEN_DERIVED_TOKENS):
            out.append(name)
    return out


def _suffix(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".parquet") or name.endswith(".pq"):
        return "parquet"
    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        return "tsv"
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return "csv"
    raise ValueError(f"Unsupported source extension: {path}")


def source_schema(path: str | Path) -> pa.Schema:
    p = Path(path)
    kind = _suffix(p)
    if kind == "parquet":
        return pq.ParquetFile(p).schema_arrow
    delim = "\t" if kind == "tsv" else ","
    reader = pacsv.open_csv(
        p,
        parse_options=pacsv.ParseOptions(delimiter=delim),
        read_options=pacsv.ReadOptions(block_size=1 << 20),
    )
    return reader.schema


def resolve_membership_mapping(path: str | Path) -> SourceMapping:
    names = set(source_schema(path).names)
    return SourceMapping(
        reaction=_choose(names, MEMBERSHIP_ALIASES["reaction"], required=True),
        template_hash=_choose(names, MEMBERSHIP_ALIASES["template_hash"], required=True),
        template_number=_choose(names, MEMBERSHIP_ALIASES["template_number"], required=False),
        reaction_id=_choose(names, MEMBERSHIP_ALIASES["reaction_id"], required=False),
        reaction_hash=_choose(names, MEMBERSHIP_ALIASES["reaction_hash"], required=False),
    )


def resolve_template_mapping(path: str | Path) -> SourceMapping:
    names = set(source_schema(path).names)
    return SourceMapping(
        template_hash=_choose(names, TEMPLATE_SOURCE_ALIASES["template_hash"], required=True),
        template_number=_choose(names, TEMPLATE_SOURCE_ALIASES["template_number"], required=False),
        template_smarts=_choose(names, TEMPLATE_SOURCE_ALIASES["template_smarts"], required=True),
    )


def iter_batches(path: str | Path, *, batch_size: int = 8192, columns: list[str] | None = None) -> Iterator[pa.RecordBatch]:
    p = Path(path)
    kind = _suffix(p)
    if kind == "parquet":
        yield from pq.ParquetFile(p).iter_batches(batch_size=batch_size, columns=columns)
        return

    delim = "\t" if kind == "tsv" else ","
    convert = pacsv.ConvertOptions(include_columns=columns) if columns else pacsv.ConvertOptions()
    reader = pacsv.open_csv(
        p,
        parse_options=pacsv.ParseOptions(delimiter=delim),
        read_options=pacsv.ReadOptions(block_size=max(1 << 20, batch_size * 256)),
        convert_options=convert,
    )
    for batch in reader:
        yield batch


def inspect_source(path: str | Path, *, kind: str = "membership") -> dict[str, Any]:
    p = Path(path)
    schema = source_schema(p)
    names = set(schema.names)
    mapping = resolve_membership_mapping(p) if kind == "membership" else resolve_template_mapping(p)
    return {
        "path": str(p),
        "format": _suffix(p),
        "columns": schema.names,
        "types": {f.name: str(f.type) for f in schema},
        "resolved_mapping": mapping.as_dict(),
        "ignored_derived_columns": derived_columns(names),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Inspect/resolve an old source file for native V5 ingestion")
    ap.add_argument("--source", required=True)
    ap.add_argument("--kind", choices=("membership", "template"), default="membership")
    args = ap.parse_args()
    print(json.dumps(inspect_source(args.source, kind=args.kind), indent=2))


if __name__ == "__main__":
    main()
