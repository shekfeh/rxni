"""Create a small Parquet pilot from an old membership source."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .source_adapter import iter_batches, source_schema


def make_pilot(source: str | Path, output: str | Path, *, rows: int = 10000) -> int:
    src = Path(source)
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise FileExistsError(out)
    schema = source_schema(src)
    writer = pq.ParquetWriter(out, schema, compression="zstd")
    written = 0
    try:
        for batch in iter_batches(src, batch_size=min(8192, rows)):
            remaining = rows - written
            if remaining <= 0:
                break
            if batch.num_rows > remaining:
                batch = batch.slice(0, remaining)
            writer.write_batch(batch)
            written += batch.num_rows
            if written >= rows:
                break
    finally:
        writer.close()
    return written


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--rows", type=int, default=10000)
    args = ap.parse_args()
    n = make_pilot(args.source, args.output, rows=args.rows)
    print(f"wrote {n} rows -> {args.output}")


if __name__ == "__main__":
    main()
