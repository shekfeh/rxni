# Getting started — RXNI v0.1.x

Use a fresh Python 3.10+ virtual environment. The independently tested `0.1.1` code is a maintenance candidate at Git commit `5e030e388118024627f98aa4de7d58fd6ee706c8`; the published software tag as audited remains `v0.1.0`.

## 1. Install and inspect

From a cloned RXNI source checkout:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -e '.[api,cli]'
rxni --version
rxni retrieve --help
rxni analyze --help
rxni conditions --help
rxni-build-library --help
```

The software installation does **not** download historical reactions, commercial datasets, ChemKN knowledge cards or RXNMapper. Optional external mapping software belongs in a separate environment.

## 2. Obtain and verify public data

Follow [Tutorial 01](tutorials/01-download-verify-library.md) to fetch the publicly available [RXNI Historical USPTO library](https://zenodo.org/records/23856755), verify the 1,444,488,348-byte archive with the published SHA256, check all internal files, and run `rxni-verify-library --smoke` on installed `0.1.1` code.

The examples below expect:

```bash
export RXNI_LIBRARY="$HOME/rxni-data/uspto_1976_2016_rxni_library"
```

This variable is used by the tutorial shell commands. The local HTTP server uses `REACTION_INTELLIGENCE_LIBRARY` instead.

## 3. Select a mapped reaction that exists in the library

```bash
REACTION=$(python - "$RXNI_LIBRARY" <<'PY'
import sys
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
pf = pq.ParquetFile(lib.representations)
print(pf.read_row_group(0, columns=['reaction'])['reaction'][0].as_py())
PY
)
```

The query must retain meaningful reactant/product atom mapping. In the audited prepatch `0.1.1` implementation, **unmapped** queries can produce the signature `NO_MAPPED_BOND_EDIT` and misleading unrelated candidates. Strict rejection of unmapped queries is a **planned code patch**, not behavior already verified in `0.1.1`.

## 4. Retrieve candidates

```bash
rxni retrieve --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --mode sem_multi --top-k 5 --candidate-budget 500 --pager never
```

This intentionally queries a stored reaction and can produce an exact self-hit. Use it only as a read-side smoke test, not as a held-out scientific benchmark. Explore [Tutorial 02](tutorials/02-retrieval.md) for JSON results and interpretation.

## 5. Continue

See the [tutorial index](index.md), [API guide](tutorials/03-analysis-api.md), [template/provenance guide](tutorials/04-template-catalogue-and-provenance.md) and [release limitations](release-limits-v0.1.md). Building a *new* library is a distinct workflow in [Tutorial 07](tutorials/07-build-library-from-sources.md).
