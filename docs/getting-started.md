# Getting started with RXNI v0.1.0

RXNI is a Python reaction-semantics and retrieval toolkit. The software release is **0.1.0**; its frozen scientific-method baseline is **Reaction Intelligence v6.7.0**. It does **not** promise chemical reaction success, independent retrieval accuracy, or an electron-pushing mechanism prediction.

## 1. Install

A fresh virtual environment is recommended. The declared Python support floor is 3.10.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
```

For API and optional rich interface:

```bash
python -m pip install '.[api,cli]'
```

For YAML-backed optional adapter functionality:

```bash
python -m pip install '.[chemkn]'
```

Do not expect RXNI to download a reaction corpus or proprietary ChemKN knowledge assets during installation.

## 2. Check the CLI

```bash
rxni --help
rxni retrieve --help
rxni analyze --help
rxni-summary --help
rxni-batch --help
```

`rxni`, `rxni-summary`, and `rxni-batch` are the tested installed entry points. The actual arguments of `analyze` and batch commands should be read from your installed version's `--help`; this document does not invent undocumented options.

## 3. Obtain a compatible library

RXNI retrieval depends on a **locally prepared, validated library** comprising compatible structural/metadata indexes, retrieval V2 representations and locators. For command examples, set:

```bash
export RXNI_LIBRARY='/absolute/path/to/validated/library'
```

`RXNI_LIBRARY` is only a shell variable used by these tutorials; the internal API's environment variable is `REACTION_INTELLIGENCE_LIBRARY` (see [API tutorial](tutorials/02-analysis-api.md)).

A Parquet file containing reactions is **not** by itself a complete retrieval-ready library. Different source datasets cannot be mixed with incompatible index/key namespaces.

## 4. Run a retrieval query

```bash
rxni retrieve \
  --library "$RXNI_LIBRARY" \
  --reaction 'CC(=O)O.OCC>>CC(=O)OCC' \
  --format json \
  --top-k 5 \
  --candidate-budget 500
```

This is an **illustrative input**, not a gold-standard result. Its rank positions, matches and conditions depend on your particular library. For a deterministic choice of a real reaction from an installed fixture, see [Tutorial 01](tutorials/01-retrieval.md).

## 5. Continue

Use the [CLI reference](cli-reference.md) for command syntax and [the API tutorial](tutorials/02-analysis-api.md) for JSON analysis. Read [release limitations](release-limits-v0.1.md) before interpreting results.
