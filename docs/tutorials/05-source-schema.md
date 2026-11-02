# Tutorial 05 — Input source-schema profiles

## Fresh-machine source-adapter acceptance test

The public source-schema reader was tested with an independently prepared **16-reaction** mapped Parquet fixture (`reaction_key`, `reaction`, `native_cleaned_reaction`, and source evidence fields). The CLI returned `SOURCE_SCHEMA_COMPATIBLE`; `iter_source` read all 16 records. To repeat on **your own prepared** aligned source Core (not an unrelated published V2 file):

```bash
export SOURCE_CORE='/absolute/path/to/reaction_core_v5.parquet'
export SOURCE_PROFILE='examples/profiles/crd25-uspto-v023a.json'
python -m rxni.retrieval.builders.source_schema inspect \
  --source-core "$SOURCE_CORE" --input-schema "$SOURCE_PROFILE"
python - "$SOURCE_CORE" "$SOURCE_PROFILE" <<'PY'
import sys
from rxni.retrieval.builders.source_schema import load_profile, iter_source
rows = list(iter_source(sys.argv[1], load_profile(sys.argv[2])))
print('Accepted mapped-source rows:', len(rows))
assert rows and all(rxn.count('>') == 2 for _, _, rxn in rows)
PY
```

A correct source profile **only** accepts input columns. It does not independently create compatible structural/metadata indexes or certify a scientific library. It does not re-key historical template identities.


Source-schema profiles control how the representation builder reads the source CSV or Parquet reaction table. They **do not** redefine the frozen retrieval V2 format or reconcile incompatible key namespaces.

## Example v0.23a source profile

`examples/profiles/crd25-uspto-v023a.json` maps `reaction_key` and `reaction` in the v0.23a exports.

```bash
python -m rxni.retrieval.builders.source_schema inspect \
  --source-core /path/to/reaction_core_v5.parquet \
  --input-schema examples/profiles/crd25-uspto-v023a.json
```

The inspection command does not write a new V2 library.

## Example representation build

Only run this on a validated source/foundation pairing, writing to a new private directory:

```bash
python -m rxni.retrieval.builders.representations \
  --library /path/to/compatible/library-foundation \
  --input-core /path/to/reaction_core_v5.parquet \
  --input-schema examples/profiles/crd25-uspto-v023a.json \
  --output-dir /path/to/new/private/build \
  --workers 2 --chunk-size 128 --limit 100
```

A built representation alone does **not** make a library retrieval-ready. Check structural-index keys, metadata lookup, representation keys, locator and manifest consistency.

For historical source fields use `examples/profiles/historical-v67.json`; for a custom UTF-8 CSV see `examples/profiles/custom-csv.json`. See the full [source schema guide](../input-schema-profiles.md).
