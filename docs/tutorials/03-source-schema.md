# Tutorial 03 — Input source-schema profiles

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
