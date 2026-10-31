# RXNI source-schema profiles (experimental v0.1)

**Separate the input schema from the retrieval storage ABI.** A JSON source-schema
profile tells the V2 *representation builder* which fields of a Parquet or CSV
source correspond to opaque `reaction_key` and forward reaction SMILES. It does
not redefine the `retrieval_representations_v2` Parquet columns, locator, ranking,
or scientific event rules.

Example (`examples/profiles/crd25-uspto-v023a.json`):

```json
{
  "schema_version": 1,
  "profile_id": "native-v023a-reaction-core",
  "source_format": "parquet",
  "columns": {
    "reaction_key": "reaction_key",
    "reaction_smiles": "reaction"
  }
}
```

The historical v6.7 source uses `forward_reaction_smiles` instead and can use
the corresponding example profile. For CSV sources, set `source_format` to
`csv`; UTF-8 CSV input is streamed without loading the entire table.

Inspect **without writing output**:

```bash
python -m rxni.retrieval.builders.source_schema inspect \
  --source-core /path/to/reaction_core_v5.parquet \
  --input-schema examples/profiles/crd25-uspto-v023a.json
```

Build V2 representations (with the upstream library/foundation supplied
separately):

```bash
python -m rxni.retrieval.builders.representations \
  --library /path/to/library-foundation \
  --input-core /path/to/reaction_core_v5.parquet \
  --input-schema examples/profiles/crd25-uspto-v023a.json \
  --output-dir /path/to/new-private-build \
  --workers 2 --chunk-size 128 --limit 100
```

## Contract boundaries

- The source MUST contain a stable reaction key, and each selected reaction
  MUST be `reactants>agents>products` SMILES. No reaction keys are invented.
- Optional `valid` and `feature_error` columns support preexisting filtering.
  Invalid Boolean encodings are errors, not silently accepted.
- Input format is either **Parquet** or **CSV**. A profile is *not* a chemistry
  parser or arbitrary structure-to-reaction converter.
- The build writes `source_binding_v1.json` privately in its output directory.
  Resuming with a different source SHA256 or profile SHA256 is rejected.
- The resulting V2 manifest may include `source_schema_adapter` provenance,
  but its existing physical V2 schema, locator generation, and identifiers are
  unchanged. Existing `--library`-only builds follow the historical path.
- **Retrieval readiness is separate.** Structural-index keys, metadata lookup
  keys, template memberships, and V2 representation keys must describe the
  *same population and key namespace*. A profile cannot repair mismatched
  indexes, missing auxiliary files, or incompatible retrieval serialization.
- Do not mix a V2 representation built from a new Core into an unrelated old
  five-file retrieval library merely because the schema columns align.

For flexible **read-side** file locations, see `docs/library-formats.md` and
`rxni-library.json`. A future full-library adapter may additionally support
verified SQLite index/metadata column mappings and explicit index rebuilding.
