# Custom Retrieval V2 physical schemas

The public **writer** always emits the frozen canonical 10-column Parquet schema
`rxni_retrieval_v2_compact10_v1`. This extension adds **read-only, explicitly
configured** compatibility for user-supplied V2 Parquet layouts.

Do not confuse this with *source reaction dataset profiles*, SQLite index
compatibility, or changes to RXNI scientific semantics. Those are separate
versioned contracts.

## Explicit custom profile

A schema profile JSON must include:

- `schema_profile_version`: `1`
- `profile_id`: unique string beginning `rxni_custom_`
- `columns`: exact 10 canonical RXNI V2 field names mapped to their physical
  column names (one-to-one; can be renamed/reordered, but values must be Arrow
  string or large_string).
- `extra_columns`: explicit physical extra-column declarations, each with an
  Arrow type and a role of `provenance` or `annotation`. Extra columns are **not**
  interpreted as additional retrieval signals.

See `examples/custom-v2-schema.json` for a full example.

The V2 manifest for a custom layout must declare:

```json
{
  "physical_schema_profile": "rxni_custom_example_v1",
  "custom_schema_profile_sha256": "<actual SHA256 of the profile JSON>"
}
```

The existing manifest requirements (build ID, record generation, source and
representation checksums, record counts) remain applicable. Replace the
example hash with a real hash; never edit or relabel an old immutable manifest
in place.

### Validate a supplied Parquet and manifest

```bash
python -m rxni.retrieval.storage.custom_profiles \
  /path/to/library/retrieval_representations_v2.parquet \
  --profile /path/to/custom-v2-schema.json \
  --manifest /path/to/library/retrieval_v2_manifest.json
```

### Use a custom profile with the stored-record provider

```python
from rxni.retrieval.storage.provider import V2StoredRecordProvider

provider = V2StoredRecordProvider(
    "/path/library/retrieval_representations_v2.parquet",
    "/path/library/representation_locator_v2.sqlite",
    "/path/library/retrieval_v2_manifest.json",
    schema_profile="/path/to/custom-v2-schema.json",
)
record = provider.get(record_id="known_reaction_key")
```

Alternatively, `RXNI_V2_SCHEMA_PROFILE=/path/to/custom-v2-schema.json` opts
in at the process level **for any code that constructs this provider**, without
changing the provider factory. It is not a guarantee that every other RXNI
entrypoint uses this provider. Use explicit provider configuration when possible.

Schema compatibility does **not** establish chemical correctness, semantic
compatibility, ranking parity, or full-library compatibility. Custom profiles
cannot alter representation generation, supply missing scientific fields,
change the canonical codec, or modify structural/SQLite indexes.

When an additional field is scientifically meaningful, build a documented
schema/semantic version and validate it separately; do not declare it a
provenance column merely to bypass checks.
