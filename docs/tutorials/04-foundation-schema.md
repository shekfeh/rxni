# Tutorial 04 — Foundation-schema profiles

The foundation schema adapter maps physical SQLite table and column **names** to RXNI's existing logical fields. It does not alter hash algorithms, reaction-key namespaces or V2 stored records.

## Configure a library descriptor

Inside a compatible `rxni-library.json`, add a path to the optional foundation profile:

```json
{
  "foundation_schema_profile": "config/foundation-schema.json"
}
```

This is a **fragment** of an existing library descriptor, not a standalone complete descriptor. See `examples/profiles/foundation-schema-example.json` for a full mapping example.

## Preflight check

```python
from rxni.retrieval.libraries import resolve_library
from rxni.retrieval.libraries.foundation_schema import (
    load_foundation_profile, validate_foundation_schema, audit_key_membership,
)

lib = resolve_library('/absolute/path/to/validated/library')
profile = load_foundation_profile(lib.foundation_profile)
print(validate_foundation_schema(lib.structural_index, lib.metadata_lookup, profile))
print(audit_key_membership(
    lib.structural_index, lib.metadata_lookup, lib.locator,
    profile, sample_size=32,
))
```

This is a **bounded sample** and does not certify full-library key integrity. Every field's semantics must still be established by the source producer. Read [foundation schema contracts](../foundation-schema-profiles.md) before substituting profiles.
