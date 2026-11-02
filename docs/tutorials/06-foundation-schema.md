# Tutorial 06 — Foundation-schema profiles

## Fresh-machine standard and renamed-schema tests

The **published** Zenodo library resolves with `foundation_profile=None`: this is valid because its standard SQLite schema does not need a custom mapping. Do **not** call `load_foundation_profile(None)`.

```python
from rxni.retrieval.libraries import resolve_library
lib = resolve_library('/absolute/path/to/uspto_1976_2016_rxni_library')
if lib.foundation_profile is None:
    print('Standard schema: custom foundation profile not required')
else:
    from rxni.retrieval.libraries.foundation_schema import (
        load_foundation_profile, validate_foundation_schema, audit_key_membership)
    profile = load_foundation_profile(lib.foundation_profile)
    print(validate_foundation_schema(lib.structural_index, lib.metadata_lookup, profile))
    print(audit_key_membership(lib.structural_index, lib.metadata_lookup,
                               lib.locator, profile, sample_size=32))
```

In the independent WSL test, a separately constructed identity profile validated the **real** published structural/metadata SQLite files and linked **32/32** sampled structural keys through metadata and V2 locator. An invalid physical column was rejected. Separate **small** SQLite fixtures with renamed tables/columns also validated with **32/32** key joins against the original V2 locator. These are **bounded** checks, not complete key-set/position audits.

A read-only QA fixture for repeating these checks is supplied separately in the documentation patch package (`qa/rxni_tutorial_audit_fixture_v011.py`). Use only the published dataset as input and a new private work directory.


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
if lib.foundation_profile is None:
    print('Standard library schema: no custom profile configured')
else:
    profile = load_foundation_profile(lib.foundation_profile)
    print(validate_foundation_schema(lib.structural_index, lib.metadata_lookup, profile))
    print(audit_key_membership(
        lib.structural_index, lib.metadata_lookup, lib.locator,
        profile, sample_size=32,
    ))
```

This is a **bounded sample** and does not certify full-library key integrity. Every field's semantics must still be established by the source producer. Read [foundation schema contracts](../foundation-schema-profiles.md) before substituting profiles.
