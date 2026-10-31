# RXNI foundation-schema profiles (R3S2)

The **source schema** profile (R3S1) maps CSV/Parquet reaction inputs to the
V2 representation builder. A **foundation schema** profile maps existing
SQLite *structural index* and *metadata lookup* tables/columns to the frozen
RXNI retrieval semantics. These are separate interfaces; neither converts an
arbitrary source Core into a complete retrieval-ready library.

A library descriptor (`rxni-library.json`, schema_version 1) may add:

```json
"foundation_schema_profile": "config/foundation-schema.json"
```

The referenced profile uses the format in `examples/profiles/foundation-schema-example.json`.
It defines explicit `structural_index.table`, `metadata_lookup.table` and
logical→physical `columns` mappings. Six **structural hash fields** and the
`reaction_key` must be mapped. The metadata profile must map `reaction_key`
and `forward_reaction_smiles`; it may additionally map `reaction_id`,
`retro_template_hash`, `retro_template_number` and
`membership_occurrence_count`. The latter fields affect available output:
map all fields carried by the source to preserve result parity.

**Only physical SQLite names vary.** The mathematical definitions of the
structural hashes, source reaction identity, representation generation and
retrieval/scoring logic do not change. Identically named columns that contain
different hash algorithms are **not** made compatible by this profile. The
V2 locator, Parquet records and manifest remain fixed and mandatory. Source
keys must agree **exactly** across index, metadata and stored records; there
is no silent rekeying, fallback matching, hash repair or automatic synthesis.

Profiles must be local to the declared library directory; path escapes are
rejected, and SQLite identifier values must be simple identifiers. Missing
physical columns and invalid schemas fail closed. Without a profile, the
existing legacy five-file read behavior remains unchanged.

For additional read-only preflight:

```python
from rxni.retrieval.libraries import resolve_library
from rxni.retrieval.libraries.foundation_schema import (
    load_foundation_profile, validate_foundation_schema, audit_key_membership,
)
lib = resolve_library('path/to/library')
p = load_foundation_profile(lib.foundation_profile)
print(validate_foundation_schema(lib.structural_index, lib.metadata_lookup, p))
print(audit_key_membership(lib.structural_index, lib.metadata_lookup, lib.locator, p, sample_size=32))
```

`audit_key_membership` is a bounded sample, **not** an exhaustive full-library
key-consistency proof. The V2 provider additionally enforces a locator/store
key match at record fetch time. For a complete schema migration, test result
parity against the same V2 library with both original and renamed SQLite
schemas, and require zero differences in ranked members and score vectors.

An incorrect schema *mapping* could silently misrepresent chemistry even when
all columns exist. RXNI verifies physical schema and key connectivity, not
whether the producer's semantic hashes were derived using the right rules.
