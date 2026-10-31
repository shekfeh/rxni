# RXNI Retrieval V2 — Physical storage schema profiles

## Canonical layout for all new public libraries

`rxni_retrieval_v2_compact10_v1` is the **existing** 10-column public RXNI
retrieval codec. The columns and their order are:

1. `record_id` — string
2. `reaction` — string
3. `semantic_profile_json` — string
4. `semantics_json` — string
5. `reaction_centers_json` — string
6. `localized_fingerprints_json` — string
7. `multichannel_fingerprint_json` — string
8. `reaction_key` — string
9. `representation_generation` — string
10. `representation_build_id` — string

Do not rename or reorder these fields under this profile. A different physical
storage layout requires an explicit new profile, not a silent writer change.

## Historical V6.7 compatibility

`rxni_retrieval_v2_inherited14_v1` describes the original private Runtime V6.7
writer, which inherited 12 V1 columns and appended generation/build identifiers.
Its additional columns are:

- `forward_edit_signature_hash` (uint64)
- `library_schema_version` (string)
- `reaction_key_version` (string)
- `system_i_representation_version` (string)

For the G2K CRD25 and USPTO pilots (128 total records), all four are null. The
common ten columns, decoded records and locator snapshots matched exactly.
Legacy layouts are accepted by the public schema validator only when their
four legacy-only fields are null; otherwise they require explicit scientific
migration review. Do not silently drop non-null historical values.

## Write/read policy

New public builds always use the compact 10-column schema. The public writer
adds the optional manifest field `physical_schema_profile`, leaving the
existing `schema_version`, V2 codec, ranking algorithms and locator format
unchanged. Old manifests without this field remain valid: readers infer the
known layout from the actual file schema. A mismatch between manifest-declared
and observed profile fails validation.

Validate a candidate library with:

```bash
python -m rxni.retrieval.storage.schema_profiles \
    /path/candidate_library/retrieval_representations_v2.parquet \
    --manifest /path/candidate_library/retrieval_v2_manifest.json
```

This schema validator does not establish semantic correctness or retrieval
ranking parity. G2L must establish functional compatibility with historical
14-column V2 data; private-vs-public V6.7 retrieval ranking requires a
separate paired runtime test.

## Provenance and immutable scientific contracts

- Semantic classifier: validated RXNI / private V6.7 lineage.
- Builder source: migrated Library Builder v0.2.3.
- Storage profile: independently declared (compact10 or inherited14).
- Preserve all previous manifests, hash evidence, source datasets, frozen V2
  decoders, and original artefacts. Do not rewrite historical files in place.
