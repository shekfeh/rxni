# RXNI documentation — learning path

**Software:** published `v0.1.0` and audited maintenance `0.1.1` source candidate. **Scientific baseline:** v6.7.0. A validated local library is needed for indexed retrieval. The historical dataset is publicly accessible on [Zenodo record 23856755](https://zenodo.org/records/23856755), but scientific validation and rights attestation remain separate from availability.

## First-time users

1. [Getting started](getting-started.md) — installation and environment sanity checks.
2. [Tutorial 01 — Download and verify a library](tutorials/01-download-verify-library.md) — published archive, complete hashes and mapped smoke query.
3. [Tutorial 02 — Reproducible retrieval](tutorials/02-retrieval.md) — ranked candidates, top templates and JSON result.
4. [Tutorial 03 — HTTP analysis API](tutorials/03-analysis-api.md) — health, schema, structured analysis and raw-output restriction.
5. [Tutorial 04 — Template SMARTS and provenance](tutorials/04-template-catalogue-and-provenance.md) — template identity and scientific limits.

## Library developers

6. [Tutorial 05 — Source-schema profiles](tutorials/05-source-schema.md) — CSV/Parquet source fields, actual inspection and read tests.
7. [Tutorial 06 — Foundation-schema profiles](tutorials/06-foundation-schema.md) — structural/metadata SQLite mapping, bounded key audit.
8. [Tutorial 07 — Build a new library](tutorials/07-build-library-from-sources.md) — rxnutils/AiZynthTrain background, strict input evidence, Core, native auxiliary and V2 candidate (not publication).

## Technical references

- [CLI reference](cli-reference.md) · [API quickstart](api-quickstart.md) · [Library formats](library-formats.md)
- [Architecture](architecture.md) · [Version contracts](version-contracts.md) · [Historical compatibility IDs](compatibility-identifiers.md)
- [Library-builder overview](library-building.md) · [Input-schema contracts](input-schema-profiles.md) · [Foundation-schema contracts](foundation-schema-profiles.md)
- [Library distribution](library-distribution.md) · [Release limitations](release-limits-v0.1.md) · [Data rights](licensing-and-data.md) · [Distribution-rights checklist](distribution-rights-checklist.md)
- [0.1.1 maintenance-candidate notes](release/RELEASE_NOTES_v0.1.1_CANDIDATE.md)

## Mapping requirement and publication gate

**All reaction query examples in this learning path use mapped reactions.** v0.1.1 was observed to accept unmapped reactions that produced `NO_MAPPED_BOND_EDIT` and scientifically weak retrieval. The planned strict-rejection code patch must be tested before these docs are deployed as the current runtime contract. A mapper integration, if added, is optional and separate from the native scientific runtime; no neural mapping dependency is required for basic RXNI operation.

Examples and installed scripts remain in `examples/`, notably `examples/scripts/retrieve_from_core.sh` and `examples/profiles/*.json`. The public Zenodo dataset is **not** a source Core; do not substitute an unrelated Core or mix index/key namespaces. Do not commit the large historical library into ordinary Git history.
