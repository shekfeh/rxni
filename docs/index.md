# RXNI documentation

RXNI v0.1.0 — public software package; scientific method baseline: Reaction Intelligence v6.7.0.

## First steps

1. [Getting started](getting-started.md) — install RXNI and locate a validated library.
2. [CLI reference](cli-reference.md) — tested entry points, options and examples.
3. [Tutorial 01: query a reaction library](tutorials/01-retrieval.md) — reproducible retrieval on a prepared reaction-core Parquet.
4. [Tutorial 02: reaction analysis API](tutorials/02-analysis-api.md) — API request, response and failure behavior.
5. [Tutorial 03: source schema profiles](tutorials/03-source-schema.md) — CSV/Parquet schemas and V2 builders.
6. [Tutorial 04: foundation schema profiles](tutorials/04-foundation-schema.md) — remap validated SQLite schema names.

## Technical references already present

- [Architecture](architecture.md) — scientific layers and retrieval integration.
- [Version and compatibility contracts](version-contracts.md).
- [Source-schema profiles](input-schema-profiles.md).
- [Foundation-schema profiles](foundation-schema-profiles.md).
- [Release limitations](release-limits-v0.1.md).
- [Distribution-rights checklist](distribution-rights-checklist.md).
- [Licensing and external data](licensing-and-data.md).

## Examples

- `examples/data/demo_reactions.smi` — tiny, original illustrative reaction strings; **not a retrieval index**.
- `examples/scripts/retrieve_from_core.sh` — choose one valid fixture reaction and query a prebuilt library.
- `examples/profiles/*.json` — four schema profiles. A schema profile does not construct a complete retrieval library by itself.

The examples require an existing validated retrieval library. They **do not** download or redistribute CRD2025, USPTO, ChemKN cards, Organic Syntheses, Carey or March datasets.
