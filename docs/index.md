# RXNI documentation

RXNI software v0.1.0; scientific method baseline v6.7.0. A compatible local reaction library is required for retrieval; source software installation alone does not download one.

## First steps

1. [Getting started](getting-started.md) — install RXNI, inspect CLI and configure a library.
2. [CLI reference](cli-reference.md) — commands, modes, options and human-facing reports.
3. [Tutorial 01: query a reaction library](tutorials/01-retrieval.md) — basic retrieval.
4. [Tutorial 02: reaction analysis API](tutorials/02-analysis-api.md) — API requests and responses.
5. [Tutorial 03: source schema profiles](tutorials/03-source-schema.md) — CSV/Parquet source mapping.
6. [Tutorial 04: foundation schema profiles](tutorials/04-foundation-schema.md) — verified SQLite layouts.
7. [Tutorial 05: download and use a library](tutorials/05-download-and-use-library.md) — standalone historical USPTO dataset (publication URL pending).
8. [Tutorial 06: complete library construction](tutorials/06-build-library-from-sources.md) — rxnutils, AiZynthTrain, aligned exports, Core, native indexes, retrieval and certification.
9. [Tutorial 07: template SMARTS and provenance](tutorials/07-template-catalogue-and-provenance.md) — template identities and scientific controls.

## Technical and publication references

- [Historical compatibility identifiers](compatibility-identifiers.md) — stable public aliases versus immutable V2 and historical scientific fields.
- [Unreleased 0.1.1 maintenance candidate](release/RELEASE_NOTES_v0.1.1_CANDIDATE.md).

- [Architecture](architecture.md)
- [Version contracts](version-contracts.md)
- [Library builder overview](library-building.md)
- [Historical USPTO release status](library-distribution.md)
- [Source-schema profiles](input-schema-profiles.md)
- [Foundation-schema profiles](foundation-schema-profiles.md)
- [Release limitations](release-limits-v0.1.md)
- [Distribution rights checklist](distribution-rights-checklist.md)
- [Licensing and external data](licensing-and-data.md)

## Examples and installed scripts

- `examples/data/demo_reactions.smi` — synthetic examples, **not** a retrieval index.
- `examples/scripts/retrieve_from_core.sh` — query example against a prebuilt, compatible library.
- `examples/profiles/*.json` — explicit input/storage mappings, *not chemistry conversions*.
- `scripts/release/finalize_uspto_historical.py` — **private packaging** of the selected historical Ring031 scientific candidate, not publication authorization.
- `scripts/release/verify_published_library.py` — extraction/checksum and read-side smoke validation.
- `scripts/release/replay_aizynth_template_hashes.py` — fingerprint hash replay in an isolated rxnutils environment.

No third-party corpus should be added to Git history or declared publicly redistributable without a separate rights review.
