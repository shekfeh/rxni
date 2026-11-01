# Repository review: prepared changes (1 November 2026)

Based on the 1 Nov independent review of `shekfeh/rxni` main (`9bacd435`).

| Finding | Action |
|---|---|
| Immutable v0.1.0 lacks new docs/scripts | Prepare untagged `0.1.1` maintenance candidate. No retagging. |
| Two missing documented examples | Add `examples/data/demo_reactions.smi` and `examples/scripts/retrieve_from_core.sh`. |
| Stale MANIFEST.in | Replace with explicit full public-doc, tutorial, script, test and example allowlist. |
| No GitHub CI | Add 3.10–3.12 unit tests, installed CLI smoke, wheel/sdist build and archive-content check. |
| Inconsistent public release status | Identify `v0.1.0` as published and `0.1.1` as pending. |
| FGI authorship wording | Use established project-authorship record; separate from third-party rights. |
| Lowe USPTO dataset licensing | Record that the identified upstream Figshare source shows CC0, without assuming all derivative ingredients are cleared. |
| Data status/namespace mismatch | Use `uspto_1976_2016_rxni_library`, include archive digest/bytes, leave DOI placeholder until obtained. |
| Legacy names in public examples | Use existing public `ReactionRetriever` alias; add `StoredRecordProvider` compatibility alias and explain frozen scientific identifiers. |
| Download verifier only as a script | Add installed `rxni-verify-library` with the same read-only purpose. |
| WSL validation absent | Remains tomorrow's separate independent test. No acceptance claim today. |

## Deliberate nonchanges / outstanding blockers

- **No scientific runtime refactoring:** older `V66Retriever`/`System-I` internal classes, scientific build identifiers and source pins are frozen compatibility/provenance. Removing them would require a separate validated migration.
- **No fabricated DOI:** publishing a draft to Zenodo does not supply an immutable published DOI link until user provides it.
- **No blanket rights signoff:** Lowe source CC0 is supported by the user's source check; the derivative supply chain still needs its own explicit provenance clearance.
- **Scientific accuracy:** the Ring031 scientific review queue and broader retrieval-quality validation remain separate from technical parity.
- **No data in Git:** the ZIP is separately hosted; no large SQLite/Parquet files in the repository.
- **Release tag:** only after CI, user review, Zenodo metadata, and tomorrow's WSL result; no push/tag performed by this kit.
