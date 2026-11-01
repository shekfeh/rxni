# RXNI v0.1.1 — candidate release notes (NOT tagged or published)

Proposed documentation/packaging maintenance release after immutable RXNI v0.1.0.

- Add both formerly missing documented examples.
- Package all public library tutorials, release scripts and documentation in source distributions.
- Add a portable `rxni-verify-library` console command; no external corpus is bundled.
- Add GitHub Actions Python 3.10, 3.11 and 3.12 unit/packaging/entrypoint checks.
- Distinguish RXNI software Apache-2.0, project-authored FGI rules, parent Lowe USPTO CC0 source claim, and proposed derivative-dataset CC BY 4.0.
- Standardize archive extraction name and document historical frozen compatibility identifiers without altering science or storage.

**No scientific/runtime ABI change**: original semantic/scoring algorithms, source reaction identities, reaction fingerprints, template hashes, SQLite schema and stored V2 layout are unchanged.

**Hold before tagging**: check CI; inspect source-distribution and wheel contents; link the actual published Zenodo dataset DOI (not just its uploaded draft); evaluate remaining scientific/rights limitations; perform the independent WSL fresh-machine evaluation. Do not mutate or move `v0.1.0`.
