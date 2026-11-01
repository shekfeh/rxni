# RXNI publication rights — required before release

## Open decisions

- The identified Lowe 1976–2016 Figshare USPTO source is marked **CC0**; preserve its citation. **CC BY 4.0** is selected for the derived RXNI dataset, conditional on ownership/rights for all other ingredients. A source CC0 marking alone does not prove that every intermediate transformation is cleared.

- **Selected:** original RXNI software uses Apache-2.0 (`LICENSE` and SPDX metadata). Determine correct copyright/attribution notices, contributor ownership and any additional required NOTICE content before publication.
- Verify original ownership and contributors of each Python source and bundled resource.
- R4H inspected 11 semantic resources: ten development-only provenance/audit files are excluded from the public package and archived in private R4I evidence. The retained 82-contract FGI registry is confirmed project-authored in `docs/release/RIGHTS_REVIEW.md`; preserve checksum and provenance; separately review any other third-party material.
- Keep private ChemKN, commercial/proprietary books, user-specific data, and non-redistributable datasets outside the public repository.
- Review optional dependency license notices and any external packaged database/GLF artifacts independently from the code.
- Confirm README scientific claims correspond exactly to validated tests.

Until this review is complete, packaging success is not permission to publish the code or any datasets.
