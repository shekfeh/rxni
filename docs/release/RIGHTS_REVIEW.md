# Continuing rights and future-release review

**RXNI v0.1.0 has already been publicly released.** Further release/dataset approval requires a separate rights and scientific review; this file does not retract the published tag.

- [x] Review derivation, authorship, attribution and redistribution rights for `src/rxni/semantic/contracts/fgi_contracts_v1.json` (retained 82-entry contract; accepted reference checksum is recorded privately).
- [ ] Review third-party code, vendored material, embedded reaction examples, dataset terms, license compatibility and any required notices.
- [ ] Check Git staged files and commit history for keys, private user data, proprietary reaction datasets, model weights and local file paths.
- [ ] Keep patent-derived data and private held-out truth outside public Git commits, release assets and package distributions unless separately cleared.
- [ ] Review public README, installation instructions, CLI/API docs and package description for accurate supported features (do not advertise a ChemKN installation extra that is not declared).
- [ ] Record clean-wheel and sdist validation, an installed-wheel API parity test and SHA256 of final distributions.
- [ ] Select GitHub owner, visibility and release license; obtain explicit final approval before any GitHub push, tag publication or GitHub release.

## Dataset provenance distinction (1 November 2026)

- Daniel Lowe's identified original Figshare USPTO dataset is marked CC0.
- The derived RXNI USPTO retrieval archive has CC BY 4.0 selected for Zenodo; approval requires confirmation of the complete extraction/derivative source chain.
- The archived dataset has passed a clean RXNI v0.1.0 installation and a real retrieval smoke, but the Ring031 review case remains disclosed.

## FGI provenance decision (31 October 2026)

- **Project authorship confirmed:** the 82 directional FGI rules and generated
  detection implementation originated in the Reaction Intelligence project.
- **Scientific development evidence:** heuristic validation and refinement used
  extracted USPTO 1960–2016 reactions/templates; these source datasets are
  independent of the distributable software and are not licensed by this
  provenance declaration.
- **Frozen contract SHA256:**
  `c5cfca057c283068b40a3dfb1e881c294c73e6196353deff2d5f5242c60667c2`.
- **Scope:** the project-author provenance item is resolved. Third-party code,
  dataset redistribution, sensitive information, and final GitHub publication
  must still pass the separate tracked-content review.
