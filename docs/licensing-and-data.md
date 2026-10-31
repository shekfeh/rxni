# License, knowledge sources and redistributed data

The original RXNI software is offered under **Apache License, Version 2.0**. Read the complete terms in the project-root `LICENSE` file. Its SPDX identifier is `Apache-2.0`.

This license declaration applies to material that the RXNI project has the rights to license. **It does not relicense third-party books, descriptions, reaction images, databases or scientific knowledge cards.** Installed Python dependencies retain their respective licenses.

## v0.1 packaging scope

- RXNI's scientific, retrieval and interface software is licensed Apache-2.0, subject to source-rights review.
- The optional ChemKN interface accepts explicitly supplied normalized JSONL cards; no textbook-specific adapters or knowledge cards are bundled. External providers retain responsibility for their own data rights.
- Any Wikipedia-derived text/figures may require separate CC BY-SA compliance, attribution and file-level checks; do not imply Apache-2.0 covers them.
- After R4I cleanup, the single distributed semantic resource is `semantic/contracts/fgi_contracts_v1.json` (82 FGI contracts), subject to authorship confirmation. Ten development-only resources were moved into private evidence. Successful packaging does not itself establish third-party rights.
- The tiny `examples/data/demo_reactions.smi` examples are original illustrative reactions, not copied corpus records.

## Release gating

Before public publication, classify every bundled data/resource file using its actual contents and provenance. If an item lacks confirmed redistribution rights, remove it from distributable artifacts or replace it with a cleared independent alternative **without silently breaking the frozen scientific interface**.

Keep source-code licensing, knowledge-content licensing, and corpus/database redistribution records separate. No private development evidence or source-corpus files belong in the public distribution.
