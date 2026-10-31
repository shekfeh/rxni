# RXNI architecture

RXNI provides an independently packaged scientific reaction-intelligence runtime:

- `rxni.semantic`: mapped-graph reaction semantics, event resolution, validated directional functional-group conversions, and canonical semantic profiles.
- `rxni.retrieval`: structured candidate generation and ranking, native stored-record codec, structural lookup, and library-format adapters.
- `rxni.integration`: query orchestration and shared retriever interface.
- `rxni.method`: method and condition-oriented views.
- `rxni.public`: typed service and response models.
- `rxni.cli` / `rxni.cli_analysis` / `rxni.cli_summary`: command-line access.
- `rxni.api`: optional HTTP application.

The public package release version is independent from both the validated semantic baseline and the persisted representation format. See [Version contracts](version-contracts.md) and [Library formats](library-formats.md).

RXNI expects explicitly configured external libraries; the Python package does not bundle the reaction datasets, prebuilt private profile indices, or historical development reports. The built library's provenance and license should be evaluated separately from the source-code license.

Retrieval rankings are evidence-oriented candidate recommendations, not proof of reaction mechanism or successful reaction conditions. Any optional knowledge enrichment should identify missing sources rather than imply that they were consulted.
