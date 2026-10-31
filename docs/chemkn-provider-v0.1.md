# Optional ChemKN provider interface (RXNI v0.1)

RXNI core reaction semantics, V2 retrieval and condition evidence operate without
any ChemKN cards or an external knowledge repository. RXNI does not ship or
search publisher-specific knowledge cards by default.

## Explicit normalized JSONL input

Set `RXNI_KNOWLEDGE_CARDS=/absolute/path/to/cards.jsonl` for the HTTP/service
interface, or pass `--knowledge-cards /absolute/path/to/cards.jsonl` to the
legacy analysis/summary CLIs. Optionally provide `RXNI_METHOD_LINKS` or
`--method-links` for *reviewed* method-to-card references. Paths are never
discovered automatically.

Each JSONL record must contain `source`, `source_id`, `canonical_name`, and
`source_kind` (nonempty strings). Optional fields follow `KnowledgeCard` in
`rxni.chemkn.knowledge_schema`. The `compatibility` object, when supplied,
uses `CompatibilitySignature` fields; unknown keys are rejected. Duplicate
`source:source_id` identifiers are rejected. `source` denotes a generic
provider identifier, not a list of hardcoded publishers.

Example (illustrative, original test data):

```json
{"source":"example_open","source_id":"demo-001","canonical_name":"Illustrative acyl substitution","source_kind":"curated_rule","reaction_family":"nucleophilic_acyl_substitution","retrieval_tags":["acyl","substitution"],"compatibility":{"required_events":[],"required_event_mode":"none"},"evidence":{"license":"CC0-1.0","status":"synthetic_example"}}
```

## Important limits

- Cards are **evidence**, not scientific or experimental outcome guarantees.
- Empty or missing configuration reports `unavailable`, not fabricated matches.
- Compatibility rules are evaluated for *all* sources, without special treatment.
- External provider integration is intentionally minimal in v0.1; loader plugin
  architecture, caching, indexing and rights enforcement require a future design.
- Knowledge source licenses and attribution are the provider's responsibility.
- The examples are fabricated for demonstrations, not extracted from textbooks.

Published RXNI does not contain publisher-specific adapter implementations or
textbook card databases. Legacy private adapters remain in private R4H backups.
