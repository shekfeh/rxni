# RXNI public names and frozen compatibility identifiers

**Public package:** `rxni` (release series `0.1.x`; next proposed `0.1.1`).
The software version is independent of the scientific rule baseline, template identities, V2 stored ABI, or the physical filename layout of a published dataset.

| Stable public interface | Preserved compatibility name / source claim | Why it remains |
|---|---|---|
| `rxni.integration.ReactionRetriever` | `V66Retriever` implementation | Validated scoring API; an alias already exists. |
| `rxni.retrieval.StoredRecordProvider` | `V2StoredRecordProvider` | Exact persisted stored-provider contract. |
| `reaction_representations.parquet` | `retrieval_representations_v2` | Physical storage generation; NOT user-facing dataset title. |
| `representation_locator.sqlite` | SQLite table `representation_locator_v2` | Unchanged index/locator ABI. |
| `rxni-library.json` | descriptor version 1 | Relocatable file mapping without data conversion. |
| Project-developed semantic rules | historical runtime-v6.x semantic source identifiers | Recorded scientific producer provenance; do not falsify or strip. |
| `retro_template_hash` | historical AiZynthTrain `hash_from_bits()` | Validated 42,128 historical fingerprint-based identities. |

Do not rename internal SQLite tables, stored Parquet fields, ABI generation strings, original hash identifiers, build IDs or provenance SHA256 identifiers for cosmetic reasons. Such changes require separate scientific and storage migrations. `V66Retriever` and legacy provider names remain importable for backwards compatibility; use stable aliases in new public examples.

**Historical dataset:** `uspto_1976_2016_rxni_library/` in the downloadable archive. The 1976–2016 label denotes the parent USPTO source range, not independently verified dates for all 812,916 indexed records. The historical candidate's exact semantic source claims and unresolved review case must remain in the archive metadata, even if public usage tutorials emphasize RXNI v0.1.x rather than internal runtime labels.
