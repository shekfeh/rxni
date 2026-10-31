# RXNI retrieval-library inputs (v0.1)

RXNI **does not require** one fixed dataset source, but a read-ready retrieval
library requires a structural index, metadata lookup, stored retrieval records,
locator and retrieval manifest. v0.1 retains the tested V2 storage codec.

## Directory compatibility

Existing five-file retrieval-v2 directories open unchanged. The required files
are `structural_index_v1.sqlite`, `runtime_lookup_v1.sqlite`,
`retrieval_representations_v2.parquet`, `representation_locator_v2.sqlite`,
and `retrieval_v2_manifest.json`.

## Flexible file layouts

Place `rxni-library.json` in the library directory:

```json
{
  "schema_version": 1,
  "library_id": "example-library",
  "representation_format": "retrieval_representations_v2",
  "semantic_engine_version": "6.7",
  "profile_policy": "profile-canonical-v1",
  "artifacts": {
    "structural_index": "index/structure.sqlite",
    "metadata_lookup": "index/metadata.sqlite",
    "representations": "store/profiles.parquet",
    "locator": "index/locator.sqlite",
    "retrieval_manifest": "store/manifest.json"
  }
}
```

Artifact paths are relative to the library root and must resolve inside it.
The retrieval manifest must contain `representation_generation` equal to
`retrieval_representations_v2` and a nonempty `build_id`.

`semantic_engine_version` and `profile_policy` are **claims from the library
producer**; RXNI does not infer chemical equivalence or scientific correctness
from these strings. The existing code continues to use its validated V2
stored-record reader and frozen scoring equations.

## Support boundaries

* **Existing libraries already carrying V2 profiles**: readable immediately.
* **Libraries with the five V2 artifacts but different filenames**: readable
  through a descriptor.
* **Historical V1-only profiles**: discovered, but explicitly require conversion
  or rebuilding; they are never silently decoded as V2.
* **Dataset Core without auxiliary indices and profiles**: requires the separate
  library-building workflow.

Use `rxni.retrieval.libraries.inspect_library(path)` for format inspection,
then `rxni.integration.retriever.ReactionRetriever(path)` for actual retrieval.

Format acceptance is not a benchmark of ranking, template correctness or
abstention quality.
