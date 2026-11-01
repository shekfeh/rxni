# Building RXNI Core, native and retrieval libraries

The current public `rxni.builder` module contains **both** the immutable Core writer and an explicit v6.7-aware native/retrieval auxiliary pipeline. It also retains separate scientific compatibility gates. For the full reproducible workflow, see [Tutorial 06](tutorials/06-build-library-from-sources.md).

## Quick navigation

- [Install/use a downloaded library](tutorials/05-download-and-use-library.md)
- [End-to-end source preparation, AiZynthTrain, Core, auxiliary and packaging](tutorials/06-build-library-from-sources.md)
- [Template SMARTS, identity namespaces and scientific provenance](tutorials/07-template-catalogue-and-provenance.md)
- [Historical USPTO distribution status](library-distribution.md)
- [Source CSV/Parquet schema adapters](input-schema-profiles.md)
- [Foundation SQLite schema adapters](foundation-schema-profiles.md)

## Actual public builder commands

```bash
rxni-build-library profiles
rxni-build-library inspect --export-dir /path/to/aligned/export --reaction-column reaction
rxni-build-library build-core --export-dir /path/to/aligned/export \
  --reaction-column reaction --dataset MY_DATASET --out-dir /new/output/core
rxni-build-library validate-core --core /new/output/core --rdkit --report /new/output/audit.json
rxni-build-library inspect-runtime --root /path/to/rxni --runtime v6.7
rxni-build-library build-auxiliary --core /new/output/core --work-dir /new/output/auxiliary \
  --root /path/to/rxni --runtime v6.7 --workers 4 --through v2
rxni-build-library assemble-candidate --core /new/output/core --work-dir /new/output/auxiliary \
  --root /path/to/rxni --runtime v6.7 \
  --lookup-db /new/output/auxiliary/native/runtime_lookup_v1.sqlite \
  --out-dir /new/output/candidate
```

**Status distinction:** `build-core` creates a 4-Parquet Core, *not* a retrieval library. `build-auxiliary` now resolves **public installed** producer modules `rxni.builder._native.build_all`, `rxni.retrieval.builders.representations`, and `rxni.builder.lookup_cache`; the older claim that the code requires private upstream runtime imports is stale. Nevertheless, the resulting V2 files are **candidate outputs** until full index, semantic and ranking parity and source rights are certified.

## Critical template-identity distinction

The current strict public Core writer verifies `retro_template_hash == SHA256(exact retro SMARTS)`. The historical AiZynthTrain `TemplateHash` is calculated using `ReactionTemplate.hash_from_bits()`. Historical 42,128-template replay was 42,128/42,128 successful. These identifiers are not interchangeable. A new Core source export must preserve the native hash separately and re-key all memberships intentionally; the historical Ring031 library must retain its original native identity namespace. See Tutorial 06 §2–3 before attempting new datasets.

## Release boundary

Do not commit 2+ GB historical library data to ordinary Git. Keep data in a separate downloadable archive with DOI and attach only scripts/docs and final immutable download links to this repository. `rxni-library.json` can map clean public filenames without altering the V2 storage ABI. **No dataset is cleared for public redistribution merely by passing RXNI format tests.**
