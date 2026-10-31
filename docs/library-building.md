# Build an RXNI Core library

The versionless Python namespace `rxni.builder` contains the verified v0.2.3 source-suite Library Builder Core implementation. The public console command is `rxni-build-library`.

The Core builder requires a **complete, lineage-attested source export** directory containing `reaction_core_v5.parquet`, `reaction_template_membership_v5.parquet`, `template_library_v5.parquet` and `manifest.json` with `reaction_policy: native_cleaned_template_extraction_input`. Its `build-core` command computes a normalized four-Parquet Core, *not* an immediately searchable V2 candidate library.

```
rxni-build-library inspect --export-dir /path/to/aligned/export --reaction-column reaction
rxni-build-library build-core --export-dir /path/to/aligned/export --reaction-column reaction --dataset MY_DATASET --out-dir /new/output/core
rxni-build-library validate-core --core /new/output/core --rdkit --report /new/output/core_audit.json
```

The G2C split-only `train_reaction_core_v5.parquet` files are **not** complete source exports. Align and filter the membership and template tables from the same parent export, verify join/key coverage and lineage, then separately compile foundation V5 and V2 representation indexes. Never use the full-corpus indexes to evaluate held-out conditions.

**Current boundary:** The archived `build-auxiliary` and `assemble-candidate` commands still depend on nonpublic upstream `runtime_v6_5` / `runtime_v6_7` producers and a same-workspace runtime lookup. They are not verified standalone public index-build commands. Do not advertise or execute them as independent public capabilities until the producer migration and ABI parity tests pass. No conditional scientific accuracy is claimed.
