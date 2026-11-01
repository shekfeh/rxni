# Tutorial 06 — Build an RXNI reaction library: source data → Core → searchable retrieval

**Audience:** advanced cheminformatics researchers and RXNI library maintainers.  **Software contract:** RXNI 0.1.0 (verify against the installed release).  **Status:** executable producer stages exist in the public `rxni.builder` package, but a new candidate is not a scientifically approved dataset until all stated gates pass.

> **Two distinct goals:** (A) build a *new*, independently sourced RXNI dataset using the current public builder; (B) reproduce the **historical USPTO v6.7 P6/Ring031** library exactly. This tutorial provides the public pipeline for (A) and the lineage/acceptance requirements for (B). Running the current CLI does **not** automatically regenerate the previously frozen historical Ring031 semantic artifact bit-for-bit.

## 1. Architecture and naming boundaries

```text
Raw source reaction records / mapped reaction SMILES
     | acquisition and normalization (external: rxnutils / AiZynthTrain)
     v
Template extraction: RetroTemplate + TemplateHash + reaction membership
     | verified native-cleaned reaction and identity mapping
     v
Aligned source export: 3 Parquet files + manifest.json
     |  rxni-build-library inspect / build-core / validate-core
     v
RXNI Core: reactions.parquet | templates.parquet | memberships.parquet | evidence.parquet | manifest.json
     |  rxni-build-library build-auxiliary
     v
Native foundation: structural index | member lookup | template SMARTS catalog | intermediate tables
     |  stored retrieval representation materialization (v6.7 query policy)
     v
V2-compatible retrieval: stored profiles | locator SQLite | manifest
     |  rxni-build-library assemble-candidate + separate certification
     v
Portable six-file public distribution + rxni-library.json + checksums
```

The versionless **library display name** is independent of the internal ABI: `retrieval_representations_v2` remains the frozen storage generation. The V5 Core and native foundation are distinct concepts. Structural hash definitions, semantic policy, scoring equations, template identity and storage layout are separate compatibility contracts. Do not silently change one by editing another.

## 2. Upstream processing: rxnutils and AiZynthTrain

### Sources and environment

Use publicly authorized source records and preserve exact source versions, identifiers, reaction mapping and cleaning procedures. Historical USPTO patent-reaction corpora have been distributed by multiple projects with different subsets, mapping and filtering; a generic “USPTO” label is insufficient provenance.

The archived [AiZynthTrain repository](https://github.com/MolecularAI/aizynthtrain) describes its USPTO training pipeline. Its README now points to [AiZynthModels](https://github.com/MolecularAI/aizynthmodels) as the successor. For **historical reproduction**, pin the original pipeline rather than silently substituting a newer model/data-cleaning package. Upstream preprocessing is provided by [reaction-utils/rxnutils](https://github.com/MolecularAI/reaction_utils), including the mapped `uspto_data_mapped.csv` described in AiZynthTrain's README.

The original AiZynthTrain flow is:

1. Prepare USPTO source reaction data using `rxnutils` to create `uspto_data_mapped.csv`.
2. Create the pipeline work folder using `python configs/uspto/setup_folder.py PATH_TO_YOUR_FOLDER` in the selected AiZynthTrain checkout.
3. Adapt and execute the repository's `configs/uspto/template_pipeline.sh` (written as SLURM workflow; adapt launcher/files/paths for local machines).
4. Retain per-reaction `RetroTemplate` / `TemplateHash`, source keys, mapped/native-cleaned forward reaction strings and all extraction errors.
5. Treat the outputs `uspto_template_library.csv` / `uspto_ringbreaker_template_library.csv` as template-selection artifacts, not automatically as a complete RXNI source export. The expansion-model pipeline is optional for building an RXNI retrieval library; training a neural retrosynthesis model is **not** an RXNI prerequisite.

Inspect the exact upstream checkout's options instead of assuming they never changed:

```bash
conda run --no-capture-output -n aizynthtrain \
  python -m aizynthtrain.utils.template_runner --help
```

In the inspected historical implementation, the template runner accepts `--input_path`, `--output_path`, required `--radius`, and required `--smiles_column`; it also exposes `--expand_ring`, `--expand_hetero`, `--ringbreaker_column` and optional `--batch` ranges. Configure these to the **real column names in the prepared mapped reaction data**, not an invented default. For example, only after verifying the input column `ReactionSmiles` exists:

```bash
conda run --no-capture-output -n aizynthtrain python \
  -m aizynthtrain.utils.template_runner \
  --input_path /data/prepared/reactions.csv \
  --output_path /data/derived/templates.csv \
  --radius 1 \
  --smiles_column ReactionSmiles
```

This is a **syntax illustration**, not the reconstructed historical USPTO extraction configuration. The original pipeline may apply additional selection, ring expansion and cleaning. Its pinned configuration and output checksums are the authority.

### Correct historical template-hash convention

AiZynthTrain calls:

```python
_, retro_template = reaction.generate_reaction_template(...)
smarts = retro_template.smarts
native_template_hash = retro_template.hash_from_bits()
```

For the historical 42,128-template catalogue, `ReactionTemplate(smarts, direction='retro').hash_from_bits()` reproduced **42,128/42,128** stored template IDs with `rxnutils` in the `aizynthtrain` environment (RDKit 2023.09.6). This hash uses sorted fingerprint-bit differences, *not* SHA256 of the raw SMARTS string.

The current public RXNI **Core source bridge instead requires**:

```python
retro_template_hash = sha256(exact_retro_smarts.encode('utf-8')).hexdigest()
```

That strict surrogate is used by `rxni.builder.bridge.prepare` and `rxni.builder.knowledge.make_core`. It **cannot** accept the historical AiZynthTrain `TemplateHash` values as the field named `retro_template_hash` unchanged. A schema-profile field rename is NOT a hash-policy conversion.

**Required adapter for newly imported AiZynthTrain data:** in a distinct, audited export step, preserve the original `TemplateHash` as `native_template_hash`, derive a new RXNI-Core `retro_template_hash=SHA256(exact SMARTS)` and re-key the reaction→template membership table using the same one-to-one verified mapping. Preserve both identities and never alter source SMARTS to force a match. If two different native template identities map to a single raw-SMARTS digest, record that relationship explicitly and reject any ambiguous or accidental many-to-many join; do not invent identity equality.

For reproduction of the historical V5/native retrieval library, retain its original fingerprint-based identifiers in the **historical path** and do not rebuild or rewrite them merely to satisfy the new Core ingestion policy. The two lineages require different, explicitly versioned identity profiles.

## 3. Build a complete, aligned export from upstream results

### Required files and minimum fields

Directory `SOURCE_EXPORT`:

| File | Necessary columns | Purpose |
|---|---|---|
| `reaction_core_v5.parquet` | Source-scoped `reaction_key`, authenticated mapped/native-cleaned `reaction` | Forward reaction identity and source records |
| `reaction_template_membership_v5.parquet` | `reaction_key`, `retro_template_hash` | Reaction-to-template membership and documentary evidence |
| `template_library_v5.parquet` | `retro_template_hash`, `retro_smarts` (or `retro_template_smarts`), optional `retro_template_number`, `native_template_hash` | Original SMARTS + *distinct* identity namespaces |
| `manifest.json` | `reaction_policy`, dataset/source lineage, source hashes | Source legitimacy and policy attestation |

Exact accepted source columns are implemented in `src/rxni/builder/bridge.py`. `rxni-build-library inspect` reports the resolved mapping. **A source `reaction_key` joining export tables is not necessarily the final RXNI Core reaction identity**: the Core writer re-derives reaction identity from SHA256 of the exact mapped reaction string and preserves source-scoped IDs as evidence.

The export manifest must explicitly contain:

```json
{
  "reaction_policy": "native_cleaned_template_extraction_input",
  "dataset": "MY_LICENSED_DATASET",
  "source_lineage": {
    "input_dataset": "record exact source release and URL",
    "cleaning_pipeline": "record pinned git commit and parameters",
    "template_hash_policy": "source native hash_from_bits retained separately; RXNI bridge uses sha256 exact SMARTS"
  }
}
```

This is a minimal **illustrative manifest**, not a complete provenance record. Keep separate checksums of source CSV, mapping output, template-selection output, final Parquet exports and Core. For historical v0.23a sources, the `raw_rebuilds/crd25_from_raw_v023_a` and `raw_rebuilds/uspto_17_25_from_raw_v023_a` directories contain rebuilt **CSV** sources; the export Parquet belongs downstream under `library_candidates/*/export_final/`. Record the actual CSV→Parquet lineage; never point a Parquet reader at a raw-rebuild directory.

**Full export invariants (fail closed):** all mapped reaction strings have reactant and product atom maps, no dangling membership keys, all template hashes resolve to one authentic SMARTS string, no template without a membership, and no reaction→template evidence is silently dropped. Preserve source observation identifiers and duplicate counts. Never use a mapper-only reaction when the template was extracted from a different cleaned reaction.

### Inspect a source export before building

```bash
export SOURCE_EXPORT='/absolute/path/to/aligned/export_final'
rxni-build-library inspect \
  --export-dir "$SOURCE_EXPORT" \
  --reaction-column reaction
```

The following direct sanity check reveals the hash convention rather than assuming it:

```bash
python - "$SOURCE_EXPORT" <<'PY'
import hashlib, sys
from pathlib import Path
import pyarrow.parquet as pq
p=Path(sys.argv[1])/'template_library_v5.parquet'
pf=pq.ParquetFile(p)
cols=pf.schema_arrow.names
scol='retro_smarts' if 'retro_smarts' in cols else 'retro_template_smarts'
rows=pf.read_row_group(0,columns=['retro_template_hash',scol]).slice(0,10).to_pylist()
for r in rows:
    k=r['retro_template_hash']; s=r[scol]
    print(k==hashlib.sha256(s.encode('utf8')).hexdigest(), k)
PY
```

The public Core builder **rejects** nonmatching hashes. If this test yields `False`, do not rerun with a misleading column mapping. Use the explicit native-ID→Core-ID rekey workflow above, while preserving the original metadata and source lineage.

## 4. Build and validate the immutable RXNI Core

```bash
export RXNI_ROOT='/path/to/local/rxni/checkout'
python -m pip install -e "$RXNI_ROOT"

export BUILD_ROOT='/path/to/private/output'
mkdir -p "$BUILD_ROOT"

rxni-build-library build-core \
  --export-dir "$SOURCE_EXPORT" \
  --reaction-column reaction \
  --dataset MY_LICENSED_DATASET \
  --out-dir "$BUILD_ROOT/core"

rxni-build-library validate-core \
  --core "$BUILD_ROOT/core" \
  --rdkit \
  --report "$BUILD_ROOT/core_audit.json"
```

The output Core contains exactly:

```text
core/
  reactions.parquet      # one reaction identity, mapped SMILES, evidence count
  templates.parquet      # hash, authentic SMARTS, native hash, source metadata
  memberships.parquet    # reaction/template membership + occurrence count
  evidence.parquet       # individual source observations
  manifest.json          # counts + SHA256 provenance, no retrieval indexes
```

Core identity is deliberately independent of structural/semantic representations. Core is **not yet searchable by RXNI's indexed retrieval**. A brute-force Core search is not part of the current public five-file retrieval ABI.

The emitted Core manifest initially records `CORE_BUILT_PENDING_AUDIT`; keep the independent `validate-core` report as proof that counted source evidence and all hashes were checked. Do not silently relabel the Core manifest as scientifically approved.

## 5. Build the native structural foundation and retrieval auxiliary

Public v0.1.0 has installed modules for the native builder, lookup cache and retrieval representation producer. Its older `docs/library-building.md` contained a stale warning about nonpublic producer dependencies; the public code now resolves `rxni.builder._native.build_all`, `rxni.builder.lookup_cache` and `rxni.retrieval.builders.representations`. However, **executable producers alone do not establish complete scientific equivalence with an accepted historical library**.

```bash
rxni-build-library profiles
rxni-build-library inspect-runtime \
  --root "$RXNI_ROOT" \
  --runtime v6.7
```

Create isolated build outputs (do not overwrite existing libraries):

```bash
rxni-build-library build-auxiliary \
  --core "$BUILD_ROOT/core" \
  --work-dir "$BUILD_ROOT/auxiliary" \
  --root "$RXNI_ROOT" \
  --runtime v6.7 \
  --workers 4 \
  --validation-sample 500 \
  --chunk-size 2048 \
  --through v2
```

The builder's stages are:

| Stage | Code module | Product and reason |
|---|---|---|
| 1. Normalize | `_native.normalize` | Source evidence, clean reactions, membership staging SQLite |
| 2. Graph keys | `_native.structural_keys` | Versioned graph-edit and local-center structural hashes |
| 3. Native representations | `_native.representations` | Historical System-I representation intermediate; not the public V2 records |
| 4. Member rows | `_native.retrieval_members` | Candidate identity, template metadata and provenance |
| 5. Template aggregation | `_native.templates` | Authentic SMARTS and template-member counts, exact-source preservation |
| 6. Structural indexes | `_native.indexes` | `structural_index_v1.sqlite` and auxiliary index files |
| 7. Validate | `_native.validate` | Native structural/template integrity checks |
| 8. Retrieval representation materialization | `rxni.retrieval.builders.representations` | Stored V2 Parquet, locator and manifest using the installed semantic engine |
| 9. Lookup | `rxni.builder.lookup_cache` | `runtime_lookup_v1.sqlite` keyed by reaction identity |

These stages may create multiple gigabytes of intermediate SQLite/Parquet data. The historical example occupied ~5 GB in its native workspace while its core five-file runtime subset occupied ~2.1 GB for the older v6.6 semantic version. Plan ample temporary disk space and record software/environment snapshots.

The CLI writes `auxiliary_build_status.json` with status **`CANDIDATE_SOURCES_READY_PENDING_PARITY`**, not a release approval. It verifies native template preservation through `rxni.builder.native_verify` but still explicitly marks semantic compatibility and ranking parity as untested.

### Safe recovery after interrupted building

The builder supports `--resume`, but only when the original locked Core, source SHA256 hashes and build parameters match. Do not use `--resume` to paper over a chemistry-policy change; start a fresh workspace instead.

```bash
rxni-build-library build-auxiliary \
  --core "$BUILD_ROOT/core" \
  --work-dir "$BUILD_ROOT/auxiliary" \
  --root "$RXNI_ROOT" \
  --runtime v6.7 --workers 4 \
  --validation-sample 500 --chunk-size 2048 \
  --through v2 --resume
```

### Known modeling limitation in the lookup cache

The present `lookup_cache.py` creates `member_lookup` with `reaction_key` as its primary key and uses `INSERT OR REPLACE`. If one reaction key legitimately has **multiple template memberships**, the lookup is only one row per reaction and cannot faithfully expose every member's template assignment. The historical baseline has one extra metadata-only entry; its reachable key-set and latest stored locator match, but this does not eliminate the general one-to-many issue in newly built corpora. Quantify duplicate memberships, define the expected lookup policy and review any provenance loss before claiming one-to-many completeness.

## 6. Assemble the five-file *candidate*

```bash
rxni-build-library assemble-candidate \
  --core "$BUILD_ROOT/core" \
  --work-dir "$BUILD_ROOT/auxiliary" \
  --root "$RXNI_ROOT" \
  --runtime v6.7 \
  --lookup-db "$BUILD_ROOT/auxiliary/native/runtime_lookup_v1.sqlite" \
  --out-dir "$BUILD_ROOT/candidate"
```

The command requires the lookup database from the **same Core/build workspace**. It creates SQLite backup snapshots, preserves the corresponding Parquet/manifest, records SHA256 checks and emits `candidate_manifest.json` with `CANDIDATE_PENDING_SEMANTIC_AND_FUNCTIONAL_VALIDATION`. It is not a publisher.

Never splice a new Core's structural or metadata SQLite files into historical representations just because row counts or column names look similar. All reaction keys, graph-hash channels, template identities, locator positions and stored build IDs must agree exactly.

### Optional source-schema profiles for nonstandard input

The **source-schema profile** described by `docs/input-schema-profiles.md` governs input CSV/Parquet mapping to a representation producer. It does **not** change the physical V2 ABI. A separate **foundation-schema profile** described by `docs/foundation-schema-profiles.md` maps validated SQLite names and columns. A third layer, `rxni-library.json`, maps distributable filenames. These are independent and are not chemical identity conversions.

For example, inspect an explicitly profiled Core-derived representation source:

```bash
python -m rxni.retrieval.builders.source_schema inspect \
  --source-core /path/to/reaction_core_v5.parquet \
  --input-schema "$RXNI_ROOT/examples/profiles/crd25-uspto-v023a.json"
```

To build V2 independently, the producer's `--library`, `--input-core`, `--input-schema`, `--output-dir`, `--workers`, `--chunk-size` and `--resume` options must be used with an actual *aligned* foundation; `--input-core` does not magically create the structural/lookup SQLite databases.

## 7. Distinguish new-library building from historical P6/Ring031 reproduction

The **selected historical library** was rebuilt and attested in an earlier scientific campaign:

```text
$ROOT/runs/retrieval_v67_reprofile/v67_retrieval_v2_full_candidate_P6_Ring031_01/
  retrieval_representations_v2.parquet
  representation_locator_v2.sqlite
  retrieval_v2_manifest.json
  V67_TRUTHFUL_SEMANTIC_PROVENANCE.json
  v67_build_input_snapshot.json
  full_profile_comparison.json
  p6_to_ring_comparison/
  ...
```

Its manifest pins detector and resolver SHA256 hashes and states v6.7 post-P6/Ring031 semantic authority. It contains 812,916 stored records and 397 Parquet row groups. Of those, 145,632 `semantic_profile_json`, 13,484 reaction centers, 31,179 localized fingerprints and 58,186 multichannel fingerprints differed from the original historical v6.6 representation. Differences in all `semantics_json` records can reflect diagnostics/provenance, not necessarily changed primary-event assignments.

Its stored-representation **build ID begins `v66_`** because the unchanged V2 builder retains that opaque prefix. Do not interpret this as stored v6.6 semantics. Its physical schema is V2 and its live semantic generation was via the accepted v6.7 implementation. The source provenance explicitly says `CANDIDATE_BUILD_COMPLETE_NOT_PROMOTED`; promotion is a separate chemical and retrieval validation decision.

For that historic dataset, an independently built Core with SHA256 SMARTS IDs could have **different template/key namespace** from the frozen AiZynthTrain native foundation. Therefore do not substitute the newly built Core for the historical V5 structural indexes or lookup without a complete, verified identity crosswalk and full regression.

## 8. Certification and key-connectivity requirements

Minimum acceptance for **every** new dataset (not only the historical example):

1. Full SHA256 verification of actual source/export/Core/auxiliary files; build ID coherence across stored Parquet, locator and manifest.
2. `PRAGMA integrity_check` on each distributed SQLite snapshot (use SQLite backup API if sources have WAL files).
3. Full exact key-set agreement between structural, metadata and representation locator, allowing only explicitly explained unreachable metadata-only records; no missing reachable candidate representation.
4. Full Parquet-key-to-locator **row group and row offset** verification, not only key-set equality.
5. Full template hash/SMARTS coverage and parseability. For the historical AiZynthTrain identity policy, use `hash_from_bits()` in a pinned environment rather than raw SMARTS SHA256.
6. Retrieval parity between the staged source and candidate distribution on diverse mapped reactions; inspect candidate keys, templates and score vectors. A fixed number of self-hit examples is only portability evidence.
7. Query-side semantic vs stored semantic compatibility on accepted examples, including FGI and non-FGI competitors. Review specific changed classifications separately.
8. Independent chemical validation and named-reaction challenges, with known remaining review cases reported honestly.
9. Dataset source/derived-data redistribution-rights audit, attribution and license separate from RXNI Apache-2.0 software.
10. Clean installation and offline library test using only the released archive, with no legacy tree mounted.

**Do not use a library containing train+test reactions to claim held-out accuracy.** Record split provenance and test information leakage in any benchmark.

## 9. Create a *portable* publication candidate

This repository provides `scripts/release/finalize_uspto_historical.py` to package the **specific validated historical P6/Ring031** dataset from its known components, and `scripts/release/verify_published_library.py` to check the final archive after extraction. The historical packaging script is **not a general dataset builder**.

```bash
export ROOT="$HOME/0_LinuxSpace/chemRAG/reaction_intelligence"
export HIST="$ROOT/artifacts/libraries/template_library_v5_native_v01"
export RING="$ROOT/runs/retrieval_v67_reprofile/v67_retrieval_v2_full_candidate_P6_Ring031_01"
export STAGED="$ROOT/runs/v680/library_distribution/uspto_historical_v67_ring031_r1"
export FINAL="$ROOT/runs/v680/library_distribution/uspto_historical_r1_FINAL"

python "$RXNI_ROOT/scripts/release/finalize_uspto_historical.py" \
  --staged "$STAGED" \
  --historical "$HIST" \
  --ring "$RING" \
  --out "$FINAL" \
  --smoke 3

python "$RXNI_ROOT/scripts/release/verify_published_library.py" \
  --library "$FINAL" --smoke
```

This finalizer:

- Creates a **new** output directory and never mutates the historical/native/Ring031 sources or previous staged results.
- Snapshots three SQLite databases safely and copies the immutable Ring031 representation Parquet.
- Includes `template_catalog.parquet` (authentic 42,128-template SMARTS).
- Emits `rxni-library.json` so all six public data files use human-friendly filenames.
- Corrects the **physical** locator SHA256 after SQLite re-snapshotting while preserving the original scientific build ID and representation SHA256. Removes only private absolute paths and the unshipped locator Parquet checksum from the distributed manifest.
- Verifies **all** Parquet keys and positional locator coordinates, cross-DB key coverage and SMARTS parse/validation; optionally compares retrieval signatures with the staged source.
- Records original Ring031 manifest/provenance SHA256 and accepted source-rule pins in `SCIENTIFIC_PROVENANCE.json`, plus publication gate status and technical checks in `RELEASE_AUDIT.json`.
- Writes `SHA256SUMS`. The finalizer **does not** authorize publication; the output is explicitly marked pending scientific and data-rights review.

Do not modify internal `representation_generation`, locator table name, stored Parquet column names, semantic fields or `build_id` to create an apparently versionless dataset; only the **public filesystem names** are renamed.

## 10. Data release and publication layout

Place docs, validation code and a release-page link in Git. Do **not** add multi-gigabyte Parquet or SQLite files directly to the source repository. Prefer a data repository such as Zenodo with DOI, full provenance, data-rights attestation, archive checksums, and a link from the RXNI GitHub README. A GitHub Release may carry appropriately sized compressed data assets, but each individual release asset must remain below the current 2 GiB limit. The full uncompressed six-file folder is larger than this; use an approved archive/split strategy if opting for GitHub Releases.

After publication, use Tutorial 05 for installation and Tutorial 07 for reading the template catalogue, evidence limitations and scientific provenance. For an API integration, consult `docs/api-quickstart.md`.

## 11. Failure modes worth distinguishing

| Failure | Interpretation / recommended action |
|---|---|
| `Template hash mismatch ... != SHA256(SMARTS)` | Historic native fingerprint hash supplied to strict Core bridge; construct an **audited identity adapter**, do not silently rehash the historical native catalogue. |
| `Orphan membership` | Reaction or template absent from the *same* source export; fix join lineage before building. |
| `Core reaction is not identical ... native cleaned` | Wrong reaction representation; preserve the template-extraction input exactly. |
| `Stored representation generation mismatch` | ABI mismatch, not a filename problem; no automatic storage conversion. |
| Source and locator keys differ | Mixed datasets, incompatible reaction-key policy, partial build or stale indexes. |
| `CANDIDATE_SOURCES_READY_PENDING_PARITY` | Builder produced files; independent ranking/scientific validation is not yet executed. |
| 35 RDKit SMARTS warnings | Historical template parse/validation warnings; do not claim template applicability has been tested. |

## 12. Primary code references

- `rxni.builder.cli`, `rxni.builder.bridge`, `rxni.builder.knowledge`, `rxni.builder.auxiliary`, `rxni.builder.native_runtime`, `rxni.builder.native_verify`.
- `rxni.builder._native.build_all` and `rxni.builder.lookup_cache`.
- `rxni.retrieval.builders.representations`, `rxni.retrieval.record`, `rxni.retrieval.storage.provider`, `rxni.retrieval.libraries`.
- [AiZynthTrain](https://github.com/MolecularAI/aizynthtrain) `aizynthtrain/utils/template_runner.py` and its README / USPTO pipeline configs.
- [reaction-utils](https://github.com/MolecularAI/reaction_utils) `rxnutils/chem/template.py` `ReactionTemplate.hash_from_bits()`.
- Existing RXNI docs: `input-schema-profiles.md`, `foundation-schema-profiles.md`, `library-formats.md`, `licensing-and-data.md`, `distribution-rights-checklist.md`.

For historical reconstruction, save the exact git commits, source licenses, upstream CSV/Parquet hashes, source cleaning/mapping parameters, template extraction configuration, and final RXNI code SHA256 before relying on reproducibility claims.
