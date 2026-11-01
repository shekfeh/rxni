# RXNI USPTO 1976–2016 library — archive construction and Zenodo deposit

**Scientific basis:** historical USPTO extraction, native V5 structural/lookup foundation, v6.7 P6/Ring v0.3.1 semantic reprofile, preserved V2 physical storage ABI.  This is a dataset deposit **separate** from RXNI's source-code GitHub/Zenodo release.

## Repository and dataset

- RXNI software and tutorials: https://github.com/shekfeh/rxni
- Scientific source: `runs/retrieval_v67_reprofile/v67_retrieval_v2_full_candidate_P6_Ring031_01`
- Local five-file staging candidate: `runs/v680/library_distribution/uspto_historical_v67_ring031_r1`
- New **final distribution directory:** `runs/v680/library_distribution/uspto_1976_2016_rxni_library`
- New **single archive:** `runs/v680/library_distribution/uspto_1976_2016_rxni_library.zip`

The finalizer creates a new directory under the requested name, rather than renaming or modifying the original Ring031/scientific source.  It copies five retrieval files plus 42,128 SMARTS, creates version-independent filenames using `rxni-library.json`, verifies complete key and Parquet-location consistency, and retains the v6.7 science provenance.

## Commands (RXNI venv; run from your workstation)

```bash
ROOT="$HOME/0_LinuxSpace/chemRAG/reaction_intelligence"
RXNI="$HOME/Projects/Github_projects/rxni"
BASE="$ROOT/runs/v680/library_distribution"
LIB="$BASE/uspto_1976_2016_rxni_library"

# 1. Finalize existing staged P6/Ring031 scientific source to NEW name.
#    The destination must not yet exist.
"$ROOT/.venv/bin/python" "$RXNI/scripts/release/finalize_uspto_historical.py" \
  --staged "$BASE/uspto_historical_v67_ring031_r1" \
  --historical "$ROOT/artifacts/libraries/template_library_v5_native_v01" \
  --ring "$ROOT/runs/retrieval_v67_reprofile/v67_retrieval_v2_full_candidate_P6_Ring031_01" \
  --out "$LIB" \
  --smoke 3

# 2. Generate a comprehensive dataset README / real schema inventory
#    and one ZIP64 archive with full CRC validation.
"$ROOT/.venv/bin/python" "$RXNI/scripts/release/prepare_zenodo_dataset.py" \
  --library "$LIB" \
  --archive "$BASE/uspto_1976_2016_rxni_library.zip" \
  --nested-json all

# 3. Independently verify finalized folder and its real retrieval behavior.
"$ROOT/.venv/bin/python" "$RXNI/scripts/release/verify_published_library.py" \
  --library "$LIB" --smoke

# 4. Verify payload digests and the archive digest.
(cd "$LIB" && sha256sum -c SHA256SUMS)
(cd "$BASE" && sha256sum -c uspto_1976_2016_rxni_library.zip.sha256)
```

**If you already ran the prior finalizer** into `uspto_historical_r1_FINAL`, and the directory contains the six version-free files, `rxni-library.json`, `RELEASE_AUDIT.json`, and `SCIENTIFIC_PROVENANCE.json`, you may move that complete folder to the desired name instead of repeating step 1.  First inspect both paths and never overwrite an existing destination:

```bash
# ONLY if step 1 has already been completed, and LIB does not exist:
# mv "$BASE/uspto_historical_r1_FINAL" "$LIB"
```

## Contents of the archive

The ZIP stores one top-level `uspto_1976_2016_rxni_library/` folder containing:

| File | Role |
|---|---|
| `structural_index.sqlite` | candidate generation for six historical graph-hash channels |
| `reaction_lookup.sqlite` | keyed reaction/template metadata; includes one known metadata-only row |
| `reaction_representations.parquet` | 812,916 P6/Ring031 retrieval records, semantic analyses and fingerprints |
| `representation_locator.sqlite` | exact `reaction_key` → Parquet row-group / row-offset lookup |
| `retrieval_manifest.json` | scientific source/build identity, V2 storage ABI, hashes |
| `template_catalog.parquet` | 42,128 original retro SMARTS, numerical counts, hash-consensus fields |
| `rxni-library.json` | maps clean filenames to five required RXNI reader components |
| `README.md` | comprehensive file descriptions, full physical column data dictionary, complete discovered metadata JSON keys, usage examples |
| `NESTED_JSON_FIELDS.md` | observed nested JSON paths in all 812,916 records when `--nested-json all` is used; dynamic numeric and long hex keys normalized; depth capped |
| `SCHEMA_INVENTORY.json` | machine-readable Parquet/SQLite/JSON schema inventory |
| `SCIENTIFIC_PROVENANCE.json` | selected scientific source and unclosed research-review queue |
| `RELEASE_AUDIT.json` | technical validation, explicit rights/science gates |
| `SHA256SUMS` | SHA256 for every included file except itself |

The `.zip.sha256` sidecar is **outside** the ZIP.  Uploading that small sidecar alongside the ZIP to Zenodo is optional and useful.

## Interpreting the year range

The label **1976–2016** is the known coverage range of the parent USPTO patent-reaction corpus used by the standard literature, not a claim that each calendar year is represented or independently audited within this filtered 812,916-reaction subset.  The minimal runtime library omits patent publication dates and patent numbers.  For a source-period assertion specific to this derivative, attach the original ingest/source-manifest evidence.

## Zenodo upload workflow

1. Visit https://zenodo.org/ and create **New upload**. Select **Dataset** as the resource type. Use the separate dataset record rather than the automated GitHub software archive.
2. Upload `uspto_1976_2016_rxni_library.zip` (and optionally `uspto_1976_2016_rxni_library.zip.sha256`). The 1,444,488,348-byte ZIP archive is below Zenodo's default 50 GB per-record limit, subject to actual archive size.
3. Use the public title **RXNI USPTO 1976–2016 Reaction Retrieval Library — Release 1 (for RXNI v0.1.x)**; retain actual technical producer identifiers in the archive provenance.  Supply the actual creator(s), institutional affiliations and funding, rather than inferred credentials.
4. In *Related identifiers*, relate the dataset to `https://github.com/shekfeh/rxni` as **software** or **is supplemented by**, as appropriate; add a DOI for the RXNI software if already available. A GitHub account connection does **not** automatically attach a manually deposited dataset to a GitHub release.
5. Describe the *scope*: 812,916 retrieval records, 42,128 authentic retro SMARTS, reproducible hash matching, explicit structural and semantic runtime lineage; distinguish technical parity from scientific accuracy.
6. Do **not** select a data license without establishing rights to redistribute the source extraction/derivatives, including database rights where applicable. RXNI code's Apache-2.0 license is separate from the data license.
7. Keep as a Zenodo **draft** until the rights/provenance review, known scientific review case, and clean-install verification are concluded. A successful archive creation does not itself authorize public publication.
8. After public publication, copy the **actual DOI and archive SHA256** into the RXNI README and `docs/library-distribution.md`, and optionally into `CITATION.cff` for the dataset citation link. Do not invent or predict a DOI.

Zenodo documentation: https://help.zenodo.org/docs/deposit/create-new-upload/ and https://help.zenodo.org/docs/github/enable-repository/.  GitHub-Zenodo integration archives **GitHub source-code releases** automatically; manual dataset deposits are independent records.

## Verified 1 November 2026 technical candidate

The final ZIP has 1,444,488,348 bytes and SHA256 `5f31b689c942359c95e43debc0b8e6a8141fbb200d890fdea544ac3cd3a2a1ca`. A clean isolated Python 3.10 installation of the published `v0.1.0` software successfully installed, verified all 12 checksummed dataset files, checked 812,916 reactions and 42,128 templates, and retrieved nonempty candidates. The parent Lowe Figshare source is marked CC0; CC BY 4.0 was selected for this derivative Zenodo draft. No public DOI is supplied yet. Do not claim tomorrow's WSL test has run.

## Publication gates

- [ ] Validate the source parent corpus's stated year span against specific historical source metadata.
- [ ] Confirm redistribution rights for source reaction records, mapped SMILES, extracted templates and fingerprints; document attribution.
- [ ] Review the outstanding Ring031 scientific case in `SCIENTIFIC_PROVENANCE.json` and avoid overstating acceptance if unresolved.
- [ ] Inspect `RELEASE_AUDIT.json`, all 812,916 locator positions, full Parquet/SQLite schema and `SHA256SUMS` after running the finalizer.
- [x] Verify archive extraction and retrieval in a clean RXNI v0.1.0 environment on this development workstation.
- [ ] Independent fresh-machine WSL installation and tutorials acceptance (tomorrow).
- [ ] Verify no user-specific or private local paths, logs or temporary SQLite WAL/SHM files are included.
- [ ] Link the dataset DOI from the public RXNI repo only **after** deposit publication.

## Relevant scientific primary sources

- Thakkar A. et al. *Datasets and Their Influence on the Development of Computer Assisted Synthesis Planning Tools in the Pharmaceutical Domain*. Chemical Science (2019): https://doi.org/10.1039/C9SC04944D
- Original `reaction_utils`/`rxnutils` project (used for fingerprint-derived template identities): https://github.com/MolecularAI/reaction_utils
- Historical `aizynthtrain` templates pipeline: https://github.com/MolecularAI/aizynthtrain
