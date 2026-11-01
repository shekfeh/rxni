# RXNI Historical USPTO library: distribution, audit and publication plan

**Current selected source:** `runs/retrieval_v67_reprofile/v67_retrieval_v2_full_candidate_P6_Ring031_01`.

**Data status (1 November 2026):** the 1,444,488,348-byte ZIP is prepared, full checksums and clean v0.1.0 installation/retrieval passed. The original Ring031 scientific source remains an unpromoted candidate; the known review case is disclosed. The dataset was being uploaded to a Zenodo draft; no verified public DOI/URL has been supplied. Full derivative-source rights review and tomorrow's independent WSL test remain open.

## Distribution specification

Release label: **RXNI Historical USPTO Reaction Library — Release 1**. The public title intentionally omits the implementation words “V2” and “System-I”. The internal physical ABI and scientific provenance do not.

| File | Purpose |
|---|---|
| `structural_index.sqlite` | Candidate-generation graph structural hashes |
| `reaction_lookup.sqlite` | Reaction-string and template/member metadata keyed by `reaction_key` |
| `reaction_representations.parquet` | 812,916 frozen v6.7 P6/Ring031 stored reaction profiles and fingerprints |
| `representation_locator.sqlite` | `reaction_key` → Parquet row group/offset, build ID and V2 generation |
| `retrieval_manifest.json` | Build identity, scientific contract, hashes, representation generation |
| `template_catalog.parquet` | 42,128 original retrospective SMARTS |
| `rxni-library.json` | Maps human-facing names to exactly five runtime artifacts |
| `SCIENTIFIC_PROVENANCE.json`, `RELEASE_AUDIT.json` | Selected source lineage, unresolved gates and audited results |
| `SHA256SUMS`, `README.md` | File integrity and usage/rights notice |

Use `rxni.retrieval.libraries.resolve_library()` to load the descriptor. Keep internal V2 Parquet fields and the SQLite table `representation_locator_v2` untouched. A filename rename is not an ABI migration.

**Selected archive:** `uspto_1976_2016_rxni_library.zip`
**SHA256:** `5f31b689c942359c95e43debc0b8e6a8141fbb200d890fdea544ac3cd3a2a1ca`
**Extracted folder:** `uspto_1976_2016_rxni_library/`
**Compatible software:** released `v0.1.0` (independent clean-venv retrieval verified).

## Confirmed historical evidence

- 812,916 Ring031 representations, 397 row groups, exact structural-key↔Ring031 locator key set.
- Historical lookup has 812,917 rows; one extra metadata-only key, unreachable from structural candidates. Do not erase it without an identity/provenance justification.
- 42,128 complete SMARTS, zero blank/missing lookup templates, 42,128 RDKit parsable, 35 template warnings, zero reported RDKit validation errors.
- Full `rxnutils` historical hash replay matched 42,128/42,128 native template hashes with RDKit 2023.09.6.
- 12/12 Ring031 stored-record/source and public live-profile agreement; 48/48 parity on separate older v6.6 candidate only.
- Stored semantic authority pins v6.7 P6/Ring031 detector/resolver SHA256; V2 builder prefix `v66_` is a storage identifier, not scientific semantics.

## Finalizer and verifier

```bash
export ROOT="$HOME/0_LinuxSpace/chemRAG/reaction_intelligence"
export RXNI="$HOME/Projects/Github_projects/rxni"

"$ROOT/.venv/bin/python" "$RXNI/scripts/release/finalize_uspto_historical.py" \
 --staged "$ROOT/runs/v680/library_distribution/uspto_historical_v67_ring031_r1" \
 --historical "$ROOT/artifacts/libraries/template_library_v5_native_v01" \
 --ring "$ROOT/runs/retrieval_v67_reprofile/v67_retrieval_v2_full_candidate_P6_Ring031_01" \
 --out "$ROOT/runs/v680/library_distribution/uspto_1976_2016_rxni_library" \
 --smoke 3

"$ROOT/.venv/bin/python" "$RXNI/scripts/release/verify_published_library.py" \
 --library "$ROOT/runs/v680/library_distribution/uspto_1976_2016_rxni_library" --smoke
```

The finalizer refuses to overwrite `--out`, captures the original scientific receipt SHA256 rather than distributing private absolute path information, verifies every key and locator coordinate and creates a new distribution. It does **not** certify source rights or independent chemical accuracy.

## Remaining gates before release

- [x] Create ZIP, check SHA256 and perform independent clean-venv technical retrieval on the development machine.
- [ ] Perform the second-machine WSL acceptance tomorrow; no result is claimed today.
- [ ] Verify final published DOI and update download instructions after Zenodo publication.
- [ ] Confirm remaining derivative inputs and attribution; parent Lowe Figshare source is marked CC0 and derivative Zenodo CC BY 4.0 was selected.

- [ ] Independent review/approval of third-party corpus and derivative SMARTS redistribution rights, original source names, attribution and data license. Do not infer rights from RXNI Apache-2.0 software license or USPTO patent publication alone.
- [ ] Adjudicate or accurately disclose known Ring031 scientific review case and any other scientific acceptance limits; preserve source-receipt status.
- [ ] Run the finalizer and inspect `RELEASE_AUDIT.json` and exact `SHA256SUMS` on **the user's real data**, then verify from an isolated fresh RXNI environment with only the extracted library.
- [ ] Document the extra metadata-only lookup row and the one-row-per-reaction lookup policy for multi-template membership.
- [ ] Decide dataset version, DOI host and immutable release citation; include exact checksum and source lineage in Zenodo metadata.
- [ ] Inspect archives for unwanted local paths, private evidence, unlicensed contents and identifiable/private data.
- [ ] Publish and test download/extraction and full checksums; then update RXNI README/docs with the actual URL/DOI.

**Hosting:** Keep only documentation and small validation scripts in Git. GitHub ordinary repository objects are limited to 100 MB; GitHub Release assets are limited to less than 2 GiB each. Zenodo accepts 50 GB across up to 100 files per record. Prefer Zenodo for these ~2+ GB data and link the immutable DOI from GitHub. Confirm current host limits at publication.

## Related docs

- `tutorials/05-download-and-use-library.md`
- `tutorials/06-build-library-from-sources.md`
- `tutorials/07-template-catalogue-and-provenance.md`
- `library-building.md`, `input-schema-profiles.md`, `foundation-schema-profiles.md`, `licensing-and-data.md`.
