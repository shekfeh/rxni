# RXNI v0.1.0 — scope of current validation

## Verified in migration and R3

- 48-case frozen semantic oracle parity.
- 384 paired reaction cases tested for both semantic and retrieval parity (192 R3A, 128 R3B, 64 corrected CRD v0.23a R3B), not expert-labeled gold standard.
- R3S1: versioned source CSV/Parquet schema profiles; 9 unit tests and 12-record two-source A/B representation comparisons.
- R3S2: foundation SQLite profiles; 10 unit tests, 12 retrieval A/B queries and 24 sampled key-membership checks per dataset.
- R3C: standalone wheel packaging and installed CLI retrieval; API five-endpoint status checks passed for each of two small candidate libraries.

## R4 validation completed after R3

- R4C: metadata and lazy API import behavior verified; frozen 48/48 semantic oracle unchanged.
- R4D/E: Python 3.10 fresh-environment dependency resolution and installed wheel/CLI/API passed on both small fixtures; one environment does not establish every advertised Python/platform combination.
- R4F: documentation/profile files correctly packaged in a source distribution. Apache-2.0 source license selected during R4G, while resource provenance remains unapproved.

## Not yet established

- Full-size V2-library integrity and exact source coverage across v0.23a exports.
- Rigorous independent reaction family and subfamily accuracy, relevance ranking metrics, and reaction-condition prediction accuracy.
- Production load/concurrency safety beyond the serialized SQLite-backed API model.
- Full clean-install matrix and license/third-party data redistribution clearance.

An exact match against the frozen historical runtime is migration parity; it is not independent chemistry evidence.
