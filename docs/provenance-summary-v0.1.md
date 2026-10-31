# RXNI v0.1 semantic resource provenance

The R4H audit classified 11 non-Python scientific files collected from the
2026-10-30 RXNI development tree. The following is a technical release record,
not a legal copyright opinion.

- Retained as runtime-required: `src/rxni/semantic/contracts/fgi_contracts_v1.json`.
  This contains 82 directional-FGI validation contracts; the runtime
  `semantic/contracts/registry.py` directly loads it. SHA256:
  `c5cfca057c283068b40a3dfb1e881c294c73e6196353deff2d5f5242c60667c2`.
- Excluded from public source and distributions: two earlier port provenance
  records, two validation specifications, one duplicate draft contract registry
  and five audit TSVs. All original bytes are preserved in the private
  R4I source-backup evidence directory.
- The omitted 509-reaction audit includes 19 known-failure examples with
  reaction SMILES, requiring independent source-data provenance review before
  any future redistribution.
- The contract JSON has traceable lineage to RXNI's private frozen scientific
  baseline; attribution and ownership remain to be explicitly approved
  before claiming the public release is legally cleared.
- No ChemKN textbook datasets or publisher-specific loaders are bundled.

The Apache-2.0 license covers only materials the authors are authorized
to license. External scientific knowledge providers are independent.
