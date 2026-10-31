# RXNI v0.1.0: scientific validation scope

This is a public-facing summary of the **bounded** RXNI v0.1.0 technical evaluation.
The underlying per-query evidence, patent-derived input records and dataset artifacts remain private.

- **Public/private production parity:** 400/400 held-out queries produced strictly equal compared API projections between private runtime V6.7 and public RXNI, using shared training-only retrieval libraries.
- **Full-index public retrieval smoke:** 400/400 held-out queries returned HTTP 200 with top-10 candidates; no foreign indexed reaction keys were observed.
- **Library sizes:** 11,953 indexed CRD25 transformations and 11,857 indexed USPTO 2017–2025 transformations; held-out transformation overlap with the indexed libraries was zero.
- **Training-only raw component benchmark:** 400 held-out reactions evaluated. Top-five rank-weighted component consensus yielded 81.25% of reactions with at least one matching recorded component, and 31.0% with every distinct recorded component recovered. Mean per-reaction recall was 57.35%, versus 30.22% for the global training-frequency baseline.

## Interpretation and limits

The benchmark evaluates retrospective recovery of raw reaction middle-field components, **not** verified solvent/reagent/catalyst roles, reaction success, feasible synthetic conditions, or a trained reaction-condition predictor. Exact component matching preserves stereochemistry, removes atom maps and does not normalize tautomers or salts. Transformation and document-level separation were audited; a patent-family-disjoint split has not been established. The pooled 400-query sample spans two separately prepared datasets and is not a measure of overall chemistry coverage.

This document is not publication or redistribution clearance for any licensed third-party dataset, retained FGI contract or source material.
