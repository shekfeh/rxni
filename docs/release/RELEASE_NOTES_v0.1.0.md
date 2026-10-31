# RXNI v0.1.0 — Initial Public Release

RXNI is a Python 3.10+ library for structural analysis, semantic classification,
retrieval and method/condition evidence for mapped chemical transformations.

## Included in v0.1.0

- Public Python API, CLI and optional FastAPI service.
- Mapping-aware reaction graph analysis, directional FGI semantics and non-FGI
  transformation classification.
- Configurable reaction-library and storage schema profiles.
- Retrieval against user-provided compatible reaction libraries, with provenance
  and compatibility validation.
- ChemKN provider interfaces without bundled proprietary knowledge assets.

## Validation summary

- Independent private v6.7/public RXNI production-API parity: **400/400**
  accepted masked queries with identical output, across CRD25 and USPTO17–25.
- Installed-wheel full-index retrieval parity against the accepted public
  baseline: **48/48**, with no response mismatches (24 per dataset).
- Separate retrospective train-only condition-component evaluation on 400
  held-out queries: **81.25%** have at least one recorded component recovered
  among the top five consensus predictions; **31%** recover all recorded
  components. These are raw-component recovery measures, not validated
  solvent/reagent/catalyst-role prediction accuracies.
- PEP 517 wheel/sdist, Python 3.10 clean-install/API/CLI checks and Twine
  distribution-metadata validation passed for the certified release candidate.

## Data and licensing boundary

The FGI rules were developed within Reaction Intelligence and refined with
heuristic tests using extracted USPTO 1960–2016 reaction/template data. The
raw reaction corpora, patent-derived training records, private held-out truth,
model weights and proprietary ChemKN knowledge assets are **not included** in
this software distribution. The underlying source datasets retain their own
terms and are not licensed or redistributed by this release.

Apache-2.0 applies to original RXNI project code, subject to the terms and
notices for any separately identified third-party components.

## Installation

```bash
python3.10 -m pip install './rxni-0.1.0-py3-none-any.whl[api,cli]'
```

Run this command in the directory containing the GitHub Release wheel.
Users supply and configure their own compatible reaction libraries; none is
downloaded by RXNI.

For scope and scientific limitations, consult the project documentation,
including `docs/release-limits-v0.1.md` and
`docs/provenance-summary-v0.1.md`.
