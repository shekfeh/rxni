# RXNI v0.1.1 — strict input mapping and evidence boundaries

## Public input acceptance

Queries **must contain complete, consistent atom mapping** on the reactant and product sides. RXNI rejects missing atom maps, partially mapped heavy atoms, duplicate map labels on either side, element changes under the same map label, and mappings with fewer than two persistent mapped heavy atoms. Agents are exempt from atom-map coverage. This policy is a conservative public-query gate, not a claim that any accepted mapping is chemically or mechanistically correct. Previously accepted unmapped inputs are **intentionally incompatible**.

For read-only retrieval, the stored library and all historical hashes, V2 field names, indexes, and scores remain unchanged. A valid mapped query that has no mapped bond edit must **not** create an `exact_edit` candidate channel from the `NO_MAPPED_BOND_EDIT` sentinel.

## Optional external mapper (off by default)

No neural dependency is added to the `rxni` core install. To explicitly map a reaction using an external interpreter where `rxnmapper` is separately installed:

```bash
rxni-map-reaction --mapper-python /absolute/path/to/mapper/python --reaction 'CCBr.N>>CCN'
```

RXNI validates the returned mapping and checks that reactant/product molecular identities and the agent field were not changed. An administrator can opt in to mapping before retrieval by setting `RXNI_MAPPER_PYTHON=/absolute/path/to/mapper/python`. Never expose the mapper executable path as untrusted HTTP request data. The optional worker reports mapping provenance; it is not scientific proof of correct atom correspondences. RXNMapper is not bundled with RXNI and may have its own weights, dependencies, and terms. The planned Pixi `rxni` and `mapper` environments are deferred to **v0.1.2**.

## Retrieval diagnostics

A `retrieval_quality` section is exposed in Python/CLI raw result and the typed HTTP diagnostic schema. It reports query and candidate primary-event consistency, channel counts, exact-edit and broad-fallback proportions, candidate-budget saturation, and the fact that **single-budget sensitivity is not evaluated**. These are descriptive diagnostics, not ranking accuracy estimates.

## Condition-evidence policy

A method-family label is *inferred*, not an experimentally validated reagent/solvent/temperature recommendation. The historical public data does not by itself supply a sufficiently curated experimental evidence contract. `condition_evidence_assessment.recommendation_status` therefore **abstains**; `experimentally_validated_recommendations` is empty. Query-supplied solvents/reagents remain query context and must not be mistaken for external experimental precedent. No existing heuristic results are silently promoted.

## RDKit record trace

```bash
rxni-audit-rdkit --library /path/to/uspto_1976_2016_rxni_library \
  --max-rows 5000 --output /private/path/rdkit_findings.jsonl
```

This opt-in read-only report records `reaction_key` and Parquet row coordinates for parse/sanitization failures and phosphorus/isolated-H structural flags. It does **not** automatically delete records, reinterpret chemistry, or certify sanitization warnings as erroneous structures. Inspect the relevant stored source records before making exclusions.

## Audit expectations

Test strict mapping rejections through CLI, service, and HTTP (422); mapped N-alkylation retrieval parity and top-five rank/score invariance; optional mapper failure/provenance; condition abstention; budget diagnostics at two budgets; RDKit record trace; historical V2 manifest unchanged; rebuilt wheel/sdist/installed entry points; and the complete Tutorials 01–07 with refreshed documentation. Full library scientific accuracy, historical reconstruction, and multi-template one-to-many lookup behavior remain independent acceptance questions.
