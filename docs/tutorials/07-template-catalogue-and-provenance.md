# Tutorial 07 — Template SMARTS, hash identities and scientific provenance

**Context:** RXNI Historical USPTO Reaction Library (v6.7 P6/Ring031 representation, internal V2 storage).

## 1. Why a sixth file matters

The RXNI retriever's five-file runtime uses structural keys, reaction metadata, precomputed stored representations and locator positions. Ranked candidates may contain `template_hash`, but retrieval itself does not load the SMARTS string. The sixth file, `template_catalog.parquet`, restores access to the exact retrospective template SMARTS associated with those hashes. It is supplementary scientific data, not an ABI requirement.

Historical audit of the *specific V5 template catalogue used by this foundation*:

| Audit | Confirmed result |
|---|---:|
| Unique template records | 42,128 |
| Nonempty SMARTS | 42,128 |
| Conflicting hash→SMARTS entries | 0 |
| Retrieval lookup template hashes covered | 42,128/42,128 |
| RDKit `ReactionFromSmarts` parsed | 42,128 |
| RDKit `Validate()` warnings | 35 |
| RDKit `Validate()` errors | 0 reported |
| rxnutils fingerprint hash replay | 42,128/42,128 |

`Validate()` warnings are a distinct measure from invalid parse; a syntactically valid SMARTS does **not** imply a chemically generalizable or applicable transformation. Test actual retrosynthetic application on the known parent reaction as a separate experimental protocol if required.

## 2. Read the template catalogue

```python
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library

root = '/data/RXNI-Historical-USPTO-r1'
lib = resolve_library(root)
catalog = pq.read_table(
    f'{root}/template_catalog.parquet',
    columns=['retro_template_hash', 'retro_template_number', 'retro_template_smarts', 'n_members']
)
print(catalog.num_rows)
row = catalog.slice(0, 1).to_pylist()[0]
print(row)
```

To join a retrieved `template_hash` to SMARTS:

```python
by_hash = {
    x['retro_template_hash']: x['retro_template_smarts']
    for x in catalog.select(['retro_template_hash', 'retro_template_smarts']).to_pylist()
}

# `result` is the dictionary returned by V66Retriever.retrieve(...)
for member in result.get('members', [])[:5]:
    key = member.get('template_hash')
    print(key, by_hash.get(key))
```

The `retro_template_hash` value is a template identity; it is not the atomic reaction key nor a scoring feature on its own. Template supports/occurrences are not guaranteed to equal independent patent counts because of duplicate-source observations.

## 3. Reconstruct the historical hash policy

The historical AiZynthTrain extraction used the result of `retro_template.hash_from_bits()`; `rxnutils` computes fingerprint-bit differences and hashes the sorted bit/count pairs. This is **not** a raw-SMARTS checksum. A full 42,128-template replay succeeded with RDKit 2023.09.6 in the existing `aizynthtrain` Mamba environment.

Reproduce with the optional helper (first export catalogue rows to JSONL if needed):

```bash
conda run --no-capture-output -n aizynthtrain python \
  scripts/release/replay_aizynth_template_hashes.py \
  --jsonl /path/to/template_hash_full_42128.jsonl \
  --report /path/to/template_hash_replay_report.json
```

The tool writes a checksum-attested audit including the reaction-utils and RDKit versions. Do not add `rxnutils` to the minimal RXNI runtime environment only to inspect stored template SMARTS; keep it separate for build/reproduction tasks.

**Do not conflate hash namespaces:** the historical native `TemplateHash` is fingerprint-based; RXNI's current strict Core builder demands `SHA256(exact SMARTS)` on its *new Core export interface*. For conversions, preserve both and re-key memberships explicitly. A generic column-renaming schema profile is insufficient.

## 4. Read scientific provenance correctly

`SCIENTIFIC_PROVENANCE.json` and `retrieval_manifest.json` have different roles:

- `retrieval_manifest.json`: historical build ID, V2 physical ABI, stored data SHA256, scientific authority statements and positional locator checksum; its locator digest must be reattested after creating a safe SQLite snapshot.
- `SCIENTIFIC_PROVENANCE.json`: original Ring031 source receipt SHA256, detector/resolver/code-source pins, known review queue, source Core SHA256 and explicit nonpromotion/rights status.
- `RELEASE_AUDIT.json`: independently evaluated physical integrity, complete locator positioning, template completeness, and staged-to-final retrieval sample parity.
- `rxni-library.json`: relocatable filenames; it intentionally does not change scientific semantics or column meanings.
- `SHA256SUMS`: exact *distributed file* checksums; verify after download.

The opaque `v66_...` build ID prefix means the existing V2 writer was reused. Scientific authority is explicit: `runtime_v6_7.semantic accepted through P6 and Ring v0.3.1 scoped patches`. The source was marked `CANDIDATE_BUILD_COMPLETE_NOT_PROMOTED`, and its scientific review queue included `0006e7559667b6841e0a7ea26957d85a1aaa6e4fcddfaa7cabda9baed0466a12`. A technical packaging pass must not silently erase those caveats.

## 5. Test a SMARTS rather than only parse it

Reaction SMARTS can be parsed with RDKit:

```python
from rdkit.Chem import rdChemReactions
rxn = rdChemReactions.ReactionFromSmarts(row['retro_template_smarts'])
if rxn is None:
    raise ValueError('Invalid reaction SMARTS')
warnings, errors = rxn.Validate(silent=True)
print('warnings:', warnings, 'errors:', errors)
```

Actual retrosynthetic applicability would additionally require running the correctly oriented template on an appropriate mapped product and evaluating mapped precursor matches. For the historical pipeline, `ReactionTemplate(smarts, direction='retro').apply(product_smiles)` uses `rdchiral` and may generate multiple or zero candidates. A validator should distinguish application failure, invalid atom mappings, omitted reagents and different but chemically equivalent reactant representations. A successful template run is not experimental yield evidence.

## 6. Avoid scientific overclaims

- `V2_READABLE` is loader compatibility, not chemical correctness.
- 12/12 Ring031 query agreement and template coverage are valuable portability evidence, not an independent generalization benchmark.
- The 48/48 historical parity cases were executed against the **older v6.6 distribution**; they do not certify the newer Ring031 scientific results.
- The 42,128/42,128 replay validates one **historical identity policy** using pinned tool versions, not that a newly extracted template has identical application semantics.
- Document dataset composition, unknown/variable patent extraction quality, duplicate provenance, mapping/valence warnings, missing conditions and intellectual-property/data-use restrictions.

See Tutorial 05 for usage and Tutorial 06 for source-to-library reconstruction.
