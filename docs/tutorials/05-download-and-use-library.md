# Tutorial 05 — Download, verify and use an RXNI reaction library

**Audience:** new RXNI user · **Software:** RXNI 0.1.0 · **Dataset:** RXNI Historical USPTO Reaction Library, Release 1 (pending rights clearance and public URL).

This tutorial describes the *version-independent distribution layout*. The physical representation format remains `retrieval_representations_v2`; names of public files are not storage-format versions.

## 1. What the library does

A compatible library allows RXNI to find chemically related recorded reactions using graph-derived candidate keys, semantic reaction profiles and precomputed fingerprints. The historical dataset contains **812,916** stored reaction representations and **42,128** unique retrospective reaction-template SMARTS. The scientific semantic authority for the selected dataset is **v6.7, after P6 and Ring v0.3.1**. Structural retrieval preserves older v6.5-compatible hash definitions and v6.6-compatible storage/scoring; these are independent version axes.

RXNI retrieval is **not** an experimental outcome predictor, a guarantee of successful synthesis or an explicit electron-pushing mechanism generator. Some recorded structures have RDKit warnings, and a small chemical review queue remains.

## 2. Install RXNI

Requirements: Python 3.10+, sufficient disk for the extracted dataset (~2.4 GB, check actual release size), RAM suitable for the selected candidate budget, and network only to obtain the package/data. A fresh environment is preferred.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
# In a cloned RXNI source checkout:
python -m pip install -e ".[api,cli]"
rxni --version
rxni retrieve --help
```

The current software repository is <https://github.com/shekfeh/rxni>. Source installation does not automatically download reaction libraries. `rxnutils` and AiZynthTrain are **not required** to run the read-only retrieval engine.

## 3. Obtain the dataset

**Release download URL/DOI: PENDING PUBLICATION.** Do not invent a release URL until the rights audit and upload are complete. Once a URL has been published, download the dataset archive and extract to a stable data directory, *outside* your Git checkout.

Expected directory:

```text
RXNI-Historical-USPTO-r1/
    rxni-library.json
    structural_index.sqlite
    reaction_lookup.sqlite
    reaction_representations.parquet
    representation_locator.sqlite
    retrieval_manifest.json
    template_catalog.parquet          # optional to the 5-file runtime, included with this dataset
    SCIENTIFIC_PROVENANCE.json
    RELEASE_AUDIT.json
    README.md
    SHA256SUMS
```

Do not rename individual files after extraction unless you also update `rxni-library.json`. The descriptor selects the five actual runtime paths. The template catalogue is auxiliary, not a sixth mandatory input to the stored-record reader.

## 4. Verify download integrity and the installed reader

```bash
export RXNI_LIBRARY="$HOME/data/RXNI-Historical-USPTO-r1"
cd "$RXNI_LIBRARY"
sha256sum -c SHA256SUMS
```

The public read-only verification helper included in this repository performs the same checks, validates the manifest and SQLite files, and optionally executes one retrieval:

```bash
python /path/to/rxni/scripts/release/verify_published_library.py \
  --library "$RXNI_LIBRARY" --smoke
```

Expected fields from Python:

```python
from rxni.retrieval.libraries import inspect_library
print(inspect_library('/path/to/RXNI-Historical-USPTO-r1'))
# status: V2_READABLE
# layout: descriptor-v1
```

`V2_READABLE` establishes that the loader can resolve the V2 ABI. It **does not** establish the chemical correctness of every reaction.

## 5. Run a structurally supported retrieval query

The reaction must be represented as `reactants>agents>products`. Atom-mapped reactions generally provide the strongest graph-edit evidence; an unmapped illustrative input may return fewer or less useful matches.

```bash
rxni retrieve \
  --library "$RXNI_LIBRARY" \
  --reaction 'CC(=O)O.OCC>>CC(=O)OCC' \
  --mode sem_multi \
  --top-k 5 \
  --candidate-budget 500 \
  --show-templates \
  --show-candidates \
  --pager never
```

The example reaction is illustrative and **is not a promised match**. For a guaranteed corpus-derived query, select a stored reaction:

```bash
python - "$RXNI_LIBRARY" <<'PY'
import sys
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
pf = pq.ParquetFile(lib.representations)
print(pf.read_row_group(0, columns=['reaction'])['reaction'][0].as_py())
PY
```

Copy the printed reaction into `--reaction` (quote it carefully). Be aware that a corpus reaction used as its own query can retrieve itself; **self-hit is not an independent generalization benchmark**.

## 6. Interpret retrieval results

A result may contain:

- **Structural generation channels:** exact graph-edit and local-center/topology candidate sets.
- **Query semantic profile:** primary/specific events, graph features, FGI and non-FGI hypotheses.
- **Candidate members:** mapped reaction strings, reaction keys, rank-relevant scores and template hashes.
- **Template aggregates:** retrieved templates and support counts.
- **Method/condition landscape:** heuristic or evidence-backed summaries, *not a validated synthetic procedure*.

For an editable JSON output, use:

```bash
rxni analyze \
  --library "$RXNI_LIBRARY" \
  --reaction 'CC(=O)O.OCC>>CC(=O)OCC' \
  --mode sem_multi --top-k 10 \
  --candidate-budget 1000 \
  --member-limit 10 \
  --review events --review fgi --review axes \
  --format json \
  --output "$HOME/rxni-example-result.json" \
  --pager never
```

`--output` persists the complete JSON result. See `rxni analyze --help` for the installed release's exact options.

## 7. Condition-aware queries and method landscape

Supply known context explicitly; do not infer reagent/solvent identity solely from reactant/product chemistry:

```bash
rxni conditions \
  --library "$RXNI_LIBRARY" \
  --reaction 'CCBr.N>>CCN' \
  --solvent DMF --base K2CO3 \
  --temperature-c 60 \
  --profile-methods \
  --top-k 5 \
  --summary full \
  --pager never
```

This illustrates argument syntax, not recommended conditions for that reaction. `--profile-methods` is optional and contextual. The historical library does not guarantee normalized experimental yield, duration, solvent or full paragraph-level procedure evidence for every record. Chemical knowledge cards require separately supplied, appropriately licensed data.

## 8. Python API

```python
from rxni.integration.retriever import V66Retriever

library = '/absolute/path/to/RXNI-Historical-USPTO-r1'
retriever = V66Retriever(library)
try:
    result = retriever.retrieve(
        'CC(=O)O.OCC>>CC(=O)OCC',
        mode='sem_multi', top_k=5, max_members=500, member_limit=10
    )
    for member in result.get('members', [])[:5]:
        print(member.get('reaction_key'), member.get('template_hash'), member.get('scores'))
finally:
    retriever.close()
```

`V66Retriever` is the current **Python class name**, retained for the frozen scoring implementation. It can read the new v6.7/Ring031 stored profiles and does not imply that all query semantics are v6.6.

### Optional local HTTP API

```bash
export REACTION_INTELLIGENCE_LIBRARY="$RXNI_LIBRARY"
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

For request payloads, consult `docs/api-quickstart.md`; endpoint contracts should not be guessed from the CLI.

## 9. Examine the included SMARTS catalogue

The optional `template_catalog.parquet` has 42,128 rows and includes `retro_template_hash`, `retro_template_number` and `retro_template_smarts`. Example:

```python
import pyarrow.parquet as pq
from rdkit.Chem import rdChemReactions

catalog = pq.read_table('/path/to/RXNI-Historical-USPTO-r1/template_catalog.parquet',
                        columns=['retro_template_hash', 'retro_template_smarts'])
row = catalog.slice(0, 1).to_pylist()[0]
print('ID:', row['retro_template_hash'])
print('Retrosynthetic SMARTS:', row['retro_template_smarts'])
assert rdChemReactions.ReactionFromSmarts(row['retro_template_smarts']) is not None
```

The template identifiers use **AiZynthTrain/rxnutils `hash_from_bits()`**, not SHA256 of the literal SMARTS string. Their complete replay (42,128/42,128) was validated in the separate AiZynthTrain environment. Do **not** recompute these IDs using a raw-SMARTS hash when joining the metadata lookup to the catalogue.

## 10. Troubleshooting

| Symptom | Likely cause / action |
|---|---|
| `ModuleNotFoundError: rxni` | Install RXNI into the Python interpreter actually running the command. |
| `NOT_RETRIEVAL_READY` | Missing files or descriptor paths; inspect `rxni-library.json`. |
| Build-ID mismatch | Mixed representation file, manifest and locator from different builds; never mix them. |
| No candidates | Input mapping, corpus coverage, budget or retrieval mode; no match is not proof of no possible reaction. |
| RDKit hydrogen/valence warnings | Some historical patent structures challenge RDKit normalization; inspect the specific reaction. |
| Templates show hashes but not SMARTS | The five-file runtime exposes IDs; inspect the separate `template_catalog.parquet` to retrieve SMARTS. |

## 11. Scientific provenance and limitations

The selected Ring031 representation build records 812,916 reaction entries and post-P6/Ring031 semantic provenance. A known scientific review case remains; the accepted staged build was previously marked `CANDIDATE_BUILD_COMPLETE_NOT_PROMOTED`. The 12/12 recorded portability/live-profile test demonstrates agreement on 12 samples; it is not whole-corpus scientific accuracy. Read `SCIENTIFIC_PROVENANCE.json`, `RELEASE_AUDIT.json`, `docs/release-limits-v0.1.md` and the accompanying publication notes before making strong scientific claims.

**Data rights are separate from software rights.** RXNI's Apache-2.0 software license does not automatically grant permission to redistribute this historical USPTO-derived reaction corpus and derived template table. Use public distribution only after rights verification.
