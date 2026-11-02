# Tutorial 02 — Reproducible retrieval against a validated library

**Goal:** select a genuine mapped reaction from the published indexed V2 representations, then run CLI and Python retrieval. **No separately downloaded V5 Core is needed for a first-time user.** A matching V5 source Core is an alternative for developers; do not substitute a Core from another dataset generation.

## Prerequisites

Follow [Tutorial 01](01-download-verify-library.md), install RXNI and verify the library archive. Set `RXNI_LIBRARY` to your **extracted** library directory.

## 1. Select a mapped reaction from the published representations

```bash
REACTION=$(python - "$RXNI_LIBRARY" <<'PY'
import sys
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
pf = pq.ParquetFile(lib.representations)
reaction = pf.read_row_group(0, columns=['reaction'])['reaction'][0].as_py()
assert reaction.count('>') == 2
print(reaction)
PY
)
```

The source here is `reaction_representations.parquet`, resolved through `rxni-library.json`; it is **not** the source dataset's `reaction_core_v5.parquet`. The stored reaction is atom-mapped. This yields an in-corpus query with a possible self-hit.

## 2. Run indexed retrieval and persist JSON

```bash
mkdir -p "$HOME/rxni-results"
rxni retrieve \
  --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --mode sem_multi --top-k 5 --candidate-budget 500 \
  --format json --pager never \
  > "$HOME/rxni-results/retrieval.json"

python - <<'PY'
import json
from pathlib import Path
r = json.loads((Path.home()/'rxni-results/retrieval.json').read_text())
print('Primary event:', r.get('primary_event'))
print('Members:', len(r.get('members', [])))
print('Templates:', len(r.get('templates', [])))
assert r.get('members'), 'No retrieved members'
PY
```

The independent WSL audit obtained `n_alkylation`, 500 candidates and nonempty candidates and templates. Repeating it against another library does not promise those counts. A score of `1.0` for an identical stored reaction is a **self-hit**, not held-out evaluation.

## 3. Use the public Python API

```bash
python - "$RXNI_LIBRARY" <<'PY'
import sys
import pyarrow.parquet as pq
from rxni.integration import ReactionRetriever
from rxni.retrieval.libraries import resolve_library

root = sys.argv[1]
lib = resolve_library(root)
reaction = pq.ParquetFile(lib.representations).read_row_group(
    0, columns=['reaction']
)['reaction'][0].as_py()
retriever = ReactionRetriever(root)
try:
    result = retriever.retrieve(
        reaction, mode='sem_multi', top_k=5,
        max_members=500, member_limit=10,
    )
    print(result.get('primary_event'), result.get('candidate_member_count'))
    for m in result.get('members', [])[:5]:
        print(m.get('reaction_key'), m.get('template_hash'), m.get('scores', {}).get('sem_multi'))
finally:
    retriever.close()
PY
```

The independent WSL test reproduced the top five CLI reaction keys and `sem_multi` scores. Python and CLI agreement supports interface consistency, not scientific accuracy.

## 4. Alternative: a matching source Core

If you **already possess** a `reaction_core_v5.parquet` from the **same indexed population**, you may run the historic helper:

```bash
export RXNI_CORE='/absolute/path/to/matching/reaction_core_v5.parquet'
bash examples/scripts/retrieve_from_core.sh
```

The shipped helper was also tested successfully with the **published representations Parquet** supplied through `RXNI_CORE`, because it only reads an eligible `reaction` string. For clarity, use the new explicit V2 query example above instead of pretending the Zenodo archive contains a separate source Core.

## 5. Interpreting results and strict mapping

Graph-edit retrieval requires meaningful, consistent reactant/product atom mapping. The audited prepatch `0.1.1` code can accept unmapped SMILES and generate `NO_MAPPED_BOND_EDIT`; do not rely on such results. The planned **strict mapping preflight** will reject unmapped/insufficiently mapped queries. Recheck runtime behavior after the code change before publishing these docs.

Candidate-generation channels include `exact_edit`, local mapped-center radii (`raw_r0`, `raw_r1`, `raw_r2`), `component_evidence` and broad `component_topology`. A large candidate budget can introduce many different-primary-event matches; inspect class agreement instead of interpreting raw counts as relevance. Recorded condition fields and method-family labels are not validated synthetic instructions.

**For scientific benchmarking,** choose reactions absent from the retrieval library, freeze a held-out split, avoid document/transformation leakage, record source lineage and independently score chemical relevance. See [Tutorial 03](03-analysis-api.md) for HTTP and [Tutorial 04](04-template-catalogue-and-provenance.md) for template identities.
