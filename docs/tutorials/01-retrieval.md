# Tutorial 01 — Reproducible retrieval against a validated library

**Goal:** select an actual reaction from a Parquet reaction Core, then run JSON retrieval against the matching compatible library. This avoids inventing a reaction expected to exist in your database.

## Prerequisites

- RXNI installed and `rxni --help` succeeds.
- A tested retrieval library directory, containing compatible indexes and V2 representations.
- The matching `reaction_core_v5.parquet` file. Do **not** substitute a Core from another generation.

```bash
export RXNI_LIBRARY='/absolute/path/to/validated/candidate_library'
export RXNI_CORE='/absolute/path/to/matching/reaction_core_v5.parquet'
```

## Choose a reaction and query it

The following script uses PyArrow, already a core RXNI dependency, and tolerates the historical and v0.23a reaction column names. It does not modify either file.

```bash
REACTION=$(python - <<'PY'
import os
import pyarrow.parquet as pq

pf = pq.ParquetFile(os.environ['RXNI_CORE'])
columns = set(pf.schema.names)
field = 'reaction' if 'reaction' in columns else 'forward_reaction_smiles'
if field not in columns:
    raise RuntimeError('No supported reaction SMILES column')
for batch in pf.iter_batches(batch_size=128, columns=[field]):
    for reaction in batch.column(0).to_pylist():
        if isinstance(reaction, str) and reaction.count('>') == 2:
            left, middle, right = reaction.split('>')
            if left and right:
                print(reaction)
                raise SystemExit(0)
raise RuntimeError('No eligible reaction in the source Core')
PY
)

rxni retrieve \
  --library "$RXNI_LIBRARY" \
  --reaction "$REACTION" \
  --format json \
  --top-k 5 \
  --candidate-budget 500
```

For a reusable script, run `examples/scripts/retrieve_from_core.sh` from the project root after setting `RXNI_LIBRARY` and `RXNI_CORE`.

## Interpret the result

Inspect its candidate records, structural/semantic provenance, and any condition precedents. A query already present in the library may retrieve itself; this is useful for a smoke test but **not** a held-out evaluation. Condition precedents do not guarantee reaction success or transferability.

## Next task

Change the query to a reaction absent from the candidate library and compare the highest-ranked analogues. For scientific benchmarking, freeze query selection and ground-truth relevance independently of RXNI.
