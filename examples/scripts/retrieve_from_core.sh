#!/usr/bin/env bash
# Retrieve an actual reaction from an aligned Parquet Core, not an invented match.
# Requires RXNI_LIBRARY and RXNI_CORE to point to the SAME source population.
if [ -z "${RXNI_LIBRARY:-}" ] || [ -z "${RXNI_CORE:-}" ]; then
  printf 'Set RXNI_LIBRARY and RXNI_CORE to a matching library/core.\n' >&2
  exit 2
fi
if [ ! -f "$RXNI_CORE" ]; then
  printf 'Core Parquet not found: %s\n' "$RXNI_CORE" >&2
  exit 2
fi
REACTION="$(python - <<'PY'
import os
import pyarrow.parquet as pq
pf = pq.ParquetFile(os.environ['RXNI_CORE'])
field = next((x for x in ('reaction','forward_reaction_smiles','mapped_reaction_smiles')
              if x in pf.schema_arrow.names), None)
if not field:
    raise SystemExit('No supported reaction column')
for batch in pf.iter_batches(columns=[field], batch_size=256):
    for value in batch.column(0).to_pylist():
        if isinstance(value,str) and value.count('>') == 2:
            a,_,b = value.split('>')
            if a and b:
                print(value)
                raise SystemExit(0)
raise SystemExit('No valid reaction found in Core')
PY
)" || exit 1
rxni retrieve --library "$RXNI_LIBRARY" --reaction "$REACTION" --format json --top-k 5 --candidate-budget 500
