# Tutorial 03 — Reaction analysis with the HTTP API

**Goal:** use RXNI's locally hosted FastAPI service with a validated library and a mapped reaction. This test covers API structure and the raw-output boundary, not independent chemical prediction accuracy.

## 1. Install and configure

Install API dependencies from the RXNI checkout:

```bash
python -m pip install -e '.[api]'
export REACTION_INTELLIGENCE_LIBRARY="$RXNI_LIBRARY"
```

Follow [Tutorial 01](01-download-verify-library.md) before starting the server. `RXNI_LIBRARY` is the tutorial's shell variable; `REACTION_INTELLIGENCE_LIBRARY` is the API startup setting.

## 2. Select a mapped query and create JSON request

```bash
python - "$RXNI_LIBRARY" <<'PY'
import json, sys
from pathlib import Path
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
rxn = pq.ParquetFile(lib.representations).read_row_group(
    0, columns=['reaction']
)['reaction'][0].as_py()
request = {'reaction': rxn, 'options': {'top_k': 5, 'candidate_budget': 500}}
Path('/tmp/rxni_mapped_request.json').write_text(json.dumps(request))
print('Mapped API request ready')
PY
```

The request uses a mapped stored reaction rather than the previously documented unmapped esterification. **This is in-corpus self-query testing**, not a new-chemistry benchmark.

## 3. Start the local server

In one terminal:

```bash
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

The library must be configured **before** the first scientific request. App import is allowed before service initialization. Do not expose a local dataset/API publicly without a separate security review.

## 4. Execute the HTTP requests

In another terminal, with the environment active:

```bash
curl -sS http://127.0.0.1:8000/healthz
curl -sS http://127.0.0.1:8000/v1/info
curl -sS http://127.0.0.1:8000/v1/schema/request

curl -sS -o /tmp/rxni_api_response.json \
  -w 'analyze HTTP %{http_code}\n' \
  http://127.0.0.1:8000/v1/analyze \
  -H 'Content-Type: application/json' \
  --data-binary @/tmp/rxni_mapped_request.json
```

In the independent WSL audit, `/healthz`, `/v1/info`, `/v1/schema/request` and `/v1/analyze` all returned HTTP **200**. The analysis response included `reaction`, `structure`, `semantics`, `retrieval`, `conditions`, `diagnostics`, `provenance`, `runtime`, `schema_version`, and `summary` (plus a `chemkn` section). These are observations for the tested code, not a promise that all sections contain scientific evidence.

## 5. Confirm raw-result restriction

```bash
curl -sS -o /tmp/rxni_raw_denial.json -w 'raw-policy HTTP %{http_code}\n' \
  http://127.0.0.1:8000/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{"reaction":"[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]","options":{"include_raw_runtime_result":true}}'
```

Expected default denial: HTTP **403**. This tests an explicit access boundary only; it is not a full API penetration/security test. The standalone mapped example here is syntactically illustrative and not assumed to have matching corpus precedents.

## 6. Programmatic app factory

```python
from pathlib import Path
from rxni.api.config import ApiSettings
from rxni.api.app import create_app

app = create_app(settings=ApiSettings(library=Path('/absolute/path/to/validated/library')))
```

## 7. Mapping and provenance caution

The prepatch `0.1.1` API accepted unmapped inputs but could yield misleading `NO_MAPPED_BOND_EDIT` retrieval. The planned strict mapping preflight must reject unmapped/insufficiently mapped inputs on **all relevant entry points** (CLI, Python and HTTP). The intended HTTP validation status/error schema must be verified when code changes land; this tutorial deliberately does **not** invent an HTTP status for that not-yet-implemented behavior.

See [Tutorial 02](02-retrieval.md) for Python parity and [API quickstart](../api-quickstart.md) for endpoint reference.
