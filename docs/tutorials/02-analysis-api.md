# Tutorial 02 — Reaction analysis with the HTTP API

**Goal:** obtain a structured reaction analysis response using RXNI's validated FastAPI interface.

## Step 1 — install API extras

```bash
python -m pip install '.[api]'
```

## Step 2 — start the server

```bash
export REACTION_INTELLIGENCE_LIBRARY='/absolute/path/to/validated/library'
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

The library path must be configured **before** the first server request. The ASGI module can be imported without starting the scientific service; initialization is deferred.

## Step 3 — confirm service and schema

In another terminal:

```bash
curl -sS http://127.0.0.1:8000/healthz
curl -sS http://127.0.0.1:8000/v1/info
curl -sS http://127.0.0.1:8000/v1/schema/request
```

## Step 4 — analyze a reaction

```bash
curl -sS http://127.0.0.1:8000/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{
    "reaction": "CC(=O)O.OCC>>CC(=O)OCC",
    "options": {"top_k": 5, "candidate_budget": 500}
  }'
```

A successful analysis may contain sections including `reaction`, `structure`, `semantics`, `retrieval`, `conditions`, `diagnostics`, and `provenance`. The scientific interpretation depends on available library records and the runtime's own limitations.

## Step 5 — inspect the raw-result policy

```bash
curl -sS -o /dev/null -w '%{http_code}\n' \
  http://127.0.0.1:8000/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{"reaction":"CC(=O)O.OCC>>CC(=O)OCC","options":{"include_raw_runtime_result":true}}'
```

The R3C and R4D–R4E validation suites observed HTTP `403` for this unauthorized option. This is a tested policy, not a substitute for a comprehensive API security audit.

## Programmatic server creation

The validated factory uses explicit settings:

```python
from pathlib import Path
from rxni.api.config import ApiSettings
from rxni.api.app import create_app

app = create_app(settings=ApiSettings(library=Path('/absolute/path/to/validated/library')))
```

For basic CLI retrieval without a server, see [Tutorial 01](01-retrieval.md).
