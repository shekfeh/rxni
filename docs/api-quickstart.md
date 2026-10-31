# RXNI HTTP API — quick reference

Install `rxni[api]`, configure a validated local library, and run the ASGI entry point:

```bash
export REACTION_INTELLIGENCE_LIBRARY=/absolute/path/to/validated/library
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

A completed R4D–R4E test ran the same app through FastAPI's `TestClient` and confirmed the endpoints below for two distinct test libraries. Actual remote-server requests should be re-tested against your installed environment.

| Endpoint | Method | Expected successful code |
| --- | --- | --- |
| `/healthz` | GET | 200 |
| `/v1/info` | GET | 200 |
| `/v1/schema/request` | GET | 200 |
| `/v1/analyze` | POST | 200 for valid supported input |

```bash
curl --fail-with-body http://127.0.0.1:8000/healthz
curl --fail-with-body http://127.0.0.1:8000/v1/info
curl --fail-with-body http://127.0.0.1:8000/v1/schema/request

curl --fail-with-body -sS http://127.0.0.1:8000/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{"reaction":"CC(=O)O.OCC>>CC(=O)OCC","options":{"top_k":5,"candidate_budget":500}}'
```

The requests above use a simple illustrative reaction. Scientific usefulness and retrieval matches are library-dependent. `curl --fail-with-body` is available on recent curl releases; omit that flag if unsupported.

A request for `include_raw_runtime_result=true` is rejected by the validated default API policy (HTTP 403). Do not expose local files, sensitive libraries or the raw runtime endpoint publicly without an independent security review.

See [Tutorial 02](tutorials/02-analysis-api.md) for a complete guided API example.
