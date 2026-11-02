# RXNI HTTP API — quick reference

Requires an installed RXNI with API dependencies and a complete compatible local library. Obtain/verify the [Zenodo historical library](https://zenodo.org/records/23856755) following [Tutorial 01](tutorials/01-download-verify-library.md).

```bash
python -m pip install -e '.[api]'
export REACTION_INTELLIGENCE_LIBRARY="$RXNI_LIBRARY"
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

The independent fresh-machine WSL test of the `0.1.1` source returned HTTP **200** for `GET /healthz`, `GET /v1/info`, `GET /v1/schema/request` and valid `POST /v1/analyze`, and HTTP **403** for a forbidden raw-runtime result request. This is interface acceptance, **not** a comprehensive security or chemical-validation audit.

In a second terminal:

```bash
curl -sS http://127.0.0.1:8000/healthz
curl -sS http://127.0.0.1:8000/v1/info
curl -sS http://127.0.0.1:8000/v1/schema/request
```

Construct a request from the published library's **mapped** stored reaction:

```bash
python - "$RXNI_LIBRARY" <<'PY'
import json, sys
from pathlib import Path
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
rxn = pq.ParquetFile(lib.representations).read_row_group(0, columns=['reaction'])['reaction'][0].as_py()
Path('/tmp/rxni_api_request.json').write_text(json.dumps({
    'reaction':rxn, 'options':{'top_k':5,'candidate_budget':500}
}))
PY

curl -sS -o /tmp/rxni_api_response.json -w 'HTTP %{http_code}\n' \
  http://127.0.0.1:8000/v1/analyze \
  -H 'Content-Type: application/json' \
  --data-binary @/tmp/rxni_api_request.json
```

For the forbidden raw result boundary:

```bash
curl -sS -o /dev/null -w '%{http_code}\n' \
  http://127.0.0.1:8000/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{"reaction":"[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]","options":{"include_raw_runtime_result":true}}'
```

The experimental sample is syntactically mapped but not promised to match this corpus. The intended strict mapping code patch will reject missing/inconsistent atom mapping at CLI, Python and HTTP boundaries; HTTP error semantics must be tested and documented **after** that code patch, not guessed now. A method label or stored evidence field is not a validated synthetic procedure.

See [Tutorial 03](tutorials/03-analysis-api.md) for the guided workflow and [release limitations](release-limits-v0.1.md) before chemical interpretation.
