# RXNI CLI reference — audited 0.1.1 source candidate

The original published software tag is `v0.1.0`; the audited `0.1.1` source was installed independently from commit `5e030e388118024627f98aa4de7d58fd6ee706c8`. Use the help of your installed version as the source of exact available flags.

## Installed entry points

```bash
rxni --version
rxni --help
rxni retrieve --help
rxni analyze --help
rxni conditions --help
rxni-verify-library --help
rxni-build-library --help
rxni-summary --help
rxni-batch --help
```

The first-time-user path is [Tutorial 01](tutorials/01-download-verify-library.md), followed by [Tutorial 02](tutorials/02-retrieval.md).

## Example retrieval (mapped query)

After downloading and verifying the published dataset, export `RXNI_LIBRARY` and select a real mapped reaction:

```bash
REACTION=$(python - "$RXNI_LIBRARY" <<'PY'
import sys
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
print(pq.ParquetFile(lib.representations).read_row_group(0, columns=['reaction'])['reaction'][0].as_py())
PY
)
rxni retrieve --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --mode sem_multi --format json --top-k 5 --candidate-budget 500 --pager never
```

A reaction needs meaningful reactant/product atom mapping for graph-edit interpretation. The tested prepatch 0.1.1 code did **not** reject missing mapping and could return irrelevant `NO_MAPPED_BOND_EDIT` candidates. A strict mapping validator is a **future code-change dependency** for these docs. `--candidate-budget` is a candidate cap, not a quality guarantee; larger budgets may broaden chemically different events.

## Structured analysis

```bash
rxni analyze --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --mode sem_multi --top-k 10 --candidate-budget 1000 \
  --member-limit 10 --review events --review fgi --review axes \
  --format json --output "$HOME/rxni-analysis.json" --pager never
```

Do not assume that all `retrieve` flags are accepted by `analyze` and `conditions`; consult the installed `--help` subcommands.

## Method-family landscape

```bash
rxni conditions --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --profile-methods --top-k 5 --summary full --pager never
```

This explores method-family evidence, not validated solvent/base/temperature/yield recommendations. Inspect available candidate semantics, support and missing conditions before drawing chemical conclusions.

## Verification and builder

```bash
rxni-verify-library --library "$RXNI_LIBRARY" --smoke
rxni-build-library profiles
rxni-build-library inspect-runtime --root /absolute/path/to/rxni/checkout --runtime v6.7
```

The installed verifier was tested on the published Zenodo library. The public builder was exercised on a small aligned 16-reaction fixture, through candidate assembly and a self-retrieval. New libraries require full provenance and separate scientific certification; see [Tutorial 07](tutorials/07-build-library-from-sources.md).

## HTTP configuration

```bash
export REACTION_INTELLIGENCE_LIBRARY="$RXNI_LIBRARY"
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

See [Tutorial 03](tutorials/03-analysis-api.md). Do not expose local datasets or raw runtime results on a public server without security review.
