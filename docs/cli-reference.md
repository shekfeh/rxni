# RXNI CLI reference — v0.1.0

The following entry points and `retrieve` arguments were exercised by independent-install regression scripts during R4D–R4E. They are intended for a **validated local library**.

## Top-level commands

| Command | Purpose | Verification |
| --- | --- | --- |
| `rxni --help` | Main interface help | Installed wheel smoke |
| `rxni retrieve --help` | Retrieval argument discovery | Main interface subcommand |
| `rxni analyze --help` | Analysis argument discovery | Command documented in public README; use local help to confirm syntax |
| `rxni-summary --help` | Summary interface | Installed wheel smoke |
| `rxni-batch --help` | Batch interface | Installed wheel smoke |

**Do not assume an option accepted by `retrieve` is also accepted by `analyze` or `batch`.** Consult their actual help in your installed RXNI version.

## Retrieving reactions

```bash
rxni retrieve \
  --library /absolute/path/to/candidate_library \
  --reaction 'CC(=O)O.OCC>>CC(=O)OCC' \
  --format json \
  --top-k 5 \
  --candidate-budget 500
```

Options exercised by clean-install tests:

- `--library`: path to compatible retrieval library.
- `--reaction`: reaction SMILES in `reactants>agents>products` notation; `>>` indicates an empty agent section.
- `--format json`: machine-readable JSON output.
- `--top-k 5`: request up to five top results.
- `--candidate-budget 500`: configure candidate consideration budget.

The schema and field names in JSON outputs depend on the versioned scientific contract. Avoid treating a matching top reaction as proof that a suggested condition is experimentally suitable.

## Environment configuration

For API startup:

```bash
export REACTION_INTELLIGENCE_LIBRARY=/absolute/path/to/validated/library
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

The name of this environment variable is inherited from the frozen compatibility layer. It does not imply a private code dependency.

## CLI troubleshooting

- `No template library configured`: configure `REACTION_INTELLIGENCE_LIBRARY` before starting the server, or use explicit `ApiSettings`.
- Missing index / locator / metadata: rebuild or validate the complete library. A source-schema profile alone cannot repair it.
- Invalid reaction SMILES: ensure both reactant and product sections are present and that input is chemically parsable.
- No relevant match: inspect the candidate library coverage and ranking rather than assuming an algorithm failure.

For a reproducible tutorial using a selected reaction from a real Parquet fixture, see [Tutorial 01](tutorials/01-retrieval.md).
