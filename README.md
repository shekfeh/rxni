# RXNI — Reaction Intelligence

**Development release candidate: 0.1.0**. RXNI provides deterministic reaction-semantic analysis, library-based structural/semantic retrieval, and exploration of recorded reaction conditions. It does **not** establish experimental reaction success or mechanistic electron-pushing predictions.

## Install from source

In a supported Python environment (Python 3.10 or later):

```bash
python -m pip install -e .
# API / legacy rich interface, only if needed:
python -m pip install -e '.[api,cli]'
```

The optional ChemKN interface accepts explicitly supplied normalized JSONL knowledge cards (`RXNI_KNOWLEDGE_CARDS` or `--knowledge-cards`) and bundles no publisher-specific adapters or datasets. No external library is downloaded with the core package. See `docs/chemkn-provider-v0.1.md`.

## Inspect the command line

```bash
rxni --version
rxni --help
rxni retrieve --help
rxni analyze --help
rxni-summary --help
rxni-batch --help
```

Running actual retrieval requires a locally prepared, validated RXNI library. Consult `docs/architecture.md`, `docs/input-schema-profiles.md`, and `docs/foundation-schema-profiles.md` for compatibility constraints. Schema profiles map supported source or SQLite layouts; they do not convert incompatible hashing algorithms or scientific representation formats.

## HTTP API

Install with `.[api]`, configure the library location, then start the ASGI app:

```bash
export REACTION_INTELLIGENCE_LIBRARY=/path/to/validated/rxni/library
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

The API configuration is resolved on the first ASGI invocation. Importing `rxni.api.config` requires neither a configured library nor FastAPI. For embedded integration, call `rxni.api.create_app(settings=...)` with explicit configuration.

## Validation and scope

The migration and integration evidence includes 384 paired reaction-case parity evaluations, schema-profile checks, and independently installed CLI/HTTP API checks. These establish migration equivalence and interface compatibility within the tested conditions, **not** independent chemical accuracy or large-library retrieval quality. Full-library integrity and separate scientific benchmarks remain release gates. See `docs/release-limits-v0.1.md`.

## License and data rights

RXNI original software is licensed under the **Apache License, Version 2.0** (SPDX `Apache-2.0`); see [LICENSE](LICENSE). This does not grant rights to third-party chemistry knowledge collections, commercial textbooks, external databases or bundled scientific resources of unverified origin. Complete the resource-by-resource provenance audit before public publication. See [licensing and data](docs/licensing-and-data.md) and [the distribution-rights checklist](docs/distribution-rights-checklist.md).

## Documentation and tutorials

Start with the [documentation index](docs/index.md), [getting started](docs/getting-started.md), [CLI reference](docs/cli-reference.md), [API reference](docs/api-quickstart.md), and [four practical tutorials](docs/index.md#first-steps). Examples assume an existing compatible reaction library and do not bundle third-party source corpora.

## Rule provenance and data boundary

RXNI's 82 directional FGI rules are project-developed heuristics, iteratively
refined against USPTO 1960–2016 reaction/template datasets. The underlying
source reaction datasets and private evaluation truth are not bundled.
See [the provenance summary](docs/provenance-summary-v0.1.md) and
[release limitations](docs/release-limits-v0.1.md).
