# RXNI examples

Start with [`../docs/getting-started.md`](../docs/getting-started.md), then the four tutorials in `../docs/tutorials/`.

`data/demo_reactions.smi` contains original toy reaction strings, **not** a reaction library or scientific benchmarking gold set.

`scripts/retrieve_from_core.sh` selects a usable reaction from a supplied Parquet Core and queries its matching indexed library. Set `RXNI_LIBRARY` and `RXNI_CORE` before execution. The script is not intended to be run as part of an automated scientific-quality benchmark.

The existing `profiles/` directory provides four JSON profiles for supported source/foundation schema forms. These must be used with compatible source and index contracts.
