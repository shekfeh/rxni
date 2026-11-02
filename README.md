# RXNI — Reaction Intelligence

**Software release status:** `v0.1.0` is the published release. This documentation refresh targets the public `0.1.1` maintenance code at commit `5e030e388118024627f98aa4de7d58fd6ee706c8`; a `v0.1.1` Git tag and release assets were not verified as part of this documentation patch. **Scientific baseline:** Reaction Intelligence v6.7.0. RXNI supports reaction-semantic analysis, indexed structural/semantic retrieval, and exploration of recorded reaction methods. It does **not** establish experimental success, validated synthetic conditions, or explicit electron-pushing mechanisms.

## Install the software

Python 3.10 or newer is required. For a cloned public RXNI checkout:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -e '.[api,cli]'
rxni --version
rxni --help
rxni-verify-library --help
```

The optional ChemKN interface accepts separately supplied normalized JSONL knowledge cards (`RXNI_KNOWLEDGE_CARDS` / `--knowledge-cards`); no publisher-specific database or historical reaction library is bundled in the software installation. See [ChemKN provider](docs/chemkn-provider-v0.1.md).

## Download the published historical reaction library

**[RXNI Historical USPTO Reaction Library — Zenodo record 23856755](https://zenodo.org/records/23856755)**. The publicly accessible archive is independent of the GitHub source checkout.

```bash
mkdir -p "$HOME/rxni-data"
cd "$HOME/rxni-data"
curl -fL --retry 5 -C - \
  -o uspto_1976_2016_rxni_library.zip \
  'https://zenodo.org/records/23856755/files/uspto_1976_2016_rxni_library.zip?download=1'
echo '5f31b689c942359c95e43debc0b8e6a8141fbb200d890fdea544ac3cd3a2a1ca  uspto_1976_2016_rxni_library.zip' | sha256sum -c -
unzip -t uspto_1976_2016_rxni_library.zip
unzip -q uspto_1976_2016_rxni_library.zip
export RXNI_LIBRARY="$HOME/rxni-data/uspto_1976_2016_rxni_library"
(cd "$RXNI_LIBRARY" && sha256sum -c SHA256SUMS)
rxni-verify-library --library "$RXNI_LIBRARY" --smoke
```

The archive is **1,444,488,348 bytes** with the SHA256 shown above. The released library has **812,916** stored representations and **42,128** template SMARTS. `rxni-verify-library` is available in the tested `0.1.1` installed code; users of the earlier `v0.1.0` release may instead use the checked-out `scripts/release/verify_published_library.py` verification helper. The standalone verification step checks technical integrity and a retrieval smoke, **not** chemical accuracy or independent redistribution-rights clearance. The actual Zenodo record should be consulted for its canonical DOI and current citation metadata; none is guessed here.

## First retrieval — mapped reactions required for scientific interpretation

The audited v0.1.1 runtime still accepts unmapped input, but it can produce the non-informative `NO_MAPPED_BOND_EDIT` signature and unrelated apparent `exact_edit` hits. **Do not use unmapped inputs.** These revised tutorials are prepared for a forthcoming **strict mapping preflight** that will reject unmapped or insufficiently mapped reaction queries; **the runtime code change must land and pass tests before these docs are published as describing enforced behavior**.

Get a real mapped query from the published library:

```bash
REACTION=$(python - "$RXNI_LIBRARY" <<'PY'
import sys
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
pf = pq.ParquetFile(lib.representations)
print(pf.read_row_group(0, columns=['reaction'])['reaction'][0].as_py())
PY
)
rxni retrieve --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --mode sem_multi --top-k 5 --candidate-budget 500 --pager never
```

This is an **in-corpus self-retrieval test**, not an independent retrieval-accuracy benchmark. Consult [Getting started](docs/getting-started.md) and the [tutorials, in learning order](docs/index.md).

## HTTP API

```bash
export REACTION_INTELLIGENCE_LIBRARY="$RXNI_LIBRARY"
uvicorn rxni.api.app:app --host 127.0.0.1 --port 8000
```

See [Tutorial 03](docs/tutorials/03-analysis-api.md) for health, schema, analysis and restricted raw-result tests. Use `rxni retrieve --help`, `rxni analyze --help`, `rxni conditions --help`, `rxni-build-library --help` for the CLI version installed locally.

## Public library builder and version boundaries

The installed public builder can create a **new** Core and V2 auxiliary from a properly aligned, licensed source export. This was independently exercised on a **16-reaction** fixture through Core validation, native producer execution, V2 assembly and self-retrieval. It does **not** recreate the historical P6/Ring031 dataset bit-for-bit or certify chemical accuracy. See [Tutorial 07: library construction](docs/tutorials/07-build-library-from-sources.md).

New RXNI Core template IDs use `SHA256(exact retrosynthetic SMARTS)` while the published historical AiZynthTrain catalogue retains `ReactionTemplate.hash_from_bits()` IDs. The namespaces are **not interchangeable**. The historical retrieval manifest retains its original opaque V2 build ID prefix (`v66_`) while declaring v6.7 semantic authority; this is not a request to modify the frozen storage ABI.

## Scientific limitations, data rights and provenance

Historical scientific provenance retains a known review case and an unpromoted source-candidate status (`CANDIDATE_BUILD_COMPLETE_NOT_PROMOTED`). Technical portability and the Zenodo publication do **not** override these scientific caveats. Read the archive's `SCIENTIFIC_PROVENANCE.json` and `RELEASE_AUDIT.json`, [release limits](docs/release-limits-v0.1.md) and [distribution status](docs/library-distribution.md).

RXNI original software is **Apache-2.0** ([LICENSE](LICENSE)); the processed historical dataset's stated license is **CC BY 4.0**, with the identified Daniel Lowe USPTO upstream source marked **CC0** at Figshare. These declarations do not themselves resolve every possible upstream contribution or independent rights review. See [licensing/data](docs/licensing-and-data.md) and the [rights checklist](docs/distribution-rights-checklist.md).

The **82 directional FGI rules** are project-authored heuristics developed using USPTO 1960–2016 reaction/template data; underlying source datasets and private evaluation truth are not included with the software. See [provenance summary](docs/provenance-summary-v0.1.md).
