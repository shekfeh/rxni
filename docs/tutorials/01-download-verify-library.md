# Tutorial 01 — Download, verify and use the published RXNI library

**Audience:** first-time RXNI users. **Software:** published `v0.1.0` or audited `0.1.1` source; `rxni-verify-library` is available in tested `0.1.1`. **Dataset:** [RXNI Historical USPTO Reaction Library, Release 1](https://zenodo.org/records/23856755).

RXNI's core Python package does not include the large reaction library. The library is independently published and has **812,916** stored reaction representations and **42,128** retrospective SMARTS templates. The selected stored semantic authority is v6.7 after P6 and Ring v0.3.1. Original V2 storage, hash, scoring and semantic version axes are distinct. A working dataset is not evidence of universal chemical accuracy, synthesis feasibility or mechanistic electron pushing.

## 1. Install RXNI

Python 3.10+ is required. In a public RXNI source checkout:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -e '.[api,cli]'
rxni --version
rxni retrieve --help
```

`rxnutils` and AiZynthTrain are not required for published-library retrieval. Optional neural mappers should be installed separately from the native RXNI runtime.

## 2. Download the released data

The verified public landing page is <https://zenodo.org/records/23856755>. Consult that page for citation/DOI metadata; the URL alone is sufficient to identify this record and no unverified DOI is invented here.

```bash
mkdir -p "$HOME/rxni-data"
cd "$HOME/rxni-data"
curl -fL --retry 5 -C - \
  -o uspto_1976_2016_rxni_library.zip \
  'https://zenodo.org/records/23856755/files/uspto_1976_2016_rxni_library.zip?download=1'
```

Archive: `uspto_1976_2016_rxni_library.zip` · **1,444,488,348 bytes** · SHA256:

```text
5f31b689c942359c95e43debc0b8e6a8141fbb200d890fdea544ac3cd3a2a1ca
```

## 3. Verify before extraction; verify again afterward

```bash
cd "$HOME/rxni-data"
echo '5f31b689c942359c95e43debc0b8e6a8141fbb200d890fdea544ac3cd3a2a1ca  uspto_1976_2016_rxni_library.zip' | sha256sum -c -
unzip -t uspto_1976_2016_rxni_library.zip
unzip -q uspto_1976_2016_rxni_library.zip
export RXNI_LIBRARY="$HOME/rxni-data/uspto_1976_2016_rxni_library"
(cd "$RXNI_LIBRARY" && sha256sum -c SHA256SUMS)
```

The archive contains `rxni-library.json`, `structural_index.sqlite`, `reaction_lookup.sqlite`, `reaction_representations.parquet`, `representation_locator.sqlite`, `retrieval_manifest.json`, the optional `template_catalog.parquet`, and README, audit, checksum and scientific-provenance documents. **Do not rename internal files** without updating the descriptor. The public names are independent of the underlying V2 physical ABI.

## 4. Run the installed validator

```bash
rxni-verify-library --library "$RXNI_LIBRARY" --smoke
```

In the independent WSL2 test of `0.1.1`, the validator returned `TECHNICAL_PASS`, with 12 verified files, 812,916 reactions, 42,128 templates and a passing real-retrieval smoke. If you are using older published `v0.1.0` without the installed validator, use the checked-out `scripts/release/verify_published_library.py --library "$RXNI_LIBRARY" --smoke` instead.

**Meaning of PASS:** checksum, SQLite, V2 layout, locator/storage and one executable retrieval were validated. This is not a held-out chemical accuracy result, experimental-condition validation or rights clearance.

## 5. Retrieve an actual mapped reaction

**Atom maps are essential for graph-edit retrieval.** The independently audited prepatch `0.1.1` runtime accepted unmapped queries, which could produce `NO_MAPPED_BOND_EDIT` and irrelevant `exact_edit` candidates. A subsequent strict-mapping runtime is intended to **reject** such inputs; these docs must be released together with its validated enforcement. Neither tutorial relies on the neural mapper to run.

Select the first reaction stored in the installed representations:

```bash
REACTION=$(python - "$RXNI_LIBRARY" <<'PY'
import sys
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
lib = resolve_library(sys.argv[1])
pf = pq.ParquetFile(lib.representations)
reaction = pf.read_row_group(0, columns=['reaction'])['reaction'][0].as_py()
assert isinstance(reaction, str) and reaction.count('>') == 2
print(reaction)
PY
)
printf 'Mapped query: %s\n' "$REACTION"

rxni retrieve \
  --library "$RXNI_LIBRARY" \
  --reaction "$REACTION" \
  --mode sem_multi --top-k 5 --candidate-budget 500 \
  --show-templates --show-candidates --pager never
```

The WSL audit's first stored reaction was classified `n_alkylation` (`C_Cl_to_N`), returned 500 same-primary-event members at a 500-candidate budget, and ranked its exact in-corpus record first. Those values describe **one tested dataset/fixture**, not universal guarantees. In-corpus self-retrieval is a portability smoke, **not generalization**. Broader candidate budgets may dilute same-event agreement.

## 6. Structured JSON analysis

```bash
rxni analyze \
  --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --mode sem_multi --top-k 10 --candidate-budget 1000 \
  --member-limit 10 --review events --review fgi --review axes \
  --format json --output "$HOME/rxni-analysis.json" --pager never
```

The audited JSON contained `primary_event`, `ranking_mode`, `members`, `templates`, structural keys, diagnostics and timing details. To assess actual chemistry, inspect individual candidate transformations and query–candidate event agreement; the presence of output alone is not a relevance benchmark.

## 7. Method landscape (not a validated reaction-condition recommender)

```bash
rxni conditions \
  --library "$RXNI_LIBRARY" --reaction "$REACTION" \
  --profile-methods --top-k 5 --summary full --pager never
```

For the audited in-library N-alkylation the method landscape included `substitution.base_assisted`; it did **not** justify a specific experimental solvent, base, temperature, success probability or yield. A mapped but unsupported out-of-corpus query can still fall back to broad `component_topology` retrieval with many unrelated candidates. **Do not treat a method-family score or aggregate support count as validated experimental conditions.** RDKit warnings observed for some historical structures are retained in logs rather than silently treated as chemistry evidence.

## 8. Use the Python interface

```python
import os
import pyarrow.parquet as pq
from rxni.retrieval.libraries import resolve_library
from rxni.integration import ReactionRetriever

library = os.environ['RXNI_LIBRARY']
lib = resolve_library(library)
reaction = pq.ParquetFile(lib.representations).read_row_group(
    0, columns=['reaction']
)['reaction'][0].as_py()
retriever = ReactionRetriever(library)
try:
    result = retriever.retrieve(
        reaction, mode='sem_multi', top_k=5,
        max_members=500, member_limit=10,
    )
    for member in result.get('members', [])[:5]:
        print(member.get('reaction_key'), member.get('template_hash'), member.get('scores'))
finally:
    retriever.close()
```

This standalone Python snippet reads its own mapped query. See [Tutorial 02](02-retrieval.md) for CLI/Python parity and [Tutorial 03](03-analysis-api.md) for HTTP.

## 9. Inspect template SMARTS and science provenance

The `template_catalog.parquet` preserves the historical AiZynthTrain/rxnutils `ReactionTemplate.hash_from_bits()` identities, **not SHA256 of literal SMARTS**. A new Core uses `SHA256(exact SMARTS)` and must explicitly re-key memberships to that new namespace. See [Tutorial 04](04-template-catalogue-and-provenance.md).

Read `SCIENTIFIC_PROVENANCE.json` and `RELEASE_AUDIT.json`: the original P6/Ring031 scientific candidate was **not automatically promoted** by its technical release; a review case persists. A 12-case portability/live-profile check is not an independent whole-library chemical validation. Upstream Lowe USPTO CC0 and processed dataset CC BY 4.0 at Zenodo are distinct from the RXNI software's Apache-2.0 license and do not eliminate the need for rights attribution/review.

## 10. Troubleshooting

| Symptom | Check |
|---|---|
| `rxni` or validator command missing | Confirm `which python`, `which rxni`, package version; `rxni-verify-library` needs installed 0.1.1 code. |
| Download checksum failure | Redownload; never ignore an incorrect SHA256. |
| `NOT_RETRIEVAL_READY` | Check `rxni-library.json` and the five referenced files. |
| Unmapped/incomplete reaction rejected | Provide authenticated atom mapping; do not infer edit semantics from unmapped SMILES. |
| No exact/local-center hits | The query may be out of coverage; inspect fallback channels and reaction-event agreement. |
| RDKit valence/hydrogen warnings | Record and trace affected structures; nonfatal warnings do not establish chemical accuracy. |
| Template hash has no visible SMARTS | Join it to `template_catalog.parquet`; do not use a raw-SMARTS SHA256 for historical hashes. |

**Next:** [Tutorial 02 — reproducible retrieval](02-retrieval.md).
