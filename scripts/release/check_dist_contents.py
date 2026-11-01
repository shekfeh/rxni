#!/usr/bin/env python3
"""Fail CI if release source archives omit documented tutorials or tools."""
import pathlib
import sys
import tarfile
import zipfile

MUST_HAVE = [
    'README.md', 'LICENSE', 'MANIFEST.in',
    'docs/library-distribution.md', 'docs/zenodo-dataset-deposit.md',
    'docs/tutorials/05-download-and-use-library.md',
    'docs/tutorials/06-build-library-from-sources.md',
    'docs/tutorials/07-template-catalogue-and-provenance.md',
    'scripts/release/finalize_uspto_historical.py',
    'scripts/release/prepare_zenodo_dataset.py',
    'scripts/release/verify_published_library.py',
    'examples/data/demo_reactions.smi',
    'examples/scripts/retrieve_from_core.sh',
]

def main(dist):
    d = pathlib.Path(dist)
    sdists = list(d.glob('rxni-*.tar.gz'))
    wheels = list(d.glob('rxni-*.whl'))
    if len(sdists) != 1 or len(wheels) != 1:
        raise SystemExit('Expected exactly one RXNI sdist and one wheel')
    with tarfile.open(sdists[0], 'r:gz') as f:
        files = {p.partition('/')[2] for p in f.getnames()}
    missing = sorted(set(MUST_HAVE) - files)
    if missing:
        raise SystemExit('sdist missing: ' + ', '.join(missing))
    with zipfile.ZipFile(wheels[0]) as f:
        names = f.namelist()
        required = {'rxni/verify_library.py', 'rxni/semantic/contracts/fgi_contracts_v1.json'}
        missing = required - set(names)
        eps = [x for x in names if x.endswith('.dist-info/entry_points.txt')]
        if len(eps) != 1:
            raise SystemExit('Wheel missing entry points')
        entry = f.read(eps[0]).decode()
        if 'rxni-verify-library' not in entry:
            raise SystemExit('Wheel missing installed rxni-verify-library command')
    if missing:
        raise SystemExit('wheel missing: ' + ', '.join(sorted(missing)))
    print('PASS: sdist documentation/examples/scripts and installed verifier present')

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv)>1 else 'dist')
