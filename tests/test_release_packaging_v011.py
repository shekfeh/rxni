"""Public release/packaging sanity checks; no USPTO corpus is required."""
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
from rxni import __version__
from rxni.integration import ReactionRetriever
from rxni.integration.retriever import V66Retriever
from rxni.retrieval import StoredRecordProvider, V2StoredRecordProvider
from rxni.verify_library import verify_checksums

ROOT = Path(__file__).resolve().parents[1]


def test_version_and_aliases():
    assert __version__ == '0.1.1'
    assert ReactionRetriever is V66Retriever
    assert StoredRecordProvider is V2StoredRecordProvider


def test_examples_are_present():
    assert (ROOT / 'examples/data/demo_reactions.smi').is_file()
    assert (ROOT / 'examples/scripts/retrieve_from_core.sh').is_file()


def test_integrity_checker_rejects_corrupted_file():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        (root / 'one.txt').write_text('one')
        h = hashlib.sha256(b'one').hexdigest()
        (root / 'SHA256SUMS').write_text(f'{h}  one.txt\n')
        assert verify_checksums(root) == 1
        (root / 'one.txt').write_text('tampered')
        with pytest.raises(ValueError, match='SHA256 mismatch'):
            verify_checksums(root)


def test_checker_rejects_path_traversal():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        (root / 'SHA256SUMS').write_text('0' * 64 + '  ../outside\n')
        with pytest.raises(ValueError, match='Unsafe|Malformed'):
            verify_checksums(root)


def test_verifier_entrypoint_help():
    p = subprocess.run([sys.executable, '-m', 'rxni.verify_library', '--help'],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert '--smoke' in p.stdout
