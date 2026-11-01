"""Public RXNI distribution identity and scientific-baseline separation."""
from __future__ import annotations
import importlib.util
import ast
from pathlib import Path
import unittest


class TestPublicVersionContract(unittest.TestCase):
    def test_public_package_version(self):
        import rxni
        self.assertEqual(rxni.__version__, '0.1.1')
        self.assertEqual(rxni.SCIENTIFIC_BASELINE_VERSION, '6.7.0')
        self.assertEqual(rxni.RUNTIME_VERSION, rxni.SCIENTIFIC_BASELINE_VERSION)

    def test_semantic_identifiers_are_frozen(self):
        from rxni.semantic.version import (RUNTIME_VERSION,
            SEMANTIC_ENGINE_VERSION, SEMANTIC_RULESET_VERSION)
        self.assertEqual(RUNTIME_VERSION, '6.7.0')
        self.assertTrue(SEMANTIC_ENGINE_VERSION.startswith('unified-semantic-'))
        self.assertTrue(SEMANTIC_RULESET_VERSION.startswith('unified-rules-'))

    def test_stored_generation_not_rebranded(self):
        # Avoid eagerly importing the optional PyArrow retrieval stack just
        # to check a persisted-format constant.
        import rxni
        path = Path(rxni.__file__).parent / 'retrieval' / 'version.py'
        tree = ast.parse(path.read_text(encoding='utf-8'))
        values = {target.id: ast.literal_eval(node.value)
                  for node in tree.body if isinstance(node, ast.Assign)
                  for target in node.targets if isinstance(target, ast.Name)}
        self.assertEqual(values['STORED_REPRESENTATION_GENERATION'],
                         'retrieval_representations_v2')

    def test_renamed_cli_modules_present(self):
        for name in ('rxni.cli_analysis','rxni.cli_summary','rxni.batch_cli'):
            with self.subTest(name=name):
                self.assertIsNotNone(importlib.util.find_spec(name))
        for name in ('rxni.cli_v67','rxni.cli_v67_summary','rxni.batch_cli_v67'):
            with self.subTest(name=name):
                self.assertIsNone(importlib.util.find_spec(name))


if __name__ == '__main__':
    unittest.main()
