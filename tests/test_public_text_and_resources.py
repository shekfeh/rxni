"""Guard cleaned public-facing names without changing serialized contracts."""
from __future__ import annotations
import unittest
from pathlib import Path

class TestPublicTextAndResources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.package=Path(__file__).resolve().parents[1]/'src'/'rxni'
    def test_help_is_rxni_branded(self):
        # Inspect the parser AST without importing PyArrow/RDKit-dependent modules.
        # Real --help execution is a separate integration gate.
        import ast
        for filename in ('cli_analysis.py', 'cli_summary.py'):
            source=(self.package/filename).read_text(encoding='utf-8')
            tree=ast.parse(source, filename=filename)
            parsers=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name=='build_parser']
            self.assertEqual(len(parsers),1,filename)
            parser_source=ast.get_source_segment(source,parsers[0])
            self.assertIn('rxni',parser_source.lower())
            self.assertNotIn('reaction-intelligence',parser_source.lower())
            self.assertNotIn('Runtime-v6.7',parser_source)
    def test_package_versions_and_storage_abi(self):
        # Static contract test: avoid importing storage modules in minimal setups.
        import ast
        def literal(path, name):
            tree=ast.parse(path.read_text(encoding='utf-8'))
            assignments=[node for node in tree.body if isinstance(node,ast.Assign)
                         and any(isinstance(t,ast.Name) and t.id==name for t in node.targets)]
            self.assertEqual(len(assignments),1,name)
            return ast.literal_eval(assignments[0].value)
        self.assertEqual(literal(self.package/'__init__.py','__version__'),'0.1.1')
        self.assertEqual(literal(self.package/'semantic/version.py','RUNTIME_VERSION'),'6.7.0')
        self.assertEqual(literal(self.package/'retrieval/version.py','STORED_REPRESENTATION_GENERATION'),
                         'retrieval_representations_v2')
    def test_archived_development_resources_absent(self):
        self.assertFalse((self.package/'docs').exists())
        self.assertFalse(list(self.package.rglob('*MIGRATION_MANIFEST.json')))
        self.assertFalse((self.package/'semantic/nonfgi/detectors.py.pre_condensation_v67').exists())
    def test_scientific_contract_asset_retained(self):
        self.assertTrue((self.package/'semantic/contracts/fgi_contracts_v1.json').is_file())
    def test_public_docs_present(self):
        root=self.package.parents[1]
        self.assertTrue((root/'docs/architecture.md').is_file())
        self.assertTrue((root/'docs/version-contracts.md').is_file())

if __name__=='__main__':unittest.main()
