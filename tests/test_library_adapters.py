"""Storage-layout unit tests; run with PYTHONPATH pointing at src/."""
import json
import tempfile
import unittest
from pathlib import Path
from rxni.retrieval.libraries import (DESCRIPTOR_FILENAME, LEGACY_V2_NAMES,
                                      LibraryFormatError, inspect_library, resolve_library)


class TestLibraryLayout(unittest.TestCase):
    def populate(self, root, names):
        for p in names.values():
            f=root/p
            f.parent.mkdir(parents=True, exist_ok=True)
            if p.endswith(".json"):
                f.write_text(json.dumps({"build_id":"test", "representation_generation":"retrieval_representations_v2"}))
            else:
                f.write_bytes(b"fixture")

    def test_legacy_v2(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            self.populate(root,LEGACY_V2_NAMES)
            a=resolve_library(root)
            self.assertEqual(a.layout,"legacy-v2-directory")
            self.assertEqual(inspect_library(root)["status"], "V2_READABLE")

    def test_descriptor_alternate_artifact_names(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            paths={key:f"data/{key}.json" if key=="retrieval_manifest" else f"data/{key}.bin" for key in LEGACY_V2_NAMES}
            self.populate(root,paths)
            desc={"schema_version":1,"representation_format":"retrieval_representations_v2",
                  "library_id":"test-library","artifacts":paths,"semantic_engine_version":"6.7"}
            (root/DESCRIPTOR_FILENAME).write_text(json.dumps(desc))
            a=resolve_library(root/DESCRIPTOR_FILENAME)
            self.assertEqual(a.layout,"descriptor-v1")
            self.assertEqual(a.library_id,"test-library")
            self.assertEqual(a.semantic_engine_claim,"6.7")

    def test_v1_only_rejected_explicitly(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            (root/"system_i_representations_v1.parquet").write_bytes(b"fixture")
            with self.assertRaisesRegex(LibraryFormatError,"V1-only"):
                resolve_library(root)

    def test_core_requires_materialization(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            (root/"reaction_core_v5.parquet").write_bytes(b"fixture")
            with self.assertRaisesRegex(LibraryFormatError,"Reaction Core"):
                resolve_library(root)

    def test_descriptor_escape_denied(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            desc={"schema_version":1,"representation_format":"retrieval_representations_v2",
                  "artifacts":{k:"../outside" for k in LEGACY_V2_NAMES}}
            (root/DESCRIPTOR_FILENAME).write_text(json.dumps(desc))
            with self.assertRaisesRegex(LibraryFormatError,"within the library"):
                resolve_library(root)

    def test_wrong_generation_denied(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            self.populate(root,LEGACY_V2_NAMES)
            (root/LEGACY_V2_NAMES["retrieval_manifest"]).write_text(json.dumps({
                "build_id":"test","representation_generation":"unrecognized"}))
            with self.assertRaisesRegex(LibraryFormatError,"does not match"):
                resolve_library(root)


if __name__ == "__main__": unittest.main()
