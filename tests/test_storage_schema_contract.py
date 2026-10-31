"""Public RXNI V2 stored-record codec and schema regression checks."""
from __future__ import annotations

import unittest

import pyarrow as pa

from rxni.retrieval.storage.codec import (
    V2_STORAGE_COLUMNS,
    retrieval_record_to_storage_row_v2,
    retrieval_record_from_storage_row_v2,
    retrieval_storage_schema,
)


class TestStorageSchemaContract(unittest.TestCase):
    def test_explicit_codec_owned_schema(self):
        expected = (
            "record_id", "reaction", "semantic_profile_json", "semantics_json",
            "reaction_centers_json", "localized_fingerprints_json",
            "multichannel_fingerprint_json", "reaction_key",
            "representation_generation", "representation_build_id",
        )
        self.assertEqual(V2_STORAGE_COLUMNS, expected)
        self.assertEqual(retrieval_storage_schema().names, list(expected))
        self.assertTrue(all(pa.types.is_string(f.type) for f in retrieval_storage_schema()))

    def test_native_encoding_and_v2_decoding(self):
        original = {
            "record_id": "reaction-example",
            "reaction": "[CH3:1][Br:2]>>[CH3:1][OH:2]",
            "semantic_profile": {"P:substitution": 3.0},
            "semantics": {"primary_events": []},
            "reaction_centers": {},
            "localized_fingerprints": {"centers": []},
            "multichannel_fingerprint": {"channels": {}},
        }
        row = retrieval_record_to_storage_row_v2(original, build_id="test-build")
        row["reaction_key"] = "reaction-example"
        table = pa.Table.from_pylist([row], schema=retrieval_storage_schema())
        self.assertEqual(table.num_rows, 1)
        decoded = retrieval_record_from_storage_row_v2(
            table.to_pylist()[0], expected_build_id="test-build"
        )
        for key in original:
            self.assertEqual(decoded[key], original[key])
        self.assertEqual(decoded["metadata"]["representation_build_id"], "test-build")

    def test_mismatched_build_rejected(self):
        record = {
            "record_id": "x", "reaction": "[CH3:1]Br>>[CH3:1]O",
            "semantic_profile": {}, "semantics": {}, "reaction_centers": {},
            "localized_fingerprints": {"centers": []},
            "multichannel_fingerprint": {"channels": {}},
        }
        row = retrieval_record_to_storage_row_v2(record, build_id="one")
        with self.assertRaisesRegex(ValueError, "build mismatch"):
            retrieval_record_from_storage_row_v2(row, expected_build_id="two")


if __name__ == "__main__":
    unittest.main()
