"""Regression checks for canonical/legacy Retrieval-V2 storage profiles."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import pyarrow as pa
import pyarrow.parquet as pq

from rxni.retrieval.storage.codec import V2_STORAGE_COLUMNS, retrieval_storage_schema
from rxni.retrieval.storage.schema_profiles import (
    CANONICAL_SCHEMA_PROFILE, LEGACY_SCHEMA_PROFILE,
    CANONICAL_COLUMNS, LEGACY_COLUMNS, LEGACY_ONLY_COLUMNS,
    classify_physical_schema, validate_parquet, sha256,
)


class TestV2SchemaProfiles(unittest.TestCase):
    def _write(self, root, columns, extra=None, tagged=True):
        row = {name: None for name, _ in columns}
        row.update({"record_id":"r1", "reaction_key":"r1", "reaction":"C>>CO",
                    "semantic_profile_json":"{}", "semantics_json":"{}",
                    "reaction_centers_json":"{}", "localized_fingerprints_json":"{}",
                    "multichannel_fingerprint_json":"{}",
                    "representation_generation":"retrieval_representations_v2",
                    "representation_build_id":"b1"})
        row.update(extra or {})
        schema=pa.schema([pa.field(name,pa.uint64() if typ=="uint64" else pa.string()) for name,typ in columns])
        f=root/"retrieval_representations_v2.parquet"
        pq.write_table(pa.Table.from_pylist([row],schema=schema),f)
        manifest={"representations_sha256":sha256(f),"rows":1,"row_groups":1,
                  "build_id":"b1","representation_generation":"retrieval_representations_v2"}
        if tagged:manifest["physical_schema_profile"]=classify_physical_schema(schema)
        m=root/"retrieval_v2_manifest.json"
        m.write_text(json.dumps(manifest))
        return f,m

    def test_existing_compact_codec_remains_unchanged(self):
        self.assertEqual(tuple(V2_STORAGE_COLUMNS),tuple(name for name,_ in CANONICAL_COLUMNS))
        self.assertEqual(classify_physical_schema(retrieval_storage_schema()),CANONICAL_SCHEMA_PROFILE)

    def test_legacy_null_compatibility(self):
        with tempfile.TemporaryDirectory() as temp:
            p,m=self._write(Path(temp),LEGACY_COLUMNS)
            report=validate_parquet(p,m)
            self.assertEqual(report["profile"],LEGACY_SCHEMA_PROFILE)
            self.assertFalse(report["canonical_for_new_writes"])
            self.assertTrue(all(v==0 for v in report["legacy_only_nonnull_counts"].values()))

    def test_public_compact_schema(self):
        with tempfile.TemporaryDirectory() as temp:
            p,m=self._write(Path(temp),CANONICAL_COLUMNS)
            report=validate_parquet(p,m)
            self.assertTrue(report["canonical_for_new_writes"])

    def test_nonnull_legacy_fields_fail_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            p,m=self._write(Path(temp),LEGACY_COLUMNS,{"library_schema_version":"v5"})
            with self.assertRaisesRegex(ValueError,"Non-null legacy-only"):
                validate_parquet(p,m)

    def test_unknown_columns_and_wrong_order_rejected(self):
        with self.assertRaisesRegex(ValueError,"Unsupported"):
            classify_physical_schema(pa.schema([pa.field("reaction",pa.string())]))
        with self.assertRaisesRegex(ValueError,"Unsupported"):
            classify_physical_schema(pa.schema([pa.field(k,pa.string()) for k in reversed(V2_STORAGE_COLUMNS)]))

    def test_manifest_declared_profile_guard(self):
        with tempfile.TemporaryDirectory() as temp:
            p,m=self._write(Path(temp),CANONICAL_COLUMNS)
            manifest=json.loads(m.read_text());manifest["physical_schema_profile"]=LEGACY_SCHEMA_PROFILE
            m.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"Manifest declares"):
                validate_parquet(p,m)


if __name__ == "__main__":
    unittest.main()
