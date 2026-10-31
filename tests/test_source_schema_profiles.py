"""Source profile contract tests; do not require PyArrow for CSV cases."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from rxni.retrieval.builders.source_schema import (
    SourceSchemaError, load_profile, describe, iter_source,
)


class TestSourceProfiles(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.addCleanup(self.t.cleanup)
        self.root = Path(self.t.name)
        self.csv = self.root / "reactions.csv"
        with self.csv.open("w", newline="", encoding="utf-8") as f:
            w=csv.DictWriter(f, fieldnames=["id", "transformation", "use_row", "problem"])
            w.writeheader()
            w.writerow({"id": "case-1", "transformation": "CCO>>CC=O", "use_row":"true", "problem":""})
            w.writerow({"id": "case-2", "transformation": "CCBr.O>>CCO", "use_row":"false", "problem":""})
            w.writerow({"id": "case-3", "transformation": "CCBr.O>>CCO", "use_row":"true", "problem":"parse_error"})
            w.writerow({"id": "case-4", "transformation": "CCBr.O>>CCO", "use_row":"1", "problem":""})
        self.schema = self.root / "schema.json"
        self.raw = {"schema_version":1,"profile_id":"csv-test","source_format":"csv",
                    "columns":{"reaction_key":"id", "reaction_smiles":"transformation",
                               "valid":"use_row","feature_error":"problem"}}
        self.write_profile()

    def write_profile(self):
        self.schema.write_text(json.dumps(self.raw), encoding="utf-8")

    def test_csv_mapping_and_preserved_keys(self):
        profile=load_profile(self.schema)
        rows=list(iter_source(self.csv, profile))
        self.assertEqual([(i,k) for i,k,_ in rows], [(0,"case-1"),(1,"case-4")])
        self.assertEqual(rows[0][-1], "CCO>>CC=O")

    def test_csv_preflight(self):
        self.assertEqual(describe(self.csv, load_profile(self.schema))["status"],"SOURCE_SCHEMA_COMPATIBLE")

    def test_mismatched_format_rejected(self):
        self.raw["source_format"] = "parquet";self.write_profile()
        with self.assertRaisesRegex(SourceSchemaError,"requires .parquet"):
            list(iter_source(self.csv,load_profile(self.schema)))

    def test_missing_key_column_rejected(self):
        self.raw["columns"]["reaction_key"]="wrong";self.write_profile()
        with self.assertRaisesRegex(SourceSchemaError,"Missing source columns"):
            list(iter_source(self.csv,load_profile(self.schema)))

    def test_missing_required_mapping_rejected(self):
        del self.raw["columns"]["reaction_key"];self.write_profile()
        with self.assertRaisesRegex(SourceSchemaError,"Required column mapping"):
            load_profile(self.schema)

    def test_unknown_profile_fields_rejected(self):
        self.raw["rename_scoring"] = True;self.write_profile()
        with self.assertRaisesRegex(SourceSchemaError,"Unknown profile fields"):
            load_profile(self.schema)

    def test_bad_validity_flag_rejected(self):
        p=self.csv
        data=p.read_text().replace("true,", "maybe,",1)
        p.write_text(data)
        with self.assertRaisesRegex(SourceSchemaError,"Unknown validity flag"):
            list(iter_source(p,load_profile(self.schema)))

    def test_invalid_three_component_reaction_rejected(self):
        p=self.csv
        p.write_text(p.read_text().replace("CCO>>CC=O", "CCO"))
        with self.assertRaisesRegex(SourceSchemaError,"reactants>agents>products"):
            list(iter_source(p,load_profile(self.schema)))

    def test_parquet_both_v023a_and_legacy(self):
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError:
            self.skipTest("PyArrow unavailable")
        for reaction_column in ("reaction", "forward_reaction_smiles"):
            with self.subTest(reaction_column=reaction_column):
                p=self.root / "reactions.parquet"
                pq.write_table(pa.table({"reaction_key":["k1", "k2"],
                                         reaction_column:["CCO>>CC=O", "CCBr.O>>CCO"]}),p)
                self.raw={"schema_version":1,"profile_id":"parquet-test","source_format":"parquet",
                          "columns":{"reaction_key":"reaction_key","reaction_smiles":reaction_column}}
                self.write_profile()
                rows=list(iter_source(p,load_profile(self.schema)))
                self.assertEqual(len(rows),2)


if __name__=="__main__": unittest.main()
