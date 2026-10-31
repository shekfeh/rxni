"""Custom schema profile and V2 stored provider regression tests."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from rxni.retrieval.storage.codec import retrieval_record_to_storage_row_v2
from rxni.retrieval.storage.custom_profiles import (
    CANONICAL_NAMES, load_custom_profile, validate_schema,
    normalize_row, validate_custom_parquet,
)
from rxni.retrieval.storage.provider import V2StoredRecordProvider


class TestCustomV2Profiles(unittest.TestCase):
    def make_fixture(self, root: Path):
        mapping = {name: name for name in CANONICAL_NAMES}
        mapping['reaction'] = 'mapped_reaction'
        custom = {'schema_profile_version': 1, 'profile_id':'rxni_custom_fixture_v1',
                  'columns':mapping,
                  'extra_columns': {'source_document_id': {'type':'string','role':'provenance'}}}
        config = root/'custom_schema.json'
        config.write_text(json.dumps(custom))
        data = {'record_id':'r1', 'reaction':'[CH3:1]Br>>[CH3:1]O',
                'semantic_profile': {'P:substitution': 1.0},
                'semantics': {'primary_events': []}, 'reaction_centers':{},
                'localized_fingerprints':{'centers':[]},
                'multichannel_fingerprint': {'channels':{}}}
        encoded = retrieval_record_to_storage_row_v2(data, build_id='b1')
        encoded['reaction_key'] = 'r1'
        encoded['mapped_reaction'] = encoded.pop('reaction')
        encoded['source_document_id'] = 'doc1'
        names = ['source_document_id', 'representation_build_id', 'reaction_key',
                 'mapped_reaction'] + [n for n in CANONICAL_NAMES if n not in
                                       ('representation_build_id','reaction_key','reaction')]
        schema = pa.schema([pa.field(name, pa.string()) for name in names])
        f = root/'retrieval_representations_v2.parquet'
        pq.write_table(pa.Table.from_pylist([encoded],schema=schema),f)
        hf = hashlib.sha256(f.read_bytes()).hexdigest()
        cp = hashlib.sha256(config.read_bytes()).hexdigest()
        manifest = root/'retrieval_v2_manifest.json'
        manifest.write_text(json.dumps({
            'build_id':'b1', 'representation_generation':'retrieval_representations_v2',
            'representations_sha256':hf, 'rows':1,'row_groups':1,
            'physical_schema_profile':'rxni_custom_fixture_v1',
            'custom_schema_profile_sha256':cp}))
        locator = root/'representation_locator_v2.sqlite'
        with sqlite3.connect(locator) as conn:
            conn.execute('''CREATE TABLE representation_locator_v2
                         (reaction_key TEXT, row_group_id INTEGER, row_in_group INTEGER,
                          build_id TEXT, representation_generation TEXT)''')
            conn.execute('INSERT INTO representation_locator_v2 VALUES(?,?,?,?,?)',
                         ('r1',0,0,'b1','retrieval_representations_v2'))
        return f,config,manifest,locator

    def test_custom_file_and_provider(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            f,config,mf,locator=self.make_fixture(root)
            summary=validate_custom_parquet(f,config,mf)
            self.assertEqual(summary['status'],'CUSTOM_SCHEMA_READ_VALIDATED')
            self.assertEqual(summary['rows'],1)
            provider=V2StoredRecordProvider(f,locator,mf,schema_profile=config)
            got=provider.get(record_id='r1')
            self.assertEqual(got['semantic_profile'],{'P:substitution':1.0})
            self.assertEqual(got['reaction'],'[CH3:1]Br>>[CH3:1]O')
            self.assertEqual(provider.snapshot()['physical_schema_profile'],'rxni_custom_fixture_v1')

    def test_unknown_column_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            f,config,mf,locator=self.make_fixture(root)
            profile=load_custom_profile(config)
            original=pq.read_table(f)
            wrong=original.append_column('new_unapproved',pa.array(['x']))
            with self.assertRaisesRegex(ValueError,'unexpected'):
                validate_schema(profile,wrong.schema)

    def test_ambiguous_mapping_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            _,config,_,_=self.make_fixture(root)
            raw=json.loads(config.read_text())
            raw['columns']['record_id']=raw['columns']['reaction_key']
            config.write_text(json.dumps(raw))
            with self.assertRaisesRegex(ValueError,'one-to-one'):
                load_custom_profile(config)

    def test_manifest_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            f,config,mf,locator=self.make_fixture(root)
            raw=json.loads(mf.read_text())
            raw['custom_schema_profile_sha256']='0'*64
            mf.write_text(json.dumps(raw))
            with self.assertRaisesRegex(ValueError,'SHA256'):
                V2StoredRecordProvider(f,locator,mf,schema_profile=config)

    def test_schema_type_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            _,config,_,_=self.make_fixture(root)
            profile=load_custom_profile(config)
            fields=[pa.field(name,pa.int64() if name=='mapped_reaction' else pa.string())
                    for name in set(profile.columns.values())|set(profile.extra_columns)]
            with self.assertRaisesRegex(ValueError,'must map to a string'):
                validate_schema(profile,pa.schema(fields))

    def test_no_implicit_mapping_or_drop(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            _,config,_,_=self.make_fixture(root)
            profile=load_custom_profile(config)
            with self.assertRaisesRegex(ValueError,'Missing source field'):
                normalize_row(profile,{'reaction_key':'r1'})


if __name__ == '__main__':
    unittest.main()
