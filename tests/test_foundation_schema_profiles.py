"""Profile validation tests using stdlib SQLite and disposable on-disk fixtures."""
from __future__ import annotations
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

from rxni.retrieval.libraries.foundation_schema import (
    CANONICAL_STRUCTURAL, FoundationSchemaError, audit_key_membership,
    load_foundation_profile, quote_ident, validate_foundation_schema,
)
from rxni.retrieval.libraries import resolve_library, LibraryFormatError


def make_profile(root, *, bad=None):
    src={'schema_version':1,'profile_id':'renamed-index-v1',
         'structural_index':{'table':'rxni_hits',
           'columns':{'reaction_key':'rxn_id',**{name:'hash_'+str(i) for i,name in enumerate(CANONICAL_STRUCTURAL)}}},
         'metadata_lookup':{'table':'rxni_meta','columns':{
             'reaction_key':'key_text','forward_reaction_smiles':'smiles_text',
             'retro_template_hash':'template_id','membership_occurrence_count':'occurrences'}}}
    if bad: bad(src)
    path=root/'foundation.json'; path.write_text(json.dumps(src))
    return path


def populate(root):
    profile=load_foundation_profile(make_profile(root))
    structural=root/'structure.sqlite'; metadata=root/'meta.sqlite'; locator=root/'locator.sqlite'
    with sqlite3.connect(structural) as c:
        cols=['rxn_id TEXT PRIMARY KEY'] + [f'hash_{i} TEXT' for i in range(6)]
        c.execute('CREATE TABLE rxni_hits('+','.join(cols)+')')
        c.execute('INSERT INTO rxni_hits VALUES (?,?,?,?,?,?,?)',('k1',*['h']*6))
    with sqlite3.connect(metadata) as c:
        c.execute('CREATE TABLE rxni_meta (key_text TEXT,smiles_text TEXT,template_id TEXT,occurrences INTEGER)')
        c.execute('INSERT INTO rxni_meta VALUES (?,?,?,?)',('k1','CC>>CO','t1',1))
    with sqlite3.connect(locator) as c:
        c.execute('CREATE TABLE representation_locator_v2 (reaction_key TEXT,build_id TEXT,representation_generation TEXT,row_group_id INTEGER,row_in_group INTEGER)')
        c.execute('INSERT INTO representation_locator_v2 VALUES (?,?,?,?,?)',('k1','b1','retrieval_representations_v2',0,0))
    return profile, structural, metadata, locator


class TestFoundationProfile(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def test_renamed_tables_and_columns(self):
        profile,s,m,l=populate(self.root)
        self.assertEqual(validate_foundation_schema(s,m,profile)['status'],'FOUNDATION_SCHEMA_VALID')
        self.assertEqual(audit_key_membership(s,m,l,profile)['status'],'PASS')

    def test_missing_physical_column_rejected(self):
        profile,s,m,l=populate(self.root)
        with sqlite3.connect(s) as c:c.execute('ALTER TABLE rxni_hits RENAME COLUMN hash_5 TO lost_hash')
        with self.assertRaisesRegex(FoundationSchemaError,'physical columns absent'):
            validate_foundation_schema(s,m,profile)

    def test_broken_key_membership_is_not_accepted(self):
        profile,s,m,l=populate(self.root)
        with sqlite3.connect(m) as c:c.execute("UPDATE rxni_meta SET key_text='wrong'")
        report=audit_key_membership(s,m,l,profile)
        self.assertEqual(report['status'],'KEY_COVERAGE_MISMATCH')
        self.assertEqual(report['missing_metadata'],['k1'])

    def test_no_locator_key_detected(self):
        profile,s,m,l=populate(self.root)
        with sqlite3.connect(l) as c:c.execute('DELETE FROM representation_locator_v2')
        self.assertEqual(audit_key_membership(s,m,l,profile)['missing_v2_locator'],['k1'])

    def test_untrusted_identifier_rejected(self):
        for val in ('table;DROP TABLE x','name.with.dot','has space','', 'a"'):
            with self.assertRaises(FoundationSchemaError):quote_ident(val)

    def test_unknown_json_fields_rejected(self):
        f=make_profile(self.root,bad=lambda x:x.update({'allow_legacy_fallback':True}))
        with self.assertRaisesRegex(FoundationSchemaError,'unknown'):
            load_foundation_profile(f)

    def test_missing_required_semantic_field_rejected(self):
        f=make_profile(self.root,bad=lambda x:x['metadata_lookup']['columns'].pop('forward_reaction_smiles'))
        with self.assertRaisesRegex(FoundationSchemaError,'missing'):
            load_foundation_profile(f)

    def test_duplicate_physical_columns_rejected(self):
        f=make_profile(self.root,bad=lambda x:x['structural_index']['columns'].update({'forward_raw_center_hash_r0':'rxn_id'}))
        with self.assertRaisesRegex(FoundationSchemaError,'duplicate'):
            load_foundation_profile(f)

    def test_descriptor_accepts_explicit_profile_and_rejects_escape(self):
        profile,s,m,l=populate(self.root)
        for item,p in [('structural_index',s),('metadata_lookup',m),('locator',l)]:
            p.rename(self.root/({'structural_index':'structure.db','metadata_lookup':'metadata.db','locator':'locator.db'}[item]))
        (self.root/'records.parquet').write_bytes(b'fixture')
        (self.root/'manifest.json').write_text(json.dumps({'representation_generation':'retrieval_representations_v2','build_id':'b1'}))
        artifacts={'structural_index':'structure.db','metadata_lookup':'metadata.db',
                   'locator':'locator.db','representations':'records.parquet',
                   'retrieval_manifest':'manifest.json'}
        desc={'schema_version':1,'representation_format':'retrieval_representations_v2',
              'foundation_schema_profile':'foundation.json','artifacts':artifacts}
        d=self.root/'rxni-library.json';d.write_text(json.dumps(desc))
        self.assertEqual(resolve_library(self.root).foundation_profile,profile.path)
        desc['foundation_schema_profile']='../escape.json';d.write_text(json.dumps(desc))
        with self.assertRaisesRegex(LibraryFormatError,'within the library'):resolve_library(self.root)

    def test_default_five_file_behavior_remains(self):
        from rxni.retrieval.libraries import LEGACY_V2_NAMES
        root=self.root/'old';root.mkdir()
        for k,v in LEGACY_V2_NAMES.items():
            p=root/v
            if k=='retrieval_manifest':p.write_text(json.dumps({'representation_generation':'retrieval_representations_v2','build_id':'b1'}))
            else:p.write_bytes(b'fixture')
        self.assertIsNone(resolve_library(root).foundation_profile)


if __name__=='__main__':unittest.main()
