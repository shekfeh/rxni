"""Small independent builder/source and score-contract regression fixtures.

This is not full-corpus historical parity or one-to-many lookup certification.
"""
import hashlib
import tempfile
import unittest
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from rxni.builder.knowledge import make_core, audit_core
from rxni.retrieval.record import build_retrieval_record
from rxni.retrieval.scoring import member_scores
from rxni.retrieval.storage.codec import (
    retrieval_record_to_storage_row_v2,retrieval_record_from_storage_row_v2)


class BuilderRegressionV011(unittest.TestCase):
    def test_multiple_template_membership_preserved_in_core(self):
        # A single reaction legitimately has >1 template linkage. The Core
        # should retain both; final reaction-keyed lookup may still be lossy.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src=root/'src';src.mkdir()
            reaction='[CH3:1][Cl:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]'
            smarts=['[C:1][Cl:2]>>[C:1].[Cl:2]',
                    '[C:1][N:2]>>[C:1].[N:2]']
            entries=[{'retro_template_hash':hashlib.sha256(s.encode()).hexdigest(),
                      'retro_smarts':s,'retro_template_number':str(i),
                      'native_template_hash':'native-'+str(i)}
                     for i,s in enumerate(smarts)]
            pq.write_table(pa.Table.from_pylist(entries),src/'template_library_v5.parquet')
            pq.write_table(pa.Table.from_pylist([{'reaction_key':'obs-01','reaction':reaction,
                'native_cleaned_reaction':reaction,'document_id':'P1'}]),src/'reaction_core_v5.parquet')
            pq.write_table(pa.Table.from_pylist([{'reaction_key':'obs-01',
                'retro_template_hash':entry['retro_template_hash'],
                'retro_template_number':entry['retro_template_number'],'source_dataset':'TEST'}
                for entry in entries]),src/'reaction_template_membership_v5.parquet')
            (src/'manifest.json').write_text('{"reaction_policy":"native_cleaned_template_extraction_input"}')
            manifest=make_core(out_dir=root/'core',export_dir=src,dataset='TEST')
            self.assertEqual(manifest['counts']['reactions'],1)
            self.assertEqual(manifest['counts']['templates'],2)
            self.assertEqual(manifest['counts']['memberships'],2)
            self.assertEqual(audit_core(root/'core')['status'],'PASS')
            self.assertEqual(pq.ParquetFile(root/'core/memberships.parquet').metadata.num_rows,2)

    def test_ranking_and_record_codec_determinism(self):
        # Small deterministic contract check; not source-vs-public scientific parity.
        rxn='[CH3:1][Cl:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]'
        rec=build_retrieval_record(rxn,record_id='test-rxn')
        row=retrieval_record_to_storage_row_v2(rec,build_id='v011-regression')
        row['reaction_key']='test-rxn'
        restored=retrieval_record_from_storage_row_v2(row,expected_build_id='v011-regression')
        a=dict(member_scores(rec,rec))
        b=dict(member_scores(rec,restored))
        self.assertEqual(a,b)

if __name__=='__main__': unittest.main()
