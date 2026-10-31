import hashlib
import tempfile
import unittest
from pathlib import Path
from rxni.builder.bridge import sha256_text,_check_forward_smiles,REACTION_CANDIDATES
from rxni.builder.profiles import select
from rxni.builder.knowledge import make_core,audit_core
from rxni.builder.auxiliary import prepare_runtime_inputs

class ContractTests(unittest.TestCase):
    def test_reaction_priority_and_map_regex(self):
        self.assertIn('reaction',REACTION_CANDIDATES)
        self.assertLess(REACTION_CANDIDATES.index('reaction'),REACTION_CANDIDATES.index('mapping_input'))
        self.assertEqual(_check_forward_smiles('[CH3:1]O>>[CH2:1]=O'),'[CH3:1]O>>[CH2:1]=O')
        with self.assertRaises(ValueError):_check_forward_smiles('c1ccccc1>>c1ccccc1')
        with self.assertRaises(ValueError):_check_forward_smiles('BrC1=C(C=C(C=C1)I)F.CO[H].[Na]>>BrC1=C(C=C(C=C1)I)OC')
    def test_profile_is_version_guarded(self):
        self.assertEqual(select('v6.7').representation_format,'retrieval_representations_v2')
        for x in ('v6.5','v6.6','v7.0'):
            with self.assertRaises(ValueError):select(x)
    def test_template_identity(self):
        s='[C:1]>>[C:1]'
        self.assertEqual(sha256_text(s),hashlib.sha256(s.encode()).hexdigest())

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:
    pa=pq=None

@unittest.skipUnless(pa is not None,'requires pyarrow (available in user aizynthtrain environment)')
class CoreRoundTrip(unittest.TestCase):
    def test_build_audit_and_rehydrate(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td);src=td/'src';src.mkdir()
            smarts='[CH3:1][OH:2]>>[CH3:1].[OH:2]';h=sha256_text(smarts)
            rxn='[CH3:1][OH:2]>>[CH2:1]=[O:2]'
            pq.write_table(pa.Table.from_pylist([dict(retro_template_hash=h,retro_smarts=smarts,
                            retro_template_number='1',native_template_hash='native-1')]),src/'template_library_v5.parquet')
            pq.write_table(pa.Table.from_pylist([
                dict(reaction_key='obs-1',reaction=rxn,native_cleaned_reaction=rxn,mapper_reaction='[CH3:1][OH:2].[N:3]>>[CH2:1]=[O:2]',mapping_input='CO>>C=O',document_id='A'),
                dict(reaction_key='obs-2',reaction=rxn,native_cleaned_reaction=rxn,mapper_reaction='[CH3:1][OH:2].[N:3]>>[CH2:1]=[O:2]',mapping_input='CO>>C=O',document_id='B')]),src/'reaction_core_v5.parquet')
            pq.write_table(pa.Table.from_pylist([
                dict(reaction_key='obs-1',retro_template_hash=h,retro_template_number='1',source_dataset='TEST'),
                dict(reaction_key='obs-2',retro_template_hash=h,retro_template_number='1',source_dataset='TEST')]),src/'reaction_template_membership_v5.parquet')
            (src/'manifest.json').write_text('{"reaction_policy":"native_cleaned_template_extraction_input"}')
            dest=td/'core'
            m=make_core(out_dir=dest,export_dir=src,dataset='TEST')
            self.assertEqual(m['counts'],dict(reactions=1,templates=1,memberships=1,evidence=2))
            self.assertEqual(audit_core(dest)['status'],'PASS')
            self.assertEqual(len(pq.read_table(dest/'templates.parquet')),1)
            prepared=td/'work_inputs';bridge=prepare_runtime_inputs(dest,prepared)
            self.assertEqual(bridge['counts']['memberships'],2)
            self.assertEqual(len(pq.read_table(prepared/'source_memberships.parquet')),2)
    def test_missing_smarts_rejected_without_output(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td);src=td/'src';src.mkdir()
            pq.write_table(pa.Table.from_pylist([{'retro_template_hash':'h','retro_smarts':None}]),src/'template_library_v5.parquet')
            pq.write_table(pa.Table.from_pylist([{'reaction_key':'obs','reaction':'[C:1]>>[C:1]'}]),src/'reaction_core_v5.parquet')
            pq.write_table(pa.Table.from_pylist([{'reaction_key':'obs','retro_template_hash':'h'}]),src/'reaction_template_membership_v5.parquet')
            (src/'manifest.json').write_text('{"reaction_policy":"native_cleaned_template_extraction_input"}')
            with self.assertRaises(ValueError):make_core(out_dir=td/'out',export_dir=src,dataset='TEST')
            self.assertFalse((td/'out').exists())
if __name__=='__main__':unittest.main(verbosity=2)
