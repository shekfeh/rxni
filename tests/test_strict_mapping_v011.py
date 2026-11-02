"""Strict public input and scientifically non-promotional evidence tests."""
import unittest
from unittest.mock import patch
from rxni.input_validation import inspect_mapping, ReactionMappingError
from rxni.optional_mapping import prepare_query
from rxni.retrieval_diagnostics import assess_conditions, summarize_retrieval

GOOD='[CH3:1][Cl:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]'

class StrictMappingTests(unittest.TestCase):
    def test_accepted_mapped(self):
        self.assertEqual(inspect_mapping(GOOD)['status'],'mapped_complete')
    def test_agents_exempt(self):
        valid='[CH3:1][Cl:2]>O>[CH3:1][Cl:2]'
        self.assertEqual(inspect_mapping(valid)['agents_mapping_required'],False)
    def test_unmapped_rejected(self):
        with self.assertRaises(ReactionMappingError) as ctx:
            inspect_mapping('CCCl.NC>>CCNC')
        self.assertEqual(ctx.exception.code,'incomplete_atom_mapping')
    def test_partial_rejected(self):
        with self.assertRaises(ReactionMappingError):
            inspect_mapping('[CH3:1]Cl>>[CH3:1][OH:2]')
    def test_duplicate_rejected(self):
        with self.assertRaises(ReactionMappingError) as ctx:
            inspect_mapping('[CH3:1][Cl:1]>>[CH3:1][OH:2]')
        self.assertEqual(ctx.exception.code,'duplicate_atom_mapping')
    def test_element_inconsistent(self):
        with self.assertRaises(ReactionMappingError):
            inspect_mapping('[CH3:1][Cl:2]>>[CH3:1][OH:2]')
    def test_no_implicit_mapper(self):
        with patch.dict('os.environ',{},clear=True):
            with self.assertRaises(ReactionMappingError):
                prepare_query('CC>>CO')
    def test_diagnostics_do_not_promote_experiment(self):
        assessment=assess_conditions(query_conditions_supplied=True)
        self.assertEqual(assessment['experimentally_validated_recommendations'],[])
        self.assertTrue(assessment['recommendation_status'].startswith('abstain'))
        quality=summarize_retrieval(primary_event='amidation',members=[{'metadata':{'candidate_primary_event':'n_alkylation'},'generation_channels':['component_topology']}],generation_counts={},selected_channel_counts={},candidate_budget=10,selected_count=1)
        self.assertEqual(quality['primary_event_agreement_fraction'],0)
        self.assertEqual(quality['broad_fallback_member_count'],1)

if __name__=='__main__':unittest.main()
