"""Non-scoring retrieval/condition evidence interpretation diagnostics."""
from __future__ import annotations
from collections import Counter


def summarize_retrieval(*, primary_event, members, generation_counts,
                        selected_channel_counts, candidate_budget, selected_count,
                        available_counts=None):
    event_counts = Counter(str((m.get('metadata') or {}).get('candidate_primary_event') or 'unassigned') for m in members)
    exact = sum('exact_edit' in (m.get('generation_channels') or []) for m in members)
    weak = sum(not set(m.get('generation_channels') or []).intersection({'exact_edit','raw_r2','raw_r1'}) for m in members)
    agreement = event_counts.get(str(primary_event), 0) if primary_event else 0
    n = len(members)
    return {
        'status':'observational_diagnostics_only',
        'primary_event':primary_event,
        'candidate_primary_event_counts':dict(event_counts),
        'scored_member_count':n,
        'primary_event_agreement_fraction':(agreement/n if n and primary_event else None),
        'exact_edit_member_count':exact,
        'broad_fallback_member_count':weak,
        'broad_fallback_fraction':weak/n if n else None,
        'generation_counts':dict(generation_counts),
        'generation_available_counts':dict(available_counts or {}),
        'selected_channel_counts':dict(selected_channel_counts),
        'candidate_budget':int(candidate_budget),
        'selected_count':int(selected_count),
        'candidate_budget_saturated':int(selected_count)>=int(candidate_budget),
        'budget_sensitivity':'not_evaluated_single_budget',
        'note':'Repeat at another candidate budget to assess ranking sensitivity; these are not chemical accuracy metrics.'
    }


def assess_conditions(*, query_conditions_supplied=False):
    # RXNI does not currently carry independently verified experimental
    # role/temperature/yield evidence in the public historical library.
    return {
        'recommendation_status':'abstain_no_verified_precedent',
        'experimentally_validated_recommendations':[],
        'inferred_method_families':'heuristic_not_experimental',
        'query_conditions_supplied':bool(query_conditions_supplied),
        'reason':'No verified experimental condition metadata sufficient to certify a recommendation is attached to this query.'
    }
