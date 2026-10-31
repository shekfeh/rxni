from __future__ import annotations

from typing import Any

from .mapped_graph import canonical_bond_class

EVENT_TOPOLOGY_VERSION = 'i_event_topology_0.2-ring-v02'
SEMANTIC_EVENT_VERSION = 'i_semantic_event_0.2-ring-v02'


def _atom_element(facts: dict[str, Any], mn: int) -> str:
    for side in ('reactant', 'product'):
        rec = facts.get(side, {}).get('atoms', {}).get(int(mn))
        if rec:
            return str(rec.get('element') or '?')
    return '?'


def _bond_class(facts: dict[str, Any], a: int, b: int) -> str:
    return canonical_bond_class(_atom_element(facts, a), _atom_element(facts, b))


def _edit_catalog(facts: dict[str, Any]) -> list[dict[str, Any]]:
    """Deterministic objective edit catalog for event support provenance."""
    edits: list[dict[str, Any]] = []

    for i, rec in enumerate(sorted(
        facts.get('formed_bonds', []) or [],
        key=lambda r: (min(int(r['a']), int(r['b'])), max(int(r['a']), int(r['b'])), float(r.get('order', 0.0)))
    )):
        a, b = int(rec['a']), int(rec['b'])
        edits.append({
            'edit_id': f'F{i}', 'kind': 'formed_bond', 'maps': sorted((a, b)),
            'bond_class': _bond_class(facts, a, b),
            'product_order': float(rec.get('order', 0.0)),
        })

    for i, rec in enumerate(sorted(
        facts.get('broken_bonds', []) or [],
        key=lambda r: (min(int(r['a']), int(r['b'])), max(int(r['a']), int(r['b'])), float(r.get('order', 0.0)))
    )):
        a, b = int(rec['a']), int(rec['b'])
        edits.append({
            'edit_id': f'B{i}', 'kind': 'broken_bond', 'maps': sorted((a, b)),
            'bond_class': _bond_class(facts, a, b),
            'reactant_order': float(rec.get('order', 0.0)),
        })

    for i, rec in enumerate(sorted(
        facts.get('bond_order_changes', []) or [],
        key=lambda r: (min(int(r['a']), int(r['b'])), max(int(r['a']), int(r['b'])),
                       float(r.get('reactant_order', 0.0)), float(r.get('product_order', 0.0)))
    )):
        a, b = int(rec['a']), int(rec['b'])
        edits.append({
            'edit_id': f'O{i}', 'kind': 'bond_order_change', 'maps': sorted((a, b)),
            'bond_class': _bond_class(facts, a, b),
            'reactant_order': float(rec.get('reactant_order', 0.0)),
            'product_order': float(rec.get('product_order', 0.0)),
        })

    for i, rec in enumerate(sorted(
        facts.get('attachment_deltas', []) or [],
        key=lambda r: int(r['center_map'])
    )):
        center = int(rec['center_map'])
        maps = {center}
        for x in list(rec.get('lost', []) or []) + list(rec.get('gained', []) or []):
            if x.get('map') is not None:
                maps.add(int(x['map']))
        edits.append({
            'edit_id': f'A{i}', 'kind': 'attachment_delta', 'maps': sorted(maps),
            'center_map': center,
            'center_element': rec.get('center_element'),
            'lost_elements': sorted(str(x.get('element') or '?') for x in rec.get('lost', []) or []),
            'gained_elements': sorted(str(x.get('element') or '?') for x in rec.get('gained', []) or []),
            'hydrogen_delta': int(rec.get('hydrogen_delta', 0)),
        })

    return edits


def _topology_from_edits(facts: dict[str, Any], edits: list[dict[str, Any]]) -> dict[str, Any]:
    formed = [x for x in edits if x['kind'] == 'formed_bond']
    broken = [x for x in edits if x['kind'] == 'broken_bond']
    order = [x for x in edits if x['kind'] == 'bond_order_change']
    attach = [x for x in edits if x['kind'] == 'attachment_delta']

    support_maps = sorted({m for e in edits for m in e.get('maps', [])})
    formed_classes = sorted({x.get('bond_class') for x in formed if x.get('bond_class')})
    broken_classes = sorted({x.get('bond_class') for x in broken if x.get('bond_class')})
    order_classes = sorted({x.get('bond_class') for x in order if x.get('bond_class')})
    attachment_center_elements = sorted({str(x.get('center_element')) for x in attach if x.get('center_element')})

    intercomponent_formed = []
    for rec in facts.get('intercomponent_formed_bonds', []) or []:
        pair = {int(rec['a']), int(rec['b'])}
        if pair <= set(support_maps):
            intercomponent_formed.append(rec)

    return {
        'version': EVENT_TOPOLOGY_VERSION,
        'support_maps': support_maps,
        'n_formed_bonds': len(formed),
        'n_broken_bonds': len(broken),
        'n_bond_order_changes': len(order),
        'n_attachment_deltas': len(attach),
        'formed_bond_classes': formed_classes,
        'broken_bond_classes': broken_classes,
        'order_change_bond_classes': order_classes,
        'attachment_center_elements': attachment_center_elements,
        'n_intercomponent_formed_bonds': len(intercomponent_formed),
        'ring_formation_count': len((facts.get('ring_lineage_v2') or {}).get('formations', []) or []),
        'ring_opening_count': len((facts.get('ring_lineage_v2') or {}).get('openings', []) or []),
        'ring_expansion_count': len((facts.get('ring_lineage_v2') or {}).get('expansions', []) or []),
        'ring_contraction_count': len((facts.get('ring_lineage_v2') or {}).get('contractions', []) or []),
        'ring_rewiring_count': len((facts.get('ring_lineage_v2') or {}).get('rewirings', []) or []),
        'legacy_ring_formation_count': len(facts.get('ring_formations', []) or []),
        'legacy_ring_opening_count': len(facts.get('ring_openings', []) or []),
        'component_merge_count': len(facts.get('component_merges', []) or []),
        'component_split_count': len(facts.get('component_splits', []) or []),
        'common_cycle_rank_delta': int(facts.get('common_cycle_rank_delta', 0)),
    }


def whole_reaction_event_topology(facts: dict[str, Any]) -> dict[str, Any]:
    edits = _edit_catalog(facts)
    out = _topology_from_edits(facts, edits)
    out['edit_ids'] = [x['edit_id'] for x in edits]
    return out


def build_semantic_events(facts: dict[str, Any], primary_events: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Expose primary semantics as explicit per-event objects.

    Batch07 intentionally does not flatten multiple primary conversions into one label.
    Each resolved primary record becomes one SemanticEvent with deterministic edit
    provenance. Exact duplicate records with identical semantic identity/support are
    collapsed, while two events on the same connected edit center remain distinct.
    """
    catalog = _edit_catalog(facts)
    common = set(int(x) for x in facts.get('mapping', {}).get('common_maps', []) or [])
    out: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()

    for rec in primary_events:
        support = {int(m) for m in (rec.get('support_maps') or []) if int(m) in common}
        if rec.get('center_map') is not None and int(rec['center_map']) in common:
            support.add(int(rec['center_map']))
        for k in ('acyl_center_maps', 'ring_maps'):
            support.update(int(m) for m in (rec.get(k) or []) if int(m) in common)

        support_edits = [e for e in catalog if support and (set(e.get('maps', [])) & support)]
        # If a detector supplied no mapped support, keep the event but do not manufacture edit support.
        topology = _topology_from_edits(facts, support_edits) if support_edits else {
            'version': EVENT_TOPOLOGY_VERSION,
            'support_maps': sorted(support),
            'n_formed_bonds': 0, 'n_broken_bonds': 0, 'n_bond_order_changes': 0,
            'n_attachment_deltas': 0, 'formed_bond_classes': [], 'broken_bond_classes': [],
            'order_change_bond_classes': [], 'attachment_center_elements': [],
            'n_intercomponent_formed_bonds': 0,
            'ring_formation_count': 0, 'ring_opening_count': 0,
            'ring_expansion_count': 0, 'ring_contraction_count': 0, 'ring_rewiring_count': 0,
            'legacy_ring_formation_count': 0, 'legacy_ring_opening_count': 0,
            'component_merge_count': 0, 'component_split_count': 0,
            'common_cycle_rank_delta': int(facts.get('common_cycle_rank_delta', 0)),
        }

        key = (
            rec.get('event'), rec.get('subtype'), rec.get('center_map'),
            tuple(sorted(support)), tuple(e['edit_id'] for e in support_edits),
        )
        if key in seen:
            continue
        seen.add(key)

        event_id = f'E{len(out)}'
        out.append({
            'event_id': event_id,
            'event': rec.get('event'),
            'subtype': rec.get('subtype'),
            'center_map': rec.get('center_map'),
            'support_maps': sorted(support),
            'support_edit_ids': [e['edit_id'] for e in support_edits],
            'event_topology': topology,
            'detector': rec.get('detector'),
            'mapping_scope': rec.get('mapping_scope'),
            'directional_conversion': rec.get('directional_conversion'),
            'redox_direction': rec.get('redox_direction'),
            'component_role': rec.get('component_role'),
            'status': 'primary_event',
            'version': SEMANTIC_EVENT_VERSION,
        })

    return out
