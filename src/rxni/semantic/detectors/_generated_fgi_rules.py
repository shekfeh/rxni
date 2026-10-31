"""AUTO-GENERATED directional-FGI rule core, locally owned by Runtime-v6.7.

Originally generated from frozen v6.5.6 source during the v6.6 migration.
Runtime-v6.7 changes only namespace ownership here; detector logic is unchanged.
No validation/resolution policy belongs here.
"""

from __future__ import annotations


from rxni.semantic.states.common import (
    atom as _atom,
    attachments as _attachments,
    attachment_map as _att_map,
    hydrogen_delta as _hydrogen_delta,
    bond_order as _bond_between,
    element_attachments as _element_attachments,
    support_maps as _conversion_support,
)
from rxni.semantic.states.carbon import (
    anhydride_partner_centers as _anhydride_partner_centers,
    carbon_state as _carbon_fg_state,
)
from rxni.semantic.states.nitrogen import (
    is_azide_n as _is_azide_n,
    nitrogen_state as _nitrogen_fg_state,
)
from rxni.semantic.states.sulfur import (
    sulfur_state as _sulfur_fg_state,
)

from typing import Any
from rxni.semantic.structure.mapped_graph import canonical_bond_class
HALOGENS = {'F', 'Cl', 'Br', 'I'}

def _has_conserved_carbonyl(facts: dict[str, Any], mn: int) -> bool:
    if _atom(facts, 'reactant', mn).get('element') != 'C' or _atom(facts, 'product', mn).get('element') != 'C':
        return False

    def has_co(side: str) -> bool:
        return any((a.get('element') == 'O' and float(a.get('bond_order', 0)) >= 1.9 for a in _attachments(facts, side, mn)))
    return has_co('reactant') and has_co('product')

def _external_gain(facts: dict[str, Any], center: int, att: dict[str, Any]) -> bool:
    gm = _att_map(att)
    if gm is None:
        return False
    return facts['reactant']['map_to_component'].get(gm) != facts['reactant']['map_to_component'].get(center)

def _support(*maps: int | None) -> list[int]:
    return sorted({int(x) for x in maps if x is not None})

def _same_reactant_origin(facts: dict[str, Any], *maps: int | None) -> bool:
    vals = [m for m in maps if m is not None]
    if not vals:
        return False
    comps = {facts['reactant']['map_to_component'].get(int(m)) for m in vals}
    return None not in comps and len(comps) == 1

def _has_external_formed_at(facts: dict[str, Any], center: int) -> bool:
    return any((center in (int(x['a']), int(x['b'])) for x in facts.get('intercomponent_formed_bonds', [])))

def _component_heavy_count(facts: dict[str, Any], side: str, component_id: int | None) -> int:
    if component_id is None:
        return 0
    comps = facts[side].get('components', [])
    if component_id < 0 or component_id >= len(comps):
        return 0
    mol = comps[component_id].get('mol')
    if mol is None:
        return 0
    return sum((1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1))

def _center_component_provenance(facts: dict[str, Any], center: int | None) -> dict[str, Any]:
    if center is None:
        return {}
    rci = facts['reactant']['map_to_component'].get(int(center))
    pci = facts['product']['map_to_component'].get(int(center))
    common = set(facts['mapping']['common_maps'])
    source_maps = sorted((m for m in common if facts['reactant']['map_to_component'].get(m) == rci))
    product_maps = sorted((m for m in common if facts['product']['map_to_component'].get(m) == pci))
    contributing_components: dict[int, int] = {}
    for m in product_maps:
        rc = facts['reactant']['map_to_component'].get(m)
        if rc is not None:
            contributing_components[rc] = contributing_components.get(rc, 0) + 1
    source_common = len(source_maps)
    dominant_other = max((n for (rc, n) in contributing_components.items() if rc != rci), default=0)
    r_heavy = _component_heavy_count(facts, 'reactant', rci)
    merged = len(contributing_components) > 1
    role = 'substrate'
    if merged and r_heavy <= 8 and (dominant_other >= max(4, source_common * 2)):
        role = 'incorporated_reagent'
    return {'source_component': rci, 'product_component': pci, 'source_component_heavy_atoms': r_heavy, 'source_component_common_maps': source_common, 'product_component_common_maps': len(product_maps), 'contributing_reactant_components': contributing_components, 'component_role': role}

def _organic_nitrile_center(facts: dict[str, Any], side: str, center: int) -> bool:
    if _carbon_fg_state(facts, side, center) != 'nitrile':
        return False
    atts = _attachments(facts, side, center)
    nitrile_n = [x for x in atts if x.get('element') == 'N' and float(x.get('bond_order', 0.0)) >= 2.8]
    other = [x for x in atts if x not in nitrile_n]
    return len(nitrile_n) == 1 and any((x.get('element') == 'C' for x in other))

def _ester_alkyl_context(facts: dict[str, Any], side: str, center: int) -> str | None:
    for oatt in _element_attachments(facts, side, center, 'O', 0.8, 1.2):
        om = _att_map(oatt)
        if om is None:
            continue
        for x in _attachments(facts, side, om):
            cm = _att_map(x)
            if cm is None or cm == center or x.get('element') != 'C':
                continue
            ca = _atom(facts, side, cm)
            if int(ca.get('total_h', 0)) >= 3 and int(ca.get('heavy_degree', 0)) <= 1:
                return 'methyl'
            if int(ca.get('total_h', 0)) >= 2:
                return 'alkyl'
    return None

def _conversion_records(name: str, center: int | None, support_maps: list[int], direction: str='neutral', family: str='functional_group_interconversion', source_state: str | None=None, target_state: str | None=None, primary_eligible: bool=True, specific_event_name: str | None=None, extra: dict[str, Any] | None=None) -> list[dict[str, Any]]:
    base = {'event': name, 'subtype': name, 'center_map': center, 'support_maps': sorted(set(support_maps)), 'source_functional_group': source_state, 'target_functional_group': target_state, 'conversion_family': family, 'redox_direction': direction, 'primary_eligible': bool(primary_eligible), 'mapping_scope': 'same_mapped_functional_center_and_validated_local_environment'}
    if extra:
        base.update(extra)
    out = [dict(base, role='directional_conversion')]
    if direction in {'oxidation', 'reduction'}:
        out.append({'event': direction, 'role': 'axis_event', 'subtype': name, 'center_map': center, 'support_maps': sorted(set(support_maps)), 'primary_eligible': bool(primary_eligible), 'mapping_scope': 'derived_from_mapped_directional_conversion'})
    if primary_eligible:
        spec = dict(base, role='specific_event', conversion_specific=True)
        if specific_event_name:
            spec['event'] = specific_event_name
            spec['subtype'] = name
            spec['directional_conversion'] = name
        out.append(spec)
    return out

def _has_other_heavy_edits(facts: dict[str, Any], excluded_maps: set[int]) -> bool:
    for key in ('formed_bonds', 'broken_bonds'):
        for b in facts.get(key, []):
            if not {int(b['a']), int(b['b'])} <= excluded_maps:
                return True
    for b in facts.get('bond_order_changes', []):
        if not {int(b['a']), int(b['b'])} <= excluded_maps:
            return True
    return False

def _lost_gained_for_center(facts: dict[str, Any], center: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    for d in facts.get('attachment_deltas', []):
        if int(d['center_map']) == int(center):
            return (list(d.get('lost', [])), list(d.get('gained', [])))
    return ([], [])

def _donor_class(facts: dict[str, Any], donor: int) -> str:
    a = _atom(facts, 'reactant', donor)
    if a.get('element') != 'C':
        return str(a.get('element', 'unknown')).lower()
    if a.get('aromatic'):
        return 'aryl_carbon'
    atts = _attachments(facts, 'reactant', donor)
    if any((x.get('element') == 'C' and float(x.get('bond_order', 0.0)) >= 2.8 for x in atts)):
        return 'alkynyl_carbon'
    if any((x.get('element') == 'C' and float(x.get('bond_order', 0.0)) >= 1.8 for x in atts)):
        return 'alkenyl_carbon'
    return 'alkyl_carbon'

def detect_functional_group_conversions(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    common = facts['mapping']['common_maps']
    seen: set[tuple[str, int | None]] = set()

    def add(name, center, direction='neutral', family='functional_group_interconversion', source=None, target=None, support=None, primary=True, specific_event_name=None, extra=None):
        key = (name, center)
        if key in seen:
            return
        seen.add(key)
        provenance = _center_component_provenance(facts, center)
        if provenance.get('component_role') == 'incorporated_reagent' and direction in {'oxidation', 'reduction'} and (center is not None) and (not _has_external_formed_at(facts, int(center))):
            primary = False
        merged_extra = dict(provenance)
        if extra:
            merged_extra.update(extra)
        out.extend(_conversion_records(name, center, support or _conversion_support(facts, center), direction, family, source, target, primary, specific_event_name=specific_event_name, extra=merged_extra))
    for mn in common:
        el = _atom(facts, 'reactant', mn).get('element')
        if el == 'C':
            r = _carbon_fg_state(facts, 'reactant', mn)
            p = _carbon_fg_state(facts, 'product', mn)
            external = _has_external_formed_at(facts, mn)
            reduction_no_external = {('aldehyde', 'primary_alcohol'): 'aldehyde_to_primary_alcohol', ('ketone', 'secondary_alcohol'): 'ketone_to_secondary_alcohol', ('ester', 'aldehyde'): 'ester_to_aldehyde', ('ester', 'primary_alcohol'): 'ester_to_primary_alcohol', ('carboxylic_acid', 'primary_alcohol'): 'carboxylic_acid_to_primary_alcohol', ('nitrile', 'amine_carbon'): 'nitrile_to_amine', ('amide', 'amine_carbon'): 'amide_to_amine'}
            if (r, p) in reduction_no_external and (not external):
                add(reduction_no_external[r, p], mn, 'reduction', 'reduction', r, p)
            unrestricted_table = {('primary_alcohol', 'aldehyde'): ('primary_alcohol_to_aldehyde', 'oxidation'), ('secondary_alcohol', 'ketone'): ('secondary_alcohol_to_ketone', 'oxidation'), ('primary_alcohol', 'carboxylic_acid'): ('primary_alcohol_to_carboxylic_acid', 'oxidation'), ('aldehyde', 'carboxylic_acid'): ('aldehyde_to_carboxylic_acid', 'oxidation'), ('aldehyde', 'nitrile'): ('aldehyde_to_nitrile', 'oxidation'), ('amide', 'nitrile'): ('amide_to_nitrile', 'neutral')}
            if (r, p) in unrestricted_table:
                (name, direction) = unrestricted_table[r, p]
                family = direction if direction in {'oxidation', 'reduction'} else 'functional_group_interconversion'
                add(name, mn, direction, family, r, p)
            if r == 'aldehyde' and p == 'ester' and external and (_hydrogen_delta(facts, mn) < 0):
                add('aldehyde_to_ester', mn, 'oxidation', 'oxidation', r, p, specific_event_name='oxidative_esterification')
            if r == '' and p == 'aldehyde':
                ra = _atom(facts, 'reactant', mn)
                if not ra.get('aromatic') and int(ra.get('total_h', 0)) >= 2 and (_hydrogen_delta(facts, mn) <= -1) and any((x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8 for x in _attachments(facts, 'product', mn))):
                    add('c_h_to_aldehyde', mn, 'oxidation', 'oxidation', 'alkyl_c_h', 'aldehyde', specific_event_name='c_h_to_aldehyde_oxidation')
            if r == 'alkyl_halide' and p == 'aldehyde':
                benzylic = any((att.get('element') == 'C' and att.get('map') is not None and bool(_atom(facts, 'reactant', int(att['map'])).get('aromatic')) for att in _attachments(facts, 'reactant', mn)))
                product_carbonyl_o = any((att.get('element') == 'O' and float(att.get('bond_order', 0.0)) >= 1.8 for att in _attachments(facts, 'product', mn)))
                if benzylic and product_carbonyl_o:
                    add('benzylic_halide_to_aldehyde', mn, 'oxidation', 'oxidation', 'benzylic_halide', 'aldehyde')
            if r == 'aryl_halide' and p in {'aldehyde', 'ketone'}:
                (lost, gained) = _lost_gained_for_center(facts, mn)
                if any((x.get('element') in HALOGENS for x in lost)) and any((x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8 for x in gained)):
                    add('heteroaryl_halide_to_oxo', mn, 'neutral', 'functional_group_interconversion', r, p, specific_event_name='heteroaryl_hydrolysis_to_oxo')
            if r == 'nitrile' and p == 'aldehyde' and _organic_nitrile_center(facts, 'reactant', mn):
                add('nitrile_to_aldehyde', mn, 'reduction', 'reduction', r, p)
            if r == 'ester' and p == 'carboxylic_acid':
                alkyl = _ester_alkyl_context(facts, 'reactant', mn)
                extra = {'ester_alkyl_group': alkyl}
                if alkyl == 'methyl':
                    extra['specific_conversion'] = 'methyl_ester_to_carboxylic_acid'
                add('ester_to_carboxylic_acid', mn, 'neutral', 'functional_group_interconversion', r, p, extra=extra)
            elif r == 'carboxylic_acid' and p == 'ester':
                (ester_o, ester_donor) = _ester_oxygen_and_donor(facts, 'product', mn, require_new_o_c=True)
                add('carboxylic_acid_to_ester', mn, 'neutral', 'functional_group_interconversion', r, p, specific_event_name='esterification', extra={'ester_alkyl_group': _ester_alkyl_context(facts, 'product', mn), 'ester_oxygen_map': ester_o, 'ester_donor_map': ester_donor, 'mapping_scope_note': 'same_acyl_center_not_mechanistic_oxygen_provenance'})
            elif r == 'nitrile' and p == 'amide' and _organic_nitrile_center(facts, 'reactant', mn):
                add('nitrile_to_amide', mn, 'neutral', 'functional_group_interconversion', r, p, specific_event_name='nitrile_hydration')
            elif r == 'nitrile' and p == 'carboxylic_acid' and _organic_nitrile_center(facts, 'reactant', mn):
                add('nitrile_to_carboxylic_acid', mn, 'neutral', 'functional_group_interconversion', r, p)
            if r == '' and p == 'ketone':
                ra = _atom(facts, 'reactant', mn)
                if not ra.get('aromatic') and int(ra.get('total_h', 0)) >= 2 and any((x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8 for x in _attachments(facts, 'product', mn))):
                    add('methylene_to_ketone', mn, 'oxidation', 'oxidation', 'methylene', 'ketone', extra={'conversion_parent': 'c_h_to_carbonyl'})
            if r in {'primary_alcohol', 'secondary_alcohol', 'tertiary_alcohol'} and p in {'alkyl_halide', 'aryl_halide'}:
                add('alcohol_to_halide', mn, 'neutral', 'functional_group_interconversion', r, p)
            elif r in {'alkyl_halide', 'aryl_halide'} and p in {'primary_alcohol', 'secondary_alcohol', 'tertiary_alcohol'}:
                add('halide_to_alcohol', mn, 'neutral', 'functional_group_interconversion', r, p)
            if r in {'alkyl_halide', 'aryl_halide'} and p == '' and (_hydrogen_delta(facts, mn) > 0) and (not external):
                name = 'aryl_halide_to_hydrogen' if _atom(facts, 'reactant', mn).get('aromatic') else 'alkyl_halide_to_hydrogen'
                add(name, mn, 'reduction', 'reduction', r, 'hydrocarbon')
            if r in {'aldehyde', 'ketone'} and external:
                formed_n = []
                for b in facts.get('intercomponent_formed_bonds', []):
                    if mn not in (b['a'], b['b']):
                        continue
                    other = int(b['b'] if int(b['a']) == mn else b['a'])
                    if _atom(facts, 'reactant', other).get('element') == 'N':
                        formed_n.append(other)
                single_n = [n for n in formed_n if 0.8 <= _bond_between(facts, 'product', mn, n) <= 1.2]
                if single_n and (not _element_attachments(facts, 'product', mn, 'O', 1.8)) and (_hydrogen_delta(facts, mn) > 0):
                    add('carbonyl_to_amine_reductive_amination', mn, 'reduction', 'reductive_amination', r, 'amine_carbon', _support(mn, *single_n), specific_event_name='reductive_amination')
            (lost, gained) = _lost_gained_for_center(facts, mn)
            azide_gains = [g for g in gained if _att_map(g) is not None and _is_azide_n(facts, 'reactant', int(_att_map(g)))]
            if azide_gains:
                azn = int(_att_map(azide_gains[0]))
                epoxide_open = any((rec.get('subfamily') == 'epoxide_ring_opening' and mn in set(rec.get('conserved_maps', [])) for rec in facts.get('ring_openings', [])))
                if epoxide_open:
                    add('epoxide_to_azido_alcohol', mn, 'neutral', 'functional_group_interconversion', 'epoxide', 'azido_alcohol', _support(mn, azn), primary=False, extra={'specific_event': 'epoxide_azidolysis'})
                else:
                    lost_hal = [x for x in lost if x.get('element') in HALOGENS]
                    lost_oh = [x for x in lost if x.get('element') == 'O' and int(x.get('degree', 0)) <= 1]
                    lost_sulf = [x for x in lost if x.get('element') == 'O' and any((str(y[0]) == 'S' for y in x.get('neighbor_shell') or []))]
                    if lost_hal:
                        aromatic = bool(_atom(facts, 'reactant', mn).get('aromatic'))
                        name = 'aryl_halide_to_azide' if aromatic else 'alkyl_halide_to_azide'
                        add(name, mn, 'neutral', 'functional_group_interconversion', 'aryl_halide' if aromatic else 'alkyl_halide', 'organic_azide', _support(mn, azn), specific_event_name='azide_substitution', extra={'leaving_group': lost_hal[0].get('element')})
                    elif lost_sulf:
                        add('sulfonate_to_azide', mn, 'neutral', 'functional_group_interconversion', 'sulfonate', 'organic_azide', _support(mn, azn), specific_event_name='azide_substitution')
                    elif lost_oh:
                        add('alcohol_to_azide', mn, 'neutral', 'functional_group_interconversion', 'alcohol', 'organic_azide', _support(mn, azn), specific_event_name='azide_substitution')
        elif el == 'N':
            r = _nitrogen_fg_state(facts, 'reactant', mn)
            p = _nitrogen_fg_state(facts, 'product', mn)
            if r == 'nitro' and p in {'amine_nitrogen', 'aza_nitrogen'}:
                carbon_neighbors = [x for x in _attachments(facts, 'reactant', mn) if x.get('element') == 'C']
                retained_carbon_maps = {int(_att_map(x)) for x in carbon_neighbors if _att_map(x) is not None} & {int(_att_map(x)) for x in _attachments(facts, 'product', mn) if x.get('element') == 'C' and _att_map(x) is not None}
                product_noncarbon_heavy = [x for x in _attachments(facts, 'product', mn) if x.get('element') not in {'C', 'H'}]
                if carbon_neighbors and retained_carbon_maps and (not product_noncarbon_heavy):
                    context = 'aryl' if any((_att_map(x) is not None and _atom(facts, 'reactant', int(_att_map(x))).get('aromatic') for x in carbon_neighbors)) else 'alkyl'
                    add('nitro_to_amine', mn, 'reduction', 'reduction', r, p, extra={'nitro_context': context})
            elif r == 'azide_n' and p in {'amine_nitrogen', 'aza_nitrogen'}:
                add('azide_to_amine', mn, 'reduction', 'reduction', r, p)
            elif r == 'aza_n_oxide' and p == 'aza_nitrogen':
                support = set(_conversion_support(facts, mn))
                primary = not _has_other_heavy_edits(facts, support)
                add('aza_n_oxide_to_aza_nitrogen', mn, 'reduction', 'reduction', r, p, primary=primary, extra={'specific_event': 'aza_n_oxide_deoxygenation'})
            elif r == 'amine_n_oxide' and p == 'amine_nitrogen':
                support = set(_conversion_support(facts, mn))
                primary = not _has_other_heavy_edits(facts, support)
                add('tertiary_amine_n_oxide_to_amine', mn, 'reduction', 'reduction', r, p, primary=primary)
            elif r == 'aza_nitrogen' and p == 'aza_n_oxide':
                add('aza_nitrogen_to_n_oxide', mn, 'oxidation', 'oxidation', r, p)
            elif r == 'amine_nitrogen' and p == 'amine_n_oxide':
                add('tertiary_amine_to_n_oxide', mn, 'oxidation', 'oxidation', r, p)
        elif el == 'S':
            r = _sulfur_fg_state(facts, 'reactant', mn)
            p = _sulfur_fg_state(facts, 'product', mn)
            table = {('sulfide', 'sulfoxide'): ('sulfide_to_sulfoxide', 'oxidation'), ('sulfide', 'sulfone'): ('sulfide_to_sulfone', 'oxidation'), ('sulfoxide', 'sulfone'): ('sulfoxide_to_sulfone', 'oxidation'), ('sulfone', 'sulfoxide'): ('sulfone_to_sulfoxide', 'reduction'), ('sulfoxide', 'sulfide'): ('sulfoxide_to_sulfide', 'reduction'), ('sulfone', 'sulfide'): ('sulfone_to_sulfide', 'reduction')}
            if (r, p) in table:
                (name, direction) = table[r, p]
                add(name, mn, direction, direction, r, p)
    for bc in facts.get('bond_order_changes', []):
        (a, b) = (int(bc['a']), int(bc['b']))
        if _atom(facts, 'reactant', a).get('element') != 'C' or _atom(facts, 'reactant', b).get('element') != 'C':
            continue
        (ro, po) = (float(bc['reactant_order']), float(bc['product_order']))
        if ro <= 1.2 and po >= 1.8 and (_hydrogen_delta(facts, a) < 0) and (_hydrogen_delta(facts, b) < 0):
            (la, ga) = _lost_gained_for_center(facts, a)
            (lb, gb) = _lost_gained_for_center(facts, b)
            heavy_loss = [x for x in la + lb if x.get('element') != 'H']
            heavy_gain = [x for x in ga + gb if x.get('element') != 'H']
            if not heavy_loss and (not heavy_gain):
                name = 'alkane_to_alkene'
                extra = {}
                carbonyl_adj = any((any((att.get('element') == 'C' and bool(att.get('carbonyl_like')) for att in _attachments(facts, 'product', x))) for x in (a, b)))
                spec_name = 'carbonyl_alpha_beta_dehydrogenation' if carbonyl_adj else 'dehydrogenation'
                add(name, a, 'oxidation', 'oxidation', 'alkane', 'alkene', _support(a, b), specific_event_name=spec_name, extra=extra)
        if ro >= 1.8 and po <= 1.2:
            oa = {x.get('map') for x in _element_attachments(facts, 'product', a, 'O', 0.8, 1.2) if x.get('map') is not None}
            ob = {x.get('map') for x in _element_attachments(facts, 'product', b, 'O', 0.8, 1.2) if x.get('map') is not None}
            shared_o = sorted(oa & ob)
            if shared_o and ro < 2.8:
                add('alkene_to_epoxide', a, 'oxidation', 'oxidation', 'alkene', 'epoxide', _support(a, b, *shared_o))
                continue

            def gained_oh(center: int) -> list[int]:
                outm = []
                for d in facts['attachment_deltas']:
                    if int(d['center_map']) != center:
                        continue
                    for x in d.get('gained', []):
                        if x.get('element') == 'O' and int(x.get('degree', 0)) <= 1 and (x.get('map') is not None):
                            outm.append(int(x['map']))
                return outm
            (ga, gb) = (gained_oh(a), gained_oh(b))
            if ga and gb and (ro < 2.8):
                add('alkene_to_diol', a, 'oxidation', 'oxidation', 'alkene', 'vicinal_diol', _support(a, b, *ga, *gb))
                continue
            if _hydrogen_delta(facts, a) > 0 and _hydrogen_delta(facts, b) > 0 and (not _has_external_formed_at(facts, a)) and (not _has_external_formed_at(facts, b)):
                if ro >= 2.8 and po <= 1.2:
                    name = 'alkyne_to_alkane'
                elif ro >= 2.8 and po >= 1.8:
                    name = 'alkyne_to_alkene'
                else:
                    name = 'alkene_to_alkane'
                spec_name = 'complete_alkyne_hydrogenation' if name == 'alkyne_to_alkane' else 'alkene_hydrogenation' if name == 'alkene_to_alkane' else 'partial_alkyne_hydrogenation'
                add(name, a, 'reduction', 'reduction', 'pi_bond', 'reduced_pi_bond', _support(a, b), specific_event_name=spec_name)
    for d in facts.get('attachment_deltas', []):
        c = int(d['center_map'])
        if _atom(facts, 'reactant', c).get('element') != 'C':
            continue
        lost_pi_c = [x for x in d.get('lost', []) if x.get('element') == 'C' and float(x.get('bond_order', 0.0)) >= 1.8]
        gained_co = [x for x in d.get('gained', []) if x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8]
        if lost_pi_c and gained_co and (_carbon_fg_state(facts, 'product', c) == 'aldehyde'):
            add('alkene_oxidative_cleavage_to_aldehyde', c, 'oxidation', 'oxidation', 'alkene_carbon', 'aldehyde')
    return out

def _ring_formation_touches(facts: dict[str, Any], *maps: int | None) -> bool:
    targets = {int(m) for m in maps if m is not None}
    if not targets:
        return False
    for rec in facts.get('ring_formations', []):
        if targets & set(rec.get('conserved_maps') or rec.get('ring_maps') or []):
            return True
    return False

def _is_weinreb_amide_center(facts: dict[str, Any], center: int) -> bool:
    if _carbon_fg_state(facts, 'reactant', center) != 'amide':
        return False
    for natt in _element_attachments(facts, 'reactant', center, 'N', 0.8, 1.6):
        shell = list(natt.get('neighbor_shell') or [])
        elems = {str(x[0]) for x in shell}
        if 'O' in elems and 'C' in elems:
            return True
    return False

def _new_acyl_attachment(facts: dict[str, Any], center: int, element: str, max_order: float=1.6) -> list[int]:
    out: list[int] = []
    for x in _attachments(facts, 'product', center):
        if x.get('element') != element or x.get('map') is None:
            continue
        bo = float(x.get('bond_order', 0.0))
        if not 0.8 <= bo <= max_order:
            continue
        m = int(x['map'])
        if _bond_between(facts, 'reactant', center, m) < 0.2:
            out.append(m)
    return sorted(set(out))

def _ester_oxygen_and_donor(facts: dict[str, Any], side: str, center: int, require_new_o_c: bool=False) -> tuple[int | None, int | None]:
    for oatt in _element_attachments(facts, side, center, 'O', 0.8, 1.2):
        om = _att_map(oatt)
        if om is None:
            continue
        for x in _attachments(facts, side, om):
            cm = _att_map(x)
            if cm is None or cm == center or x.get('element') != 'C':
                continue
            if require_new_o_c and _bond_between(facts, 'reactant', om, int(cm)) > 0.2:
                continue
            return (int(om), int(cm))
    return (None, None)

def _new_ester_oxygen_and_donor(facts: dict[str, Any], center: int) -> tuple[int | None, int | None]:
    for om in _new_acyl_attachment(facts, center, 'O', 1.2):
        for x in _attachments(facts, 'product', om):
            cm = _att_map(x)
            if cm is None or cm == center or x.get('element') != 'C':
                continue
            return (om, int(cm))
    return (None, None)

def detect_acyl_conversions(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, int | None]] = set()

    def add_conversion(name: str, center: int, source: str, target: str, event: str | None=None, support: list[int] | None=None, extra: dict[str, Any] | None=None):
        key = (name, center)
        if key in seen:
            return
        seen.add(key)
        merged = _center_component_provenance(facts, center)
        if extra:
            merged.update(extra)
        out.extend(_conversion_records(name, center, support or _conversion_support(facts, center), 'neutral', 'acyl_conversion', source, target, True, specific_event_name=event, extra=merged))
    anhydride_pairs: set[tuple[int, int]] = set()
    for c in facts['mapping']['common_maps']:
        if _carbon_fg_state(facts, 'reactant', c) != 'carboxylic_acid':
            continue
        if _carbon_fg_state(facts, 'product', c) != 'anhydride':
            continue
        for other in _anhydride_partner_centers(facts, 'product', c):
            if other not in facts['mapping']['common_maps']:
                continue
            if _carbon_fg_state(facts, 'reactant', other) != 'carboxylic_acid':
                continue
            pair = tuple(sorted((int(c), int(other))))
            if pair in anhydride_pairs:
                continue
            anhydride_pairs.add(pair)
            intra = _same_reactant_origin(facts, *pair)
            out.append({'event': 'intramolecular_anhydride_formation' if intra else 'anhydride_formation', 'role': 'specific_event', 'subtype': 'carboxylic_acids_to_anhydride', 'center_map': pair[0], 'acyl_center_maps': list(pair), 'support_maps': list(pair), 'mapping_scope': 'paired_mapped_acyl_centers_joined_by_product_anhydride_oxygen'})
    for c in facts['mapping']['common_maps']:
        if _atom(facts, 'reactant', c).get('element') != 'C' or not _has_conserved_carbonyl(facts, c):
            continue
        r = _carbon_fg_state(facts, 'reactant', c)
        p = _carbon_fg_state(facts, 'product', c)
        (lost, gained) = _lost_gained_for_center(facts, c)
        lost_elements = [x.get('element') for x in lost]
        gained_elements = [x.get('element') for x in gained]
        if r == 'ester' and p == 'ester' and ('O' in lost_elements) and ('O' in gained_elements):
            gained_o = next((x for x in gained if x.get('element') == 'O'), None)
            if gained_o is not None and _external_gain(facts, c, gained_o):
                add_conversion('transesterification', c, 'ester', 'ester', event='transesterification')
        if r == 'amide' and p == 'amide' and ('N' in lost_elements) and ('N' in gained_elements):
            gained_n = next((x for x in gained if x.get('element') == 'N'), None)
            if gained_n is not None and _external_gain(facts, c, gained_n):
                add_conversion('transamidation', c, 'amide', 'amide', event='transamidation')
        new_ns = _new_acyl_attachment(facts, c, 'N', 1.6)
        if r in {'carboxylic_acid', 'ester', 'acyl_halide', 'haloformate'} and new_ns:
            nmap = new_ns[0]
            cyclic = int(facts.get('common_cycle_rank_delta', 0)) > 0 and _ring_formation_touches(facts, c, nmap)
            if r == 'carboxylic_acid':
                (direction, event) = ('carboxylic_acid_to_amide', 'lactamization' if cyclic else 'amidation')
            elif r == 'ester':
                (direction, event) = ('ester_to_amide', 'lactamization' if cyclic else 'amidation')
            elif r == 'acyl_halide':
                (direction, event) = ('acyl_halide_to_amide', 'lactamization' if cyclic else 'amidation')
            else:
                (direction, event) = ('haloformate_to_carbamate', 'carbamate_formation')
            add_conversion(direction, c, r, 'amide' if r != 'haloformate' else 'carbamate', event=event, support=_support(c, nmap), extra={'nucleophile_map': nmap, 'cyclic': cyclic})
        if r == 'carboxylic_acid' and p == 'acyl_halide':
            gained_x = next((x for x in gained if x.get('element') in HALOGENS), None)
            if gained_x is not None:
                add_conversion('carboxylic_acid_to_acyl_halide', c, r, p, event='acyl_halide_formation', support=_support(c, _att_map(gained_x)), extra={'halogen': gained_x.get('element')})
        if r == 'carboxylic_acid' and p == 'ester':
            (om, donor) = _new_ester_oxygen_and_donor(facts, c)
            if om is not None and int(facts.get('common_cycle_rank_delta', 0)) > 0 and _ring_formation_touches(facts, c, om):
                out.append({'event': 'lactonization', 'role': 'specific_event', 'subtype': 'intramolecular_carboxylic_acid_to_ester', 'center_map': c, 'oxygen_map': om, 'donor_map': donor, 'support_maps': _support(c, om, donor), 'mapping_scope': 'mapped_acid_to_ester_plus_validated_ring_formation'})
        if _is_weinreb_amide_center(facts, c) and p == 'ketone':
            new_cs = _new_acyl_attachment(facts, c, 'C', 1.2)
            if new_cs:
                add_conversion('weinreb_amide_to_ketone', c, 'weinreb_amide', 'ketone', event='weinreb_amide_to_ketone', support=_support(c, new_cs[0]), extra={'incoming_carbon_map': new_cs[0], 'incoming_carbon_class': _donor_class(facts, new_cs[0])})
        if r == 'acyl_halide' and p == 'ketone':
            new_cs = _new_acyl_attachment(facts, c, 'C', 1.2)
            for cm in new_cs:
                if _atom(facts, 'reactant', cm).get('aromatic'):
                    add_conversion('acyl_halide_to_aryl_ketone', c, 'acyl_halide', 'ketone', event='aromatic_acylation', support=_support(c, cm), extra={'incoming_aromatic_carbon_map': cm})
                    break
    return out

def detect_nitrile_thioamide_and_ester(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for c in facts['mapping']['common_maps']:
        if not _organic_nitrile_center(facts, 'reactant', int(c)):
            continue
        p_atts = _attachments(facts, 'product', int(c))
        has_s_double = any((x.get('element') == 'S' and float(x.get('bond_order', 0.0)) >= 1.8 for x in p_atts))
        has_n_single = any((x.get('element') == 'N' and float(x.get('bond_order', 0.0)) <= 1.2 for x in p_atts))
        if has_s_double and has_n_single:
            out.extend(_conversion_records('nitrile_to_thioamide', int(c), _conversion_support(facts, int(c)), 'neutral', 'functional_group_interconversion', 'nitrile', 'thioamide', True, specific_event_name='nitrile_to_thioamide'))
            continue
        if _carbon_fg_state(facts, 'product', int(c)) == 'ester':
            out.extend(_conversion_records('nitrile_to_ester', int(c), _conversion_support(facts, int(c)), 'neutral', 'functional_group_interconversion', 'nitrile', 'ester', True, specific_event_name='nitrile_to_ester'))
    return out

def detect_structural_semantic_extensions(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, int | None]] = set()

    def atom(side: str, mn: int) -> dict[str, Any]:
        try:
            return _atom(facts, side, int(mn)) or {}
        except Exception:
            return {}

    def atts(side: str, mn: int) -> list[dict[str, Any]]:
        try:
            return list(_attachments(facts, side, int(mn)) or [])
        except Exception:
            return []

    def bond(side: str, a: int, b: int) -> float:
        try:
            return float(_bond_between(facts, side, int(a), int(b)) or 0.0)
        except Exception:
            return 0.0

    def fg(side: str, mn: int) -> str:
        try:
            return str(_carbon_fg_state(facts, side, int(mn)) or '')
        except Exception:
            return ''

    def emit_specific(event: str, center: int | None, subtype: str | None=None, support=None, redox: str | None=None, extra: dict[str, Any] | None=None) -> None:
        key = (event, center)
        if key in seen:
            return
        seen.add(key)
        rec = {'event': event, 'role': 'specific_event', 'subtype': subtype or event, 'center_map': center, 'support_maps': list(support or ([] if center is None else _support(center))), 'mapping_scope': 's1_mapped_structural_semantic_closure'}
        if redox:
            rec['redox_direction'] = redox
        if extra:
            rec.update(extra)
        out.append(rec)

    def emit_conversion(name: str, center: int, direction: str, source: str, target: str, specific: str | None=None, support=None) -> None:
        key = (specific or name, center)
        if key in seen:
            return
        seen.add(key)
        family = direction if direction in {'oxidation', 'reduction'} else 'functional_group_interconversion'
        out.extend(_conversion_records(name, center, support or _conversion_support(facts, center), direction, family, source, target, True, specific_event_name=specific, extra={'mapping_scope': 's1_mapped_structural_semantic_closure'}))
    deltas: dict[int, dict[str, Any]] = {}
    maps: set[int] = set()
    for d in facts.get('attachment_deltas', []) or []:
        cm = d.get('center_map')
        if cm is None:
            continue
        c = int(cm)
        maps.add(c)
        deltas[c] = d
        for x in list(d.get('lost', []) or []) + list(d.get('gained', []) or []):
            m = _att_map(x)
            if m is not None:
                maps.add(int(m))
    for fb in facts.get('intercomponent_formed_bonds', []) or []:
        for k in ('a', 'b'):
            if fb.get(k) is not None:
                maps.add(int(fb[k]))
    for key in ('bond_order_changes', 'bond_order_deltas', 'changed_bonds', 'bond_changes'):
        for rec in facts.get(key, []) or []:
            if not isinstance(rec, dict):
                continue
            for k in ('a', 'b', 'atom1', 'atom2', 'map1', 'map2'):
                v = rec.get(k)
                if isinstance(v, int):
                    maps.add(int(v))

    def elems(xs):
        return [str(x.get('element') or '') for x in xs]

    def lost(c):
        return list(deltas.get(c, {}).get('lost', []) or [])

    def gained(c):
        return list(deltas.get(c, {}).get('gained', []) or [])

    def hdelta(c):
        try:
            return int(deltas.get(c, {}).get('hydrogen_delta', _hydrogen_delta(facts, c)))
        except Exception:
            return 0

    def mapped_heavy_neighbors(side: str, center: int, *, exclude_maps=(), exclude_elements=()) -> set[int]:
        excluded_maps = {int(x) for x in exclude_maps if x is not None}
        excluded_elements = set(exclude_elements)
        out = set()
        for x in atts(side, center):
            el = str(x.get('element') or '')
            if not el or el == 'H' or el in excluded_elements:
                continue
            m = _att_map(x)
            if m is None:
                continue
            m = int(m)
            if m not in excluded_maps:
                out.add(m)
        return out

    def mapped_element_neighbors(side: str, center: int, element: str) -> set[int]:
        out = set()
        for x in atts(side, center):
            if x.get('element') != element:
                continue
            m = _att_map(x)
            if m is not None:
                out.add(int(m))
        return out

    def intercomponent_bond_was_formed(a: int, b: int) -> bool:
        target = {int(a), int(b)}
        for fb in facts.get('intercomponent_formed_bonds', []) or []:
            if {int(fb['a']), int(fb['b'])} == target:
                return True
        return False
    for c in sorted(maps):
        if atom('reactant', c).get('element') != 'C':
            continue
        gx = [x for x in gained(c) if x.get('element') in HALOGENS]
        if not gx or hdelta(c) >= 0:
            continue
        external_x = []
        for x in gx:
            xm = _att_map(x)
            if xm is None:
                continue
            xm = int(xm)
            existed_on_reactant = bool(atom('reactant', xm).get('element'))
            if not existed_on_reactant or intercomponent_bond_was_formed(c, xm):
                external_x.append(x)
        if not external_x:
            continue
        lost_nonhalogen_heavy = [x for x in lost(c) if x.get('element') not in {'', 'H'} and x.get('element') not in HALOGENS]
        if lost_nonhalogen_heavy:
            continue
        r_skeleton = mapped_heavy_neighbors('reactant', c, exclude_elements=HALOGENS)
        p_skeleton = mapped_heavy_neighbors('product', c, exclude_elements=HALOGENS)
        if r_skeleton != p_skeleton:
            continue
        subtype = 'c_h_halogenation'
        alpha_acid = False
        for a in atts('reactant', c):
            m = _att_map(a)
            if m is not None and fg('reactant', int(m)) == 'carboxylic_acid':
                alpha_acid = True
                break
        if alpha_acid:
            subtype = 'alpha_carboxylic_acid_halogenation'
        emit_specific('c_h_halogenation', c, subtype=subtype, support=_support(c, *[_att_map(x) for x in external_x]))

    def oxygen_attachment_signature(side: str, center: int) -> tuple[tuple[int | None, float], ...]:
        sig = []
        for x in atts(side, center):
            if x.get('element') != 'O':
                continue
            m = _att_map(x)
            sig.append((None if m is None else int(m), round(float(x.get('bond_order', 0.0)), 2)))
        return tuple(sorted(sig, key=lambda z: (-1 if z[0] is None else z[0], z[1])))
    for c in sorted(maps):
        ra = atom('reactant', c)
        if ra.get('element') != 'C' or ra.get('aromatic'):
            continue
        lo = [x for x in lost(c) if x.get('element') == 'O']
        reactant_o_single = [x for x in atts('reactant', c) if x.get('element') == 'O' and 0.8 <= float(x.get('bond_order', 0.0)) <= 1.2]
        reactant_o_double = [x for x in atts('reactant', c) if x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8]
        product_o = [x for x in atts('product', c) if x.get('element') == 'O']
        heavy_gained = [x for x in gained(c) if x.get('element') not in {'', 'H'}]
        lost_non_o_heavy = [x for x in lost(c) if x.get('element') not in {'', 'H', 'O'}]
        r_skeleton = mapped_heavy_neighbors('reactant', c, exclude_elements={'O'})
        p_skeleton = mapped_heavy_neighbors('product', c, exclude_elements={'O'})
        lost_o_maps = {int(m) for x in lo if (m := _att_map(x)) is not None}
        oxygen_migrated = False
        for om in lost_o_maps:
            po = atom('product', om)
            if po.get('element') != 'O':
                continue
            if bond('product', c, om) <= 0.1:
                oxygen_migrated = True
                break
        if oxygen_migrated:
            continue
        other_oxygen_reorganization = False
        for d in sorted(maps):
            if d == c:
                continue
            if atom('reactant', d).get('element') != 'C' and atom('product', d).get('element') != 'C':
                continue
            rsig = oxygen_attachment_signature('reactant', d)
            psig = oxygen_attachment_signature('product', d)
            if rsig != psig:
                other_oxygen_reorganization = True
                break
        if other_oxygen_reorganization:
            continue
        if lo and reactant_o_single and (not reactant_o_double) and (not product_o) and (not heavy_gained) and (not lost_non_o_heavy) and (hdelta(c) > 0) and (r_skeleton == p_skeleton):
            emit_conversion('carbon_oxygen_to_hydrogen_deoxygenation', c, 'reduction', 'c_o_substituted_carbon', 'hydrocarbon')
    for c in sorted(maps):
        if atom('reactant', c).get('element') != 'C':
            continue
        lsi = [x for x in lost(c) if x.get('element') == 'Si']
        go = [x for x in gained(c) if x.get('element') == 'O']
        if lsi and go:
            emit_conversion('organosilane_to_alcohol', c, 'oxidation', 'organosilane_carbon', 'alcohol')
    for c in sorted(maps):
        if atom('reactant', c).get('element') != 'C':
            continue
        r = fg('reactant', c)
        pfg = fg('product', c)
        r_carbons = mapped_element_neighbors('reactant', c, 'C')
        p_carbons = mapped_element_neighbors('product', c, 'C')
        carbon_skeleton_conserved = r_carbons == p_carbons
        if r == 'carboxylic_acid' and pfg == 'nitrile' and carbon_skeleton_conserved:
            emit_conversion('carboxylic_acid_to_nitrile', c, 'neutral', r, pfg)
        if r == 'amide' and pfg == 'aldehyde' and carbon_skeleton_conserved and (hdelta(c) > 0):
            emit_conversion('amide_to_aldehyde', c, 'reduction', r, pfg)
        if r == 'carboxylic_acid' and pfg == 'aldehyde' and carbon_skeleton_conserved and (hdelta(c) > 0):
            emit_conversion('carboxylic_acid_to_aldehyde', c, 'reduction', r, pfg)
        if r in {'aldehyde', 'ketone'}:
            pc = [x for x in atts('product', c) if x.get('element') == 'C' and float(x.get('bond_order', 0.0)) >= 2.8]
            ro = [x for x in atts('reactant', c) if x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8]
            if pc and ro:
                emit_specific('carbonyl_to_alkyne', c, subtype='deoxygenative_carbonyl_to_alkyne', support=_support(c, *[_att_map(x) for x in pc]))
        if r == 'acyl_halide':
            pc_double = [x for x in atts('product', c) if x.get('element') == 'C' and float(x.get('bond_order', 0.0)) >= 1.8]
            po_double = [x for x in atts('product', c) if x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8]
            if pc_double and po_double:
                emit_conversion('acyl_halide_to_ketene', c, 'neutral', 'acyl_halide', 'ketene')
    for c in sorted(maps):
        if not _organic_nitrile_center(facts, 'reactant', c):
            continue
        if _ring_formation_touches(facts, c):
            continue
        if atom('product', c).get('aromatic'):
            continue
        r_carbons = mapped_element_neighbors('reactant', c, 'C')
        p_carbons = mapped_element_neighbors('product', c, 'C')
        if r_carbons != p_carbons:
            continue
        patts = atts('product', c)
        n_atts = [x for x in patts if x.get('element') == 'N']
        o_atts = [x for x in patts if x.get('element') == 'O']
        if len(n_atts) == 2:
            n_orders = sorted((float(x.get('bond_order', 0.0)) for x in n_atts))
            if n_orders[0] <= 1.2 and n_orders[1] >= 1.5:
                emit_conversion('nitrile_to_amidine', c, 'neutral', 'nitrile', 'amidine')
        elif len(o_atts) >= 3 and (not n_atts):
            emit_conversion('nitrile_to_orthoester', c, 'neutral', 'nitrile', 'orthoester')
        elif n_atts and o_atts:
            if any((float(x.get('bond_order', 0.0)) >= 1.5 for x in n_atts)):
                emit_conversion('nitrile_to_imidate', c, 'neutral', 'nitrile', 'imidate')
    for n in sorted(maps):
        if atom('reactant', n).get('element') != 'N':
            continue
        r_c = [x for x in atts('reactant', n) if x.get('element') == 'C']
        p_n = [x for x in atts('product', n) if x.get('element') == 'N' and float(x.get('bond_order', 0.0)) >= 2.5]
        if r_c and p_n:
            aromatic_attachment = False
            for x in r_c:
                m = _att_map(x)
                if m is not None and atom('reactant', int(m)).get('aromatic'):
                    aromatic_attachment = True
                    break
            if aromatic_attachment:
                emit_specific('amine_to_diazonium', n, subtype='aryl_amine_to_diazonium', support=_support(n, *[_att_map(x) for x in p_n]))
    for c in sorted(maps):
        if atom('reactant', c).get('element') != 'C':
            continue
        r_hal = [x for x in atts('reactant', c) if x.get('element') in HALOGENS]
        p_n = [x for x in atts('product', c) if x.get('element') == 'N' and float(x.get('bond_order', 0.0)) >= 2.8]
        if len(r_hal) >= 2 and p_n:
            emit_specific('primary_amine_to_isocyanide', c, subtype='isocyanide_formation', support=_support(c, *[_att_map(x) for x in p_n]))
    for n in sorted(maps):
        if atom('reactant', n).get('element') != 'N':
            continue
        rc = [x for x in atts('reactant', n) if x.get('element') == 'C']
        pc = [x for x in atts('product', n) if x.get('element') == 'C']
        lc = [x for x in lost(n) if x.get('element') == 'C']
        if len(rc) >= 3 and len(pc) < len(rc) and lc:
            emit_specific('n_dealkylation', n, subtype='mapped_n_c_cleavage_dealkylation', support=_support(n, *[_att_map(x) for x in lc]), redox='reduction')
    for c in sorted(maps):
        if atom('reactant', c).get('element') != 'C':
            continue
        if fg('product', c) != 'aldehyde':
            continue
        reactant_double_o = any((x.get('element') == 'O' and float(x.get('bond_order', 0.0)) >= 1.8 for x in atts('reactant', c)))
        if reactant_double_o:
            continue
        reactant_n = any((x.get('element') == 'N' and float(x.get('bond_order', 0.0)) <= 1.2 for x in atts('reactant', c)))
        le = elems(lost(c))
        if reactant_n and 'N' in le and ('C' in le):
            emit_conversion('alpha_amino_acid_to_aldehyde', c, 'neutral', 'alpha_amino_acid_carbon', 'aldehyde')
    carbon_maps = sorted((m for m in maps if atom('reactant', m).get('element') == 'C' and atom('product', m).get('element') == 'C'))
    for (i, a) in enumerate(carbon_maps):
        for b in carbon_maps[i + 1:]:
            rb = bond('reactant', a, b)
            pb = bond('product', a, b)
            ga = gained(a)
            gb = gained(b)
            ga_o = [x for x in ga if x.get('element') == 'O']
            gb_o = [x for x in gb if x.get('element') == 'O']
            ga_n = [x for x in ga if x.get('element') == 'N']
            gb_n = [x for x in gb if x.get('element') == 'N']
            if rb >= 2.8 and 1.7 <= pb <= 2.2:

                def external_heavy_delta(center: int, partner: int) -> bool:
                    for x in lost(center) + gained(center):
                        el = str(x.get('element') or '')
                        if not el or el == 'H':
                            continue
                        m = _att_map(x)
                        if m is not None and int(m) == int(partner):
                            continue
                        return True
                    return False
                r_a = mapped_heavy_neighbors('reactant', a, exclude_maps={b})
                p_a = mapped_heavy_neighbors('product', a, exclude_maps={b})
                r_b = mapped_heavy_neighbors('reactant', b, exclude_maps={a})
                p_b = mapped_heavy_neighbors('product', b, exclude_maps={a})
                if hdelta(a) > 0 and hdelta(b) > 0 and (not external_heavy_delta(a, b)) and (not external_heavy_delta(b, a)) and (r_a == p_a) and (r_b == p_b):
                    emit_specific('partial_alkyne_hydrogenation', a, subtype='alkyne_to_alkene', support=_support(a, b), redox='reduction')
            if rb >= 2.8 and 0.8 <= pb <= 1.2:
                if fg('product', a) == 'ketone' or fg('product', b) == 'ketone':
                    kc = a if fg('product', a) == 'ketone' else b
                    emit_specific('alkyne_hydration_to_ketone', kc, subtype='alkyne_to_ketone', support=_support(a, b))
            if 1.7 <= rb <= 2.2 and 0.8 <= pb <= 1.2:
                if fg('product', a) == 'ketone' or fg('product', b) == 'ketone':
                    kc = a if fg('product', a) == 'ketone' else b
                    emit_specific('alkene_to_ketone_oxidation', kc, subtype='alkene_to_ketone', support=_support(a, b), redox='oxidation')
                if ga_n and gb_o or (gb_n and ga_o):
                    emit_specific('alkene_aminohydroxylation', a, subtype='vicinal_aminohydroxylation', support=_support(a, b, *[_att_map(x) for x in ga_n + gb_n + ga_o + gb_o]), redox='oxidation')
                total_o = ga_o + gb_o
                total_n = ga_n + gb_n
                single_o = [x for x in total_o if 0.8 <= float(x.get('bond_order', x.get('order', 0.0)) or 0.0) <= 1.2]
                if len(single_o) == 1 and (not total_n):
                    ox = single_o[0]
                    om = _att_map(ox)
                    hydration_oxygen = False
                    if om is None:
                        hydration_oxygen = True
                    else:
                        om = int(om)
                        r_ox = atom('reactant', om)
                        if not r_ox.get('element'):
                            hydration_oxygen = True
                        else:
                            r_ox_heavy = [y for y in atts('reactant', om) if y.get('element') not in {'', 'H'}]
                            hydration_oxygen = len(r_ox_heavy) == 0
                    if hydration_oxygen:
                        emit_specific('alkene_hydration', a if ga_o else b, subtype='alkene_to_alcohol', support=_support(a, b, *[_att_map(x) for x in single_o]))
    dearomatized = [m for m in maps if atom('reactant', m).get('element') == 'C' and atom('reactant', m).get('aromatic') and atom('product', m) and (not atom('product', m).get('aromatic'))]
    if len(dearomatized) >= 4:
        emit_specific('arene_partial_hydrogenation', dearomatized[0], subtype='partial_dearomatizing_hydrogenation', support=_support(*dearomatized), redox='reduction')
    for c in sorted(maps):
        if fg('reactant', c) not in {'aldehyde', 'ketone'}:
            continue
        ge = set(elems(gained(c)))
        if {'N', 'P'}.issubset(ge):
            emit_specific('aminophosphonate_multicomponent_addition', c, subtype='carbonyl_aminophosphonate_addition', support=_support(c, *[_att_map(x) for x in gained(c)]))
    for c in sorted(maps):
        if atom('reactant', c).get('element') != 'C':
            continue
        for x in atts('product', c):
            if x.get('element') != 'N' or float(x.get('bond_order', 0.0)) < 1.8:
                continue
            nm = _att_map(x)
            if nm is None or bond('reactant', c, int(nm)) > 0.1:
                continue
            nm = int(nm)
            r_o = [y for y in atts('reactant', nm) if y.get('element') == 'O']
            p_o = [y for y in atts('product', nm) if y.get('element') == 'O']
            genuine_nitroso = len(r_o) == 1 and float(r_o[0].get('bond_order', 0.0)) >= 1.8 and (not p_o)
            if genuine_nitroso:
                emit_specific('condensation', c, subtype='nitroso_condensation', support=_support(c, nm))

    def product_n_neighbors(nmap: int) -> set[int]:
        out_n = set()
        for x in atts('product', nmap):
            if x.get('element') != 'N':
                continue
            m = _att_map(x)
            if m is not None:
                out_n.add(int(m))
        return out_n
    for fb in facts.get('intercomponent_formed_bonds', []) or []:
        a = int(fb['a'])
        b = int(fb['b'])
        if atom('reactant', a).get('element') != 'N' or atom('reactant', b).get('element') != 'N':
            continue
        if atom('product', a).get('aromatic') or atom('product', b).get('aromatic'):
            continue
        nset = {a, b}
        nset.update(product_n_neighbors(a))
        nset.update(product_n_neighbors(b))
        if len(nset) != 3:
            continue
        if any((atom('product', n).get('element') != 'N' or atom('product', n).get('aromatic') for n in nset)):
            continue
        n_degrees = sorted((len(product_n_neighbors(n) & nset) for n in nset))
        if n_degrees != [1, 1, 2]:
            continue
        if not any((_is_azide_n(facts, 'product', n) for n in nset)):
            continue
        emit_specific('n_n_bond_formation', a, subtype='azide_n_n_formation', support=_support(*sorted(nset)))
    for c in sorted(maps):
        ra = atom('reactant', c)
        if ra.get('element') != 'C':
            continue
        if ra.get('aromatic'):
            continue
        if int(ra.get('total_h', 0) or 0) <= 0:
            continue
        if any((x.get('element') not in {None, '', 'H'} and float(x.get('bond_order', 0.0) or 0.0) >= 1.5 for x in atts('reactant', c))):
            continue
        go = [x for x in gained(c) if x.get('element') == 'O']
        if not go or hdelta(c) >= 0:
            continue
        allylic = False
        for n in atts('reactant', c):
            nm = _att_map(n)
            if nm is None or n.get('element') != 'C':
                continue
            for nn in atts('reactant', int(nm)):
                mm = _att_map(nn)
                if mm is None or int(mm) == c or nn.get('element') != 'C':
                    continue
                if bond('reactant', int(nm), int(mm)) >= 1.8:
                    allylic = True
                    break
            if allylic:
                break
        if allylic:
            emit_specific('allylic_oxidative_acyloxylation', c, subtype='allylic_c_o_formation', support=_support(c, *[_att_map(x) for x in go]), redox='oxidation')
    return out

def detect_directional_functional_group_extensions(facts: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    common = {int(x) for x in facts.get('mapping', {}).get('common_maps', []) or []}
    seen = set()

    def atom(side: str, mn: int) -> dict[str, Any]:
        return _atom(facts, side, int(mn))

    def atts(side: str, mn: int) -> list[dict[str, Any]]:
        return _attachments(facts, side, int(mn)) or []

    def bo(side: str, a: int, b: int) -> float:
        return float(_bond_between(facts, side, int(a), int(b)) or 0.0)

    def heavy(rows):
        return [x for x in rows if str(x.get('element') or '') not in {'', 'H'}]

    def mapped_heavy_neighbors(side: str, center: int, exclude=None) -> set[int]:
        exclude = {int(x) for x in exclude or []}
        ans = set()
        for x in atts(side, center):
            if x.get('element') == 'H':
                continue
            m = _att_map(x)
            if m is None:
                continue
            m = int(m)
            if m in exclude:
                continue
            ans.add(m)
        return ans

    def emit(directional_name: str, specific_event: str, center: int, support, source: str, target: str, *, direction: str='neutral', family: str='functional_group_interconversion', note: str):
        support = _support(*support)
        key = (directional_name, specific_event, int(center), tuple(sorted(support)))
        if key in seen:
            return
        seen.add(key)
        out.extend(_conversion_records(directional_name, int(center), support, direction, family, source, target, True, specific_event_name=specific_event, extra={'s2_5_directional_fgi': True, 'mapping_scope_note': note}))
    for bc in facts.get('bond_order_changes', []) or []:
        a = int(bc['a'])
        b = int(bc['b'])
        if atom('reactant', a).get('element') != 'C' or atom('reactant', b).get('element') != 'C':
            continue
        ro = float(bc.get('reactant_order', 0.0) or 0.0)
        po = float(bc.get('product_order', 0.0) or 0.0)
        if not (1.7 <= ro <= 2.2 and 0.8 <= po <= 1.2):
            continue
        for (center, partner) in ((a, b), (b, a)):
            (lost_c, gained_c) = _lost_gained_for_center(facts, center)
            (lost_p, gained_p) = _lost_gained_for_center(facts, partner)
            gained_o = [x for x in gained_c if x.get('element') == 'O' and 0.8 <= float(x.get('bond_order', 0.0) or 0.0) <= 1.2 and (_att_map(x) is not None)]
            if len(gained_o) != 1:
                continue
            om = int(_att_map(gained_o[0]))
            if om not in common:
                continue
            if not 0.8 <= bo('product', center, om) <= 1.2:
                continue
            p_o_heavy = heavy(atts('product', om))
            if len(p_o_heavy) != 1:
                continue
            p_o_map = _att_map(p_o_heavy[0])
            if p_o_map is None or int(p_o_map) != center:
                continue
            if atom('reactant', om).get('element') != 'O':
                continue
            r_o_heavy = heavy(atts('reactant', om))
            if not r_o_heavy:
                continue
            other_center_gain = [x for x in heavy(gained_c) if _att_map(x) != om]
            if heavy(lost_c) or heavy(lost_p) or other_center_gain or heavy(gained_p):
                continue
            if _hydrogen_delta(facts, partner) <= 0:
                continue
            pfg = _carbon_fg_state(facts, 'product', center)
            if pfg not in {'primary_alcohol', 'secondary_alcohol', 'tertiary_alcohol'}:
                continue
            emit('alkene_to_alcohol', 'alkene_to_alcohol_oxygenation', center, [a, b, om], 'alkene', pfg, note='mapped_c_c_double_to_single_plus_product_hydroxyl_oxygen_installation;oxygen_provenance_not_interpreted_mechanistically')
    for d in facts.get('attachment_deltas', []) or []:
        c = int(d['center_map'])
        if d.get('center_element') != 'C':
            continue
        ra = atom('reactant', c)
        pa = atom('product', c)
        if not (ra.get('aromatic') and pa.get('aromatic')):
            continue
        lost_h = heavy(d.get('lost', []) or [])
        gained_h = heavy(d.get('gained', []) or [])
        if len(lost_h) != 1 or lost_h[0].get('element') != 'N' or gained_h:
            continue
        if int(d.get('hydrogen_delta', 0) or 0) <= 0:
            continue
        nm = _att_map(lost_h[0])
        emit('aryl_amine_to_hydrogen', 'aryl_amine_deamination', c, [c, nm], 'aryl_or_heteroaryl_amine', 'aryl_or_heteroaryl_hydrogen', note='conserved_aromatic_carbon_loses_single_n_attachment_and_gains_h')
    for c in sorted(common):
        ca = atom('reactant', c)
        if ca.get('element') != 'C':
            continue
        if int(ca.get('total_h', 0) or 0) < 1:
            continue
        rn = [x for x in atts('reactant', c) if x.get('element') == 'N' and 1.7 <= float(x.get('bond_order', 0.0) or 0.0) <= 2.2 and (_att_map(x) is not None)]
        if len(rn) != 1:
            continue
        n = int(_att_map(rn[0]))
        if not bo('product', c, n) >= 2.8:
            continue
        if _carbon_fg_state(facts, 'product', c) != 'nitrile':
            continue
        r_n_o = [x for x in atts('reactant', n) if x.get('element') == 'O']
        p_n_o = [x for x in atts('product', n) if x.get('element') == 'O']
        if not r_n_o or p_n_o:
            continue
        r_skeleton = mapped_heavy_neighbors('reactant', c, exclude={n})
        p_skeleton = mapped_heavy_neighbors('product', c, exclude={n})
        if r_skeleton != p_skeleton:
            continue
        emit('oxime_to_nitrile', 'oxime_to_nitrile', c, [c, n, *[_att_map(x) for x in r_n_o]], 'aldoxime', 'nitrile', note='same_mapped_c_n_pair_changes_c_n_double_to_triple_while_n_loses_o')
    for bc in facts.get('bond_order_changes', []) or []:
        a = int(bc['a'])
        b = int(bc['b'])
        if atom('reactant', a).get('element') != 'C' or atom('reactant', b).get('element') != 'C':
            continue
        ro = float(bc.get('reactant_order', 0.0) or 0.0)
        po = float(bc.get('product_order', 0.0) or 0.0)
        if not (1.7 <= ro <= 2.2 and 0.8 <= po <= 1.2):
            continue
        for (c, partner) in ((a, b), (b, a)):
            pfg = _carbon_fg_state(facts, 'product', c)
            if pfg not in {'aldehyde', 'ketone'}:
                continue
            r_oxy = [x for x in atts('reactant', c) if x.get('element') == 'O' and 0.8 <= float(x.get('bond_order', 0.0) or 0.0) <= 1.2 and (_att_map(x) is not None)]
            for ox in r_oxy:
                om = int(_att_map(ox))
                if not bo('product', c, om) >= 1.8:
                    continue
                other_r_c = [x for x in atts('reactant', om) if x.get('element') == 'C' and (_att_map(x) is None or int(_att_map(x)) != c)]
                if not other_r_c:
                    continue
                p_o_heavy = heavy(atts('product', om))
                if len(p_o_heavy) != 1:
                    continue
                pm = _att_map(p_o_heavy[0])
                if pm is None or int(pm) != c:
                    continue
                emit('enol_ether_to_carbonyl', 'enol_ether_to_carbonyl', c, [c, partner, om], 'o_alkyl_enol_ether', pfg, note='mapped_enol_c_c_double_reduces_while_same_enol_oxygen_becomes_carbonyl_oxygen_and_loses_alkyl_context')
    for c in sorted(common):
        if _carbon_fg_state(facts, 'reactant', c) != 'ester' or _carbon_fg_state(facts, 'product', c) != 'ester':
            continue
        r_os = [x for x in _element_attachments(facts, 'reactant', c, 'O', 0.8, 1.2) if _att_map(x) is not None]
        for ro in r_os:
            om = int(_att_map(ro))
            if not 0.8 <= bo('product', c, om) <= 1.2:
                continue
            r_carbons = {int(_att_map(x)) for x in atts('reactant', om) if x.get('element') == 'C' and _att_map(x) is not None and (int(_att_map(x)) != c)}
            p_carbons = {int(_att_map(x)) for x in atts('product', om) if x.get('element') == 'C' and _att_map(x) is not None and (int(_att_map(x)) != c)}
            for alkyl_c in sorted(r_carbons & p_carbons):
                ra = atom('reactant', alkyl_c)
                pa = atom('product', alkyl_c)
                if ra.get('element') != 'C' or pa.get('element') != 'C':
                    continue
                if not (int(ra.get('total_h', 0) or 0) == 2 and int(pa.get('total_h', 0) or 0) >= 3):
                    continue
                (lost_k, gained_k) = _lost_gained_for_center(facts, alkyl_c)
                lost_carbons = [x for x in heavy(lost_k) if x.get('element') == 'C']
                if not lost_carbons or heavy(gained_k):
                    continue
                emit('ethyl_ester_to_methyl_ester', 'transesterification', c, [c, om, alkyl_c, *[_att_map(x) for x in lost_carbons]], 'ethyl_ester', 'methyl_ester', note='mapping_limited_same_ester_o_retained_while_alkoxy_ch2_becomes_ch3')
    return out

def detect_nitro_to_hydroxylamine(facts: dict[str, Any]) -> list[dict[str, Any]]:
    formed = facts.get('formed_bonds', []) or []
    broken = facts.get('broken_bonds', []) or []
    order_changed = facts.get('bond_order_changes', []) or []
    if formed or broken or len(order_changed) != 1:
        return []
    out: list[dict[str, Any]] = []
    for n in facts.get('mapping', {}).get('common_maps', []) or []:
        n = int(n)
        rn = _atom(facts, 'reactant', n)
        pn = _atom(facts, 'product', n)
        if rn.get('element') != 'N' or pn.get('element') != 'N':
            continue
        if int(rn.get('charge', 0) or 0) != 1:
            continue
        if int(pn.get('charge', 0) or 0) != 0:
            continue
        if int(rn.get('total_h', 0) or 0) != 0:
            continue
        if int(pn.get('total_h', 0) or 0) < 1:
            continue
        r_atts = list(_attachments(facts, 'reactant', n) or [])
        p_atts = list(_attachments(facts, 'product', n) or [])
        r_oxygens = [x for x in r_atts if x.get('element') == 'O']
        if len(r_oxygens) < 2:
            continue
        r_double_o = [x for x in r_oxygens if float(x.get('bond_order', 0.0) or 0.0) >= 1.8 and _att_map(x) is not None]
        if len(r_double_o) != 1:
            continue
        o = int(_att_map(r_double_o[0]))
        p_oxygens = [x for x in p_atts if x.get('element') == 'O' and _att_map(x) is not None and (int(_att_map(x)) == o) and (float(x.get('bond_order', 0.0) or 0.0) <= 1.2)]
        if len(p_oxygens) != 1:
            continue
        po = _atom(facts, 'product', o)
        if int(po.get('total_h', 0) or 0) < 1:
            continue
        r_c = {int(_att_map(x)) for x in r_atts if x.get('element') == 'C' and _att_map(x) is not None}
        p_c = {int(_att_map(x)) for x in p_atts if x.get('element') == 'C' and _att_map(x) is not None}
        retained_c = sorted(r_c & p_c)
        if not retained_c:
            continue
        ch = order_changed[0]
        try:
            pair = {int(ch['a']), int(ch['b'])}
        except Exception:
            continue
        if pair != {n, o}:
            continue
        if not (float(ch.get('reactant_order', 0.0) or 0.0) >= 1.8 and float(ch.get('product_order', 0.0) or 0.0) <= 1.2):
            continue
        c = retained_c[0]
        out.extend(_conversion_records('nitro_to_hydroxylamine', n, _support(n, o, c), 'reduction', 'reduction', 'nitro', 'hydroxylamine', True, specific_event_name='nitro_to_hydroxylamine', extra={'mapping_scope': 'mapped_nitro_n_o_to_n_hydroxylamine'}))
    return out
