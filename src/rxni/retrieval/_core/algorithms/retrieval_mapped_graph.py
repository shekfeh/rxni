from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from typing import Any, Iterable

from rdkit import Chem


def split_reaction(reaction: str) -> tuple[str, str, str]:
    parts = str(reaction).split('>')
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    if '>>' in reaction:
        lhs, rhs = reaction.split('>>', 1)
        return lhs, '', rhs
    raise ValueError("Reaction must be reactants>>products or reactants>agents>products")




_BOND_ELEMENT_PRIORITY = {
    'C': 0, 'N': 1, 'O': 2, 'S': 3, 'P': 4,
    'B': 5, 'Si': 6,
    'F': 7, 'Cl': 8, 'Br': 9, 'I': 10,
}


def canonical_bond_class(a: str, b: str) -> str:
    """Stable chemistry-oriented bond-class label independent of atom ordering.

    Carbon is placed first when present; common heteroatoms follow; halogens are
    placed after the atom they functionalize. Unknown elements fall back to
    lexical ordering after the known element set.
    """
    def key(e: str) -> tuple[int, str]:
        return (_BOND_ELEMENT_PRIORITY.get(e, 100), e)
    x, y = sorted((str(a), str(b)), key=key)
    return f"{x}-{y}"

def bond_order(bond: Chem.Bond) -> float:
    if bond.GetIsAromatic():
        return 1.5
    return float(bond.GetBondTypeAsDouble())


def _atom_state(atom: Chem.Atom) -> dict[str, Any]:
    return {
        'map': int(atom.GetAtomMapNum()),
        'element': atom.GetSymbol(),
        'atomic_num': int(atom.GetAtomicNum()),
        'charge': int(atom.GetFormalCharge()),
        'aromatic': bool(atom.GetIsAromatic()),
        'hybridization': str(atom.GetHybridization()),
        'degree': int(atom.GetDegree()),
        'heavy_degree': int(sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)),
        'total_h': int(atom.GetTotalNumHs(includeNeighbors=True)),
        'in_ring': bool(atom.IsInRing()),
    }


def _branch_signature(anchor: Chem.Atom, nbr: Chem.Atom, bond: Chem.Bond) -> dict[str, Any]:
    """One-hop signature for an attachment neighbor.

    The signature is used only as supporting evidence for a conserved mapped center.
    Unmapped atoms never become independent event centers.
    """
    others = []
    for b in nbr.GetBonds():
        other = b.GetOtherAtom(nbr)
        if other.GetIdx() == anchor.GetIdx():
            continue
        others.append((other.GetSymbol(), round(bond_order(b), 2), bool(other.GetIsAromatic())))
    others.sort()
    return {
        'map': int(nbr.GetAtomMapNum()) if nbr.GetAtomMapNum() > 0 else None,
        'element': nbr.GetSymbol(),
        'bond_order': round(bond_order(bond), 2),
        'aromatic': bool(nbr.GetIsAromatic()),
        'charge': int(nbr.GetFormalCharge()),
        'degree': int(nbr.GetDegree()),
        'in_ring': bool(nbr.IsInRing()),
        'neighbor_shell': others,
        'carbonyl_like': bool(
            nbr.GetSymbol() == 'C'
            and any(
                b.GetOtherAtom(nbr).GetSymbol() == 'O' and bond_order(b) >= 1.9
                for b in nbr.GetBonds()
                if b.GetOtherAtom(nbr).GetIdx() != anchor.GetIdx()
            )
        ),
        'carboxyl_like': bool(
            nbr.GetSymbol() == 'C'
            and sum(
                1
                for b in nbr.GetBonds()
                if b.GetOtherAtom(nbr).GetSymbol() == 'O'
                and b.GetOtherAtom(nbr).GetIdx() != anchor.GetIdx()
            ) >= 2
        ),
        'sulfonyl_like': bool(
            nbr.GetSymbol() == 'S'
            and sum(
                1
                for b in nbr.GetBonds()
                if b.GetOtherAtom(nbr).GetSymbol() == 'O' and bond_order(b) >= 1.9
            ) >= 2
        ),
    }


def _attachment_identity(rec: dict[str, Any]) -> tuple[Any, ...]:
    # Mapped neighbor identity is its map number. Bond-order changes are tracked
    # separately and therefore do not masquerade as replacement.
    if rec['map'] is not None:
        return ('mapped', int(rec['map']), rec['element'])
    # Unmapped neighbors are discriminated by a compact local signature.
    return (
        'unmapped',
        rec['element'],
        rec['bond_order'],
        rec['aromatic'],
        rec['charge'],
        rec['degree'],
        rec['carbonyl_like'],
        rec['carboxyl_like'],
        rec['sulfonyl_like'],
        tuple(tuple(x) for x in rec['neighbor_shell']),
    )


def _parse_side(side: str) -> list[dict[str, Any]]:
    components: list[dict[str, Any]] = []
    for comp_id, smiles in enumerate(x for x in str(side).split('.') if x):
        mol = Chem.MolFromSmiles(smiles)
        components.append({'component_id': comp_id, 'smiles': smiles, 'mol': mol})
    return components


def _side_record(side: str) -> dict[str, Any]:
    components = _parse_side(side)
    atoms: dict[int, dict[str, Any]] = {}
    map_to_component: dict[int, int] = {}
    bonds: dict[tuple[int, int], float] = {}
    attachments: dict[int, list[dict[str, Any]]] = defaultdict(list)
    rings: list[dict[str, Any]] = []
    duplicate_maps: list[int] = []
    parse_failures: list[str] = []

    for comp in components:
        ci, smiles, mol = comp['component_id'], comp['smiles'], comp['mol']
        if mol is None:
            parse_failures.append(smiles)
            continue
        for atom in mol.GetAtoms():
            mn = int(atom.GetAtomMapNum())
            if mn <= 0:
                continue
            if mn in atoms:
                duplicate_maps.append(mn)
            atoms[mn] = _atom_state(atom)
            map_to_component[mn] = ci
        for bond in mol.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()
            ma, mb = int(a.GetAtomMapNum()), int(b.GetAtomMapNum())
            bo = bond_order(bond)
            if ma > 0 and mb > 0:
                bonds[tuple(sorted((ma, mb)))] = bo
            if ma > 0:
                attachments[ma].append(_branch_signature(a, b, bond))
            if mb > 0:
                attachments[mb].append(_branch_signature(b, a, bond))

        for atom_ring in mol.GetRingInfo().AtomRings():
            ring_atoms = []
            for idx in atom_ring:
                atom = mol.GetAtomWithIdx(idx)
                ring_atoms.append({
                    'idx': int(idx),
                    'map': int(atom.GetAtomMapNum()) if atom.GetAtomMapNum() > 0 else None,
                    'element': atom.GetSymbol(),
                })
            rings.append({'component_id': ci, 'atoms': ring_atoms})

    return {
        'components': components,
        'atoms': atoms,
        'map_to_component': map_to_component,
        'bonds': bonds,
        'attachments': dict(attachments),
        'rings': rings,
        'duplicate_maps': sorted(set(duplicate_maps)),
        'parse_failures': parse_failures,
    }


def _graph_components(nodes: set[int], bonds: dict[tuple[int, int], float]) -> int:
    if not nodes:
        return 0
    adj = {n: set() for n in nodes}
    for a, b in bonds:
        if a in nodes and b in nodes:
            adj[a].add(b)
            adj[b].add(a)
    seen: set[int] = set()
    n_components = 0
    for node in nodes:
        if node in seen:
            continue
        n_components += 1
        stack = [node]
        seen.add(node)
        while stack:
            cur = stack.pop()
            for nbr in adj[cur]:
                if nbr not in seen:
                    seen.add(nbr)
                    stack.append(nbr)
    return n_components


def cycle_rank(nodes: Iterable[int], bonds: dict[tuple[int, int], float]) -> int:
    nodes = set(nodes)
    if not nodes:
        return 0
    edge_count = sum(1 for a, b in bonds if a in nodes and b in nodes)
    return edge_count - len(nodes) + _graph_components(nodes, bonds)


def _ring_subfamily(ring: dict[str, Any], side: dict[str, Any]) -> str:
    elems = [a['element'] for a in ring['atoms']]
    n = len(elems)
    comp = side['components'][ring['component_id']]
    mol = comp['mol']
    idxs = {a['idx'] for a in ring['atoms']}
    if n == 3:
        if 'O' in elems:
            return 'epoxide_ring_opening'
        if 'N' in elems:
            return 'aziridine_ring_opening'
        if 'S' in elems:
            return 'thiirane_ring_opening'
        return 'carbocycle_ring_opening'

    if mol is not None:
        for idx in idxs:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != 'C':
                continue
            has_carbonyl_o = any(
                b.GetOtherAtom(atom).GetSymbol() == 'O' and bond_order(b) >= 1.9
                for b in atom.GetBonds()
            )
            if not has_carbonyl_o:
                continue
            for b in atom.GetBonds():
                nbr = b.GetOtherAtom(atom)
                if nbr.GetIdx() not in idxs or bond_order(b) >= 1.9:
                    continue
                if nbr.GetSymbol() == 'O':
                    return 'lactone_ring_opening'
                if nbr.GetSymbol() == 'N':
                    return 'lactam_ring_opening'

    if 'O' in elems:
        return 'ether_ring_opening'
    if 'N' in elems:
        return 'amine_ring_opening'
    return 'carbocycle_ring_opening'


def _ring_lineage(reactant: dict[str, Any], product: dict[str, Any], common_maps: set[int]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    opened: list[dict[str, Any]] = []
    formed: list[dict[str, Any]] = []

    for ring in reactant['rings']:
        mapped = [a['map'] for a in ring['atoms'] if a['map'] is not None]
        if len(mapped) < 3:
            continue
        common = [m for m in mapped if m in common_maps]
        if len(common) < 3:
            continue
        conserved_fraction = len(common) / len(mapped)
        if conserved_fraction < 0.60:
            continue
        rr = cycle_rank(common, reactant['bonds'])
        pr = cycle_rank(common, product['bonds'])
        if rr > 0 and pr == 0:
            opened.append({
                'ring_maps': mapped,
                'conserved_maps': common,
                'conserved_fraction': conserved_fraction,
                'subfamily': _ring_subfamily(ring, reactant),
                'reactant_cycle_rank': rr,
                'product_cycle_rank': pr,
            })

    for ring in product['rings']:
        mapped = [a['map'] for a in ring['atoms'] if a['map'] is not None]
        if len(mapped) < 3:
            continue
        common = [m for m in mapped if m in common_maps]
        if len(common) < 3:
            continue
        conserved_fraction = len(common) / len(mapped)
        if conserved_fraction < 0.60:
            continue
        rr = cycle_rank(common, reactant['bonds'])
        pr = cycle_rank(common, product['bonds'])
        if rr == 0 and pr > 0:
            formed.append({
                'ring_maps': mapped,
                'conserved_maps': common,
                'conserved_fraction': conserved_fraction,
                'reactant_cycle_rank': rr,
                'product_cycle_rank': pr,
            })
    return opened, formed


def _attachment_deltas(reactant: dict[str, Any], product: dict[str, Any], common_maps: set[int]) -> list[dict[str, Any]]:
    out = []
    for mn in sorted(common_maps):
        ra = reactant['attachments'].get(mn, [])
        pa = product['attachments'].get(mn, [])
        rc = Counter(_attachment_identity(x) for x in ra)
        pc = Counter(_attachment_identity(x) for x in pa)
        lost_keys = list((rc - pc).elements())
        gained_keys = list((pc - rc).elements())
        if not lost_keys and not gained_keys:
            continue
        r_by = defaultdict(list)
        p_by = defaultdict(list)
        for x in ra:
            r_by[_attachment_identity(x)].append(x)
        for x in pa:
            p_by[_attachment_identity(x)].append(x)
        lost = [r_by[k].pop() for k in lost_keys]
        gained = [p_by[k].pop() for k in gained_keys]
        out.append({
            'center_map': mn,
            'center_element': reactant['atoms'][mn]['element'],
            'lost': lost,
            'gained': gained,
            'hydrogen_delta': product['atoms'][mn]['total_h'] - reactant['atoms'][mn]['total_h'],
        })
    return out


def build_mapped_graph_facts(reaction: str) -> dict[str, Any]:
    lhs, agents, rhs = split_reaction(reaction)
    reactant = _side_record(lhs)
    product = _side_record(rhs)

    r_maps = set(reactant['atoms'])
    p_maps = set(product['atoms'])
    common = r_maps & p_maps
    r_only = r_maps - p_maps
    p_only = p_maps - r_maps

    formed = []
    broken = []
    order_changes = []
    for key in sorted(set(product['bonds']) - set(reactant['bonds'])):
        if key[0] in common and key[1] in common:
            formed.append({'a': key[0], 'b': key[1], 'order': product['bonds'][key]})
    for key in sorted(set(reactant['bonds']) - set(product['bonds'])):
        if key[0] in common and key[1] in common:
            broken.append({'a': key[0], 'b': key[1], 'order': reactant['bonds'][key]})
    for key in sorted(set(reactant['bonds']) & set(product['bonds'])):
        if key[0] not in common or key[1] not in common:
            continue
        rbo, pbo = reactant['bonds'][key], product['bonds'][key]
        if abs(rbo - pbo) > 0.20:
            order_changes.append({'a': key[0], 'b': key[1], 'reactant_order': rbo, 'product_order': pbo})

    inter_formed = []
    intra_formed = []
    for rec in formed:
        a, b = rec['a'], rec['b']
        same_origin = reactant['map_to_component'].get(a) == reactant['map_to_component'].get(b)
        x = dict(rec)
        x['reactant_component_a'] = reactant['map_to_component'].get(a)
        x['reactant_component_b'] = reactant['map_to_component'].get(b)
        if same_origin:
            intra_formed.append(x)
        else:
            inter_formed.append(x)

    component_splits = []
    for rci in sorted(set(reactant['map_to_component'].values())):
        maps = [m for m in common if reactant['map_to_component'].get(m) == rci]
        product_components = sorted({product['map_to_component'].get(m) for m in maps})
        if len(product_components) > 1:
            component_splits.append({
                'reactant_component': rci,
                'product_components': product_components,
                'maps': sorted(maps),
            })

    component_merges = []
    for pci in sorted(set(product['map_to_component'].values())):
        maps = [m for m in common if product['map_to_component'].get(m) == pci]
        reactant_components = sorted({reactant['map_to_component'].get(m) for m in maps})
        if len(reactant_components) > 1:
            component_merges.append({
                'product_component': pci,
                'reactant_components': reactant_components,
                'maps': sorted(maps),
            })

    ring_openings, ring_formations = _ring_lineage(reactant, product, common)
    attachment_deltas = _attachment_deltas(reactant, product, common)

    cr_r = cycle_rank(common, reactant['bonds'])
    cr_p = cycle_rank(common, product['bonds'])

    formed_bond_classes = []
    for rec in formed:
        ea = reactant['atoms'][rec['a']]['element']
        eb = reactant['atoms'][rec['b']]['element']
        formed_bond_classes.append(canonical_bond_class(ea, eb))

    mapping_valid = not (
        reactant['duplicate_maps'] or product['duplicate_maps'] or reactant['parse_failures'] or product['parse_failures']
    )
    element_mismatches = []
    for mn in sorted(common):
        if reactant['atoms'][mn]['element'] != product['atoms'][mn]['element']:
            element_mismatches.append({
                'map': mn,
                'reactant_element': reactant['atoms'][mn]['element'],
                'product_element': product['atoms'][mn]['element'],
            })
    if element_mismatches:
        mapping_valid = False

    return {
        'reaction': reaction,
        'mapping': {
            'valid': mapping_valid,
            'common_maps': sorted(common),
            'reactant_only_maps': sorted(r_only),
            'product_only_maps': sorted(p_only),
            'duplicate_maps_reactant': reactant['duplicate_maps'],
            'duplicate_maps_product': product['duplicate_maps'],
            'element_mismatches': element_mismatches,
            'parse_failures_reactant': reactant['parse_failures'],
            'parse_failures_product': product['parse_failures'],
        },
        'reactant': reactant,
        'product': product,
        'formed_bonds': formed,
        'broken_bonds': broken,
        'bond_order_changes': order_changes,
        'intercomponent_formed_bonds': inter_formed,
        'intracomponent_formed_bonds': intra_formed,
        'formed_bond_classes': sorted(set(formed_bond_classes)),
        'attachment_deltas': attachment_deltas,
        'component_splits': component_splits,
        'component_merges': component_merges,
        'common_cycle_rank_reactant': cr_r,
        'common_cycle_rank_product': cr_p,
        'common_cycle_rank_delta': cr_p - cr_r,
        'ring_openings': ring_openings,
        'ring_formations': ring_formations,
    }
