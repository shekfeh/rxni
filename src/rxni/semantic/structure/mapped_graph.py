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


RING_LINEAGE_VERSION = "ring_lineage_v2.0"


def _ring_kind(ring: dict[str, Any], side: dict[str, Any]) -> str:
    """Return a structural ring identity independent of reaction name or conditions."""
    elems = [a['element'] for a in ring['atoms']]
    n = len(elems)
    counts = Counter(elems)

    if n == 3:
        if counts == Counter({'C': 2, 'O': 1}):
            return 'epoxide'
        if counts == Counter({'C': 2, 'N': 1}):
            return 'aziridine'
        if counts == Counter({'C': 2, 'S': 1}):
            return 'thiirane'
        if counts == Counter({'C': 3}):
            return 'cyclopropane'

    if n == 4:
        if counts == Counter({'C': 3, 'O': 1}):
            return 'oxetane'
        if counts == Counter({'C': 3, 'N': 1}):
            comp = side['components'][ring['component_id']]
            mol = comp['mol']
            idxs = {int(a['idx']) for a in ring['atoms']}
            if mol is not None:
                for idx in idxs:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() != 'C':
                        continue
                    if any(
                        b.GetOtherAtom(atom).GetSymbol() == 'O' and bond_order(b) >= 1.9
                        for b in atom.GetBonds()
                    ):
                        return 'beta_lactam'
            return 'azetidine'
        if counts == Counter({'C': 4}):
            return 'cyclobutane'

    if any(e not in {'C', 'H'} for e in elems):
        return 'heterocycle'
    return 'carbocycle'


def _ring_descriptor(
    ring: dict[str, Any], side: dict[str, Any], common_maps: set[int]
) -> dict[str, Any]:
    atoms = list(ring['atoms'])
    maps = [int(a['map']) for a in atoms if a.get('map') is not None]
    common = sorted(set(maps) & set(common_maps))
    edges: set[tuple[int, int]] = set()
    for i, atom in enumerate(atoms):
        nxt = atoms[(i + 1) % len(atoms)]
        ma, mb = atom.get('map'), nxt.get('map')
        if ma is not None and mb is not None:
            edges.add(tuple(sorted((int(ma), int(mb)))))
    return {
        'component_id': int(ring['component_id']),
        'size': len(atoms),
        'maps': sorted(set(maps)),
        'common_maps': common,
        'elements': [str(a['element']) for a in atoms],
        'element_counts': dict(sorted(Counter(str(a['element']) for a in atoms).items())),
        'mapped_edges': [list(x) for x in sorted(edges)],
        'kind': _ring_kind(ring, side),
        'subfamily': _ring_subfamily(ring, side),
    }


def _whole_side_cycle_rank(side: dict[str, Any]) -> int:
    """Cycle-space rank over full parsed molecular graphs, including unmapped atoms."""
    total = 0
    for comp in side.get('components', []) or []:
        mol = comp.get('mol')
        if mol is None:
            continue
        total += int(mol.GetNumBonds()) - int(mol.GetNumAtoms()) + 1
    return total


def _ring_lineage_v2(
    reactant: dict[str, Any],
    product: dict[str, Any],
    common_maps: set[int],
    formed_bonds: list[dict[str, Any]],
    broken_bonds: list[dict[str, Any]],
    legacy_openings: list[dict[str, Any]],
    legacy_formations: list[dict[str, Any]],
) -> dict[str, Any]:
    """Insertion/deletion-aware ring lineage and local ring-architecture facts.

    The legacy ring lineage operates on the conserved mapped induced graph.  This
    complementary lineage compares full reactant/product RDKit rings while using mapped
    atoms only for cross-side correspondence.  It therefore recognizes atom insertion,
    atom deletion, same-rank ring replacement, and small-ring formation involving a
    product-only atom.  Highly polycyclic systems suppress global SSSR birth/death
    interpretation; specific local small-ring construction remains available.
    """
    rr = [_ring_descriptor(r, reactant, common_maps) for r in reactant.get('rings', []) or []]
    pp = [_ring_descriptor(r, product, common_maps) for r in product.get('rings', []) or []]
    whole_r = _whole_side_cycle_rank(reactant)
    whole_p = _whole_side_cycle_rank(product)
    whole_delta = whole_p - whole_r
    common_r_rank = cycle_rank(common_maps, reactant['bonds'])
    common_p_rank = cycle_rank(common_maps, product['bonds'])
    common_delta = common_p_rank - common_r_rank
    polycyclic_ambiguous = len(rr) > 12 or len(pp) > 12

    candidates: list[tuple[float, int, float, int, int]] = []
    if not polycyclic_ambiguous:
        for i, r in enumerate(rr):
            rs = set(r['common_maps'])
            if len(rs) < 2:
                continue
            for j, p in enumerate(pp):
                ps = set(p['common_maps'])
                if len(ps) < 2:
                    continue
                inter = rs & ps
                n = len(inter)
                if n < 2:
                    continue
                contain = n / max(1, min(len(rs), len(ps)))
                jacc = n / max(1, len(rs | ps))
                # Size changes can legitimately insert/delete a ring atom. Same-size
                # rewiring needs stronger overlap to avoid interpreting an annulation
                # that merely shares an edge with a retained ring as ring replacement.
                size_change_ok = (
                    whole_delta == 0
                    and r['size'] != p['size']
                    and n >= 3
                    and contain >= 0.80
                    and jacc >= 0.75
                )
                same_size_ok = r['size'] == p['size'] and n >= 4 and jacc >= 0.65
                exact_lineage = set(r['maps']) == set(p['maps']) and len(r['maps']) >= 3
                if not (size_change_ok or same_size_ok or exact_lineage):
                    continue
                score = 3.0 * contain + 2.0 * jacc + 0.10 * n - 0.05 * abs(r['size'] - p['size'])
                candidates.append((score, n, jacc, i, j))

    candidates.sort(reverse=True)
    used_r: set[int] = set()
    used_p: set[int] = set()
    matches: list[dict[str, Any]] = []
    expansions: list[dict[str, Any]] = []
    contractions: list[dict[str, Any]] = []
    rewirings: list[dict[str, Any]] = []

    for score, n, jacc, i, j in candidates:
        if i in used_r or j in used_p:
            continue
        r, p = rr[i], pp[j]
        rs, ps = set(r['common_maps']), set(p['common_maps'])
        inter = sorted(rs & ps)
        contain = len(inter) / max(1, min(len(rs), len(ps)))
        used_r.add(i)
        used_p.add(j)
        rec = {
            'reactant_ring_index': i,
            'product_ring_index': j,
            'reactant_ring_size': int(r['size']),
            'product_ring_size': int(p['size']),
            'reactant_ring_maps': list(r['maps']),
            'product_ring_maps': list(p['maps']),
            'conserved_maps': inter,
            'inserted_maps': sorted(set(p['maps']) - set(r['maps'])),
            'lost_maps': sorted(set(r['maps']) - set(p['maps'])),
            'reactant_kind': r['kind'],
            'product_kind': p['kind'],
            'overlap_fraction': round(float(contain), 6),
            'jaccard_fraction': round(float(jacc), 6),
            'support_maps': sorted(set(r['maps']) | set(p['maps'])),
        }
        matches.append(rec)
        if p['size'] > r['size']:
            expansions.append(dict(rec, event='ring_expansion'))
        elif p['size'] < r['size']:
            contractions.append(dict(rec, event='ring_contraction'))
        else:
            r_maps, p_maps = set(r['maps']), set(p['maps'])
            # Same atom membership with only mapped adjacency/bond-order permutation is
            # handled by rearrangement/pi-reorganization, not by ring rewiring.
            if (
                r_maps != p_maps
                and (r.get('kind') != p.get('kind') or r.get('element_counts') != p.get('element_counts'))
            ):
                rewirings.append(dict(rec, event='ring_rewiring'))

    # A validated rewiring is a zero-net ring-lineage replacement, not a concurrent
    # formation/opening or size-change architecture. Filter after all lineages exist.

    legacy_form_sets = [set(x.get('ring_maps', []) or []) for x in legacy_formations]
    legacy_open_sets = [set(x.get('ring_maps', []) or []) for x in legacy_openings]
    common_set = set(common_maps)

    retained_product_indices = {
        int(m['product_ring_index'])
        for m in matches
        if m['reactant_ring_size'] == m['product_ring_size']
        and set(m['reactant_ring_maps']) == set(m['product_ring_maps'])
    }

    raw_formations: list[dict[str, Any]] = []
    for j, p in enumerate(pp):
        if j in used_p:
            continue
        p_maps = set(p['maps'])
        p_common = p_maps & common_set
        formed_inside = sorted(
            tuple(sorted((int(b['a']), int(b['b']))))
            for b in formed_bonds
            if int(b['a']) in p_maps and int(b['b']) in p_maps
        )
        product_only = sorted(p_maps - common_set)
        legacy = any(p_maps == x for x in legacy_form_sets)
        specialized_small = p['kind'] in {
            'epoxide', 'aziridine', 'thiirane', 'cyclopropane', 'oxetane',
            'azetidine', 'beta_lactam', 'cyclobutane',
        }
        # v0.3 closure: a single product-only mapped atom can be the missing
        # ring atom itself.  Because it belongs to an RDKit ring, it necessarily
        # bridges two ring neighbours; requiring a positive full-graph cycle-rank
        # delta and >=2 conserved ring atoms prevents ordinary product-only
        # substituents from becoming ring construction.  This closes mapped
        # Paal-Knorr/Bogert-Cook type insertions without reintroducing the old
        # conserved-map ring false positives.
        single_inserted_ring_atom_closure = bool(
            whole_delta > 0
            and len(product_only) == 1
            and len(p_common) >= 2
        )
        # A second mapping-recovery regime occurs in fused/polycyclic products:
        # one reactant ring is retained exactly, while another newly independent
        # carbocycle contains several product-only mapped atoms because those atoms
        # were present but left unmapped on the reactant side.  Positive full-graph
        # cycle rank plus an exact retained-ring anchor makes this a bounded
        # retained-scaffold polycycle recovery rather than a generic product-only
        # ring heuristic (Bally-Scholl-like mapping).
        retained_scaffold_polycyclic_closure = bool(
            whole_delta > 0
            and p['kind'] == 'carbocycle'
            and len(product_only) >= 2
            and len(p_common) >= 2
            and bool(retained_product_indices)
        )
        validated = bool(
            legacy
            or (specialized_small and len(p_common) >= 2)
            or (formed_inside and len(p_common) >= 2)
            or single_inserted_ring_atom_closure
            or retained_scaffold_polycyclic_closure
        )
        if not validated:
            continue
        if polycyclic_ambiguous and not specialized_small:
            continue

        fused = False
        spiro = False
        p_edges = {tuple(x) for x in p['mapped_edges']}
        for k in retained_product_indices:
            kept = pp[k]
            overlap = p_maps & set(kept['maps'])
            if len(overlap) == 1:
                spiro = True
            kept_edges = {tuple(x) for x in kept['mapped_edges']}
            if p_edges & kept_edges:
                fused = True

        single_inserted_ring_atom_recovery = bool(
            single_inserted_ring_atom_closure
            and not legacy
            and not specialized_small
            and not formed_inside
        )
        retained_scaffold_polycyclic_recovery = bool(
            retained_scaffold_polycyclic_closure
            and not legacy
            and not specialized_small
            and not formed_inside
            and not single_inserted_ring_atom_recovery
        )
        formation_rec = {
            'product_ring_index': j,
            'ring_maps': list(p['maps']),
            'common_maps': sorted(p_common),
            'product_only_maps': product_only,
            'ring_size': int(p['size']),
            'kind': p['kind'],
            'subfamily': p['subfamily'],
            'elements': list(p['elements']),
            'formed_ring_bonds': [list(x) for x in formed_inside],
            'fused_to_retained_ring': bool(fused),
            'spiro_to_retained_ring': bool(spiro and not fused),
            'support_maps': list(p['maps']),
            'mapping_scope': (
                'full_ring_lineage_product_ring_single_inserted_atom_validation'
                if single_inserted_ring_atom_recovery
                else (
                    'full_ring_lineage_retained_scaffold_polycyclic_recovery'
                    if retained_scaffold_polycyclic_recovery
                    else 'full_ring_lineage_product_ring_validation'
                )
            ),
        }
        if single_inserted_ring_atom_recovery:
            formation_rec['single_inserted_ring_atom_closure'] = True
        if retained_scaffold_polycyclic_recovery:
            formation_rec['retained_scaffold_polycyclic_closure'] = True
        raw_formations.append(formation_rec)

    # One new closure can generate several alternative SSSR rings. Collapse those
    # descriptions to the smallest ring for the same local closure/new-atom cause.
    form_groups: dict[tuple[Any, ...], dict[str, Any]] = {}
    for rec in raw_formations:
        if rec['formed_ring_bonds']:
            key = ('formed_bonds', tuple(tuple(x) for x in rec['formed_ring_bonds']))
        elif rec['product_only_maps']:
            key = ('product_only', tuple(rec['product_only_maps']))
        else:
            key = ('ring_maps', tuple(rec['ring_maps']))
        old = form_groups.get(key)
        if old is None or (rec['ring_size'], -len(rec['common_maps'])) < (old['ring_size'], -len(old['common_maps'])):
            form_groups[key] = rec
    formations = list(form_groups.values())

    raw_openings: list[dict[str, Any]] = []
    for i, r in enumerate(rr):
        if i in used_r:
            continue
        r_maps = set(r['maps'])
        r_common = r_maps & common_set
        broken_inside = sorted(
            tuple(sorted((int(b['a']), int(b['b']))))
            for b in broken_bonds
            if int(b['a']) in r_maps and int(b['b']) in r_maps
        )
        reactant_only = sorted(r_maps - common_set)
        legacy = any(r_maps == x for x in legacy_open_sets)
        specialized_small = r['kind'] in {
            'epoxide', 'aziridine', 'thiirane', 'cyclopropane', 'oxetane',
            'azetidine', 'beta_lactam', 'cyclobutane',
        }
        validated = bool(
            legacy
            or (specialized_small and len(r_common) >= 2)
            or (broken_inside and len(r_common) >= 2)
        )
        if not validated:
            continue
        if polycyclic_ambiguous and not specialized_small:
            continue
        raw_openings.append({
            'reactant_ring_index': i,
            'ring_maps': list(r['maps']),
            'common_maps': sorted(r_common),
            'reactant_only_maps': reactant_only,
            'ring_size': int(r['size']),
            'kind': r['kind'],
            'subfamily': r['subfamily'],
            'elements': list(r['elements']),
            'broken_ring_bonds': [list(x) for x in broken_inside],
            'support_maps': list(r['maps']),
            'mapping_scope': 'full_ring_lineage_reactant_ring_validation',
        })

    open_groups: dict[tuple[Any, ...], dict[str, Any]] = {}
    for rec in raw_openings:
        if rec['broken_ring_bonds']:
            key = ('broken_bonds', tuple(tuple(x) for x in rec['broken_ring_bonds']))
        elif rec['reactant_only_maps']:
            key = ('reactant_only', tuple(rec['reactant_only_maps']))
        else:
            key = ('ring_maps', tuple(rec['ring_maps']))
        old = open_groups.get(key)
        if old is None or (rec['ring_size'], -len(rec['common_maps'])) < (old['ring_size'], -len(old['common_maps'])):
            open_groups[key] = rec
    openings = list(open_groups.values())

    # v0.3 mapped-core annulation recovery.  Some benchmark reactions retain all
    # mapped atoms but contain unmapped atoms in a reactant aromatic ring.  Fusion
    # into a fully mapped product ring can therefore increase the cycle rank of the
    # conserved mapped core while leaving the full molecular cycle rank unchanged.
    # Treat this only as annulation evidence under a deliberately narrow contract:
    #   * zero net full-graph cycle-rank change;
    #   * positive conserved-core cycle-rank change;
    #   * no mapped atom insertion/deletion across the whole reaction;
    #   * at least two formed mapped bonds;
    #   * same-size ring lineage with an incompletely mapped reactant ring, a fully
    #     mapped product ring, no lost mapped ring atoms, and >=2 newly recruited
    #     mapped ring atoms.
    # This is intended for Pschorr/Scholl-like fused-ring mapping recovery and is
    # kept separate from generic ring formation.
    mapped_core_annulations: list[dict[str, Any]] = []
    mapped_sets_conserved = (
        set(reactant.get('atoms', {})) == set(common_maps)
        and set(product.get('atoms', {})) == set(common_maps)
    )
    if (
        whole_delta == 0
        and common_delta > 0
        and mapped_sets_conserved
        and len(formed_bonds) >= 2
        and not polycyclic_ambiguous
    ):
        for rec in matches:
            if rec['reactant_ring_size'] != rec['product_ring_size']:
                continue
            if rec.get('reactant_kind') != rec.get('product_kind'):
                continue
            if rec.get('lost_maps'):
                continue
            if len(rec.get('inserted_maps') or []) < 2:
                continue
            if len(rec.get('reactant_ring_maps') or []) >= int(rec['reactant_ring_size']):
                continue
            if len(rec.get('product_ring_maps') or []) != int(rec['product_ring_size']):
                continue
            product_ring_maps = set(rec.get('product_ring_maps') or [])
            formed_inside = sorted(
                tuple(sorted((int(b['a']), int(b['b']))))
                for b in formed_bonds
                if int(b['a']) in product_ring_maps and int(b['b']) in product_ring_maps
            )
            if len(formed_inside) < 2:
                continue
            mapped_core_annulations.append({
                'reactant_ring_index': rec['reactant_ring_index'],
                'product_ring_index': rec['product_ring_index'],
                'reactant_ring_maps': list(rec['reactant_ring_maps']),
                'product_ring_maps': list(rec['product_ring_maps']),
                'ring_maps': list(rec['product_ring_maps']),
                'ring_size': int(rec['product_ring_size']),
                'kind': rec['product_kind'],
                'inserted_maps': list(rec.get('inserted_maps') or []),
                'formed_ring_bonds': [list(x) for x in formed_inside],
                'support_maps': sorted(
                    set(rec.get('reactant_ring_maps') or [])
                    | set(rec.get('product_ring_maps') or [])
                ),
                'mapping_scope': 'mapped_core_fusion_annulation_zero_full_rank_change',
            })

    if whole_delta != 0 or expansions or contractions or formations or openings:
        rewirings = []

    result = {
        'version': RING_LINEAGE_VERSION,
        'whole_cycle_rank_reactant': int(whole_r),
        'whole_cycle_rank_product': int(whole_p),
        'whole_cycle_rank_delta': int(whole_delta),
        'polycyclic_sssr_ambiguous': bool(polycyclic_ambiguous),
        'matches': matches,
        'formations': formations,
        'openings': openings,
        'expansions': expansions,
        'contractions': contractions,
        'rewirings': rewirings,
        'reactant_rings': rr,
        'product_rings': pp,
    }
    if mapped_core_annulations:
        result['mapped_core_annulations'] = mapped_core_annulations
    return result


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
    ring_lineage_v2 = _ring_lineage_v2(
        reactant, product, common, formed, broken, ring_openings, ring_formations
    )
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
        'ring_lineage_v2': ring_lineage_v2,
    }


def build_ring_lineage_v2(reaction: str) -> dict[str, Any]:
    """Build only the v6.7 full-ring lineage extension for a mapped reaction.

    This helper exists so the non-FGI compatibility boundary can preserve the
    frozen v6.6 base graph facts while adding v6.7 ring-lineage semantics.
    """
    return build_mapped_graph_facts(reaction).get("ring_lineage_v2", {})
