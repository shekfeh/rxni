from __future__ import annotations

import hashlib
from collections import Counter
from typing import Any, Iterable

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from .retrieval_mapped_graph import build_mapped_graph_facts, split_reaction, canonical_bond_class

FINGERPRINT_VERSION = 'i_fp_0.1-batch05'
REACTION_CENTER_VERSION = 'i_center_0.1-batch05'


def _mol_from_side(side: str) -> Chem.Mol | None:
    if not str(side).strip():
        return None
    return Chem.MolFromSmiles(str(side))


def _map_to_idx(mol: Chem.Mol | None) -> dict[int, int]:
    if mol is None:
        return {}
    out: dict[int, int] = {}
    for atom in mol.GetAtoms():
        mn = int(atom.GetAtomMapNum())
        if mn > 0:
            out[mn] = int(atom.GetIdx())
    return out


def _count_payload(fp) -> dict[str, list[int]]:
    nz = fp.GetNonzeroElements()
    idx = sorted(int(i) for i in nz)
    return {'indices': idx, 'counts': [int(nz[i]) for i in idx]}


def _morgan_count(mol: Chem.Mol | None, *, radius: int, fp_size: int,
                  root_maps: Iterable[int] | None = None) -> dict[str, list[int]]:
    if mol is None:
        return {'indices': [], 'counts': []}
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=int(radius), fpSize=int(fp_size), includeChirality=True)
    kwargs: dict[str, Any] = {}
    if root_maps is not None:
        m2i = _map_to_idx(mol)
        roots = sorted({m2i[int(m)] for m in root_maps if int(m) in m2i})
        if not roots:
            return {'indices': [], 'counts': []}
        kwargs['fromAtoms'] = roots
    return _count_payload(gen.GetCountFingerprint(mol, **kwargs))


def _stable_bucket(feature: str, fp_size: int) -> int:
    digest = hashlib.sha256(feature.encode('utf-8')).digest()
    return int.from_bytes(digest[:8], 'big') % int(fp_size)


def _edit_features(facts: dict[str, Any], allowed_maps: set[int] | None = None) -> list[str]:
    feats: list[str] = []

    def allowed(*maps: int) -> bool:
        return allowed_maps is None or bool(set(int(m) for m in maps) & allowed_maps)

    for rec in facts.get('formed_bonds', []):
        a, b = int(rec['a']), int(rec['b'])
        if not allowed(a, b):
            continue
        ea = facts['reactant']['atoms'].get(a, {}).get('element', '?')
        eb = facts['reactant']['atoms'].get(b, {}).get('element', '?')
        feats.append(f"FORM:{canonical_bond_class(ea, eb)}:{float(rec['order']):.1f}")

    for rec in facts.get('broken_bonds', []):
        a, b = int(rec['a']), int(rec['b'])
        if not allowed(a, b):
            continue
        ea = facts['reactant']['atoms'].get(a, {}).get('element', '?')
        eb = facts['reactant']['atoms'].get(b, {}).get('element', '?')
        feats.append(f"BREAK:{canonical_bond_class(ea, eb)}:{float(rec['order']):.1f}")

    for rec in facts.get('bond_order_changes', []):
        a, b = int(rec['a']), int(rec['b'])
        if not allowed(a, b):
            continue
        ea = facts['reactant']['atoms'].get(a, {}).get('element', '?')
        eb = facts['reactant']['atoms'].get(b, {}).get('element', '?')
        feats.append(
            f"ORDER:{canonical_bond_class(ea, eb)}:{float(rec['reactant_order']):.1f}>{float(rec['product_order']):.1f}"
        )

    for rec in facts.get('attachment_deltas', []):
        c = int(rec['center_map'])
        if not allowed(c):
            continue
        cel = rec.get('center_element', '?')
        lost = sorted(f"{x.get('element','?')}:{float(x.get('bond_order',0.0)):.1f}:{'M' if x.get('map') is not None else 'U'}"
                      for x in rec.get('lost', []))
        gained = sorted(f"{x.get('element','?')}:{float(x.get('bond_order',0.0)):.1f}:{'M' if x.get('map') is not None else 'U'}"
                        for x in rec.get('gained', []))
        feats.append(f"ATTACH:{cel}:L[{','.join(lost)}]:G[{','.join(gained)}]:H{int(rec.get('hydrogen_delta',0)):+d}")

    if allowed_maps is None and int(facts.get('common_cycle_rank_delta', 0)):
        feats.append(f"CYCLE_RANK_DELTA:{int(facts['common_cycle_rank_delta']):+d}")
    return feats


def _hashed_feature_counts(features: Iterable[str], fp_size: int) -> dict[str, list[int]]:
    counts = Counter(_stable_bucket(f, fp_size) for f in features)
    idx = sorted(counts)
    return {'indices': idx, 'counts': [int(counts[i]) for i in idx]}


class _DSU:
    def __init__(self, nodes: Iterable[int]):
        self.parent = {int(n): int(n) for n in nodes}

    def add(self, n: int):
        n = int(n)
        self.parent.setdefault(n, n)

    def find(self, n: int) -> int:
        n = int(n)
        self.add(n)
        while self.parent[n] != n:
            self.parent[n] = self.parent[self.parent[n]]
            n = self.parent[n]
        return n

    def union(self, a: int, b: int):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def extract_reaction_centers(reaction: str, classification: dict[str, Any] | None = None) -> dict[str, Any]:
    """Return a variable-length set of disconnected mapped reaction centers.

    A center is a connected *edit complex*, not a single atom. There is intentionally
    no hard maximum: multi-event/one-pot records may contain two or more disconnected
    centers. Raw edit endpoints establish seeds; validated semantic support maps may
    enlarge those centers without creating chemistry from semantics alone.
    """
    facts = build_mapped_graph_facts(reaction)
    raw_seed: set[int] = set()
    semantic_maps: set[int] = set()
    groups: list[set[int]] = []

    for key in ('formed_bonds', 'broken_bonds', 'bond_order_changes'):
        for rec in facts.get(key, []):
            g = {int(rec['a']), int(rec['b'])}
            raw_seed.update(g)
            groups.append(g)

    for rec in facts.get('attachment_deltas', []):
        g = {int(rec['center_map'])}
        raw_seed.add(int(rec['center_map']))
        for x in list(rec.get('lost', [])) + list(rec.get('gained', [])):
            if x.get('map') is not None:
                g.add(int(x['map']))
                raw_seed.add(int(x['map']))
        groups.append(g)

    if classification is None:
        # Lazy import avoids a resolver -> fingerprints import cycle.
        from .semantic_adapter import classify_for_retrieval
        classification = classify_for_retrieval(reaction)

    for key in ('primary_events', 'specific_events', 'directional_conversions'):
        for rec in classification.get(key, []):
            g = set(int(m) for m in (rec.get('support_maps') or []) if int(m) in facts['mapping']['common_maps'])
            if rec.get('center_map') is not None and int(rec['center_map']) in facts['mapping']['common_maps']:
                g.add(int(rec['center_map']))
            g.update(int(m) for m in (rec.get('acyl_center_maps') or []) if int(m) in facts['mapping']['common_maps'])
            if g:
                semantic_maps.update(g)
                groups.append(g)

    nodes = raw_seed | semantic_maps
    if not nodes:
        return {
            'version': REACTION_CENTER_VERSION,
            'n_centers': 0,
            'centers': [],
            'raw_seed_maps': [],
            'mapping_valid': bool(facts['mapping']['valid']),
        }

    dsu = _DSU(nodes)
    for g in groups:
        g = sorted(g & nodes)
        if len(g) > 1:
            first = g[0]
            for m in g[1:]:
                dsu.union(first, m)

    # Direct adjacency in either side joins local edit seeds into one center complex.
    for side in ('reactant', 'product'):
        for a, b in facts[side]['bonds']:
            if a in nodes and b in nodes:
                dsu.union(a, b)

    comps: dict[int, set[int]] = {}
    for m in sorted(nodes):
        comps.setdefault(dsu.find(m), set()).add(m)

    centers = []
    for i, maps in enumerate(sorted(comps.values(), key=lambda x: min(x))):
        maps_sorted = sorted(maps)
        centers.append({
            'center_id': f'rc{i}',
            'maps': maps_sorted,
            'raw_seed_maps': sorted(set(maps_sorted) & raw_seed),
            'semantic_support_maps': sorted(set(maps_sorted) & semantic_maps),
            'reactant_components': sorted({facts['reactant']['map_to_component'].get(m) for m in maps_sorted
                                           if facts['reactant']['map_to_component'].get(m) is not None}),
            'product_components': sorted({facts['product']['map_to_component'].get(m) for m in maps_sorted
                                          if facts['product']['map_to_component'].get(m) is not None}),
        })

    return {
        'version': REACTION_CENTER_VERSION,
        'n_centers': len(centers),
        'centers': centers,
        'raw_seed_maps': sorted(raw_seed),
        'mapping_valid': bool(facts['mapping']['valid']),
    }


def generate_localized_fingerprints(reaction: str, *, classification: dict[str, Any] | None = None,
                                    radius: int = 2, fp_size: int = 2048) -> dict[str, Any]:
    """Per-reaction-center rooted Morgan + directional edit fingerprints."""
    facts = build_mapped_graph_facts(reaction)
    centers = extract_reaction_centers(reaction, classification=classification)
    lhs, _agents, rhs = split_reaction(reaction)
    rmol, pmol = _mol_from_side(lhs), _mol_from_side(rhs)

    out = []
    for center in centers['centers']:
        maps = set(center['maps'])
        out.append({
            **center,
            'channels': {
                'reactant_local': _morgan_count(rmol, radius=radius, fp_size=fp_size, root_maps=maps),
                'product_local': _morgan_count(pmol, radius=radius, fp_size=fp_size, root_maps=maps),
                'edit_local': _hashed_feature_counts(_edit_features(facts, maps), fp_size),
            },
        })
    return {
        'fingerprint_version': FINGERPRINT_VERSION,
        'reaction_center_version': REACTION_CENTER_VERSION,
        'radius': int(radius),
        'fp_size': int(fp_size),
        'encoding': 'folded_sparse_count',
        'n_centers': len(out),
        'centers': out,
    }


def generate_multichannel_fingerprint(reaction: str, *, classification: dict[str, Any] | None = None,
                                      radius: int = 2, fp_size: int = 2048) -> dict[str, Any]:
    """First System-I structural multi-channel reaction fingerprint.

    Channels intentionally remain separate for later ablation/weighting:
      reactant_global, product_global, reactant_local, product_local, edit.
    No semantic labels are hashed into the structural fingerprint; semantics stay an
    independent retrieval/routing axis.
    """
    facts = build_mapped_graph_facts(reaction)
    if classification is None:
        from .semantic_adapter import classify_for_retrieval
        classification = classify_for_retrieval(reaction)
    center_info = extract_reaction_centers(reaction, classification=classification)
    center_maps = {m for c in center_info['centers'] for m in c['maps']}
    lhs, _agents, rhs = split_reaction(reaction)
    rmol, pmol = _mol_from_side(lhs), _mol_from_side(rhs)

    return {
        'fingerprint_version': FINGERPRINT_VERSION,
        'reaction_center_version': REACTION_CENTER_VERSION,
        'radius': int(radius),
        'fp_size': int(fp_size),
        'encoding': 'folded_sparse_count',
        'n_centers': int(center_info['n_centers']),
        'center_maps': sorted(center_maps),
        'channels': {
            'reactant_global': _morgan_count(rmol, radius=radius, fp_size=fp_size),
            'product_global': _morgan_count(pmol, radius=radius, fp_size=fp_size),
            'reactant_local': _morgan_count(rmol, radius=radius, fp_size=fp_size, root_maps=center_maps),
            'product_local': _morgan_count(pmol, radius=radius, fp_size=fp_size, root_maps=center_maps),
            'edit': _hashed_feature_counts(_edit_features(facts), fp_size),
        },
    }
