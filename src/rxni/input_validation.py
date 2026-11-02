"""Strict public query mapping preflight; does not alter stored scientific contracts."""
from __future__ import annotations

from collections import Counter
from typing import Any
from rdkit import Chem


class ReactionMappingError(ValueError):
    def __init__(self, code: str, message: str, diagnostics: dict[str, Any]):
        self.code=code
        self.message=message
        self.diagnostics=diagnostics
        super().__init__(f"{code}: {message}")


def _reject(code: str, message: str, **details: Any) -> None:
    raise ReactionMappingError(code, message, details)


def _parse_side(raw: str, side: str):
    if not raw:
        _reject('invalid_reaction_syntax', f'{side} side is empty', side=side)
    molecules = []
    for index, fragment in enumerate(raw.split('.')):
        if not fragment:
            _reject('invalid_reaction_syntax', f'empty molecule in {side}', side=side, index=index)
        mol = Chem.MolFromSmiles(fragment, sanitize=False)
        if mol is None:
            _reject('invalid_reaction_smiles', f'unparseable molecule in {side}', side=side, index=index)
        molecules.append(mol)
    return molecules


def inspect_mapping(reaction: str) -> dict[str, Any]:
    """Validate map coverage on reactants/products; agents are deliberately exempt.

    All heavy atoms on each reaction side must have positive, unique map numbers.
    Lost/gained atoms are permitted, as are unmapped agents, but at least two
    map numbers must be shared by reactants and products. Map numbers may not
    change atomic element across the reaction. Atom-map quality does not prove
    chemical correctness of the proposed transformation.
    """
    if not isinstance(reaction, str) or reaction.count('>') != 2:
        _reject('invalid_reaction_syntax', 'expected reactants>agents>products')
    reactants, agents, products = reaction.split('>')
    sides = {}
    maps = {}
    for name, field in [('reactants', reactants), ('products', products)]:
        molecules = _parse_side(field, name)
        atoms = [a for m in molecules for a in m.GetAtoms() if a.GetAtomicNum() != 1]
        nums = [a.GetAtomMapNum() for a in atoms]
        missing = sum(n <= 0 for n in nums)
        if missing:
            _reject('incomplete_atom_mapping', f'{missing} unmapped heavy atom(s) in {name}', side=name, missing_count=missing, heavy_atoms=len(atoms))
        dupl = sorted(k for k,v in Counter(nums).items() if v > 1)
        if dupl:
            _reject('duplicate_atom_mapping', f'duplicate map labels in {name}', side=name, duplicate_labels=dupl[:20])
        sides[name] = len(atoms)
        maps[name] = {a.GetAtomMapNum(): a.GetAtomicNum() for a in atoms}
    # Agent SMILES may be absent or unmapped, but syntax must remain parseable.
    if agents:
        for i, part in enumerate(agents.split('.')):
            if not part or Chem.MolFromSmiles(part, sanitize=False) is None:
                _reject('invalid_reaction_smiles', 'invalid agent component', index=i)
    shared = set(maps['reactants']) & set(maps['products'])
    if len(shared) < 2:
        _reject('insufficient_shared_mapping', 'at least two mapped atoms must persist across sides', shared_count=len(shared))
    changed = sorted(k for k in shared if maps['reactants'][k] != maps['products'][k])
    if changed:
        _reject('inconsistent_atom_mapping', 'mapped atom changes chemical element', labels=changed[:20])
    return {'status':'mapped_complete', 'policy':'strict_heavy_atom_mapping_v1',
            'reactant_heavy_atoms':sides['reactants'], 'product_heavy_atoms':sides['products'],
            'shared_mapped_atoms':len(shared), 'agents_mapping_required':False}


def require_mapped_reaction(reaction: str) -> dict[str, Any]:
    return inspect_mapping(reaction)
