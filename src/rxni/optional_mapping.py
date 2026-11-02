"""Opt-in, externally executed RXNMapper adapter. Never an implicit dependency."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import sys

from rdkit import Chem
from rxni.input_validation import require_mapped_reaction, ReactionMappingError


class MappingBackendError(RuntimeError):
    pass


def _canonical_side(field: str) -> tuple[str, ...]:
    out = []
    for part in field.split('.'):
        if not part:
            raise MappingBackendError('empty component')
        mol = Chem.MolFromSmiles(part)
        if mol is None:
            raise MappingBackendError('invalid component')
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        out.append(Chem.MolToSmiles(mol, isomericSmiles=True))
    return tuple(sorted(out))


def map_external(reaction: str, *, mapper_python: str, timeout: int = 90) -> tuple[str, dict]:
    """Run opt-in worker; verify chemistry unchanged and all product/reactant atoms mapped.

    mapper_python is administrator-controlled configuration, never HTTP request data.
    """
    executable = Path(mapper_python).expanduser()
    if not executable.is_file():
        raise MappingBackendError(f'configured mapper interpreter not found: {executable}')
    worker = Path(__file__).with_name('_rxnmapper_worker.py')
    try:
        proc = subprocess.run([str(executable), str(worker)],
                              input=json.dumps({'reaction':reaction}),text=True,
                              capture_output=True,timeout=timeout,check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise MappingBackendError(f'external mapper process failed: {type(exc).__name__}') from exc
    if proc.returncode:
        raise MappingBackendError(f'external mapper returned status {proc.returncode}: {proc.stderr[-400:]}')
    try:
        payload = json.loads(proc.stdout)
        mapped = payload['mapped_reaction']
        if not isinstance(mapped, str): raise TypeError('mapped_reaction not text')
        original_sides = reaction.split('>')
        mapped_sides = mapped.split('>')
        if len(original_sides) != 3 or len(mapped_sides) != 3:
            raise MappingBackendError('invalid mapped reaction shape')
        for k in (0, 2):
            if _canonical_side(original_sides[k]) != _canonical_side(mapped_sides[k]):
                raise MappingBackendError('mapping changed reactant/product molecular identity')
        # Agents must remain unchanged; they are not a mapping requirement.
        if original_sides[1] != mapped_sides[1]:
            raise MappingBackendError('mapping changed the reaction agent field')
        qc = require_mapped_reaction(mapped)
    except (ValueError, TypeError, KeyError, ReactionMappingError) as exc:
        raise MappingBackendError(f'external mapper output rejected: {exc}') from exc
    return mapped, {'backend':'external_rxnmapper', 'automapped':True,
                    'mapper_python':executable.name, 'qc':qc,
                    'mapped_input_verified_against_original':True}


def prepare_query(reaction: str) -> tuple[str, dict]:
    """By default reject invalid mapping. An explicit administrator setting enables mapping."""
    try:
        qc = require_mapped_reaction(reaction)
        return reaction, {'backend':'provided_by_user','automapped':False,'qc':qc}
    except ReactionMappingError:
        interpreter = os.environ.get('RXNI_MAPPER_PYTHON')
        if not interpreter:
            raise
        return map_external(reaction, mapper_python=interpreter)


def main(argv=None):
    import argparse
    parser=argparse.ArgumentParser(description='Explicit opt-in external RXNMapper mapping')
    parser.add_argument('--reaction',required=True)
    parser.add_argument('--mapper-python',required=True,help='Python with rxnmapper installed')
    options=parser.parse_args(argv)
    mapped,provenance=map_external(options.reaction,mapper_python=options.mapper_python)
    print(json.dumps({'mapped_reaction':mapped,'mapping_provenance':provenance},indent=2))

if __name__=='__main__': main()
