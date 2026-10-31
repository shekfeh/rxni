"""Explicit, strictly verified SQLite foundation schema profiles for V2 libraries.

Only physical identifiers may vary.  Canonical structural hashes, reaction keys,
metadata semantics, representation format and ranking algorithms are frozen.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
import sqlite3
from typing import Any

CANONICAL_STRUCTURAL = (
    'forward_edit_signature_hash', 'forward_raw_center_hash_r0',
    'forward_raw_center_hash_r1', 'forward_raw_center_hash_r2',
    'forward_component_topology_hash', 'forward_component_evidence_hash',
)
CANONICAL_METADATA_REQUIRED = ('reaction_key', 'forward_reaction_smiles')
CANONICAL_METADATA_OPTIONAL = (
    'reaction_id', 'retro_template_hash', 'retro_template_number',
    'membership_occurrence_count',
)
STRUCTURAL_REQUIRED = ('reaction_key', *CANONICAL_STRUCTURAL)
IDENTIFIER = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')

class FoundationSchemaError(ValueError):
    """Source has an unsupported physical index or metadata schema."""


def quote_ident(value: str) -> str:
    if not isinstance(value, str) or not IDENTIFIER.fullmatch(value):
        raise FoundationSchemaError(f'Invalid SQLite identifier: {value!r}')
    return '"' + value + '"'


def _exact_keys(mapping: dict, required: set, allowed: set, label: str) -> None:
    missing, extra = required - set(mapping), set(mapping) - allowed
    if missing or extra:
        raise FoundationSchemaError(f'{label}: missing={sorted(missing)} unknown={sorted(extra)}')


@dataclass(frozen=True)
class FoundationProfile:
    path: Path
    profile_id: str
    structural_table: str
    structural_columns: dict[str, str]
    metadata_table: str
    metadata_columns: dict[str, str]


def load_foundation_profile(path: str | Path) -> FoundationProfile:
    path = Path(path).expanduser().resolve()
    try:
        obj = json.loads(path.read_text(encoding='utf-8'))
    except (ValueError, OSError) as exc:
        raise FoundationSchemaError(f'Invalid foundation profile {path}: {exc}') from exc
    if not isinstance(obj, dict):
        raise FoundationSchemaError('Foundation profile must be a JSON object')
    _exact_keys(obj, {'schema_version','profile_id','structural_index','metadata_lookup'},
                {'schema_version','profile_id','structural_index','metadata_lookup'}, 'profile')
    if type(obj['schema_version']) is not int or obj['schema_version'] != 1:
        raise FoundationSchemaError('Expected foundation schema_version=1')
    profile_id = obj['profile_id']
    if not isinstance(profile_id, str) or not profile_id.strip():
        raise FoundationSchemaError('Nonempty profile_id required')
    parsed = {}
    for role, required, allowed in (
        ('structural_index', set(STRUCTURAL_REQUIRED), set(STRUCTURAL_REQUIRED)),
        ('metadata_lookup', set(CANONICAL_METADATA_REQUIRED),
         set(CANONICAL_METADATA_REQUIRED + CANONICAL_METADATA_OPTIONAL)),
    ):
        block = obj[role]
        if not isinstance(block, dict):
            raise FoundationSchemaError(f'{role} must be an object')
        _exact_keys(block, {'table','columns'}, {'table','columns'}, role)
        quote_ident(block['table'])
        mapping = block['columns']
        if not isinstance(mapping, dict):
            raise FoundationSchemaError(f'{role}.columns must be an object')
        _exact_keys(mapping, required, allowed, role+'.columns')
        for canonical, physical in mapping.items():
            quote_ident(physical)
        if len(set(mapping.values())) != len(mapping):
            raise FoundationSchemaError(f'{role}: duplicate physical column names')
        parsed[role] = (block['table'], mapping)
    return FoundationProfile(path, profile_id,
                             parsed['structural_index'][0], parsed['structural_index'][1],
                             parsed['metadata_lookup'][0], parsed['metadata_lookup'][1])


def _columns(conn: sqlite3.Connection, table: str) -> dict[str, str]:
    q = quote_ident(table)
    names = [x[0] for x in conn.execute("SELECT name FROM sqlite_master WHERE type IN ('table','view')")]
    if table not in names:
        raise FoundationSchemaError(f'SQLite table/view absent: {table}')
    return {str(row[1]):str(row[2]) for row in conn.execute('PRAGMA table_info('+q+')')}


def _check_table(path: Path, table: str, mapping: dict[str,str], role: str) -> dict[str,Any]:
    if not path.is_file():
        raise FoundationSchemaError(f'{role} database missing: {path}')
    conn = sqlite3.connect(f'file:{path}?mode=ro', uri=True)
    try:
        cols = _columns(conn, table)
        missing = set(mapping.values()) - set(cols)
        if missing:
            raise FoundationSchemaError(f'{role}.{table}: physical columns absent: {sorted(missing)}')
        # Queries are prepared/compiled; this also guards badly quoted or partial views.
        select = ','.join(quote_ident(x) for x in mapping.values())
        conn.execute('SELECT '+select+' FROM '+quote_ident(table)+' LIMIT 0').fetchall()
        return {'table':table,'columns':dict(mapping),'physical_types':{k:cols[v] for k,v in mapping.items()}}
    finally:
        conn.close()


def validate_foundation_schema(structural_db: str|Path, metadata_db: str|Path,
                               profile: FoundationProfile) -> dict[str, Any]:
    """Fail closed on schema mismatch; do not infer chemical equivalence."""
    structural = _check_table(Path(structural_db), profile.structural_table,
                              profile.structural_columns, 'structural_index')
    metadata = _check_table(Path(metadata_db), profile.metadata_table,
                            profile.metadata_columns, 'metadata_lookup')
    return {'status':'FOUNDATION_SCHEMA_VALID','profile_id':profile.profile_id,
            'structural_index':structural,'metadata_lookup':metadata,
            'key_semantics':'EXACT_TEXT_IDENTITY_ONLY',
            'scientific_equivalence':'NOT_ASSERTED'}


def audit_key_membership(structural_db: str|Path, metadata_db: str|Path,
                         locator_db: str|Path, profile: FoundationProfile, *,
                         sample_size: int = 32) -> dict[str, Any]:
    """A bounded structural→metadata→V2-locator probe, not full population proof."""
    if sample_size < 1 or sample_size > 10000:
        raise FoundationSchemaError('sample_size must be 1..10000')
    validate_foundation_schema(structural_db, metadata_db, profile)
    s = sqlite3.connect(f'file:{Path(structural_db)}?mode=ro', uri=True)
    m = sqlite3.connect(f'file:{Path(metadata_db)}?mode=ro', uri=True)
    l = sqlite3.connect(f'file:{Path(locator_db)}?mode=ro', uri=True)
    try:
        # Refuse to pretend the fixed locator ABI is configurable.
        loc_cols = _columns(l, 'representation_locator_v2')
        required = {'reaction_key','build_id','representation_generation','row_group_id','row_in_group'}
        if not required.issubset(loc_cols):
            raise FoundationSchemaError('V2 locator ABI missing: '+str(sorted(required-set(loc_cols))))
        keys = [str(r[0]) for r in s.execute(
            'SELECT DISTINCT '+quote_ident(profile.structural_columns['reaction_key'])+
            ' FROM '+quote_ident(profile.structural_table)+' WHERE '+
            quote_ident(profile.structural_columns['reaction_key'])+
            ' IS NOT NULL LIMIT ?', (int(sample_size),)).fetchall()]
        if not keys:
            raise FoundationSchemaError('Structural index has no non-null reaction keys to validate')
        missing_meta=[]; missing_locator=[]
        metadata_sql = ('SELECT 1 FROM '+quote_ident(profile.metadata_table)+' WHERE '+
                        quote_ident(profile.metadata_columns['reaction_key'])+'=? LIMIT 1')
        for key in keys:
            if m.execute(metadata_sql, (key,)).fetchone() is None:missing_meta.append(key)
            if l.execute('SELECT 1 FROM representation_locator_v2 WHERE reaction_key=? LIMIT 1',
                         (key,)).fetchone() is None:missing_locator.append(key)
        return {'status':'PASS' if not missing_meta and not missing_locator else 'KEY_COVERAGE_MISMATCH',
                'sample_size':len(keys),'missing_metadata':missing_meta,
                'missing_v2_locator':missing_locator,
                'scope':'BOUNDED_SAMPLE_NOT_FULL_MEMBERSHIP_PROOF'}
    finally:
        s.close();m.close();l.close()
