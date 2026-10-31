"""Explicit runtime profile with independent semantic/schema/scoring lineage."""
from dataclasses import dataclass, asdict

@dataclass(frozen=True)
class Profile:
    name: str
    representation_format: str
    semantic_policy: str
    scoring_policy: str
    v2_module: str
    expected_runtime_files: tuple[str, ...]

V67=Profile(
    name='v6.7',
    representation_format='retrieval_representations_v2',
    semantic_policy='runtime_v6_7 semantic query profile (library-side compatibility MUST be verified)',
    scoring_policy='frozen v6.6 ranking algorithm',
    v2_module='rxni.retrieval.builders.representations',
    expected_runtime_files=(
        'structural_index_v1.sqlite', 'runtime_lookup_v1.sqlite',
        'retrieval_representations_v2.parquet','representation_locator_v2.sqlite',
        'retrieval_v2_manifest.json'
    )
)
def select(name):
    if name!='v6.7':
        raise ValueError('Only v6.7 auxiliary is implemented; no silent downgrade to another runtime')
    return V67
