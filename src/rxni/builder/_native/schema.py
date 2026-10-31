"""Schema/constants for the native Runtime-v6.5 Template Library v5 builder.

This builder intentionally treats pre-v4 source files only as sources of
reaction/template identity and provenance. All structural/System-I derived
fields are recomputed under Runtime-v6.5.
"""

from __future__ import annotations

BUILDER_VERSION = "builder_v5_native_0.1"
LIBRARY_SCHEMA_VERSION = "template_library_v5_native_0.1"
REACTION_KEY_VERSION = "reaction_key_sha256_reaction_text_v1"
STRUCTURAL_KEY_VERSION = "runtime_v6_5_query_structural_keys_v1"
SYSTEM_I_REPRESENTATION_VERSION = "system_i_retrieval_record_v1_tool92_93"
LOCATOR_VERSION = "representation_locator_v1"
EDIT_INDEX_VERSION = "edit_index_v1"
STRUCTURAL_INDEX_VERSION = "structural_index_v1"

STAGING_DB = "builder_staging_v5.sqlite"
REACTION_CORE_BASE_FILENAME = "reaction_core_base_v5.parquet"
REACTION_CORE_FILENAME = "reaction_core_v5.parquet"
STRUCTURAL_KEYS_FILENAME = "reaction_structural_keys_v5.parquet"
MEMBERSHIP_REL_FILENAME = "reaction_template_membership_v5.parquet"
PROVENANCE_FILENAME = "reaction_provenance_v5.parquet"
RETRIEVAL_MEMBERS_FILENAME = "retrieval_members_v5.parquet"
REPRESENTATIONS_FILENAME = "system_i_representations_v1.parquet"
REPRESENTATION_LOCATOR_DB = "representation_locator_v1.sqlite"
REPRESENTATION_LOCATOR_PARQUET = "representation_locator_v1.parquet"
REPRESENTATION_ERRORS_FILENAME = "system_i_representation_errors_v1.jsonl"
TEMPLATE_LIBRARY_FILENAME = "template_library_v5.parquet"
EDIT_INDEX_FILENAME = "edit_index_v1.parquet"
EDIT_INDEX_DB = "edit_index_v1.sqlite"
STRUCTURAL_INDEX_DB = "structural_index_v1.sqlite"
VALIDATION_FILENAME = "validation_v5.json"
BUILD_MANIFEST_FILENAME = "build_manifest_v5.json"
SOURCE_NORMALIZE_STATS = "normalize_stats.json"
STRUCTURAL_STATS = "structural_key_stats.json"
REPRESENTATION_STATS = "representation_stats.json"
TEMPLATE_STATS = "template_stats.json"
EDIT_INDEX_STATS = "edit_index_stats.json"

STRUCTURAL_HASH_FIELDS = (
    "forward_edit_signature_hash",
    "forward_raw_center_hash_r0",
    "forward_raw_center_hash_r1",
    "forward_raw_center_hash_r2",
    "forward_component_topology_hash",
    "forward_component_evidence_hash",
)

RAW_CENTER_FIELDS = (
    "forward_raw_center_hash_r0",
    "forward_raw_center_hash_r1",
    "forward_raw_center_hash_r2",
)

# These aliases are deliberately broad for old pre-v4 sources.
MEMBERSHIP_ALIASES = {
    "reaction": (
        "forward_reaction_smiles",
        "reaction_smiles",
        "mapped_reaction_smiles",
        "reaction",
        "rxn_smiles",
    ),
    "template_hash": (
        "retro_template_hash",
        "template_hash",
        "reaction_template_hash",
    ),
    "template_number": (
        "retro_template_number",
        "template_number",
        "template_id",
    ),
    "reaction_id": (
        "reaction_id",
        "member_id",
        "source_reaction_id",
        "id",
    ),
    "reaction_hash": (
        "forward_reaction_hash",
        "reaction_hash",
    ),
}

TEMPLATE_SOURCE_ALIASES = {
    "template_hash": (
        "retro_template_hash",
        "template_hash",
        "reaction_template_hash",
    ),
    "template_number": (
        "retro_template_number",
        "template_number",
        "template_id",
    ),
    "template_smarts": (
        "retro_template_smarts",
        "retro_template",
        "template_smarts",
        "reaction_smarts",
        "smarts",
    ),
}

# Input columns matching these patterns are never accepted as authoritative
# derived chemistry. They may remain in provenance, but they are ignored for
# construction of V5 derived data.
FORBIDDEN_DERIVED_TOKENS = (
    "raw_center",
    "edit_signature",
    "component_topology_hash",
    "component_evidence_hash",
    "fingerprint",
    "semantic_",
    "i_forward_",
    "i_raw_hash",
    "i_edit_hash",
    "i_component_hash",
    "probability_class",
    "reaction_class",
    "reaction_family",
)
