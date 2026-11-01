"""RXNI retrieval public API; preserves historically frozen ABI and code aliases."""

from .library import LibraryLayout
from .providers import (
    LiveV66RecordProvider,
    RetrievalRecordProvider,
    V5CompatibilityProvider,
)
from .record import build_retrieval_record
from .scoring import member_scores, template_score
from .storage.provider import V2StoredRecordProvider
from .structural_keys import derive_structural_keys
from .version import (
    LIVE_RECORD_VERSION,
    RETRIEVAL_API_VERSION,
    STORED_REPRESENTATION_GENERATION,
    V5_COMPAT_REPRESENTATION,
)

__all__ = [
    "LibraryLayout",
    "LiveV66RecordProvider",
    "RetrievalRecordProvider",
    "V5CompatibilityProvider",
    "V2StoredRecordProvider",
    "StoredRecordProvider",
    "build_retrieval_record",
    "member_scores",
    "template_score",
    "derive_structural_keys",
    "LIVE_RECORD_VERSION",
    "RETRIEVAL_API_VERSION",
    "STORED_REPRESENTATION_GENERATION",
    "V5_COMPAT_REPRESENTATION",
]

# M7E3B2_STRUCTURAL_KEY_API_ALIAS
# Compatibility API only. Scientific implementation remains derive_structural_keys.
from .structural_keys import derive_structural_keys
derive_query_structural_keys = derive_structural_keys

# Stable user-facing alias. The V2 stored format remains unchanged.
StoredRecordProvider = V2StoredRecordProvider
