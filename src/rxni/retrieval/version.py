"""Serialized retrieval-ABI and record-format contract identifiers.

Historic v66/v2 and system_i strings are persisted format identifiers.
Keep them unchanged for existing library compatibility; the RXNI package
version is maintained separately in `rxni.__version__`.
"""

RETRIEVAL_API_VERSION = "v66-retrieval-api-v1"
LIVE_RECORD_VERSION = "v66-live-record-v2"
STORED_REPRESENTATION_GENERATION = "retrieval_representations_v2"
V5_COMPAT_REPRESENTATION = "system_i_representations_v1"

# M4 established:
# - exact structural keys stable
# - candidate envelope stable
# - canonical semantics can change reaction centers/localized/multichannel payload
# Therefore v1 and v2 stored representations must never be silently mixed.
