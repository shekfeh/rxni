"""RXNI: reaction semantics, structural retrieval and condition intelligence.

Package-release versions are independent of the validated scientific
baseline and persisted retrieval-representation ABI.
"""

from .semantic.version import RUNTIME_VERSION as SCIENTIFIC_BASELINE_VERSION

# Compatibility export: do not reinterpret as the public package version.
RUNTIME_VERSION = SCIENTIFIC_BASELINE_VERSION
__version__ = "0.1.1"

__all__ = ["__version__", "SCIENTIFIC_BASELINE_VERSION", "RUNTIME_VERSION"]
