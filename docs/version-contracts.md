# RXNI version and compatibility contracts

RXNI uses independent version domains:

- **Package distribution:** `rxni.__version__` (`0.1.0` for the initial public release).
- **Validated semantic baseline:** `rxni.SCIENTIFIC_BASELINE_VERSION` (`6.7.0`). `rxni.semantic.version` records the scientific engine, ruleset, FGI contract, and profile-policy identifiers that produced semantic representations.
- **Stored retrieval representation:** `rxni.retrieval.version.STORED_REPRESENTATION_GENERATION` (`retrieval_representations_v2`). `RETRIEVAL_API_VERSION`, `LIVE_RECORD_VERSION`, and `V5_COMPAT_REPRESENTATION` retain historical string values because they are part of compatibility metadata.
- **Public typed response contract:** `SCHEMA_VERSION` in `rxni.public.models` is independent from the package version. The service and HTTP API also retain their historical response-version values until a versioned schema migration is tested.

Do not rewrite existing library manifests or stored records merely to change display names. A reader must reject incompatible representation generations instead of silently combining them.

The V2 `semantics` field is JSON-based. JSON object keys are strings; integer keys in a live Python dictionary become string keys during storage and decoding. Compare decoded semantics against JSON-normalized live results, while checking exact physical cell encoding separately.

Legacy environment variables such as `RUNTIME_V67_MAPPING_CONTROL` and `REACTION_INTELLIGENCE_LIBRARY` remain accepted for compatibility. They are not package-version identifiers.
