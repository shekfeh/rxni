"""Environment-backed RXNI HTTP API configuration.

Legacy environment variables are accepted for library compatibility.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


TRUE_VALUES = {"1", "true", "yes", "on"}


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in TRUE_VALUES


def _env_csv(name: str) -> tuple[str, ...]:
    raw = os.environ.get(name, "")
    return tuple(x.strip() for x in raw.split(",") if x.strip())


@dataclass(frozen=True)
class ApiSettings:
    library: Path
    mapping_control: str = "conversion_validated"
    retrieval_representation_semantic_policy: str | None = None
    cors_origins: tuple[str, ...] = ()
    allow_raw_runtime_result: bool = False
    debug_errors: bool = False
    serialize_analysis: bool = True
    api_prefix: str = "/v1"

    @classmethod
    def from_env(cls) -> "ApiSettings":
        library = (
            os.environ.get("REACTION_INTELLIGENCE_LIBRARY")
            or os.environ.get("V5")
        )
        if not library:
            raise RuntimeError(
                "No template library configured. Set "
                "REACTION_INTELLIGENCE_LIBRARY or V5."
            )

        mapping_control = os.environ.get(
            "RUNTIME_V67_MAPPING_CONTROL",
            "conversion_validated",
        )

        if mapping_control not in {"conversion_validated", "strict"}:
            raise RuntimeError(
                "RUNTIME_V67_MAPPING_CONTROL must be "
                "'conversion_validated' or 'strict'."
            )

        semantic_policy = os.environ.get(
            "RI_RETRIEVAL_REPRESENTATION_SEMANTIC_POLICY"
        )

        prefix = os.environ.get("RI_API_PREFIX", "/v1").strip()
        if not prefix.startswith("/"):
            prefix = "/" + prefix
        prefix = prefix.rstrip("/") or "/v1"

        return cls(
            library=Path(library),
            mapping_control=mapping_control,
            retrieval_representation_semantic_policy=semantic_policy,
            cors_origins=_env_csv("RI_CORS_ORIGINS"),
            allow_raw_runtime_result=_env_bool(
                "RI_ALLOW_RAW_RUNTIME_RESULT",
                False,
            ),
            debug_errors=_env_bool("RI_API_DEBUG_ERRORS", False),
            serialize_analysis=_env_bool(
                "RI_SERIALIZE_ANALYSIS",
                True,
            ),
            api_prefix=prefix,
        )
