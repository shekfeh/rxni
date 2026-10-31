"""Public Runtime v6.6 unified semantic pipeline.

Stage E makes validated FGI enforcement the canonical/default backend.

Non-FGI semantics still come from the frozen v6.5.6 reference bridge during
this migration step.  Conversion-derived truth is replaced by the native
82-event hypothesis -> validation path before redox/profile derivation.

This is intentionally conservative:
a rejected/ambiguous exact conversion is removed; no generic fallback primary
is invented in its place.
"""

from __future__ import annotations

from .derived import build_semantic_profile, derive_redox_axis
from .enforcement import enforce_validated_directional_semantics
from .types import SemanticResult
from .version import PROFILE_POLICY_CANONICAL, PROFILE_POLICY_PARITY


def classify_reaction(
    reaction: str,
    *,
    backend: str = "validated",
    profile_policy: str | None = None,
) -> SemanticResult:
    """Return the canonical Runtime-v6.6 SemanticResult.

    backends
    --------
    validated / native
        Canonical Stage-E behavior.  Frozen non-FGI semantics are overlaid with
        native validated FGI, then redox/profile are rebuilt from canonical data.
    reference
        Exact migration/reference behavior for audit only.
    """
    # Lazy import preserves validation-module independence from runtime_v6_5.
    from .reference_backend import classify_reference

    if backend == "reference":
        result = classify_reference(reaction)
        policy = profile_policy or PROFILE_POLICY_PARITY
        result.semantic_profile = build_semantic_profile(result, policy=policy)
        result.diagnostics["profile_policy"] = policy
        result.diagnostics["canonical_fgi_enforcement"] = False
        return result

    if backend not in {"validated", "native"}:
        raise ValueError(f"Unknown semantic backend: {backend!r}")

    reference = classify_reference(reaction)
    result = enforce_validated_directional_semantics(reference, reaction)

    # Single ownership path: redox only after accepted FGI validation.
    result.redox_axis = derive_redox_axis(result.directional_conversions)

    policy = profile_policy or PROFILE_POLICY_CANONICAL
    result.semantic_profile = build_semantic_profile(result, policy=policy)
    result.diagnostics["profile_policy"] = policy
    result.diagnostics["canonical_redox_source"] = (
        "accepted_directional_conversions_only"
    )
    return result
