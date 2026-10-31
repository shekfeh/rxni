"""HTTP-only response models.

Scientific/public reaction models remain owned by ``runtime_v6_6.public``.
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class HealthResponse(BaseModel):
    status: str
    service: str
    schema_version: str


class ApiInfoResponse(BaseModel):
    service_version: str
    schema_version: str
    runtime_version: str | None = None
    semantic_policy: str
    library_id: str
    retrieval_representation_semantic_policy: str | None = None
    semantic_representation_consistent: bool | None = None


class ErrorResponse(BaseModel):
    error: str
    message: str
    request_id: str | None = None
    details: dict[str, Any] | None = None
