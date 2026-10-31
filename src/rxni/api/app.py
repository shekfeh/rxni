"""HTTP API for RXNI reaction analysis and retrieval."""

from __future__ import annotations

import asyncio
import uuid
from concurrent.futures import ThreadPoolExecutor
from contextlib import asynccontextmanager
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from rxni.public import (
    SCHEMA_VERSION,
    SERVICE_VERSION,
    AnalysisRequest,
    AnalysisResult,
    ReactionIntelligenceService,
    model_dump_compat,
)

from .config import ApiSettings
from .models import ApiInfoResponse, ErrorResponse, HealthResponse


API_VERSION = "reaction-intelligence-fastapi-v0.1.2"


def _service_info(service: ReactionIntelligenceService) -> ApiInfoResponse:
    manifest = getattr(service, "retrieval_manifest", {}) or {}
    runtime_version = None

    runtime_meta = manifest.get("runtime")
    if isinstance(runtime_meta, dict):
        runtime_version = (
            runtime_meta.get("runtime_version")
            or runtime_meta.get("version")
        )

    rep_policy = getattr(
        service,
        "retrieval_representation_semantic_policy",
        None,
    )
    semantic_policy = getattr(service, "mapping_control", "unknown")

    consistent = None
    if rep_policy is not None:
        consistent = str(rep_policy) == str(semantic_policy)

    return ApiInfoResponse(
        service_version=SERVICE_VERSION,
        schema_version=SCHEMA_VERSION,
        runtime_version=(
            str(runtime_version) if runtime_version is not None else None
        ),
        semantic_policy=str(semantic_policy),
        library_id=str(getattr(service, "library_id", "unknown")),
        retrieval_representation_semantic_policy=(
            str(rep_policy) if rep_policy is not None else None
        ),
        semantic_representation_consistent=consistent,
    )


def create_app(
    *,
    settings: ApiSettings | None = None,
    service: ReactionIntelligenceService | None = None,
) -> FastAPI:
    """Create the HTTP application.

    A service may be injected for tests.  In production, one service instance
    is created at application startup and retained for the process lifetime.
    """

    settings = settings or ApiSettings.from_env()
    injected_service = service

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        app.state.settings = settings

        # Runtime-v6.6 currently owns SQLite-backed lookup objects whose
        # connections are thread-affine.  FastAPI executes ordinary ``def``
        # endpoints in a worker thread, so constructing the service during
        # lifespan and later calling it from an arbitrary worker violates
        # SQLite's same-thread contract.
        #
        # The API therefore gives the scientific runtime one dedicated worker
        # thread.  The service (and therefore its SQLite connections) is
        # created in that worker and every analysis call is executed in the
        # same worker.  This also gives us the conservative serialized
        # execution model intended for the first HTTP release.
        if not settings.serialize_analysis:
            raise RuntimeError(
                "RI_SERIALIZE_ANALYSIS=0 is not supported by FastAPI v0.1.2. "
                "Runtime-v6.6 lookup SQLite connections are thread-affine; "
                "keep RI_SERIALIZE_ANALYSIS=1 until a separately validated "
                "concurrency/connection-pool design is introduced."
            )

        executor = ThreadPoolExecutor(
            max_workers=1,
            thread_name_prefix="reaction-intelligence-runtime",
        )
        app.state.runtime_executor = executor

        loop = asyncio.get_running_loop()

        if injected_service is not None:
            # Test/injected services may have been constructed elsewhere, but
            # all analysis calls are still serialized onto the runtime worker.
            app.state.ri_service = injected_service
        else:
            def build_service() -> ReactionIntelligenceService:
                return ReactionIntelligenceService(
                    library=settings.library,
                    mapping_control=settings.mapping_control,
                    retrieval_representation_semantic_policy=(
                        settings.retrieval_representation_semantic_policy
                    ),
                )

            # Critical thread-affinity invariant: production service creation
            # happens inside the same dedicated thread used for analysis.
            app.state.ri_service = await loop.run_in_executor(
                executor,
                build_service,
            )

        try:
            yield
        finally:
            executor.shutdown(wait=True, cancel_futures=True)

    app = FastAPI(
        title="RXNI API",
        version=API_VERSION,
        description=(
            "Typed HTTP interface to RXNI with a preserved Retrieval-V2 "
            "ReactionIntelligenceService."
        ),
        lifespan=lifespan,
    )

    if settings.cors_origins:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=list(settings.cors_origins),
            allow_credentials=True,
            allow_methods=["GET", "POST"],
            allow_headers=["*"],
        )

    @app.middleware("http")
    async def request_id_middleware(request: Request, call_next):
        request_id = request.headers.get("X-Request-ID") or str(uuid.uuid4())
        request.state.request_id = request_id
        response = await call_next(request)
        response.headers["X-Request-ID"] = request_id
        return response

    @app.exception_handler(Exception)
    async def unhandled_exception_handler(request: Request, exc: Exception):
        req_id = getattr(request.state, "request_id", None)
        details = None
        if settings.debug_errors:
            details = {
                "exception_type": type(exc).__name__,
                "exception": str(exc),
            }
        payload = ErrorResponse(
            error="analysis_failed",
            message=(
                str(exc)
                if settings.debug_errors
                else "Reaction analysis failed."
            ),
            request_id=req_id,
            details=details,
        )
        return JSONResponse(
            status_code=500,
            content=model_dump_compat(payload),
        )

    @app.get(
        "/healthz",
        response_model=HealthResponse,
        tags=["system"],
    )
    def healthz() -> HealthResponse:
        return HealthResponse(
            status="ok",
            service=SERVICE_VERSION,
            schema_version=SCHEMA_VERSION,
        )

    @app.get(
        f"{settings.api_prefix}/info",
        response_model=ApiInfoResponse,
        tags=["system"],
    )
    def info(request: Request) -> ApiInfoResponse:
        return _service_info(request.app.state.ri_service)

    @app.get(
        f"{settings.api_prefix}/schema/request",
        response_model=dict[str, Any],
        tags=["schema"],
    )
    def request_schema(request: Request) -> dict[str, Any]:
        return request.app.state.ri_service.request_json_schema()

    @app.get(
        f"{settings.api_prefix}/schema/result",
        response_model=dict[str, Any],
        tags=["schema"],
    )
    def result_schema(request: Request) -> dict[str, Any]:
        return request.app.state.ri_service.result_json_schema()

    @app.post(
        f"{settings.api_prefix}/analyze",
        response_model=AnalysisResult,
        responses={
            403: {"model": ErrorResponse},
            500: {"model": ErrorResponse},
        },
        tags=["analysis"],
    )
    async def analyze(
        payload: AnalysisRequest,
        request: Request,
    ) -> AnalysisResult:
        # RXNI_G1F_API_BOUNDARY_V01: input syntax validation outside runtime exception handling.
        from fastapi import HTTPException as _rxni_boundary_http_exception
        from rdkit import Chem as _rxni_boundary_chem
        _rxni_boundary_reaction = getattr(payload, "reaction", None)
        if isinstance(_rxni_boundary_reaction, str):
            _rxni_boundary_parts = _rxni_boundary_reaction.split(">")
            _rxni_boundary_invalid = (
                len(_rxni_boundary_parts) != 3
                or not _rxni_boundary_parts[0]
                or not _rxni_boundary_parts[2]
            )
            if not _rxni_boundary_invalid:
                for _rxni_boundary_part in _rxni_boundary_parts:
                    if not _rxni_boundary_part:
                        continue  # Agents may be empty: reactants>>products.
                    for _rxni_boundary_component in _rxni_boundary_part.split("."):
                        if not _rxni_boundary_component:
                            _rxni_boundary_invalid = True
                            break
                        try:
                            _rxni_boundary_mol = _rxni_boundary_chem.MolFromSmiles(
                                _rxni_boundary_component, sanitize=False
                            )
                        except Exception:
                            _rxni_boundary_mol = None
                        if _rxni_boundary_mol is None:
                            _rxni_boundary_invalid = True
                            break
                    if _rxni_boundary_invalid:
                        break
            if _rxni_boundary_invalid:
                raise _rxni_boundary_http_exception(
                    status_code=422,
                    detail={
                        "error": "invalid_reaction_smiles",
                        "message": "Expected parseable reactants>agents>products reaction SMILES.",
                    },
                )
        if (
            payload.options.include_raw_runtime_result
            and not settings.allow_raw_runtime_result
        ):
            raise HTTPException(
                status_code=403,
                detail=(
                    "Raw runtime result exposure is disabled for this API "
                    "process."
                ),
            )

        service_obj = request.app.state.ri_service
        executor = request.app.state.runtime_executor
        loop = asyncio.get_running_loop()

        # Do not let FastAPI/AnyIO choose an arbitrary worker thread here.
        # The same dedicated runtime thread owns both service construction and
        # SQLite-backed retrieval execution.
        return await loop.run_in_executor(
            executor,
            service_obj.analyze_request,
            payload,
        )

    return app


# Importable ASGI application for:
#   uvicorn runtime_v6_6.api.app:app
#
# Configuration is intentionally environment-based at process startup.
class _DeferredASGIApp:
    """Keep the ``:app`` ASGI entry point without requiring env at import time.

    Configuration and the scientific service are initialized on the first ASGI
    invocation (normally the lifespan scope), not on importing the module.
    """

    def __init__(self):
        self._impl = None

    async def __call__(self, scope, receive, send):
        if self._impl is None:
            self._impl = create_app()
        await self._impl(scope, receive, send)


app = _DeferredASGIApp()