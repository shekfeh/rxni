"""Optional RXNI HTTP interface; importing configuration has no API side effects.

HTTP dependencies are optional. Use ``from rxni.api import create_app`` or
``uvicorn rxni.api.app:app`` after installing ``rxni[api]``.
"""
from __future__ import annotations
from importlib import import_module
from .config import ApiSettings

__all__ = ["API_VERSION", "ApiSettings", "app", "create_app"]


def create_app(*, settings=None, service=None):
    """Create the HTTP application on demand; no work at package import."""
    factory = import_module(".app", __name__).create_app
    return factory(settings=settings, service=service)


def __getattr__(name):
    if name in {"API_VERSION", "app"}:
        return getattr(import_module(".app", __name__), name)
    raise AttributeError(name)
