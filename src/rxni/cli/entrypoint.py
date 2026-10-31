"""Public RXNI console entry point.

The package entry point is `rxni`; invoke `rxni.cli.entrypoint:main`
from the distribution configuration once packaging is finalized."""

from __future__ import annotations

from .app import app


def main() -> None:
    app()


if __name__ == "__main__":
    main()
