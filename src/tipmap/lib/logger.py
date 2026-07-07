"""Logging configuration helpers for TIPMap command-line tools."""

from __future__ import annotations

import logging


def configure_logging(level: int = logging.INFO) -> None:
    """Configure a concise default logger for command-line entry points."""

    logging.basicConfig(level=level, format="%(levelname)s:%(name)s:%(message)s")
