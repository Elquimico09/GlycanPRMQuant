"""Logging setup helpers for command-line and GUI workflows."""

import logging
import sys


def configure_logging(level: int = logging.INFO, force: bool = False) -> None:
    """Configure package logging with a concise console format."""
    if force or not logging.getLogger().handlers:
        logging.basicConfig(
            level=level,
            format="%(levelname)s: %(message)s",
            stream=sys.stderr,
            force=force,
        )
    else:
        logging.getLogger().setLevel(level)

    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
