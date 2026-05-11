"""Logging setup helpers for command-line and GUI workflows."""

import logging
import sys


class QueueLogHandler(logging.Handler):
    """Logging handler that writes formatted records to a multiprocessing queue."""

    def __init__(self, queue):
        super().__init__()
        self.queue = queue

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self.queue.put(self.format(record) + "\n")
        except Exception:
            self.handleError(record)


def configure_logging(level: int = logging.INFO, force: bool = False, log_queue=None) -> None:
    """Configure package logging with a concise console format."""
    root = logging.getLogger()
    formatter = logging.Formatter("%(levelname)s: %(message)s")

    if force or not root.handlers:
        for handler in list(root.handlers):
            root.removeHandler(handler)

        if log_queue is None:
            handler = logging.StreamHandler(sys.stderr)
        else:
            handler = QueueLogHandler(log_queue)
        handler.setFormatter(formatter)
        root.addHandler(handler)
        root.setLevel(level)
    else:
        root.setLevel(level)

    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
