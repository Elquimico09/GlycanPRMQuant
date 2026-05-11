"""Helpers for locating bundled package resources."""

from pathlib import Path
import sys


def resource_path(relative_path: str) -> str:
    """
    Return an absolute path for a bundled resource.

    Works from the source tree and from PyInstaller builds where bundled files
    are unpacked under ``sys._MEIPASS``.
    """
    rel = Path(relative_path)
    if rel.is_absolute():
        return str(rel)

    bundle_root = getattr(sys, "_MEIPASS", None)
    if bundle_root:
        candidate = Path(bundle_root) / rel
        if candidate.exists():
            return str(candidate)

    repo_root = Path(__file__).resolve().parent.parent
    candidate = repo_root / rel
    if candidate.exists():
        return str(candidate)

    return str(Path.cwd() / rel)
