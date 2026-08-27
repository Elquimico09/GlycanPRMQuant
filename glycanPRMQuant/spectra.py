"""Mass-spectrometry input detection and MS2 extraction dispatch."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from glycanPRMQuant.msfileReader import extractMS2 as extract_mzml_ms2


SUPPORTED_SUFFIXES = {".mzml", ".raw"}


def detect_input_type(path: str | Path) -> str:
    suffix = Path(path).suffix.lower()
    if suffix == ".mzml":
        return "mzml"
    if suffix == ".raw":
        return "thermo_raw"
    raise ValueError(
        f"Unsupported mass-spectrometry file '{path}'. Expected a Thermo .raw or .mzML file."
    )


def validate_input_file_types(paths: Iterable[str | Path]) -> str:
    input_types = {detect_input_type(path) for path in paths}
    if not input_types:
        raise ValueError("At least one input file is required")
    if len(input_types) != 1:
        raise ValueError("A batch cannot mix Thermo .raw and .mzML input files")
    return input_types.pop()


def extract_ms2(path: str | Path, min_intensity: float = 1):
    input_type = detect_input_type(path)
    if input_type == "mzml":
        return extract_mzml_ms2(str(path), min_intensity=min_intensity)

    from glycanPRMQuant.thermo_raw import extract_thermo_ms2

    return extract_thermo_ms2(str(path), min_intensity=min_intensity)
