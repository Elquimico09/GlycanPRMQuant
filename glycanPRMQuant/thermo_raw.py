"""Thermo RAW MS2 extraction through AlphaRaw."""

from __future__ import annotations

import numpy as np
import pandas as pd


MS2_COLUMNS = [
    "scan_number",
    "rt",
    "precursor_mz",
    "fragment_mz",
    "precursor_intensity",
    "fragment_intensity",
]


def _alpharaw_reader_class():
    """Import AlphaRaw lazily so mzML-only workflows do not load pythonnet."""
    try:
        from alpharaw.thermo import ThermoRawData
    except (ImportError, RuntimeError) as exc:
        raise RuntimeError(
            "Thermo RAW support requires AlphaRaw and its Windows runtime dependencies. "
            "Install the project dependencies and run this workflow on Windows."
        ) from exc
    return ThermoRawData


def _normalize_alpharaw_ms2(
    spectrum_df: pd.DataFrame,
    peak_df: pd.DataFrame,
    min_intensity: float = 1,
) -> pd.DataFrame:
    """Convert AlphaRaw's spectrum/peak tables to the pipeline's MS2 schema."""
    spectrum_columns = {
        "spec_idx",
        "rt",
        "ms_level",
        "precursor_mz",
        "peak_start_idx",
        "peak_stop_idx",
    }
    peak_columns = {"mz", "intensity"}

    missing_spectrum = spectrum_columns.difference(spectrum_df.columns)
    missing_peak = peak_columns.difference(peak_df.columns)
    if missing_spectrum or missing_peak:
        details = []
        if missing_spectrum:
            details.append(f"spectrum columns: {sorted(missing_spectrum)}")
        if missing_peak:
            details.append(f"peak columns: {sorted(missing_peak)}")
        raise ValueError("AlphaRaw returned an unexpected table schema; missing " + ", ".join(details))

    rows = []
    for spectrum in spectrum_df.loc[spectrum_df["ms_level"] == 2].itertuples(index=False):
        precursor_mz = float(spectrum.precursor_mz)
        if not np.isfinite(precursor_mz) or precursor_mz <= 0:
            continue

        start = int(spectrum.peak_start_idx)
        stop = int(spectrum.peak_stop_idx)
        if start < 0 or stop <= start or stop > len(peak_df):
            continue

        peaks = peak_df.iloc[start:stop]
        mz_values = peaks["mz"].to_numpy(dtype=float, copy=False)
        intensity_values = peaks["intensity"].to_numpy(dtype=float, copy=False)
        mask = (
            np.isfinite(mz_values)
            & np.isfinite(intensity_values)
            & (intensity_values > min_intensity)
        )

        scan_num = getattr(spectrum, "scan_num", None)
        if scan_num is None or not np.isfinite(scan_num):
            scan_num = int(spectrum.spec_idx) + 1

        for fragment_mz, fragment_intensity in zip(mz_values[mask], intensity_values[mask]):
            rows.append(
                {
                    "scan_number": int(scan_num),
                    "rt": float(spectrum.rt),
                    "precursor_mz": precursor_mz,
                    "fragment_mz": float(fragment_mz),
                    "precursor_intensity": 0.0,
                    "fragment_intensity": float(fragment_intensity),
                }
            )

    return pd.DataFrame(rows, columns=MS2_COLUMNS)


def extract_thermo_ms2(raw_file: str, min_intensity: float = 1) -> pd.DataFrame:
    """Read centroided PRM MS2 spectra directly from a Thermo ``.raw`` file."""
    reader_class = _alpharaw_reader_class()
    reader = reader_class(
        centroided=True,
        dda=False,
        process_count=1,
        save_as_hdf=False,
    )
    reader.import_raw(raw_file)
    return _normalize_alpharaw_ms2(
        reader.spectrum_df,
        reader.peak_df,
        min_intensity=min_intensity,
    )
