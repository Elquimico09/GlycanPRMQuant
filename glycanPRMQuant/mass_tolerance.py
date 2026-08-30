"""Shared fragment-mass tolerance calculations."""

from __future__ import annotations

import numpy as np


TOLERANCE_UNITS = ("Da", "ppm")


def normalize_tolerance_unit(unit: str) -> str:
    """Return the canonical ``Da`` or ``ppm`` spelling."""
    normalized = str(unit or "").strip().lower()
    if normalized in {"da", "dalton", "daltons"}:
        return "Da"
    if normalized == "ppm":
        return "ppm"
    raise ValueError("Fragment tolerance unit must be 'Da' or 'ppm'")


def validate_tolerance(value: float, unit: str) -> tuple[float, str]:
    """Validate and normalize a fragment tolerance value/unit pair."""
    value = float(value)
    unit = normalize_tolerance_unit(unit)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("Fragment tolerance value must be positive")
    if unit == "ppm" and value >= 1e6:
        raise ValueError("Fragment ppm tolerance must be less than 1,000,000")
    return value, unit


def mass_error_ppm(observed_mz, theoretical_mz):
    """Return absolute mass error in ppm using theoretical m/z as denominator."""
    observed = np.asarray(observed_mz, dtype=float)
    theoretical = np.asarray(theoretical_mz, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs(observed - theoretical) / theoretical * 1e6
    return np.where(theoretical > 0, error, np.nan)


def within_tolerance(
    observed_mz,
    theoretical_mz,
    tolerance_value: float,
    tolerance_unit: str,
):
    """Return whether observed/theoretical pairs satisfy the selected tolerance."""
    tolerance_value, tolerance_unit = validate_tolerance(
        tolerance_value, tolerance_unit
    )
    observed = np.asarray(observed_mz, dtype=float)
    theoretical = np.asarray(theoretical_mz, dtype=float)
    if tolerance_unit == "Da":
        return np.abs(observed - theoretical) <= tolerance_value
    return mass_error_ppm(observed, theoretical) <= tolerance_value


def conservative_query_radius_da(
    observed_mz,
    tolerance_value: float,
    tolerance_unit: str,
):
    """Return a Da radius guaranteed to contain every exact tolerance match."""
    tolerance_value, tolerance_unit = validate_tolerance(
        tolerance_value, tolerance_unit
    )
    observed = np.asarray(observed_mz, dtype=float)
    if tolerance_unit == "Da":
        return np.full(observed.shape, tolerance_value, dtype=float)
    fraction = tolerance_value / 1e6
    return np.abs(observed) * fraction / (1.0 - fraction)


def equivalent_tolerance_ppm(
    tolerance_value: float,
    tolerance_unit: str,
    representative_mz: float,
) -> float:
    """Express a selected tolerance as ppm at a representative fragment m/z."""
    tolerance_value, tolerance_unit = validate_tolerance(
        tolerance_value, tolerance_unit
    )
    if tolerance_unit == "ppm":
        return tolerance_value
    representative_mz = float(representative_mz)
    if not np.isfinite(representative_mz) or representative_mz <= 0:
        return float("nan")
    return tolerance_value / representative_mz * 1e6


def mass_accuracy_reliability(
    tolerance_value: float,
    tolerance_unit: str,
    representative_mz: float,
    reference_ppm: float = 20.0,
) -> tuple[float, float]:
    """Return equivalent ppm and the tolerance-dependent accuracy reliability."""
    if not np.isfinite(reference_ppm) or reference_ppm <= 0:
        raise ValueError("Fragment mass-accuracy reference ppm must be positive")
    equivalent_ppm = equivalent_tolerance_ppm(
        tolerance_value, tolerance_unit, representative_mz
    )
    if not np.isfinite(equivalent_ppm) or equivalent_ppm <= 0:
        return equivalent_ppm, 1.0
    reliability = min(1.0, float(np.sqrt(reference_ppm / equivalent_ppm)))
    return equivalent_ppm, reliability
