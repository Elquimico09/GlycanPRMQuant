"""Automatic and user-controlled filtering of centroided MS2 noise peaks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


NOISE_AUDIT_COLUMNS = [
    "scan_number",
    "rt",
    "precursor_mz",
    "noise_filter_mode",
    "noise_estimation_source",
    "noise_estimation_peak_count",
    "raw_peak_count",
    "retained_peak_count",
    "rejected_peak_count",
    "raw_ms2_tic",
    "denoised_ms2_tic",
    "noise_median_intensity",
    "noise_mad_intensity",
    "noise_sigma_intensity",
    "estimated_noise_floor",
    "noise_lower_fraction",
    "noise_sigma_multiplier",
    "minimum_estimation_peaks",
    "neighboring_scans",
    "precursor_pool_ppm",
]


@dataclass(frozen=True)
class NoiseFilterResult:
    """Filtered peak rows and one audit record per original MS2 scan."""

    peaks: pd.DataFrame
    audit: pd.DataFrame


def normalize_noise_filter_mode(mode: str) -> str:
    """Return the canonical noise-filter mode name."""
    normalized = str(mode).strip().lower()
    aliases = {"automatic": "auto", "none": "off", "disabled": "off"}
    normalized = aliases.get(normalized, normalized)
    if normalized not in {"auto", "off", "manual"}:
        raise ValueError("MS2 noise filter mode must be auto, off, or manual")
    return normalized


def _precursor_clusters(precursor_mz: pd.Series, ppm_tolerance: float) -> pd.Series:
    """Cluster scan-level precursor m/z values for neighboring-scan pooling."""
    values = pd.to_numeric(precursor_mz, errors="coerce").to_numpy(dtype=float)
    labels = np.full(values.size, -1, dtype=int)
    finite_indices = np.flatnonzero(np.isfinite(values) & (values > 0))
    if finite_indices.size:
        ordered = finite_indices[np.argsort(values[finite_indices], kind="stable")]
        cluster = 0
        labels[ordered[0]] = cluster
        previous = values[ordered[0]]
        for index in ordered[1:]:
            current = values[index]
            delta_ppm = abs(current - previous) / max(current, previous) * 1e6
            if delta_ppm > ppm_tolerance:
                cluster += 1
            labels[index] = cluster
            previous = current

    next_cluster = int(labels.max()) + 1
    for index in np.flatnonzero(labels < 0):
        labels[index] = next_cluster
        next_cluster += 1
    return pd.Series(labels, index=precursor_mz.index, dtype=int)


def _estimate_noise_floor(
    intensities: np.ndarray,
    lower_fraction: float,
    sigma_multiplier: float,
) -> tuple[float, float, float, float]:
    """Estimate a robust floor from the lower-intensity centroid population."""
    values = np.asarray(intensities, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return 0.0, np.nan, np.nan, np.nan

    ordered = np.sort(values)
    lower_count = max(3, int(np.ceil(ordered.size * lower_fraction)))
    lower = ordered[: min(lower_count, ordered.size)]
    median = float(np.median(lower))
    mad = float(np.median(np.abs(lower - median)))
    sigma = float(1.4826 * mad)

    if sigma > 0:
        floor = median + sigma_multiplier * sigma
    elif float(np.max(values)) > median:
        # A constant low-intensity population plus larger peaks is a common
        # centroided-noise pattern. The population intensity itself is the
        # most conservative useful floor when its MAD is zero.
        floor = median
    else:
        # Equal-intensity sparse spectra do not contain enough information to
        # distinguish noise from signal; retain them rather than guessing.
        floor = 0.0
    return float(max(floor, 0.0)), median, mad, sigma


def filter_ms2_noise(
    peaks: pd.DataFrame,
    mode: str = "auto",
    manual_threshold: float = 100.0,
    *,
    minimum_estimation_peaks: int = 20,
    neighboring_scans: int = 2,
    precursor_pool_ppm: float = 20.0,
    lower_fraction: float = 0.60,
    sigma_multiplier: float = 3.0,
) -> NoiseFilterResult:
    """Filter centroided MS2 peaks and audit each scan's estimated noise floor.

    Automatic estimates use the lower 60% of positive centroid intensities and
    a robust ``median + 3 * 1.4826 * MAD`` floor. Sparse scans borrow peaks from
    up to two adjacent scans on either side with the same precursor cluster. If
    the pooled population is still too small, the scan is left unfiltered.
    """
    normalized_mode = normalize_noise_filter_mode(mode)
    if normalized_mode == "manual" and manual_threshold < 0:
        raise ValueError("Manual MS2 intensity threshold cannot be negative")
    if minimum_estimation_peaks < 3:
        raise ValueError("minimum_estimation_peaks must be at least 3")
    if neighboring_scans < 0:
        raise ValueError("neighboring_scans cannot be negative")
    if precursor_pool_ppm <= 0:
        raise ValueError("precursor_pool_ppm must be positive")
    if not 0 < lower_fraction <= 1:
        raise ValueError("lower_fraction must be in (0, 1]")
    if sigma_multiplier < 0:
        raise ValueError("sigma_multiplier cannot be negative")

    required = {"scan_number", "rt", "precursor_mz", "fragment_intensity"}
    missing = required.difference(peaks.columns)
    if missing:
        raise ValueError(f"MS2 noise filtering requires columns: {sorted(missing)}")
    if peaks.empty:
        empty_peaks = peaks.copy()
        for column in (
            "raw_ms2_tic",
            "denoised_ms2_tic",
            "estimated_noise_floor",
            "noise_filter_mode",
        ):
            empty_peaks[column] = pd.Series(
                dtype=float if column != "noise_filter_mode" else object
            )
        return NoiseFilterResult(empty_peaks, pd.DataFrame(columns=NOISE_AUDIT_COLUMNS))

    work = peaks.copy().reset_index(drop=True)
    work["fragment_intensity"] = pd.to_numeric(
        work["fragment_intensity"], errors="coerce"
    )
    valid_positive = np.isfinite(work["fragment_intensity"]) & (
        work["fragment_intensity"] > 0
    )
    work = work.loc[valid_positive].copy()
    if work.empty:
        for column in (
            "raw_ms2_tic",
            "denoised_ms2_tic",
            "estimated_noise_floor",
            "noise_filter_mode",
        ):
            work[column] = pd.Series(
                dtype=float if column != "noise_filter_mode" else object
            )
        return NoiseFilterResult(work, pd.DataFrame(columns=NOISE_AUDIT_COLUMNS))

    scan_groups = {
        scan_number: group.index.to_numpy()
        for scan_number, group in work.groupby("scan_number", sort=False)
    }
    scan_table = (
        work.groupby("scan_number", sort=False)
        .agg(rt=("rt", "median"), precursor_mz=("precursor_mz", "median"))
        .reset_index()
    )
    scan_table["_precursor_cluster"] = _precursor_clusters(
        scan_table["precursor_mz"], precursor_pool_ppm
    )

    neighbor_pools: dict[object, np.ndarray] = {}
    for _, cluster in scan_table.groupby("_precursor_cluster", sort=False):
        ordered = cluster.sort_values(["rt", "scan_number"]).reset_index(drop=True)
        scan_numbers = ordered["scan_number"].tolist()
        for position, scan_number in enumerate(scan_numbers):
            start = max(0, position - neighboring_scans)
            stop = min(len(scan_numbers), position + neighboring_scans + 1)
            indices = np.concatenate(
                [scan_groups[number] for number in scan_numbers[start:stop]]
            )
            neighbor_pools[scan_number] = work.loc[
                indices, "fragment_intensity"
            ].to_numpy(dtype=float)

    keep = pd.Series(False, index=work.index, dtype=bool)
    audit_records = []
    floors: dict[object, float] = {}
    raw_tics: dict[object, float] = {}
    denoised_tics: dict[object, float] = {}

    scan_metadata = scan_table.set_index("scan_number")
    for scan_number, indices in scan_groups.items():
        intensities = work.loc[indices, "fragment_intensity"].to_numpy(dtype=float)
        raw_tic = float(np.sum(intensities))
        median = mad = sigma = np.nan
        estimation_count = 0

        if normalized_mode == "off":
            floor = 0.0
            source = "off"
        elif normalized_mode == "manual":
            floor = float(manual_threshold)
            source = "manual"
        else:
            estimation_values = intensities
            source = "scan"
            if estimation_values.size < minimum_estimation_peaks:
                estimation_values = neighbor_pools[scan_number]
                source = "neighbor_pool"
            estimation_count = int(estimation_values.size)
            if estimation_values.size < minimum_estimation_peaks:
                floor = 0.0
                source = "insufficient_peaks_unfiltered"
            else:
                floor, median, mad, sigma = _estimate_noise_floor(
                    estimation_values, lower_fraction, sigma_multiplier
                )

        scan_keep = intensities > floor
        keep.loc[indices] = scan_keep
        retained = intensities[scan_keep]
        denoised_tic = float(np.sum(retained))
        floors[scan_number] = float(floor)
        raw_tics[scan_number] = raw_tic
        denoised_tics[scan_number] = denoised_tic
        metadata = scan_metadata.loc[scan_number]
        audit_records.append(
            {
                "scan_number": scan_number,
                "rt": float(metadata["rt"]),
                "precursor_mz": float(metadata["precursor_mz"]),
                "noise_filter_mode": normalized_mode,
                "noise_estimation_source": source,
                "noise_estimation_peak_count": estimation_count,
                "raw_peak_count": int(intensities.size),
                "retained_peak_count": int(scan_keep.sum()),
                "rejected_peak_count": int((~scan_keep).sum()),
                "raw_ms2_tic": raw_tic,
                "denoised_ms2_tic": denoised_tic,
                "noise_median_intensity": median,
                "noise_mad_intensity": mad,
                "noise_sigma_intensity": sigma,
                "estimated_noise_floor": float(floor),
                "noise_lower_fraction": float(lower_fraction),
                "noise_sigma_multiplier": float(sigma_multiplier),
                "minimum_estimation_peaks": int(minimum_estimation_peaks),
                "neighboring_scans": int(neighboring_scans),
                "precursor_pool_ppm": float(precursor_pool_ppm),
            }
        )

    filtered = work.loc[keep].copy()
    filtered["raw_ms2_tic"] = filtered["scan_number"].map(raw_tics)
    filtered["denoised_ms2_tic"] = filtered["scan_number"].map(denoised_tics)
    filtered["estimated_noise_floor"] = filtered["scan_number"].map(floors)
    filtered["noise_filter_mode"] = normalized_mode
    audit = pd.DataFrame(audit_records, columns=NOISE_AUDIT_COLUMNS)
    return NoiseFilterResult(filtered.reset_index(drop=True), audit)
