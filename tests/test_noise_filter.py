import numpy as np
import pandas as pd
import pytest

from glycanPRMQuant.noise_filter import (
    filter_ms2_noise,
    normalize_noise_filter_mode,
)


def _peak_table(scan_number, intensities, rt=1.0, precursor_mz=700.0):
    count = len(intensities)
    return pd.DataFrame(
        {
            "scan_number": [scan_number] * count,
            "rt": [rt] * count,
            "precursor_mz": [precursor_mz] * count,
            "fragment_mz": np.arange(count, dtype=float) + 100.0,
            "fragment_intensity": intensities,
        }
    )


@pytest.mark.parametrize(
    ("mode", "expected"),
    [("automatic", "auto"), ("auto", "auto"), ("off", "off"), ("manual", "manual")],
)
def test_normalizes_noise_filter_modes(mode, expected):
    assert normalize_noise_filter_mode(mode) == expected


def test_rejects_unknown_noise_filter_mode():
    with pytest.raises(ValueError, match="auto, off, or manual"):
        normalize_noise_filter_mode("dynamic")


def test_automatic_filter_removes_robust_low_intensity_population():
    peaks = _peak_table(1, [10.0] * 18 + [100.0, 200.0])

    result = filter_ms2_noise(peaks, mode="auto")

    assert result.peaks["fragment_intensity"].tolist() == [100.0, 200.0]
    audit = result.audit.iloc[0]
    assert audit["noise_estimation_source"] == "scan"
    assert audit["estimated_noise_floor"] == 10.0
    assert audit["raw_peak_count"] == 20
    assert audit["retained_peak_count"] == 2
    assert audit["raw_ms2_tic"] == 480.0
    assert audit["denoised_ms2_tic"] == 300.0
    assert result.peaks["raw_ms2_tic"].eq(480.0).all()
    assert result.peaks["denoised_ms2_tic"].eq(300.0).all()


def test_sparse_scan_uses_adjacent_same_precursor_scans():
    peaks = pd.concat(
        [
            _peak_table(1, [10.0] * 9 + [100.0], rt=1.0),
            _peak_table(2, [10.0] * 9 + [200.0], rt=1.1),
        ],
        ignore_index=True,
    )

    result = filter_ms2_noise(peaks, mode="auto")

    assert result.audit["noise_estimation_source"].eq("neighbor_pool").all()
    assert result.audit["noise_estimation_peak_count"].eq(20).all()
    assert result.peaks.groupby("scan_number")["fragment_intensity"].max().to_dict() == {
        1: 100.0,
        2: 200.0,
    }


def test_sparse_scan_without_neighbors_is_left_unfiltered():
    peaks = _peak_table(1, [10.0, 20.0, 30.0])

    result = filter_ms2_noise(peaks, mode="auto")

    assert len(result.peaks) == 3
    assert result.audit.loc[0, "noise_estimation_source"] == "insufficient_peaks_unfiltered"
    assert result.audit.loc[0, "estimated_noise_floor"] == 0.0


def test_off_mode_retains_all_positive_centroids():
    peaks = _peak_table(1, [1.0, 2.0, 3.0])

    result = filter_ms2_noise(peaks, mode="off")

    assert result.peaks["fragment_intensity"].tolist() == [1.0, 2.0, 3.0]
    assert result.audit.loc[0, "denoised_ms2_tic"] == 6.0


def test_manual_mode_uses_absolute_threshold():
    peaks = _peak_table(1, [99.0, 100.0, 101.0])

    result = filter_ms2_noise(peaks, mode="manual", manual_threshold=100.0)

    assert result.peaks["fragment_intensity"].tolist() == [101.0]
    assert result.audit.loc[0, "estimated_noise_floor"] == 100.0


def test_manual_mode_rejects_negative_threshold():
    peaks = _peak_table(1, [100.0])
    with pytest.raises(ValueError, match="cannot be negative"):
        filter_ms2_noise(peaks, mode="manual", manual_threshold=-1.0)
