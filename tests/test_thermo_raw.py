import pandas as pd
import pytest

from glycanPRMQuant import thermo_raw
from glycanPRMQuant.thermo_raw import MS2_COLUMNS, _normalize_alpharaw_ms2


def test_normalizes_alpharaw_tables_for_prm_and_filters_peaks():
    spectra = pd.DataFrame(
        {
            "spec_idx": [0, 1, 2, 3],
            "scan_num": [100, 101, 102, 103],
            "rt": [1.0, 2.5, 3.0, 4.0],
            "ms_level": [1, 2, 2, 2],
            "precursor_mz": [-1.0, 800.4, 900.5, float("nan")],
            "peak_start_idx": [0, 2, 5, 6],
            "peak_stop_idx": [2, 5, 6, 7],
        }
    )
    peaks = pd.DataFrame(
        {
            "mz": [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0],
            "intensity": [1000.0, 2000.0, 99.0, 100.0, 250.0, 500.0, 1000.0],
        }
    )

    result = _normalize_alpharaw_ms2(spectra, peaks, min_intensity=100.0)

    assert result.columns.tolist() == MS2_COLUMNS
    assert result["scan_number"].tolist() == [101, 102]
    assert result["rt"].tolist() == [2.5, 3.0]
    assert result["precursor_mz"].tolist() == [800.4, 900.5]
    assert result["fragment_mz"].tolist() == [500.0, 600.0]
    assert result["fragment_intensity"].tolist() == [250.0, 500.0]
    assert result["precursor_intensity"].tolist() == [0.0, 0.0]


def test_falls_back_to_one_based_spectrum_index_for_scan_number():
    spectra = pd.DataFrame(
        {
            "spec_idx": [9],
            "rt": [10.0],
            "ms_level": [2],
            "precursor_mz": [1000.0],
            "peak_start_idx": [0],
            "peak_stop_idx": [1],
        }
    )
    peaks = pd.DataFrame({"mz": [500.0], "intensity": [1000.0]})

    result = _normalize_alpharaw_ms2(spectra, peaks, min_intensity=0.0)

    assert result.loc[0, "scan_number"] == 10


def test_rejects_unexpected_alpharaw_schema():
    with pytest.raises(ValueError, match="unexpected table schema"):
        _normalize_alpharaw_ms2(pd.DataFrame(), pd.DataFrame())


def test_extract_configures_alpharaw_for_prm_without_nested_workers(monkeypatch):
    observed = {}

    class FakeThermoRawData:
        def __init__(self, **options):
            observed["options"] = options
            self.spectrum_df = pd.DataFrame(
                {
                    "spec_idx": [0],
                    "rt": [1.0],
                    "ms_level": [2],
                    "precursor_mz": [700.0],
                    "peak_start_idx": [0],
                    "peak_stop_idx": [1],
                }
            )
            self.peak_df = pd.DataFrame({"mz": [300.0], "intensity": [500.0]})

        def import_raw(self, path):
            observed["path"] = path

    monkeypatch.setattr(thermo_raw, "_alpharaw_reader_class", lambda: FakeThermoRawData)

    result = thermo_raw.extract_thermo_ms2("sample.raw", min_intensity=100.0)

    assert observed["path"].endswith("sample.raw")
    assert observed["options"] == {
        "centroided": True,
        "dda": False,
        "process_count": 1,
        "save_as_hdf": False,
    }
    assert result["fragment_mz"].tolist() == [300.0]
