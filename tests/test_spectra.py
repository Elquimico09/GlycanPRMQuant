from pathlib import Path

import pandas as pd
import pytest

from glycanPRMQuant import spectra


def test_detects_supported_input_types_case_insensitively():
    assert spectra.detect_input_type("sample.mzML") == "mzml"
    assert spectra.detect_input_type("sample.RAW") == "thermo_raw"


def test_rejects_unsupported_input_type():
    with pytest.raises(ValueError, match="Unsupported mass-spectrometry file"):
        spectra.detect_input_type("sample.mgf")


def test_rejects_mixed_input_batch():
    with pytest.raises(ValueError, match="cannot mix"):
        spectra.validate_input_file_types(["one.raw", "two.mzML"])


def test_dispatches_mzml(monkeypatch):
    expected = pd.DataFrame({"source": ["mzml"]})
    monkeypatch.setattr(spectra, "extract_mzml_ms2", lambda path, min_intensity: expected)

    observed = spectra.extract_ms2(Path("sample.mzML"), min_intensity=12.0)

    assert observed is expected


def test_dispatches_raw_lazily(monkeypatch):
    expected = pd.DataFrame({"source": ["raw"]})
    monkeypatch.setattr(
        "glycanPRMQuant.thermo_raw.extract_thermo_ms2",
        lambda path, min_intensity: expected,
    )

    observed = spectra.extract_ms2(Path("sample.raw"), min_intensity=12.0)

    assert observed is expected
