import pandas as pd
import pytest

import glycanPRMQuant.processmzML as process_module

from glycanPRMQuant.processmzML import (
    _count_precursor_conflicts,
    _filter_ms1_to_resolved_assignments,
    _resolve_precursor_conflicts,
    normalize_figure_filetype,
    process_mzml_pipeline,
)


def test_process_mzml_pipeline_is_importable():
    assert callable(process_mzml_pipeline)


def test_pipeline_filters_fragments_but_matches_precursors_from_raw_scans(
    tmp_path, monkeypatch
):
    raw = pd.DataFrame(
        {
            "scan_number": [1] * 20,
            "rt": [5.0] * 20,
            "precursor_mz": [700.0] * 20,
            "fragment_mz": list(range(100, 120)),
            "precursor_intensity": [1000.0] * 20,
            "fragment_intensity": [10.0] * 18 + [100.0, 200.0],
        }
    )
    observed = {}

    monkeypatch.setattr(process_module, "extract_ms2", lambda *args, **kwargs: raw)

    def fake_match_ms1(data, **kwargs):
        observed["match_ms1_rows"] = len(data)
        return pd.DataFrame(columns=["precursor_mz", "Glycan", "Adduct"])

    monkeypatch.setattr(process_module, "matchMS1", fake_match_ms1)

    process_mzml_pipeline(
        "sample.mzML", str(tmp_path), ms2_noise_filter_mode="auto"
    )

    assert observed["match_ms1_rows"] == 20
    audit = pd.read_csv(tmp_path / "ms2_noise_filter_audit.csv")
    assert audit.loc[0, "raw_peak_count"] == 20
    assert audit.loc[0, "retained_peak_count"] == 2


@pytest.mark.parametrize(
    ("value", "expected"),
    [("png", "png"), ("PDF", "pdf"), (".svg", "svg")],
)
def test_normalize_figure_filetype_accepts_supported_formats(value, expected):
    assert normalize_figure_filetype(value) == expected


def test_normalize_figure_filetype_rejects_unsupported_format():
    with pytest.raises(ValueError, match="png, pdf, svg"):
        normalize_figure_filetype("jpg")


def test_resolve_precursor_conflicts_keeps_best_fragment_evidence():
    df = pd.DataFrame(
        {
            "precursor_mz": [500.0, 500.0, 500.0],
            "Glycan": ["25000", "25000", "26000"],
            "fragment_mz": [101.0, 201.0, 301.0],
            "fragment_intensity": [100.0, 200.0, 1000.0],
        }
    )

    resolved = _resolve_precursor_conflicts(df)

    assert resolved["Glycan"].unique().tolist() == ["25000"]


def test_resolve_precursor_conflicts_uses_20_ppm_clusters():
    df = pd.DataFrame(
        {
            "precursor_mz": [500.0, 500.005, 500.02],
            "Glycan": ["25000", "26000", "27000"],
            "fragment_mz": [101.0, 201.0, 301.0],
            "fragment_intensity": [100.0, 1000.0, 2000.0],
        }
    )

    resolved = _resolve_precursor_conflicts(df)

    assert _count_precursor_conflicts(df) == 1
    assert resolved["Glycan"].tolist() == ["26000", "27000"]


def test_filter_ms1_to_resolved_assignments_uses_resolved_ms2_clusters():
    ms1 = pd.DataFrame(
        {
            "precursor_mz": [500.0, 500.005, 500.02],
            "Glycan": ["25000", "26000", "27000"],
            "Adduct": ["2H", "2H", "2H"],
        }
    )
    resolved_ms2 = pd.DataFrame(
        {
            "precursor_mz": [500.005, 500.02],
            "Glycan": ["26000", "27000"],
        }
    )

    resolved_ms1 = _filter_ms1_to_resolved_assignments(ms1, resolved_ms2)

    assert resolved_ms1["Glycan"].tolist() == ["26000", "27000"]
