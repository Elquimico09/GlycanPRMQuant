import pandas as pd
import pytest

from glycanPRMQuant.processmzML import (
    _count_precursor_conflicts,
    _filter_ms1_to_resolved_assignments,
    _resolve_precursor_conflicts,
    normalize_figure_filetype,
    process_mzml_pipeline,
)


def test_process_mzml_pipeline_is_importable():
    assert callable(process_mzml_pipeline)


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
