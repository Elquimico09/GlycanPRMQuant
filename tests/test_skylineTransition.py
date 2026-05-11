import pandas as pd

from glycanPRMQuant.skylineTransition import dedupe_ms2_fragments


def test_dedupe_ms2_fragments_clusters_by_precursor(tmp_path):
    input_csv = tmp_path / "ms2.csv"
    output_csv = tmp_path / "deduped.csv"
    pd.DataFrame(
        {
            "precursor_mz": [500.0, 500.0, 600.0],
            "fragment_mz": [101.000, 101.005, 201.0],
            "Charge": [1, 1, 2],
        }
    ).to_csv(input_csv, index=False)

    out = dedupe_ms2_fragments(str(input_csv), str(output_csv), mz_tol=0.02)

    assert output_csv.exists()
    assert len(out) == 2
    assert set(out["Charge"]) == {1, 2}
