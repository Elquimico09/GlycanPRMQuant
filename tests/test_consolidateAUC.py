import pandas as pd

from glycanPRMQuant.consolidateAUC import consolidate_auc_results


def test_consolidate_auc_results_merges_sample_tables(tmp_path):
    sample_a = tmp_path / "sample_a"
    sample_b = tmp_path / "sample_b"
    sample_a.mkdir()
    sample_b.mkdir()
    pd.DataFrame({"Glycan": ["25000", "26000"], "AUC": [10.0, 20.0]}).to_csv(
        sample_a / "sample_a_auc_values.csv", index=False
    )
    pd.DataFrame({"Glycan": ["25000"], "AUC": [30.0]}).to_csv(
        sample_b / "sample_b_auc_values.csv", index=False
    )
    output_csv = tmp_path / "combined.csv"

    consolidate_auc_results(str(tmp_path), str(output_csv))

    out = pd.read_csv(output_csv)
    assert out.columns.tolist() == ["Glycan", "sample_a", "sample_b"]
    assert out.loc[out["Glycan"] == 25000, "sample_b"].iloc[0] == 30.0
