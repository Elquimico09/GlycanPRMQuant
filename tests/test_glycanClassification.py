import pandas as pd

from glycanPRMQuant.glycanClassification import classifyGlycan


def test_classify_glycan_adds_class_and_type_columns(tmp_path):
    csv_path = tmp_path / "auc.csv"
    pd.DataFrame({"Glycan": [25000, 43110, 43001], "AUC": [1.0, 2.0, 3.0]}).to_csv(
        csv_path, index=False
    )

    out = classifyGlycan(str(csv_path))

    assert {"Class", "Type", "pos1", "pos5"}.issubset(out.columns)
    assert out.loc[out["Glycan"] == 25000, "Class"].iloc[0] == "high mannose"
    assert out.loc[out["Glycan"] == 43110, "Class"].iloc[0] == "fucosylated"
    assert out.loc[out["Glycan"] == 43001, "Class"].iloc[0] == "sialylated"
