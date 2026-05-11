import pandas as pd

from glycanPRMQuant.calculateAUC import calculateAUC


def test_calculate_auc_returns_per_adduct_and_total_tables():
    df = pd.DataFrame(
        {
            "Glycan": ["25000"] * 10,
            "Adduct": ["2H"] * 5 + ["3H"] * 5,
            "scan_number": [1, 2, 3, 4, 5] * 2,
            "rt": [0.0, 1.0, 2.0, 3.0, 4.0] * 2,
            "fragment_intensity": [0.0, 10.0, 100.0, 10.0, 0.0, 0.0, 5.0, 50.0, 5.0, 0.0],
        }
    )

    per_adduct, total = calculateAUC(
        df,
        rel_height=0.5,
        rel_height_mode="height",
        smoothing_window=0,
        plot=False,
    )

    assert set(per_adduct["Adduct"]) == {"2H", "3H"}
    assert total.loc[0, "Glycan"] == "25000"
    assert total.loc[0, "AUC"] > 0
