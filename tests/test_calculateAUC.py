import pandas as pd

from glycanPRMQuant.calculateAUC import calculateAUC, calculate_feature_auc


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


def test_calculate_feature_auc_keeps_separate_rt_features():
    rows = []
    for feature_id, rt_offset, score in [(1, 0.0, 80.0), (2, 10.0, 70.0)]:
        for scan, (rt, intensity) in enumerate(
            zip([0, 1, 2, 3, 4], [0, 10, 100, 10, 0]), start=1
        ):
            rows.append(
                {
                    "Glycan": "25000",
                    "precursor_cluster": 3,
                    "rt_feature_id": feature_id,
                    "PrecursorAdduct": "2H",
                    "scan_number": scan + feature_id * 10,
                    "rt": rt + rt_offset,
                    "fragment_intensity": intensity,
                    "candidate_score": score,
                    "target_decoy_pass": True,
                }
            )

    result = calculate_feature_auc(
        pd.DataFrame(rows),
        smoothing_window=0,
        rel_height=0.5,
        rel_height_mode="height",
    )

    assert result["rt_feature_id"].tolist() == [1, 2]
    assert result["peak_rt"].tolist() == [2.0, 12.0]
    assert result["candidate_score"].tolist() == [80.0, 70.0]
    assert (result["AUC"] > 0).all()


def test_calculate_feature_auc_reports_reextraction_provenance():
    df = pd.DataFrame(
        {
            "Glycan": ["25000"] * 6,
            "precursor_cluster": [1] * 6,
            "rt_feature_id": [2] * 6,
            "PrecursorAdduct": ["2H"] * 6,
            "scan_number": [1, 1, 2, 2, 3, 3],
            "rt": [1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
            "Fragment": ["Y1", "Y2"] * 3,
            "fragment_mz": [100.0, 200.0] * 3,
            "fragment_intensity": [0.0, 0.0, 10.0, 20.0, 0.0, 0.0],
            "quantification_peak_detected": [False, False, True, True, False, False],
            "quantification_source": [
                "unfiltered_accepted_transition_reextraction"
            ]
            * 6,
        }
    )

    result = calculate_feature_auc(
        df, smoothing_window=0, rel_height=0.5, rel_height_mode="height"
    )

    assert result.loc[0, "accepted_transition_count"] == 2
    assert result.loc[0, "quantification_scan_count"] == 3
    assert result.loc[0, "quantification_trace_point_count"] == 6
    assert result.loc[0, "detected_trace_point_count"] == 2
    assert result.loc[0, "zero_trace_point_count"] == 4
