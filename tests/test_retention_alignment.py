import pandas as pd

from glycanPRMQuant.retention_alignment import (
    build_consensus_peak_result,
    consolidate_consensus_peak_results,
)


def _feature_row(sample, glycan, feature, rt, score, auc):
    return {
        "sample": sample,
        "Glycan": glycan,
        "precursor_cluster": feature,
        "rt_feature_id": feature,
        "PrecursorAdduct": "2H",
        "peak_rt": rt,
        "start_rt": rt - 0.2,
        "end_rt": rt + 0.2,
        "AUC": auc,
        "candidate_score": score,
        "assignment_q_value": 0.01,
        "target_decoy_pass": True,
    }


def test_consensus_peak_selection_prefers_reproducible_score():
    samples = ["s1", "s2", "s3", "s4", "s5"]
    drifts = [-0.2, -0.1, 0.0, 0.1, 0.2]
    rows = []
    feature = 0
    for sample, drift in zip(samples, drifts):
        # Stable single-peak glycans provide RT-alignment landmarks.
        for glycan, reference_rt in [
            ("10000", 5.0),
            ("20000", 20.0),
            ("30000", 30.0),
        ]:
            feature += 1
            rows.append(
                _feature_row(
                    sample, glycan, feature, reference_rt + drift, 90.0, 50.0
                )
            )

    stable_scores = [82.0, 80.0, 81.0, 83.0, 79.0]
    unstable_scores = [95.0, 45.0, 50.0, 90.0, 40.0]
    for index, (sample, drift) in enumerate(zip(samples, drifts)):
        feature += 1
        rows.append(
            _feature_row(
                sample,
                "25000",
                feature,
                10.0 + drift,
                stable_scores[index],
                100.0 + index,
            )
        )
        feature += 1
        rows.append(
            _feature_row(
                sample,
                "25000",
                feature,
                15.0 + drift,
                unstable_scores[index],
                500.0 + index,
            )
        )

    result = build_consensus_peak_result(
        pd.DataFrame(rows),
        rt_tolerance_minutes=0.3,
        minimum_replicate_fraction=0.8,
    )

    glycan_groups = result.peak_groups.loc[
        result.peak_groups["Glycan"] == "25000"
    ]
    assert len(glycan_groups) == 2
    selected = glycan_groups.loc[glycan_groups["consensus_selected"]].iloc[0]
    assert abs(selected["median_aligned_rt"] - 10.0) < 0.05
    assert selected["replicate_coverage"] == 1.0
    assert selected["consensus_score"] > glycan_groups.loc[
        ~glycan_groups["consensus_selected"], "consensus_score"
    ].iloc[0]

    combined = result.selected_auc.set_index("Glycan")
    assert combined.loc["25000", "s1"] == 100.0
    assert combined.loc["25000", "s5"] == 104.0
    assert (result.alignment_models["alignment_anchor_count"] >= 3).all()


def test_consensus_peak_requires_minimum_run_coverage():
    rows = [
        _feature_row("s1", "25000", 1, 10.0, 90.0, 100.0),
        _feature_row("s2", "30000", 2, 20.0, 90.0, 100.0),
    ]

    result = build_consensus_peak_result(
        pd.DataFrame(rows),
        rt_tolerance_minutes=0.3,
        minimum_replicate_fraction=0.8,
    )

    assert not result.peak_groups["consensus_eligible"].any()
    assert result.selected_auc.empty


def test_expected_empty_run_counts_against_peak_coverage():
    rows = [
        _feature_row("s1", "25000", 1, 10.0, 90.0, 100.0),
        _feature_row("s2", "25000", 2, 10.1, 89.0, 110.0),
    ]

    result = build_consensus_peak_result(
        pd.DataFrame(rows),
        rt_tolerance_minutes=0.3,
        minimum_replicate_fraction=0.8,
        expected_samples=["s1", "s2", "s3"],
    )

    peak = result.peak_groups.iloc[0]
    assert peak["replicate_coverage"] == 2 / 3
    assert not peak["consensus_eligible"]
    assert result.selected_auc.empty


def test_consolidation_writes_selected_and_audit_tables(tmp_path):
    for sample, rt, auc in [("s1", 10.0, 100.0), ("s2", 10.1, 110.0)]:
        sample_dir = tmp_path / sample
        sample_dir.mkdir()
        pd.DataFrame([_feature_row(sample, "25000", 1, rt, 90.0, auc)]).drop(
            columns="sample"
        ).to_csv(sample_dir / f"{sample}_feature_auc_values.csv", index=False)

    output = tmp_path / "combined_auc_values.csv"
    result = consolidate_consensus_peak_results(
        str(tmp_path),
        str(output),
        rt_tolerance_minutes=0.3,
        minimum_replicate_fraction=1.0,
    )

    assert output.is_file()
    combined = pd.read_csv(output, dtype={"Glycan": str})
    assert combined.loc[0, "s1"] == 100.0
    assert combined.loc[0, "s2"] == 110.0
    assert result.peak_groups["consensus_selected"].sum() == 1
    for filename in [
        "aligned_feature_auc_values.csv",
        "consensus_peak_groups.csv",
        "combined_all_feature_auc_values.csv",
        "retention_time_alignment.csv",
    ]:
        assert (tmp_path / filename).is_file()
