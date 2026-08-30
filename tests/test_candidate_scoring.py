import numpy as np
import pandas as pd
import pytest

from glycanPRMQuant.candidate_scoring import (
    CandidateScoringConfig,
    _apply_target_decoy_statistics,
    _finalize_candidate_scores,
    score_and_resolve_candidates,
)


def _scan_table(profile):
    return pd.DataFrame(
        {
            "scan_number": np.arange(1, len(profile) + 1),
            "rt": np.arange(len(profile), dtype=float),
            "precursor_mz": 500.0,
            "precursor_intensity": profile,
            "fragment_mz": 900.0,
            "fragment_intensity": np.maximum(np.asarray(profile, dtype=float) * 20.0, 1.0),
        }
    )


def _candidate_rows(
    glycan,
    fragment_mzs,
    profiles,
    mass_error=0.001,
    scan_offset=0,
    precursor_error=1.0,
):
    rows = []
    for fragment_index, (fragment_mz, profile) in enumerate(zip(fragment_mzs, profiles)):
        for position, intensity in enumerate(profile):
            if intensity <= 0:
                continue
            scan = position + 1 + scan_offset
            rows.append(
                {
                    "scan_number": scan,
                    "rt": float(scan - 1),
                    "precursor_mz": 500.0,
                    "precursor_intensity": 0.0,
                    "fragment_intensity": float(intensity),
                    "Glycan": glycan,
                    "PrecursorAdduct": "2H",
                    "Fragment": f"frag-{glycan}-{fragment_index}",
                    "fragment_annotation": f"frag-{glycan}-{fragment_index}",
                    "Charge": 1,
                    "Adduct": "H",
                    "fragment_mz": fragment_mz,
                    "observed_fragment_mz": fragment_mz + mass_error,
                    "mz_diff": mass_error,
                    "precursor_ppm_error": precursor_error,
                }
            )
    return rows


def _permissive_config(**changes):
    values = {
        "minimum_candidate_score": 0.0,
        "minimum_distinct_fragments": 2,
        "minimum_explained_intensity": 0.0,
        "minimum_discriminative_evidence_difference": 3.0,
    }
    values.update(changes)
    return CandidateScoringConfig(**values)


def test_default_minimum_explained_intensity_is_one_percent():
    assert CandidateScoringConfig().minimum_explained_intensity == 0.01


def test_repeated_rows_do_not_inflate_candidate_score():
    profile = [10, 50, 100, 50, 10]
    rows = _candidate_rows("25000", [101.0, 201.0], [profile, profile])
    rows += _candidate_rows("26000", [301.0, 401.0], [profile, profile], mass_error=0.01)
    matched = pd.DataFrame(rows)

    baseline = score_and_resolve_candidates(
        matched, _scan_table(profile), _permissive_config()
    ).candidate_scores
    duplicated = score_and_resolve_candidates(
        pd.concat([matched, matched.iloc[[0]]], ignore_index=True),
        _scan_table(profile),
        _permissive_config(),
    ).candidate_scores

    base_score = baseline.set_index("Glycan").loc["25000", "candidate_score"]
    duplicate_score = duplicated.set_index("Glycan").loc["25000", "candidate_score"]
    assert duplicate_score == base_score
    assert duplicated.set_index("Glycan").loc["25000", "distinct_fragment_count"] == 2


def test_score_uses_specificity_explained_intensity_and_mass_accuracy():
    profile = [10, 50, 100, 50, 10]
    strong = _candidate_rows("25000", [101.0, 201.0], [profile, profile], mass_error=0.001)
    weak_profile = [50, 20, 5, 0, 0]
    weak = _candidate_rows(
        "26000", [101.0, 301.0], [profile, weak_profile], mass_error=0.019
    )
    for row in weak:
        if row["fragment_mz"] == 101.0:
            row["observed_fragment_mz"] = 101.001
    scans = _scan_table(profile)
    scans["precursor_intensity"] = 0.0
    result = score_and_resolve_candidates(
        pd.DataFrame(strong + weak),
        scans,
        _permissive_config(coisolation_score=1.01),
    )
    scores = result.candidate_scores.set_index("Glycan")

    assert scores.loc["25000", "candidate_specific_fragment_count"] == 1
    assert 0 < scores.loc["25000", "explained_intensity_fraction"] <= 1
    assert (
        scores.loc["25000", "fragment_mass_accuracy_score"]
        > scores.loc["26000", "fragment_mass_accuracy_score"]
    )
    assert scores.loc["25000", "candidate_score"] > scores.loc["26000", "candidate_score"]
    assert scores.loc["25000", "coelution_reference"] == "common_fragments"
    assert set(result.candidate_scores["resolution_status"]) == {"resolved"}
    assert set(result.resolved_rows["Glycan"]) == {"25000"}


def test_discriminative_evidence_difference_preserves_ambiguous_candidates():
    profile = [10, 50, 100, 50, 10]
    rows = _candidate_rows("25000", [101.0, 201.0], [profile, profile])
    rows += _candidate_rows("26000", [301.0, 401.0], [profile, profile])
    result = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(profile),
        _permissive_config(
            minimum_discriminative_evidence_difference=10.0,
            coisolation_score=1.01,
        ),
    )

    assert set(result.candidate_scores["resolution_status"]) == {"ambiguous"}
    assert result.resolved_rows.empty
    assert set(result.reported_rows["Glycan"]) == {"25000", "26000"}
    assert not result.candidate_scores["selected"].any()
    assert result.candidate_scores["quantification_weight"].eq(0.0).all()
    assert result.candidate_scores["runner_up_score"].notna().all()
    assert result.candidate_scores["discriminative_evidence_difference"].eq(0.0).all()
    assert "score_ratio" not in result.candidate_scores


def test_chromatographic_coherence_penalizes_shifted_fragments():
    precursor_profile = [5, 30, 100, 30, 5, 1, 1]
    aligned = [5, 30, 100, 30, 5, 0, 0]
    shifted = [0, 0, 5, 30, 100, 30, 5]
    rows = _candidate_rows("25000", [101.0, 201.0], [aligned, aligned])
    rows += _candidate_rows("26000", [301.0, 401.0], [shifted, shifted])
    scores = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(precursor_profile),
        _permissive_config(coisolation_score=1.01),
    ).candidate_scores.set_index("Glycan")

    assert scores.loc["25000", "coelution_evaluable"]
    assert scores.loc["26000", "coelution_evaluable"]
    assert scores.loc["25000", "coelution_score"] > scores.loc["26000", "coelution_score"]


def test_unresolved_candidates_are_reported_but_not_selected_for_quantification():
    profile = [10, 50, 100, 50, 10]
    rows = _candidate_rows("25000", [101.0, 201.0], [profile, profile])
    rows += _candidate_rows("26000", [301.0, 401.0], [profile, profile])

    insufficient = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(profile),
        _permissive_config(minimum_candidate_score=101.0),
    )
    disabled = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(profile),
        _permissive_config(),
        resolve=False,
    )

    assert set(insufficient.candidate_scores["resolution_status"]) == {
        "insufficient_evidence"
    }
    assert set(disabled.candidate_scores["resolution_status"]) == {
        "resolution_disabled"
    }
    for result in (insufficient, disabled):
        assert result.resolved_rows.empty
        assert not result.candidate_scores["selected"].any()
        assert result.reported_rows["reported"].all()


def test_coeluting_candidate_specific_fragments_are_marked_as_possible_coisolation():
    profile = [10, 50, 100, 50, 10]
    rows = _candidate_rows("25000", [101.0, 201.0], [profile, profile])
    rows += _candidate_rows("26000", [301.0, 401.0], [profile, profile])
    result = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(profile),
        _permissive_config(),
    )

    assert set(result.candidate_scores["resolution_status"]) == {"possible_coisolation"}
    assert result.resolved_rows.empty
    assert set(result.reported_rows["Glycan"]) == {"25000", "26000"}
    assert not result.candidate_scores["selected"].any()


def test_rt_features_keep_isobaric_candidates_at_different_elution_times():
    precursor_profile = [1, 10, 100, 10, 1, 1, 10, 80, 10, 1]
    first_peak = [1, 10, 100, 10, 1]
    second_peak = [1, 10, 80, 10, 1]
    rows = _candidate_rows("25000", [101.0, 201.0], [first_peak, first_peak])
    rows += _candidate_rows(
        "26000", [301.0, 401.0], [second_peak, second_peak], scan_offset=5
    )
    result = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(precursor_profile),
        _permissive_config(),
    )

    assert result.candidate_scores["rt_feature_id"].nunique() == 2
    assert "uncontested" in set(result.candidate_scores["resolution_status"])
    assert set(result.resolved_rows["Glycan"]) == {"25000", "26000"}


def test_uncontested_candidate_must_pass_minimum_evidence_requirements():
    profile = [10, 50, 100, 50, 10]
    rows = _candidate_rows("25000", [101.0], [profile])
    result = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(profile),
        CandidateScoringConfig(minimum_distinct_fragments=2),
    )

    assert set(result.candidate_scores["resolution_status"]) == {
        "insufficient_evidence"
    }
    assert not result.candidate_scores["selected"].any()
    assert result.resolved_rows.empty
    assert not result.reported_rows.empty


def test_no_candidate_specific_fragments_get_a_structural_status():
    profile = [10, 50, 100, 50, 10]
    rows = _candidate_rows("25000", [101.0, 201.0], [profile, profile])
    rows += _candidate_rows("26000", [101.0, 201.0], [profile, profile])

    result = score_and_resolve_candidates(
        pd.DataFrame(rows),
        _scan_table(profile),
        _permissive_config(coisolation_score=1.01),
    )

    assert result.candidate_scores["candidate_specific_fragment_count"].eq(0).all()
    assert set(result.candidate_scores["resolution_status"]) == {
        "no_discriminating_fragment_evidence"
    }
    assert not result.candidate_scores["selected"].any()
    assert result.resolved_rows.empty


def test_zero_specific_mass_outlier_is_pruned_and_audited():
    profile = [10, 50, 100, 50, 10]
    supported = _candidate_rows(
        "25000",
        [101.0, 201.0],
        [profile, profile],
        precursor_error=0.4,
    )
    unsupported = _candidate_rows(
        "26000",
        [101.0],
        [profile],
        precursor_error=6.4,
    )
    result = score_and_resolve_candidates(
        pd.DataFrame(supported + unsupported),
        _scan_table(profile),
        _permissive_config(coisolation_score=1.01),
    )
    scores = result.candidate_scores.set_index("Glycan")

    assert not scores.loc["25000", "mass_outlier_pruned"]
    assert scores.loc["26000", "mass_outlier_pruned"]
    assert not scores.loc["26000", "eligible_for_resolution"]
    assert scores.loc["26000", "candidate_rejection_reason"] == (
        "no_specific_evidence_and_precursor_mass_outlier"
    )
    assert scores.loc["26000", "precursor_error_delta_ppm"] > scores.loc[
        "26000", "mass_outlier_threshold_ppm"
    ]
    assert set(scores["resolution_status"]) == {"resolved_after_mass_pruning"}
    assert set(result.resolved_rows["Glycan"]) == {"25000"}


def _component_row(cluster, glycan, error, coelution_score=0.8, evaluable=True):
    return {
        "_precursor_cluster": cluster,
        "_rt_feature": cluster,
        "Glycan": glycan,
        "median_signed_precursor_ppm_error": error,
        "median_precursor_ppm_error": abs(error),
        "candidate_specific_fragment_count": 1,
        "explained_intensity_score": 0.5,
        "fragment_support_score": 0.5,
        "fragment_mass_accuracy_score": 0.9,
        "coelution_score": coelution_score,
        "coelution_evaluable": evaluable,
    }


def test_precursor_accuracy_is_a_calibrated_within_feature_likelihood():
    rows = [
        _component_row(cluster, f"cal-{cluster}", error)
        for cluster, error in enumerate(
            [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
        )
    ]
    rows += [
        _component_row(10, "near", 0.6),
        _component_row(10, "far", 6.0),
    ]
    finalized = _finalize_candidate_scores(
        pd.DataFrame(rows), CandidateScoringConfig()
    ).set_index("Glycan")

    assert finalized["precursor_calibration_source"].eq("empirical_iqr").all()
    assert finalized.loc["near", "precursor_relative_likelihood"] == 1.0
    assert finalized.loc["far", "precursor_relative_likelihood"] < 1e-6
    assert finalized.loc["far", "precursor_log_likelihood"] < 0
    assert 0.4 < finalized.loc["near", "precursor_likelihood_center_ppm"] < 0.9
    assert finalized.loc["near", "precursor_likelihood_sigma_ppm"] >= 0.25


def test_common_components_cancel_from_discriminative_difference():
    rows = [
        _component_row(0, "A0", 0.5, coelution_score=0.95),
        _component_row(0, "B0", 0.5, coelution_score=0.95),
        _component_row(1, "A1", 0.5, coelution_score=0.10),
        _component_row(1, "B1", 0.5, coelution_score=0.10),
    ]
    rows[0]["explained_intensity_score"] = 0.6
    rows[2]["explained_intensity_score"] = 0.6
    finalized = _finalize_candidate_scores(
        pd.DataFrame(rows), CandidateScoringConfig()
    )
    differences = []
    for _, feature in finalized.groupby("_precursor_cluster"):
        differences.append(float(feature["discriminative_score"].max() - feature["discriminative_score"].min()))

    assert differences[0] == differences[1]


def test_missing_coelution_is_neutral_for_every_candidate_in_feature():
    rows = [
        _component_row(0, "evaluable", 0.5, coelution_score=0.95, evaluable=True),
        _component_row(0, "missing", 0.5, coelution_score=np.nan, evaluable=False),
    ]
    finalized = _finalize_candidate_scores(
        pd.DataFrame(rows), CandidateScoringConfig()
    )

    assert not finalized["coelution_comparable"].any()
    assert finalized["coelution_component_used"].eq(0.5).all()


def test_fragment_mass_accuracy_weight_decreases_as_tolerance_widens():
    rows = [
        _component_row(0, "accurate", 0.5),
        _component_row(0, "inaccurate", 0.5),
    ]
    rows[0]["fragment_mass_accuracy_score"] = 1.0
    rows[1]["fragment_mass_accuracy_score"] = 0.0

    orbitrap = _finalize_candidate_scores(
        pd.DataFrame(rows),
        CandidateScoringConfig(
            fragment_mass_tolerance=20.0,
            fragment_mass_tolerance_unit="ppm",
        ),
    ).set_index("Glycan")
    ion_trap = _finalize_candidate_scores(
        pd.DataFrame(rows),
        CandidateScoringConfig(
            fragment_mass_tolerance=0.5,
            fragment_mass_tolerance_unit="Da",
        ),
    ).set_index("Glycan")

    orbitrap_gap = (
        orbitrap.loc["accurate", "candidate_score"]
        - orbitrap.loc["inaccurate", "candidate_score"]
    )
    ion_trap_gap = (
        ion_trap.loc["accurate", "candidate_score"]
        - ion_trap.loc["inaccurate", "candidate_score"]
    )

    assert orbitrap["fragment_mass_accuracy_reliability"].eq(1.0).all()
    assert ion_trap["fragment_tolerance_equivalent_ppm"].eq(1000.0).all()
    assert ion_trap["fragment_mass_accuracy_reliability"].iloc[0] == pytest.approx(
        (20.0 / 1000.0) ** 0.5
    )
    assert ion_trap["effective_fragment_mass_accuracy_weight"].iloc[0] < 0.03
    assert ion_trap_gap < orbitrap_gap


def test_target_decoy_competition_rejects_a_target_that_is_outscored():
    profile = [10, 50, 100, 50, 10]
    weak = [1, 1, 1, 1, 1]
    target = _candidate_rows(
        "25000", [101.0, 201.0], [weak, weak], mass_error=0.019
    )
    decoy = _candidate_rows(
        "25000", [301.0, 401.0], [profile, profile], mass_error=0.001
    )
    result = score_and_resolve_candidates(
        pd.DataFrame(target),
        _scan_table(profile),
        _permissive_config(maximum_assignment_q_value=1.0),
        decoy_matched_data=pd.DataFrame(decoy),
    )
    scores = result.candidate_scores

    assert scores["target_decoy_evaluated"].all()
    assert not scores["target_decoy_winner"].any()
    assert set(scores["resolution_status"]) == {"decoy_outcompeted"}
    assert not scores["selected"].any()
    assert scores["best_decoy_score"].gt(scores["target_assignment_score"]).all()
    assert not result.decoy_scores.empty
    assert set(result.target_decoy_competitions["competition_winner"]) == {"decoy"}
    assert result.decoy_scores["precursor_likelihood_center_ppm"].iloc[0] == (
        scores["precursor_likelihood_center_ppm"].iloc[0]
    )
    assert result.decoy_scores["precursor_likelihood_sigma_ppm"].iloc[0] == (
        scores["precursor_likelihood_sigma_ppm"].iloc[0]
    )
    assert result.decoy_scores["fragment_tolerance_equivalent_ppm"].iloc[0] == (
        scores["fragment_tolerance_equivalent_ppm"].iloc[0]
    )
    assert result.decoy_scores["fragment_mass_accuracy_reliability"].iloc[0] == (
        scores["fragment_mass_accuracy_reliability"].iloc[0]
    )


def test_target_decoy_q_values_use_the_global_competition_distribution():
    target_scores = pd.DataFrame(
        {
            "_precursor_cluster": np.arange(20),
            "_rt_feature": np.arange(20),
            "candidate_rank": 1,
            "candidate_score": np.linspace(80.0, 60.0, 20),
            "selected": True,
            "resolution_status": "uncontested",
        }
    )
    scored, competitions = _apply_target_decoy_statistics(
        target_scores,
        pd.DataFrame(),
        CandidateScoringConfig(maximum_assignment_q_value=0.05),
    )

    assert scored["assignment_q_value"].eq(0.05).all()
    assert scored["target_decoy_pass"].all()
    assert scored["selected"].all()
    assert competitions["competition_winner"].eq("target").all()
