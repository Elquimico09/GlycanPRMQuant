import pandas as pd

from glycanPRMQuant.quantification import reextract_accepted_transitions


def _accepted(transitions=(100.0,), start=1.0, end=3.0):
    rows = []
    for index, mz in enumerate(transitions):
        rows.append(
            {
                "Glycan": "25000",
                "precursor_cluster": 2,
                "rt_feature_id": 7,
                "PrecursorAdduct": "2H",
                "precursor_mz": 700.0,
                "feature_start_rt": start,
                "feature_apex_rt": 2.0,
                "feature_end_rt": end,
                "theoretical_fragment_mz": mz,
                "fragment_mz": mz,
                "Fragment": f"Y{index + 1}",
                "Charge": 1,
                "Adduct": "+H",
                "candidate_score": 90.0,
                "selected": True,
            }
        )
    return pd.DataFrame(rows)


def _raw(rows):
    return pd.DataFrame(
        rows,
        columns=[
            "scan_number",
            "rt",
            "precursor_mz",
            "fragment_mz",
            "fragment_intensity",
            "precursor_intensity",
        ],
    )


def test_reextracts_raw_intensity_for_only_accepted_transitions():
    accepted = _accepted((100.0,))
    raw = _raw(
        [
            (1, 1.0, 700.0, 99.999, 5.0, 1000.0),
            (1, 1.0, 700.0, 100.001, 50.0, 1000.0),
            (1, 1.0, 700.0, 200.0, 10000.0, 1000.0),
            (2, 2.0, 700.0, 150.0, 200.0, 1100.0),
            (3, 3.0, 700.0, 100.002, 7.0, 1200.0),
            (4, 4.0, 700.0, 100.0, 999.0, 1300.0),
        ]
    )

    result = reextract_accepted_transitions(
        accepted,
        raw,
        fragment_mass_tolerance=0.01,
        fragment_mass_tolerance_unit="Da",
        precursor_ppm_tolerance=10.0,
    )

    assert result.peaks["scan_number"].tolist() == [1, 2, 3]
    # Equal mass errors in scan 1 are resolved in favor of the more intense peak.
    assert result.peaks["fragment_intensity"].tolist() == [50.0, 0.0, 7.0]
    assert result.peaks["observed_fragment_mz"].iloc[0] == 100.001
    assert pd.isna(result.peaks["observed_fragment_mz"].iloc[1])
    assert result.peaks["quantification_peak_detected"].tolist() == [True, False, True]
    assert result.peaks["quantification_source"].eq(
        "unfiltered_accepted_transition_reextraction"
    ).all()
    audit = result.audit.iloc[0]
    assert audit["accepted_transition_count"] == 1
    assert audit["quantification_scan_count"] == 3
    assert audit["detected_trace_point_count"] == 2
    assert audit["zero_trace_point_count"] == 1


def test_reextraction_assigns_each_centroid_to_at_most_one_transition():
    accepted = _accepted((100.000, 100.006))
    raw = _raw([(1, 2.0, 700.0, 100.003, 500.0, 1000.0)])

    result = reextract_accepted_transitions(
        accepted,
        raw,
        fragment_mass_tolerance=0.01,
        fragment_mass_tolerance_unit="Da",
        precursor_ppm_tolerance=10.0,
    )

    assert len(result.peaks) == 2
    assert result.peaks["quantification_peak_detected"].sum() == 1
    assert result.peaks["fragment_intensity"].sum() == 500.0


def test_reextraction_rejects_other_precursors_and_outside_feature_bounds():
    accepted = _accepted((100.0,), start=1.0, end=3.0)
    raw = _raw(
        [
            (1, 2.0, 700.0, 100.0, 20.0, 1000.0),
            (2, 2.0, 710.0, 100.0, 200.0, 1000.0),
            (3, 4.0, 700.0, 100.0, 300.0, 1000.0),
        ]
    )

    result = reextract_accepted_transitions(
        accepted,
        raw,
        fragment_mass_tolerance=20.0,
        fragment_mass_tolerance_unit="ppm",
        precursor_ppm_tolerance=10.0,
    )

    assert result.peaks["scan_number"].tolist() == [1]
    assert result.peaks["fragment_intensity"].tolist() == [20.0]


def test_reextraction_retains_candidate_and_feature_metadata():
    accepted = _accepted((100.0,))
    accepted["IUPAC"] = "Gal(b1-4)GlcNAc"
    raw = _raw([(1, 2.0, 700.0, 100.0, 20.0, 1000.0)])

    result = reextract_accepted_transitions(
        accepted,
        raw,
        fragment_mass_tolerance=0.01,
        fragment_mass_tolerance_unit="Da",
        precursor_ppm_tolerance=10.0,
    )

    row = result.peaks.iloc[0]
    assert row["candidate_score"] == 90.0
    assert row["precursor_cluster"] == 2
    assert row["rt_feature_id"] == 7
    assert row["feature_start_rt"] == 1.0
    assert row["feature_end_rt"] == 3.0
    assert "IUPAC" not in result.peaks.columns

    transition = result.transitions.iloc[0]
    assert transition["IUPAC"] == "Gal(b1-4)GlcNAc"
    assert transition["quantification_transition_id"] == row[
        "quantification_transition_id"
    ]
    assert transition["quantification_source"] == (
        "unfiltered_accepted_transition_reextraction"
    )


def test_transition_ids_are_distinct_and_stable_within_feature():
    accepted = _accepted((100.0, 200.0))
    raw = _raw(
        [
            (1, 2.0, 700.0, 100.0, 20.0, 1000.0),
            (1, 2.0, 700.0, 200.0, 30.0, 1000.0),
        ]
    )

    result = reextract_accepted_transitions(
        accepted,
        raw,
        fragment_mass_tolerance=0.01,
        fragment_mass_tolerance_unit="Da",
        precursor_ppm_tolerance=10.0,
    )

    ids = result.transitions["quantification_transition_id"].tolist()
    assert ids == ["25000|2|7|2H|1", "25000|2|7|2H|2"]
    assert result.peaks["quantification_transition_id"].tolist() == ids
