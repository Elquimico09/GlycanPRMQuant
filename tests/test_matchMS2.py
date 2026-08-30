import pandas as pd
import pytest

from glycanPRMQuant import matchMS2 as ms2_module


def _one_fragment_database():
    return pd.DataFrame(
        {
            "Fragment": ["frag-a"],
            "FragmentType": ["B"],
            "IUPAC": ["test-iupac"],
            "NumericalComposition": ["15000"],
            "Composition": ["{'Hex': 5, 'HexNAc': 1}"],
            "[M+H]+": [101.0],
            "[M+2H]2+": [51.0],
            "[M+NH4]+": [118.0],
            "[M+NH4+H]2+": [59.5],
            "[M+2NH4]2+": [68.0],
        }
    )


def test_fragment_adducts_follow_the_observed_precursor_adduct():
    fragment_db = _one_fragment_database()

    protonated = ms2_module._build_adduct_table(fragment_db, ["2H"])
    mixed = ms2_module._build_adduct_table(fragment_db, ["H+NH4"])
    ammonium = ms2_module._build_adduct_table(fragment_db, ["2NH4"])

    assert set(protonated["Adduct"]) == {"+H", "+2H"}
    assert set(mixed["Adduct"]) == {"+H", "+2H", "+NH4", "+H+NH4"}
    assert set(ammonium["Adduct"]) == {"+H", "+2H", "+NH4", "+2NH4"}
    assert protonated["_PrecursorAdduct"].eq("2H").all()
    assert mixed["_PrecursorAdduct"].eq("H+NH4").all()
    assert ammonium["_PrecursorAdduct"].eq("2NH4").all()


def test_shifted_fragment_decoys_are_deterministic_and_do_not_overlap_targets():
    targets = ms2_module._build_adduct_table(_one_fragment_database(), ["2H"])
    first = ms2_module._build_decoy_adduct_table(targets, 0.02, 1234)
    repeated = ms2_module._build_decoy_adduct_table(targets, 0.02, 1234)
    different = ms2_module._build_decoy_adduct_table(targets, 0.02, 5678)

    pd.testing.assert_frame_equal(first, repeated)
    assert not first["Theo_mz"].equals(different["Theo_mz"])
    assert first["is_decoy"].all()
    assert first["decoy_neutral_mass_shift"].between(1.0, 30.0).all()
    assert len(first) == len(targets)
    for mz in first["Theo_mz"]:
        assert (targets["Theo_mz"] - mz).abs().min() > 0.02


def test_fragment_matching_supports_ppm_tolerance():
    targets = ms2_module._build_adduct_table(_one_fragment_database(), ["2H"])
    observed_mz = 101.0 * (1.0 + 15.0 / 1e6)
    observed = pd.DataFrame(
        {
            "scan_number": [1],
            "rt": [5.0],
            "precursor_mz": [500.0],
            "precursor_intensity": [0.0],
            "fragment_mz": [observed_mz],
            "fragment_intensity": [1000.0],
            "Glycan": ["15000"],
            "PrecursorAdduct": ["2H"],
            "database_mz": [500.0],
            "precursor_ppm_error": [0.0],
        }
    )

    matched = ms2_module._match_adduct_library(
        observed,
        targets,
        20.0,
        is_decoy=False,
        fragment_mass_tolerance_unit="ppm",
    )
    rejected = ms2_module._match_adduct_library(
        observed,
        targets,
        10.0,
        is_decoy=False,
        fragment_mass_tolerance_unit="ppm",
    )

    assert len(matched) == 1
    assert matched["ppm_error"].iloc[0] == pytest.approx(15.0)
    assert matched["fragment_mass_tolerance_value"].iloc[0] == 20.0
    assert matched["fragment_mass_tolerance_unit"].iloc[0] == "ppm"
    assert rejected.empty


def test_ppm_peak_clustering_uses_the_selected_unit():
    mzs = pd.Series([100.0, 100.0015]).to_numpy()

    assert ms2_module.cluster_1d(mzs, 20.0, "ppm").tolist() == [0, 0]
    assert ms2_module.cluster_1d(mzs, 10.0, "ppm").tolist() == [0, 1]


def test_matchms2_returns_target_and_decoy_matches_when_enabled(
    tmp_path, monkeypatch
):
    ms2_module._FRAG_DB_CACHE.clear()
    db_path = tmp_path / "n_glycan_db.csv"
    pd.DataFrame(
        {
            "Numerical Composition": ["15000"],
            "Composition": ["{'Hex': 5, 'HexNAc': 1}"],
            "Condensed IUPAC": ["test-iupac"],
        }
    ).to_csv(db_path, index=False)

    def fake_fragment_glycan(iupac, ion_series="ABCXYZ", max_cleavages=2):
        return pd.DataFrame(
            {
                "name": ["frag-a"],
                "series": ["B"],
                "[M+H]+": [101.0],
                "[M+2H]2+": [51.0],
            }
        )

    monkeypatch.setattr(ms2_module, "fragment_glycan", fake_fragment_glycan)
    targets = ms2_module._build_adduct_table(_one_fragment_database(), ["2H"])
    decoys = ms2_module._build_decoy_adduct_table(targets, 0.001, 1234)
    decoy_mz = decoys.loc[decoys["Adduct"] == "+H", "Theo_mz"].iloc[0]

    matched = ms2_module.matchMS2(
        pd.DataFrame(
            {
                "scan_number": [1, 1],
                "rt": [5.0, 5.0],
                "precursor_mz": [500.0, 500.0],
                "fragment_mz": [101.0, decoy_mz],
                "fragment_intensity": [1000.0, 900.0],
            }
        ),
        pd.DataFrame(
            {
                "precursor_mz": [500.0],
                "Glycan": ["15000"],
                "Adduct": ["2H"],
                "ppm_error": [0.0],
            }
        ),
        precursor_composition="15000",
        fragment_mass_tol=0.001,
        intensity_threshold=100.0,
        ppm_tol=10,
        db_path=str(db_path),
        ion_series="BY",
        max_cleavages=1,
        generate_decoys=True,
        decoy_seed=1234,
    )

    assert set(matched["is_decoy"]) == {False, True}
    assert matched.loc[~matched["is_decoy"], "observed_fragment_mz"].eq(101.0).all()
    assert matched.loc[matched["is_decoy"], "observed_fragment_mz"].eq(decoy_mz).all()


@pytest.mark.parametrize(
    ("precursor_adduct", "fragment_mz", "fragment_adduct"),
    [("H+NH4", 118.0, "+NH4"), ("2NH4", 68.0, "+2NH4")],
)
def test_matchms2_matches_ammonium_products_for_ammonium_precursors(
    tmp_path,
    monkeypatch,
    precursor_adduct,
    fragment_mz,
    fragment_adduct,
):
    ms2_module._FRAG_DB_CACHE.clear()
    db_path = tmp_path / "n_glycan_db.csv"
    pd.DataFrame(
        {
            "Numerical Composition": ["15000"],
            "Composition": ["{'Hex': 5, 'HexNAc': 1}"],
            "Condensed IUPAC": ["test-iupac"],
        }
    ).to_csv(db_path, index=False)

    def fake_fragment_glycan(iupac, ion_series="ABCXYZ", max_cleavages=2):
        return pd.DataFrame(
            {
                "name": ["frag-a"],
                "series": ["B"],
                "[M+H]+": [101.0],
                "[M+2H]2+": [51.0],
                "[M+NH4]+": [118.0],
                "[M+NH4+H]2+": [59.5],
                "[M+2NH4]2+": [68.0],
            }
        )

    monkeypatch.setattr(ms2_module, "fragment_glycan", fake_fragment_glycan)
    matched = ms2_module.matchMS2(
        pd.DataFrame(
            {
                "scan_number": [1],
                "rt": [5.0],
                "precursor_mz": [500.0],
                "fragment_mz": [fragment_mz],
                "fragment_intensity": [1000.0],
            }
        ),
        pd.DataFrame(
            {
                "precursor_mz": [500.0],
                "Glycan": ["15000"],
                "Adduct": [precursor_adduct],
                "ppm_error": [0.0],
            }
        ),
        precursor_composition="15000",
        fragment_mass_tol=0.001,
        intensity_threshold=100.0,
        ppm_tol=10,
        db_path=str(db_path),
        ion_series="BY",
        max_cleavages=1,
    )

    assert not matched.empty
    assert set(matched["PrecursorAdduct"]) == {precursor_adduct}
    assert set(matched["Adduct"]) == {fragment_adduct}


def test_matchms2_generates_fragments_and_selects_best_iupac(tmp_path, monkeypatch):
    ms2_module._FRAG_DB_CACHE.clear()
    db_path = tmp_path / "n_glycan_db.csv"
    pd.DataFrame(
        {
            "Numerical Composition": ["25000", "25000"],
            "Composition": ["{'Hex': 5, 'HexNAc': 2}", "{'Hex': 5, 'HexNAc': 2}"],
            "Condensed IUPAC": ["best-iupac", "other-iupac"],
        }
    ).to_csv(db_path, index=False)

    def fake_fragment_glycan(iupac, ion_series="ABCXYZ", max_cleavages=2):
        if iupac == "best-iupac":
            return pd.DataFrame(
                {
                    "name": ["frag-a", "frag-b"],
                    "series": ["B", "Y"],
                    "[M+H]+": [101.0, 201.0],
                    "[M+2H]2+": [51.0, 101.0],
                    "[M+NH4]+": [118.0, 218.0],
                    "[M+NH4+H]2+": [59.5, 109.5],
                }
            )
        return pd.DataFrame(
            {
                "name": ["frag-c"],
                "series": ["B"],
                "[M+H]+": [301.0],
                "[M+2H]2+": [151.0],
                "[M+NH4]+": [318.0],
                "[M+NH4+H]2+": [159.5],
            }
        )

    monkeypatch.setattr(ms2_module, "fragment_glycan", fake_fragment_glycan)

    ms2_data = pd.DataFrame(
        {
            "scan_number": [1, 1, 1],
            "rt": [5.0, 5.0, 5.0],
            "precursor_mz": [500.0, 500.0, 500.0],
            "fragment_mz": [101.0, 201.0, 999.0],
            "fragment_intensity": [1000.0, 900.0, 500.0],
        }
    )
    precursor_matches = pd.DataFrame(
        {
            "precursor_mz": [500.0, 500.0],
            "Glycan": ["25000", "25000"],
            "Adduct": ["2H", "2H"],
        }
    )

    matched = ms2_module.matchMS2(
        ms2_data,
        precursor_matches,
        precursor_composition="25000",
        fragment_mass_tol=0.01,
        intensity_threshold=100.0,
        ppm_tol=10,
        db_path=str(db_path),
        ion_series="BY",
        max_cleavages=1,
    )

    assert not matched.empty
    assert set(matched["IUPAC"]) == {"best-iupac"}
    assert matched["IUPAC_match_count"].iloc[0] == 2
    assert matched["NumericalComposition"].iloc[0] == "25000"
    assert matched["fragment_intensity"].max() == 1000.0
