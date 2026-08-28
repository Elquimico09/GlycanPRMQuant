import pandas as pd

from glycanPRMQuant import matchMS2 as ms2_module


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
