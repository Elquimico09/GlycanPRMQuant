import pandas as pd

from glycanPRMQuant import matchMS1 as ms1_module
from glycanPRMQuant.constants import PROTON_MASS


def test_calculate_precursor_masses_groups_by_composition(monkeypatch):
    raw = pd.DataFrame(
        {
            "Composition": ["{'Hex': 1}", "{'Hex': 1}", "{'Hex': 2}"],
            "Condensed IUPAC": ["iupac-a", "iupac-b", "iupac-c"],
            "Numerical Composition": ["11000", "11000", "12000"],
        }
    )
    calls = []

    def fake_calculate_mass(iupac, verbose=False):
        calls.append(iupac)
        return {"iupac-a": 100.0, "iupac-b": 999.0, "iupac-c": 200.0}[iupac]

    monkeypatch.setattr(ms1_module, "calculate_mass", fake_calculate_mass)

    db = ms1_module._calculate_precursor_masses(raw)

    assert calls == ["iupac-a", "iupac-c"]
    assert db["Composition"].tolist() == ["11000", "12000"]
    assert db["[M]"].tolist() == [100.0, 200.0]


def test_matchms1_matches_precomputed_database(tmp_path):
    ms1_module._DB_CACHE.clear()
    neutral_mass = 500.0
    precursor_mz = (neutral_mass + 2 * PROTON_MASS) / 2
    db_path = tmp_path / "precursors.csv"
    pd.DataFrame({"Composition": ["25000"], "[M]": [neutral_mass]}).to_csv(db_path, index=False)
    ms1_data = pd.DataFrame({"precursor_mz": [precursor_mz, 1500.0]})

    matched = ms1_module.matchMS1(
        ms1_data,
        ppm_tol=10,
        db_path=str(db_path),
    )

    assert len(matched) == 1
    assert matched.iloc[0]["Glycan"] == "25000"
    assert matched.iloc[0]["Adduct"] == "2H"
    assert matched.iloc[0]["precursor_mz"] < 400
