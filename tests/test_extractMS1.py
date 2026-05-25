from pathlib import Path
import ast

import pandas as pd

from glycanPRMQuant.constants import DEFAULT_PRECURSOR_DB


def test_default_precursor_database_is_packaged():
    db_path = Path(DEFAULT_PRECURSOR_DB)

    assert db_path.exists()
    db = pd.read_csv(db_path, nrows=1)
    assert {"Condensed IUPAC", "Composition", "Numerical Composition"}.issubset(db.columns)


def test_n_glycan_numerical_composition_matches_composition():
    db = pd.read_csv(DEFAULT_PRECURSOR_DB, dtype={"Numerical Composition": str})

    def expected_code(composition):
        comp = ast.literal_eval(composition)
        return "".join(
            str(value)
            for value in (
                comp.get("HexNAc", 0),
                comp.get("Hex", 0),
                comp.get("Fuc", comp.get("dHex", 0)),
                comp.get("NeuAc", comp.get("Neu5Ac", 0)),
                comp.get("NeuGc", comp.get("Neu5Gc", 0)),
            )
        )

    expected = db["Composition"].apply(expected_code)
    assert db["Numerical Composition"].equals(expected)
