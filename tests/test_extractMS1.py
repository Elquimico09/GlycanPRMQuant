from pathlib import Path

import pandas as pd

from glycanPRMQuant.constants import DEFAULT_PRECURSOR_DB


def test_default_precursor_database_is_packaged():
    db_path = Path(DEFAULT_PRECURSOR_DB)

    assert db_path.exists()
    db = pd.read_csv(db_path, nrows=1)
    assert {"Condensed IUPAC", "Composition", "Numerical Composition"}.issubset(db.columns)
