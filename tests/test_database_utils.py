import pandas as pd
import pytest

from glycanPRMQuant.constants import DEFAULT_PRECURSOR_DB
from glycanPRMQuant.database_utils import validate_glycan_database


def test_bundled_glycan_database_is_valid():
    database = validate_glycan_database(DEFAULT_PRECURSOR_DB)

    assert not database.empty


def test_glycan_database_requires_bundled_schema(tmp_path):
    db_path = tmp_path / "database.csv"
    pd.DataFrame(
        {
            "Composition": ["{'Hex': 5, 'HexNAc': 2}"],
            "Numerical Composition": ["25000"],
        }
    ).to_csv(db_path, index=False)

    with pytest.raises(ValueError, match="Condensed IUPAC"):
        validate_glycan_database(str(db_path))


def test_glycan_database_rejects_blank_condensed_iupac(tmp_path):
    db_path = tmp_path / "database.csv"
    pd.DataFrame(
        {
            "Condensed IUPAC": [""],
            "Composition": ["{'Hex': 5, 'HexNAc': 2}"],
            "Numerical Composition": ["25000"],
        }
    ).to_csv(db_path, index=False)

    with pytest.raises(ValueError, match="Condensed IUPAC.*missing or blank"):
        validate_glycan_database(str(db_path))
