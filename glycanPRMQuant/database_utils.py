"""Validation helpers for glycan structure databases."""

import os

import pandas as pd


REQUIRED_GLYCAN_DATABASE_COLUMNS = (
    "Condensed IUPAC",
    "Composition",
    "Numerical Composition",
)


def validate_glycan_database(db_path: str) -> pd.DataFrame:
    """Load and validate a database used for precursor and fragment matching."""
    if not db_path:
        raise ValueError("A glycan database is required")
    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Glycan database does not exist: {db_path}")

    ext = os.path.splitext(db_path)[1].lower()
    if ext == ".csv":
        database = pd.read_csv(db_path, dtype={"Numerical Composition": str})
    elif ext in (".xlsx", ".xls"):
        database = pd.read_excel(db_path, dtype={"Numerical Composition": str})
    else:
        raise ValueError(
            "Unsupported glycan database file type. Use CSV, XLSX, or XLS."
        )

    missing_columns = [
        column
        for column in REQUIRED_GLYCAN_DATABASE_COLUMNS
        if column not in database.columns
    ]
    if missing_columns:
        raise ValueError(
            "Glycan database is missing required column(s): "
            + ", ".join(missing_columns)
        )
    if database.empty:
        raise ValueError("Glycan database contains no rows")

    for column in REQUIRED_GLYCAN_DATABASE_COLUMNS:
        values = database[column]
        invalid = values.isna() | values.astype(str).str.strip().eq("")
        if invalid.any():
            raise ValueError(
                f"Glycan database column '{column}' contains "
                f"{int(invalid.sum())} missing or blank value(s)"
            )

    return database
