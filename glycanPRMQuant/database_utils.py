"""Validation helpers for glycan structure databases."""

import os
import re

import pandas as pd


REQUIRED_GLYCAN_DATABASE_COLUMNS = (
    "Condensed IUPAC",
    "Composition",
    "Numerical Composition",
)

_LINKED_RESIDUE_PATTERN = re.compile(r"([A-Za-z][A-Za-z0-9]*)(?=\()")
_REDUCING_END_PATTERN = re.compile(r"([A-Za-z][A-Za-z0-9]*)$")
_RESIDUE_SLOT = {
    "GlcNAc": 0,
    "GlcNAc6S": 0,
    "GlcNAcOS": 0,
    "GalNAc": 0,
    "GalNAc4S": 0,
    "GalNAcOS": 0,
    "Man": 1,
    "Man6P": 1,
    "ManOP": 1,
    "ManOS": 1,
    "Gal": 1,
    "Gal3S": 1,
    "GalOS": 1,
    "Glc": 1,
    "Fuc": 2,
    "Neu5Ac": 3,
    "Neu5AcOAc": 3,
    "Neu5Gc": 4,
}


def numerical_composition_from_iupac(iupac_string: str) -> str:
    """Derive HexNAc/Hex/Fuc/Neu5Ac/Neu5Gc counts from condensed IUPAC."""
    text = str(iupac_string).strip()
    if not text:
        raise ValueError("Condensed IUPAC cannot be blank")

    residues = _LINKED_RESIDUE_PATTERN.findall(text)
    reducing_end = _REDUCING_END_PATTERN.search(text)
    if reducing_end is None:
        raise ValueError(f"Could not identify reducing-end residue in IUPAC: {text}")
    residues.append(reducing_end.group(1))

    counts = [0, 0, 0, 0, 0]
    for residue in residues:
        try:
            counts[_RESIDUE_SLOT[residue]] += 1
        except KeyError as exc:
            raise ValueError(
                f"Unsupported residue '{residue}' in condensed IUPAC: {text}"
            ) from exc
    return "".join(str(count) for count in counts)


def _normalize_numerical_composition(value: object) -> str:
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0", text):
        return text[:-2]
    return text


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

    derived_compositions = []
    for row_number, iupac_string in enumerate(database["Condensed IUPAC"], start=2):
        try:
            derived_compositions.append(
                numerical_composition_from_iupac(iupac_string)
            )
        except ValueError as exc:
            raise ValueError(f"Invalid glycan database row {row_number}: {exc}") from exc

    supplied_compositions = database["Numerical Composition"].map(
        _normalize_numerical_composition
    )
    derived_compositions = pd.Series(derived_compositions, index=database.index)
    mismatches = supplied_compositions.ne(derived_compositions)
    if mismatches.any():
        index = mismatches[mismatches].index[0]
        row_number = database.index.get_loc(index) + 2
        raise ValueError(
            f"Glycan database row {row_number} has Numerical Composition "
            f"'{supplied_compositions.loc[index]}', but its Condensed IUPAC "
            f"requires '{derived_compositions.loc[index]}'"
        )

    return database
