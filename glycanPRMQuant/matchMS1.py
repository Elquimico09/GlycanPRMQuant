import logging

import pandas as pd
import os
from .constants import PROTON_MASS, NH4_MASS, DEFAULT_PRECURSOR_DB
from .calculate_mass import calculate_mass

logger = logging.getLogger(__name__)

_DB_CACHE = {}

def _normalize_glycan(val):
    if pd.isna(val):
        return ""
    s = str(val).strip()
    # Normalize numeric-like strings (e.g., "25000.0" -> "25000")
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except ValueError:
        pass
    if s.endswith(".0"):
        return s[:-2]
    return s

def _load_precursor_db(db_path: str, mass_offset: float = 0.0, mz_offset: float = 0.0) -> pd.DataFrame:
    key = (db_path, mass_offset, mz_offset)
    if key in _DB_CACHE:
        return _DB_CACHE[key].copy()
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database file not found: {db_path}")

    ext = os.path.splitext(db_path)[1].lower()
    if ext == ".csv":
        raw_db = pd.read_csv(db_path)
    elif ext in (".xlsx", ".xls"):
        raw_db = pd.read_excel(db_path)
    else:
        raise ValueError(f"Unsupported precursor database file type: {ext}")

    if '[M]' in raw_db.columns:
        db = raw_db.copy()
    else:
        db = _calculate_precursor_masses(raw_db)

    required_cols = ['Composition', '[M]']
    for col in required_cols:
        if col not in db.columns:
            raise ValueError(f"Missing required column '{col}' in the database")
    if mass_offset != 0.0:
        logger.info(f"Applying mass offset of {mass_offset} Da to all Glycans")
        db = db.copy()
        db['[M]'] = db['[M]'] + mass_offset
    db = db.copy()
    db['mz_2H']     = (db['[M]'] + 2*PROTON_MASS)    / 2
    db['mz_3H']     = (db['[M]'] + 3*PROTON_MASS)    / 3
    db['mz_4H']     = (db['[M]'] + 4*PROTON_MASS)    / 4
    db['mz_H_NH4']  = (db['[M]'] + PROTON_MASS + NH4_MASS) / 2
    db['mz_2NH4']   = (db['[M]'] + 2*NH4_MASS)        / 2
    if mz_offset:
        logger.info(f"Applying m/z offset of {mz_offset} Da to all adducts")
        db[['mz_2H','mz_3H','mz_4H','mz_H_NH4','mz_2NH4']] = db[['mz_2H','mz_3H','mz_4H','mz_H_NH4','mz_2NH4']] + mz_offset
    _DB_CACHE[key] = db
    return db.copy()


def _calculate_precursor_masses(raw_db: pd.DataFrame) -> pd.DataFrame:
    """
    Build a one-row-per-composition precursor database from an IUPAC database.

    The N-glycan database may contain multiple isomeric IUPAC structures for the
    same composition. Isomers have the same neutral mass, so calculate the mass
    only once for each unique Composition value and keep Numerical Composition
    as the downstream glycan identifier when available.
    """
    required_cols = ['Composition', 'Condensed IUPAC']
    for col in required_cols:
        if col not in raw_db.columns:
            raise ValueError(
                f"Missing required column '{col}' in precursor database. "
                "Expected either a precomputed '[M]' column or an IUPAC database."
            )

    id_col = 'Numerical Composition' if 'Numerical Composition' in raw_db.columns else 'Composition'
    rows = []
    failed = []

    for composition, group in raw_db.groupby('Composition', sort=False, dropna=False):
        mass = None
        source_iupac = None
        for iupac_str in group['Condensed IUPAC'].dropna():
            mass = calculate_mass(iupac_str, verbose=False)
            if mass is not None:
                source_iupac = iupac_str
                break

        if mass is None:
            failed.append(composition)
            continue

        rows.append({
            'Composition': _normalize_glycan(group[id_col].iloc[0]),
            'Composition Formula': composition,
            'Condensed IUPAC': source_iupac,
            '[M]': mass,
        })

    if not rows:
        raise ValueError("No precursor masses could be calculated from the IUPAC database")

    if failed:
        logger.warning(f"Skipped {len(failed)} compositions because mass calculation failed")

    db = pd.DataFrame(rows)
    logger.info(
        f"Calculated neutral masses for {len(db)} unique compositions "
        f"from {len(raw_db)} IUPAC database rows"
    )
    return db


def matchMS1(ms1_data, ppm_tol=10, mz_min=400, mz_max=2000, mz_offset=0.0, mass_offset=0.0, db_path=None):
    """
    Match MS1 data with glycan database by computing adduct m/z values
    from calculated or precomputed neutral masses. Supports +2H, +3H, +4H,
    +H+NH4+, and +2NH4 adducts.

    :param ms1_data: DataFrame containing MS1 data with 'precursor_mz' column.
    :param ppm_tol: Mass tolerance in parts per million.
    :param mz_min: Minimum m/z value to consider.
    :param mz_max: Maximum m/z value to consider.
    :param mz_offset: Amount (in Da) to add to all calculated adduct m/z values.
    :param mass_offset: Amount (in Da) to add to neutral masses in database.
    :param db_path: Path to glycan database CSV/Excel file. If None, uses default.
    :return: DataFrame with columns ['precursor_mz','Glycan'].
    """
    logger.info(f"Starting MS1 matching with ±{ppm_tol} ppm tolerance")
    logger.info("Loading glycan database...")

    if db_path is None:
        db_path = DEFAULT_PRECURSOR_DB
    db = _load_precursor_db(db_path, mass_offset=mass_offset, mz_offset=mz_offset)
    logger.info(f"Loaded {len(db)} glycan entries")

    adduct_cols = ['mz_2H','mz_3H','mz_4H','mz_H_NH4','mz_2NH4']
    adduct_labels = {
        'mz_2H':    '2H',
        'mz_3H':    '3H',
        'mz_4H':    '4H',
        'mz_H_NH4': 'H+NH4',
        'mz_2NH4':  '2NH4'
    }

    # Filter MS1 data to m/z range
    ms1 = ms1_data[
        (ms1_data['precursor_mz'] >= mz_min) &
        (ms1_data['precursor_mz'] <= mz_max)
    ]
    logger.info(f"Filtered MS1 to {len(ms1)} scans in {mz_min}-{mz_max} m/z")

    unique_prec = ms1['precursor_mz'].unique()
    logger.info(f"{len(unique_prec)} unique precursor m/z values to match")

    matches = []
    total_matches = 0

    for i, prec in enumerate(unique_prec, 1):
        if i % 100 == 0:
            total_matches = len(matches)
            logger.info(f"Processed {i}/{len(unique_prec)} precursors, {total_matches} matches so far")
        tol = prec * ppm_tol / 1e6

        for col in adduct_cols:
            # vectorized filter
            hits = db[
                (db[col] >= prec - tol) &
                (db[col] <= prec + tol)
            ]
            for _, row in hits.iterrows():
                matches.append({
                    'precursor_mz': prec,
                    'Glycan': _normalize_glycan(row['Composition']),
                    'Adduct': adduct_labels[col],
                    'database_mz': row[col],
                    'ppm_error': (prec - row[col]) / row[col] * 1e6
                })
    out_cols = ['precursor_mz', 'Glycan', 'Adduct', 'database_mz', 'ppm_error']
    df_out = pd.DataFrame(matches, columns=out_cols)
    total_matches = len(df_out)
    matched_precursors = df_out['precursor_mz'].nunique() if not df_out.empty else 0
    logger.info(f"Matching complete: {matched_precursors} precursors with matches, {total_matches} total matches")
    return df_out
