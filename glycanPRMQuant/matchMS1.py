import pandas as pd
import os
from .constants import PROTON_MASS, NH4_MASS, DEFAULT_PRECURSOR_DB

def matchMS1(ms1_data, ppm_tol=10, mz_min=400, mz_max=2000, mz_offset=0.0, mass_offset=0.0, db_path=None):
    """
    Match MS1 data with glycan database by computing adduct m/z values
    from the neutral mass column 'M'. Supports +2H, +3H, +4H, +H+NH4+, and +2NH4 adducts.

    :param ms1_data: DataFrame containing MS1 data with 'precursor_mz' column.
    :param ppm_tol: Mass tolerance in parts per million.
    :param mz_min: Minimum m/z value to consider.
    :param mz_max: Maximum m/z value to consider.
    :param mz_offset: Amount (in Da) to add to all calculated adduct m/z values.
    :param mass_offset: Amount (in Da) to add to neutral masses in database.
    :param db_path: Path to glycan database Excel file. If None, uses default.
    :return: DataFrame with columns ['precursor_mz','Glycan'].
    """
    print(f"Starting MS1 matching with ±{ppm_tol} ppm tolerance")
    print("Loading glycan database...")

    if db_path is None:
        db_path = DEFAULT_PRECURSOR_DB

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database file not found: {db_path}")

    db = pd.read_excel(db_path)
    print(f"Loaded {len(db)} glycan entries")

    required_cols = ['Composition', '[M]']
    for col in required_cols:
        if col not in db.columns:
            raise ValueError(f"Missing required column '{col}' in the database")

    if mass_offset != 0.0:
        print(f"Applying mass offset of {mass_offset} Da to all Glycans")
        db['[M]'] = db['[M]'] + mass_offset
        

    # Compute adduct m/z values from neutral mass 'M'
    db = db.copy()
    db['mz_2H']     = (db['[M]'] + 2*PROTON_MASS)    / 2
    db['mz_3H']     = (db['[M]'] + 3*PROTON_MASS)    / 3
    db['mz_4H']     = (db['[M]'] + 4*PROTON_MASS)    / 4
    db['mz_H_NH4']  = (db['[M]'] + PROTON_MASS + NH4_MASS) / 2
    db['mz_2NH4']   = (db['[M]'] + 2*NH4_MASS)        / 2

    adduct_cols = ['mz_2H','mz_3H','mz_4H','mz_H_NH4','mz_2NH4']
    adduct_labels = {
        'mz_2H':    '2H',
        'mz_3H':    '3H',
        'mz_4H':    '4H',
        'mz_H_NH4': 'H+NH4',
        'mz_2NH4':  '2NH4'
    }

    # Apply global offset if requested
    if mz_offset:
        print(f"Applying m/z offset of {mz_offset} Da to all adducts")
        db[adduct_cols] = db[adduct_cols] + mz_offset

    # Filter MS1 data to m/z range
    ms1 = ms1_data[
        (ms1_data['precursor_mz'] >= mz_min) &
        (ms1_data['precursor_mz'] <= mz_max)
    ]
    print(f"Filtered MS1 to {len(ms1)} scans in {mz_min}-{mz_max} m/z")

    unique_prec = ms1['precursor_mz'].unique()
    print(f"{len(unique_prec)} unique precursor m/z values to match")

    matches = {}
    total_matches = 0

    for i, prec in enumerate(unique_prec, 1):
        if i % 100 == 0:
            print(f"  Processed {i}/{len(unique_prec)} precursors, {total_matches} matches so far")
        tol = prec * ppm_tol / 1e6

        for col in adduct_cols:
            # vectorized filter
            hits = db[
                (db[col] >= prec - tol) &
                (db[col] <= prec + tol)
            ]
            for _, row in hits.iterrows():
                matches.setdefault(prec, []).append({
                    'matched_glycan': str(row['Composition']),
                    'matched_adduct': adduct_labels[col],
                    'database_mz': row[col],
                    'ppm_error': (prec - row[col]) / row[col] * 1e6
                })
                total_matches += 1

    print(f"Matching complete: {len(matches)} precursors with matches, {total_matches} total matches")

    # Build result DataFrame
    results = []
    for prec, lst in matches.items():
        glycans = '; '.join(m['matched_glycan'] for m in lst)
        results.append({'precursor_mz': prec, 'Glycan': glycans})

    df_out = pd.DataFrame(results)
    print(f"Returning {len(df_out)} matched MS1 entries")
    return df_out
