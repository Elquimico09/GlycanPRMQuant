import pandas as pd

def matchMS1(ms1_data, ppm_tol=10, mz_min=400, mz_max=2000, mz_offset=0.02):
    """
    Match MS1 data with glycan database, optionally offsetting all database m/z values.

    :param ms1_data: DataFrame containing MS1 data with 'precursor_mz' column.
    :param ppm_tol: Mass tolerance in parts per million.
    :param mz_min: Minimum m/z value to consider (default: 400).
    :param mz_max: Maximum m/z value to consider (default: 2000).
    :param mz_offset: Amount (in Da) to add to all database m/z values before matching.
    :return: DataFrame with matched MS1 data and corresponding glycan information.
    """
    print(f"Starting MS1 matching with ppm tolerance: {ppm_tol}")
    print("Loading glycan database...")
    glycan_database = pd.read_excel("database/glycan_precursor_mz_list.xlsx")
    print(f"Loaded glycan database with {len(glycan_database)} entries")

    # Identify adduct columns
    adduct_columns = [
        col for col in glycan_database.columns
        if col.startswith('[M+') and col.endswith('+')
    ]
    print(f"Found {len(adduct_columns)} adduct columns: {adduct_columns}")

    # Apply mz_offset if specified
    if mz_offset != 0:
        print(f"Applying m/z offset of {mz_offset} Da to all database values")
        glycan_database[adduct_columns] = glycan_database[adduct_columns] + mz_offset

    # Filter MS1 data to only include precursors within the specified m/z range
    print(f"Filtering MS1 data with m/z range {mz_min}-{mz_max}...")
    filtered_ms1 = ms1_data[
        (ms1_data['precursor_mz'] >= mz_min) &
        (ms1_data['precursor_mz'] <= mz_max)
    ]
    print(f"Filtered MS1 data: {len(filtered_ms1)} entries (from {len(ms1_data)} total)")

    # Get unique precursor m/z values
    unique_precursors = filtered_ms1['precursor_mz'].unique()
    print(f"Found {len(unique_precursors)} unique precursor m/z values")

    precursor_matches = {}
    match_count = 0

    # Matching loop
    print("Starting matching process...")
    for i, prec in enumerate(unique_precursors):
        if i and i % 100 == 0:
            print(f"  Processed {i}/{len(unique_precursors)} precursors, {match_count} matches so far")
        tol = prec * ppm_tol / 1e6

        for adduct_col in adduct_columns:
            # find rows where (database_mz ± tol) contains prec
            hits = glycan_database[
                (glycan_database[adduct_col] >= prec - tol) &
                (glycan_database[adduct_col] <= prec + tol)
            ]
            for _, hit in hits.iterrows():
                precursor_matches.setdefault(prec, []).append({
                    'matched_glycan': str(hit.get('Composition', '')),
                    'matched_adduct': adduct_col,
                    'database_mz': hit[adduct_col],
                    'ppm_error': ((prec - hit[adduct_col]) / hit[adduct_col]) * 1e6
                })
                match_count += 1

    print(f"Matching completed. {len(precursor_matches)} precursors had matches, total {match_count} matches.")

    # Build final DataFrame
    unique_results = []
    for prec, matches in precursor_matches.items():
        comps = '; '.join(m['matched_glycan'] for m in matches)
        unique_results.append({'precursor_mz': prec, 'Glycan': comps})

    result_df = pd.DataFrame(unique_results)
    print(f"Final matched results: {len(result_df)} entries")
    return result_df
