import pandas as pd

def matchMS2(ms2_extracted_data, precursor_matched_data, precursor_composition,
             mz_tol=0.02, intensity_threshold=1e2, ppm_tol=10):
    """
    Match MS2 data with glycan database, including:
      - M+H+ (charge 1)
      - [M+2H]2+ (charge 2)
      - M+NH4+ and [M+H+NH4]2+ only if glycan composition starts with '2'
    """
    PROTON_MASS = 1.007276
    AMMONIUM_MASS = 18.033825

    print(f"Starting MS2 matching with m/z tolerance: {mz_tol} Da")
    print("Loading MS2 database...")
    db = pd.read_csv("database/fragment_database.csv")
    db['Glycan'] = db['Glycan'].astype(str)

    # filter for the target glycan composition
    dbf = db[db['Glycan'] == precursor_composition].copy()
    print(f"Loaded {len(dbf)} fragments for glycan {precursor_composition}")

    # compute adduct m/zs
    dbf['mz_H']     = dbf['Fragment Mass'] + PROTON_MASS
    dbf['mz_2H']    = (dbf['Fragment Mass'] + 2 * PROTON_MASS) / 2
    dbf['mz_NH4']   = dbf['Fragment Mass'] + AMMONIUM_MASS
    dbf['mz_H_NH4'] = (dbf['Fragment Mass'] + PROTON_MASS + AMMONIUM_MASS) / 2

    # get matched MS1 precursors for this glycan
    glycan_matches = precursor_matched_data[
        precursor_matched_data['Glycan']
        .str.split(';')
        .apply(lambda comps: precursor_composition in [c.strip() for c in comps])
    ]
    if glycan_matches.empty:
        print("No MS1 matches for", precursor_composition)
        return pd.DataFrame()
    precs = glycan_matches['precursor_mz'].values

    # filter ms2_extracted_data by precursor m/z ± ppm and intensity
    filtered = []
    for prec in precs:
        tol = prec * ppm_tol / 1e6
        m = ms2_extracted_data[
            (ms2_extracted_data['precursor_mz'].between(prec - tol, prec + tol)) &
            (ms2_extracted_data['fragment_intensity'] >= intensity_threshold)
        ]
        if not m.empty:
            filtered.append(m.assign(Glycan=precursor_composition))
    if not filtered:
        print("No MS2 data found for precursors.")
        return pd.DataFrame()
    ms2f = pd.concat(filtered, ignore_index=True)
    print(f"Found {len(ms2f)} MS2 points matching precursor")

    # build adduct list; always include H and 2H
    adducts = [
        ('mz_H', 1),
        ('mz_2H', 2)
    ]
    # only include ammonium adducts if composition starts with '2'
    if str(precursor_composition).startswith('2'):
        adducts += [
            ('mz_NH4', 1),
            ('mz_H_NH4', 2)
        ]

    # match each fragment to the closest adduct
    matches = []
    for _, row in ms2f.iterrows():
        mz_obs = row['fragment_mz']
        best = None
        best_diff = mz_tol + 1
        for col, charge in adducts:
            hits = dbf[
                dbf[col].between(mz_obs - mz_tol, mz_obs + mz_tol)
            ]
            for _, hit in hits.iterrows():
                diff = abs(hit[col] - mz_obs)
                if diff < best_diff:
                    best_diff = diff
                    best = {
                        **row.to_dict(),
                        'Charge': charge,
                        'Fragment_mz': hit[col],
                        'mz_diff': diff,
                        'ppm_error': diff / hit[col] * 1e6
                    }
        if best:
            matches.append(best)

    if not matches:
        print("No fragment matches in database")
        return pd.DataFrame()

    matched_df = pd.DataFrame(matches)
    print(f"Matched {len(matched_df)} fragments.")
    for ch, cnt in matched_df['Charge'].value_counts().items():
        print(f" - charge {ch}: {cnt} hits")
    return matched_df
