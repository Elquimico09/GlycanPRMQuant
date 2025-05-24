import os
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN

def preprocess_ms2_data(ms2_extracted_data: pd.DataFrame, mz_tol: float = 0.1) -> pd.DataFrame:
    """
    Cluster fragment m/z values within mz_tol for each spectrum,
    average the m/z in each cluster, and sum the intensities.
    """
    df = ms2_extracted_data.copy()
    result = []

    # process spectrum by spectrum
    for (scan, precursor), grp in df.groupby(['scan_number', 'precursor_mz']):
        mzs = grp['fragment_mz'].values.reshape(-1, 1)
        labels = DBSCAN(eps=mz_tol, min_samples=1).fit_predict(mzs)
        grp = grp.assign(cluster=labels)

        for cid, sub in grp.groupby('cluster'):
            row = sub.iloc[0].copy()
            row['fragment_mz'] = sub['fragment_mz'].mean()
            row['fragment_intensity'] = sub['fragment_intensity'].sum()
            result.append(row)

    out = pd.DataFrame(result).reset_index(drop=True)
    print(f"Preprocessed MS2 data: reduced from {len(df)} to {len(out)} fragment ions")
    return out

def matchMS2(
    ms2_extracted_data: pd.DataFrame,
    precursor_matched_data: pd.DataFrame,
    precursor_composition: str,
    mz_tol: float = 0.02,
    intensity_threshold: float = 1e2,
    ppm_tol: float = 10
) -> pd.DataFrame:
    """
    Match MS2 fragments to a glycan's theoretical fragments (M+H, [M+2H]2+,
    M+NH4, and [M+H+NH4]2+).  The NH4 adducts are only considered if
    the glycan composition string starts with '2'.
    """

    PROTON_MASS   = 1.007276
    AMMONIUM_MASS = 18.033826

    print(f"Starting MS2 matching with ±{mz_tol} Da tolerance")
    # 1) load and filter MS2 fragment database
    db = pd.read_csv("database/fragment_database.csv")
    db['Glycan'] = db['Glycan'].astype(str)
    dbf = db[db['Glycan'] == precursor_composition].copy()
    print(f"Loaded {len(dbf)} entries for Glycan {precursor_composition}")

    # 2) compute adduct m/z columns
    dbf['mz_H']       = dbf['Fragment Mass'] + PROTON_MASS
    dbf['mz_2H']      = (dbf['Fragment Mass'] + 2*PROTON_MASS) / 2
    dbf['mz_NH4']     = dbf['Fragment Mass'] + AMMONIUM_MASS
    dbf['mz_H_NH4']   = (dbf['Fragment Mass'] + PROTON_MASS + AMMONIUM_MASS) / 2

    # 3) filter MS2 data to only those precursors matched in MS1
    glycan_precs = precursor_matched_data[
        precursor_matched_data['Glycan']
        .str.split(';')
        .apply(lambda comps: precursor_composition in [c.strip() for c in comps])
    ]['precursor_mz'].unique()
    if len(glycan_precs) == 0:
        print(f"No MS1 precursor matches for {precursor_composition}")
        return pd.DataFrame()

    filt = []
    for prec in glycan_precs:
        ppm_window = prec * ppm_tol / 1e6
        sel = ms2_extracted_data[
            (ms2_extracted_data['precursor_mz'].between(prec-ppm_window, prec+ppm_window)) &
            (ms2_extracted_data['fragment_intensity'] >= intensity_threshold)
        ]
        if not sel.empty:
            filt.append(sel.assign(Glycan=precursor_composition))
    if not filt:
        print(f"No MS2 spectra for {precursor_composition} after filtering")
        return pd.DataFrame()
    ms2f = pd.concat(filt, ignore_index=True)
    print(f"Found {len(ms2f)} MS2 points matching precursors")

    # 4) cluster within each spectrum
    ms2f = preprocess_ms2_data(ms2f, mz_tol=0.1)

    # 5) build adduct list
    adducts = [
        ('mz_H',     1),
        ('mz_2H',    2)
    ]
    if str(precursor_composition).startswith('2'):
        adducts += [
            ('mz_NH4',   1),
            ('mz_H_NH4', 2)
        ]

    # 6) match each observed fragment to closest adduct entry
    matched = []
    for _, frag in ms2f.iterrows():
        mz_obs = frag['fragment_mz']
        best = None
        best_diff = mz_tol + 1e-6
        for col, charge in adducts:
            hits = dbf[
                (dbf[col] >= mz_obs - mz_tol) &
                (dbf[col] <= mz_obs + mz_tol)
            ]
            for _, hit in hits.iterrows():
                diff = abs(hit[col] - mz_obs)
                if diff < best_diff:
                    best_diff = diff
                    best = {
                        **frag.to_dict(),
                        'Fragment':     hit['Fragment'],   # preserve DB's Fragment name
                        'Charge':       charge,
                        'Fragment_mz':  hit[col],
                        'mz_diff':      diff,
                        'ppm_error':    diff / hit[col] * 1e6
                    }
        if best:
            matched.append(best)

    if not matched:
        print(f"No fragment matches for glycan {precursor_composition}")
        return pd.DataFrame()

    matched_df = pd.DataFrame(matched)
    print(f"Successfully matched {len(matched_df)} fragments:")
    for ch, cnt in matched_df['Charge'].value_counts().items():
        print(f"  - {cnt} hits at charge {ch}+")

    return matched_df
