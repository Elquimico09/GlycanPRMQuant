import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from .constants import PROTON_MASS, NH4_MASS, DEFAULT_FRAGMENT_DB
import os

_FRAG_DB_CACHE = {}

def _normalize_glycan(val):
    if pd.isna(val):
        return ""
    s = str(val).strip()
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except ValueError:
        pass
    if s.endswith(".0"):
        return s[:-2]
    return s

def _load_fragment_db(db_path: str) -> pd.DataFrame:
    if db_path in _FRAG_DB_CACHE:
        return _FRAG_DB_CACHE[db_path]
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Fragment database file not found: {db_path}")
    db = pd.read_csv(db_path)
    _FRAG_DB_CACHE[db_path] = db
    return db

def cluster_1d(mzs: np.ndarray, tol: float) -> np.ndarray:
    """
    Fast 1D clustering: sort mzs, start a new cluster whenever the gap > tol.
    Returns an array of cluster labels in the original order.
    """
    if mzs.size == 0:
        return np.array([], dtype=int)
    idx = np.argsort(mzs)
    sorted_mzs = mzs[idx]
    labels = np.zeros_like(sorted_mzs, dtype=int)
    cl = 0
    for i in range(1, len(sorted_mzs)):
        if sorted_mzs[i] - sorted_mzs[i-1] > tol:
            cl += 1
        labels[i] = cl
    # invert the sort
    inv = np.empty_like(idx)
    inv[idx] = np.arange(len(idx))
    return labels[inv]


def preprocess_ms2_data(ms2: pd.DataFrame, mz_tol: float = 0.1) -> pd.DataFrame:
    """
    For each (scan_number, precursor_mz) group, cluster fragment_mz via cluster_1d
    and collapse to one row per cluster (mean m/z, summed intensity).
    """
    out = []
    grouped = ms2.groupby(['scan_number', 'precursor_mz'])
    for (scan, prec), grp in grouped:
        mzs = grp['fragment_mz'].to_numpy()
        labels = cluster_1d(mzs, mz_tol)
        for cl in np.unique(labels):
            sub = grp.iloc[labels == cl]
            row = sub.iloc[0].copy()
            row['fragment_mz'] = sub['fragment_mz'].mean()
            row['fragment_intensity'] = sub['fragment_intensity'].sum()
            out.append(row)
    result = pd.DataFrame(out).reset_index(drop=True)
    print(f"Preprocessed MS2 data: {len(ms2)} → {len(result)} fragments")
    return result


def matchMS2(
    ms2_extracted_data: pd.DataFrame,
    precursor_matched_data: pd.DataFrame,
    precursor_composition: str,
    mz_tol: float = 0.02,
    intensity_threshold: float = 1e2,
    ppm_tol: float = 10,
    db_path: str = None
) -> pd.DataFrame:
    """
    Vectorized MS2 matching using cKDTree on theoretical adduct m/z's.
    NH4 adducts are only included if glycan starts with '2'.

    :param db_path: Path to fragment database CSV file. If None, uses default.
    """
    # 1) load & filter fragment database
    if db_path is None:
        db_path = DEFAULT_FRAGMENT_DB

    db = _load_fragment_db(db_path)
    db['Glycan'] = db['Glycan'].astype(str)
    precursor_composition = _normalize_glycan(precursor_composition)
    db['Glycan_norm'] = db['Glycan'].apply(_normalize_glycan)
    dbf = db.loc[db['Glycan_norm'] == precursor_composition].copy()
    if dbf.empty:
        print(f"No database entries for {precursor_composition}")
        return pd.DataFrame()

    # 2) build a single adduct table
    #    columns: Theo_mz, Fragment, Charge
    table = []
    # +H (1+)
    table.append(pd.DataFrame({
        'Theo_mz': dbf['Fragment Mass'] + PROTON_MASS,
        'Fragment': dbf['Fragment'],
        'Charge': 1,
        'Adduct': '+H'
    }))
    # +2H (2+)
    table.append(pd.DataFrame({
        'Theo_mz': (dbf['Fragment Mass'] + 2*PROTON_MASS) / 2,
        'Fragment': dbf['Fragment'],
        'Charge': 2,
        'Adduct': '+2H'
    }))
    # conditional NH4 adducts
    if str(precursor_composition).startswith('2'):
        # +NH4 (1+)
        table.append(pd.DataFrame({
            'Theo_mz': dbf['Fragment Mass'] + NH4_MASS,
            'Fragment': dbf['Fragment'],
            'Charge': 1,
            'Adduct': '+NH4'
        }))
        # +H+NH4 (2+)
        table.append(pd.DataFrame({
            'Theo_mz': (dbf['Fragment Mass'] + PROTON_MASS + NH4_MASS) / 2,
            'Fragment': dbf['Fragment'],
            'Charge': 2,
            'Adduct': '+H+NH4'
        }))

    adduct_df = pd.concat(table, ignore_index=True)
    theo_mzs   = adduct_df['Theo_mz'].to_numpy()
    tree       = cKDTree(theo_mzs.reshape(-1,1))

    # 3) select MS2 rows whose precursor_mz matches an MS1 precursor (track adduct per precursor)
    precursor_matched_data = precursor_matched_data.copy()
    precursor_matched_data['Glycan'] = (
        precursor_matched_data['Glycan']
        .astype(str)
        .apply(_normalize_glycan)
    )
    prec_rows = precursor_matched_data[
        precursor_matched_data['Glycan'] == precursor_composition
    ][['precursor_mz', 'Adduct']]
    if prec_rows.empty:
        print(f"No MS1 precursor matches for {precursor_composition}")
        return pd.DataFrame()

    filt = []
    for _, row in prec_rows.iterrows():
        p = row['precursor_mz']
        adduct_label = row['Adduct']
        tol = p * ppm_tol / 1e6
        sel = ms2_extracted_data[
            ms2_extracted_data['precursor_mz'].between(p-tol, p+tol) &
            (ms2_extracted_data['fragment_intensity'] >= intensity_threshold)
        ]
        if not sel.empty:
            filt.append(sel.assign(Glycan=precursor_composition,
                                   PrecursorAdduct=adduct_label))
    if not filt:
        print(f"No MS2 data after filtering for {precursor_composition}")
        return pd.DataFrame()
    ms2f = pd.concat(filt, ignore_index=True)

    # 4) cluster fragments within each spectrum
    ms2f = preprocess_ms2_data(ms2f, mz_tol=0.1)

    # 5) perform vectorized adduct matching
    obs_mzs = ms2f['fragment_mz'].to_numpy()
    neighbors = tree.query_ball_point(obs_mzs.reshape(-1,1), r=mz_tol)

    matched = []
    base_cols = ['scan_number', 'rt', 'precursor_mz', 'fragment_intensity', 'Glycan', 'PrecursorAdduct']
    base_data = ms2f[base_cols].to_dict('records')

    for i, nbrs in enumerate(neighbors):
        if not nbrs:
            continue
        diffs = np.abs(theo_mzs[nbrs] - obs_mzs[i])
        j = nbrs[np.argmin(diffs)]
        hit = adduct_df.iloc[j]
        diff = diffs.min()
        row = dict(base_data[i])
        row.update({
            'Fragment':    hit['Fragment'],
            'Charge':      int(hit['Charge']),
            'Fragment_mz': float(hit['Theo_mz']),
            'mz_diff':     float(diff),
            'ppm_error':   float(diff / hit['Theo_mz'] * 1e6),
            'Adduct':      hit['Adduct']
        })
        matched.append(row)

    if not matched:
        print(f"No fragments matched for {precursor_composition}")
        return pd.DataFrame()

    out = pd.DataFrame(matched)
    print(f"Matched {len(out)} fragments for {precursor_composition}:")
    print(out['Charge'].value_counts().rename_axis('Charge').to_string())
    return out
