import logging

import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from .constants import DEFAULT_PRECURSOR_DB
from .fragment_structure import fragment_glycan
import os

logger = logging.getLogger(__name__)

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

def _read_structure_db(db_path: str) -> pd.DataFrame:
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"N-glycan database file not found: {db_path}")

    ext = os.path.splitext(db_path)[1].lower()
    if ext == ".csv":
        return pd.read_csv(db_path)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(db_path)
    raise ValueError(f"Unsupported N-glycan database file type: {ext}")


def _build_fragment_db(
    db_path: str,
    precursor_composition: str,
    ion_series: str = "ABCXYZ",
    max_cleavages: int = 2
) -> pd.DataFrame:
    ion_series = _normalize_ion_series(ion_series)
    max_cleavages = _normalize_max_cleavages(max_cleavages)
    key = (db_path, precursor_composition, ion_series, max_cleavages)
    if key in _FRAG_DB_CACHE:
        return _FRAG_DB_CACHE[key].copy()

    db = _read_structure_db(db_path)
    required_cols = ['Numerical Composition', 'Composition', 'Condensed IUPAC']
    for col in required_cols:
        if col not in db.columns:
            raise ValueError(f"Missing required column '{col}' in the N-glycan database")

    db = db.copy()
    db['NumericalComposition_norm'] = db['Numerical Composition'].apply(_normalize_glycan)
    candidates = db.loc[db['NumericalComposition_norm'] == precursor_composition].copy()
    if candidates.empty:
        _FRAG_DB_CACHE[key] = pd.DataFrame()
        return pd.DataFrame()

    rows = []
    skipped = 0
    for _, candidate in candidates.iterrows():
        iupac_str = candidate['Condensed IUPAC']
        try:
            fragments = fragment_glycan(
                iupac_str,
                ion_series=ion_series,
                max_cleavages=max_cleavages
            )
        except Exception:
            skipped += 1
            continue

        if fragments.empty:
            continue

        fragments = fragments.copy()
        fragments['Fragment'] = fragments['name']
        fragments['FragmentType'] = fragments['series']
        fragments['IUPAC'] = iupac_str
        fragments['NumericalComposition'] = precursor_composition
        fragments['Composition'] = candidate['Composition']
        if 'glytoucan_id' in candidate.index:
            fragments['GlyTouCanID'] = candidate['glytoucan_id']
        rows.append(fragments)

    if not rows:
        frag_db = pd.DataFrame()
    else:
        frag_db = pd.concat(rows, ignore_index=True)

    if skipped:
        logger.warning(f"Skipped {skipped} IUPAC structures for {precursor_composition} because fragmentation failed")
    logger.info(
        f"Generated {len(frag_db)} theoretical fragments for "
        f"{candidates.shape[0] - skipped}/{candidates.shape[0]} IUPAC structures of {precursor_composition}"
    )
    _FRAG_DB_CACHE[key] = frag_db
    return frag_db.copy()


def _build_adduct_table(dbf: pd.DataFrame, precursor_composition: str) -> pd.DataFrame:
    table = [
        pd.DataFrame({
            'Theo_mz': dbf['[M+H]+'],
            'Fragment': dbf['Fragment'],
            'FragmentType': dbf['FragmentType'],
            'Charge': 1,
            'Adduct': '+H',
            'IUPAC': dbf['IUPAC'],
            'NumericalComposition': dbf['NumericalComposition'],
            'Composition': dbf['Composition'],
        }),
        pd.DataFrame({
            'Theo_mz': dbf['[M+2H]2+'],
            'Fragment': dbf['Fragment'],
            'FragmentType': dbf['FragmentType'],
            'Charge': 2,
            'Adduct': '+2H',
            'IUPAC': dbf['IUPAC'],
            'NumericalComposition': dbf['NumericalComposition'],
            'Composition': dbf['Composition'],
        }),
    ]

    if str(precursor_composition).startswith('2'):
        table.extend([
            pd.DataFrame({
                'Theo_mz': dbf['[M+NH4]+'],
                'Fragment': dbf['Fragment'],
                'FragmentType': dbf['FragmentType'],
                'Charge': 1,
                'Adduct': '+NH4',
                'IUPAC': dbf['IUPAC'],
                'NumericalComposition': dbf['NumericalComposition'],
                'Composition': dbf['Composition'],
            }),
            pd.DataFrame({
                'Theo_mz': dbf['[M+NH4+H]2+'],
                'Fragment': dbf['Fragment'],
                'FragmentType': dbf['FragmentType'],
                'Charge': 2,
                'Adduct': '+H+NH4',
                'IUPAC': dbf['IUPAC'],
                'NumericalComposition': dbf['NumericalComposition'],
                'Composition': dbf['Composition'],
            }),
        ])

    adduct_df = pd.concat(table, ignore_index=True)
    adduct_df = adduct_df.dropna(subset=['Theo_mz'])
    return adduct_df


def _normalize_ion_series(ion_series: str) -> str:
    ion_series = (ion_series or "ABCXYZ").upper().replace(",", "").replace(" ", "")
    allowed = set("ABCXYZ")
    invalid = sorted(set(ion_series) - allowed)
    if invalid:
        raise ValueError(
            f"Unsupported fragment ion series character(s): {''.join(invalid)}. "
            "Use any combination of A, B, C, X, Y, Z."
        )
    if not ion_series:
        raise ValueError("Fragment ion series cannot be empty")
    return ion_series


def _normalize_max_cleavages(max_cleavages: int) -> int:
    try:
        value = int(max_cleavages)
    except (TypeError, ValueError):
        raise ValueError("Max cleavages must be an integer")
    if value < 1:
        raise ValueError("Max cleavages must be >= 1")
    return value

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


def preprocess_ms2_data(
    ms2: pd.DataFrame,
    fragment_mass_tol: float = 0.1
) -> pd.DataFrame:
    """
    For each (scan_number, precursor_mz) group, cluster fragment_mz via cluster_1d
    and collapse to one row per cluster (mean m/z, summed intensity).
    """
    out = []
    grouped = ms2.groupby(['scan_number', 'precursor_mz'])
    for (scan, prec), grp in grouped:
        mzs = grp['fragment_mz'].to_numpy()
        labels = cluster_1d(mzs, fragment_mass_tol)
        for cl in np.unique(labels):
            sub = grp.iloc[labels == cl]
            row = sub.iloc[0].copy()
            row['fragment_mz'] = sub['fragment_mz'].mean()
            row['fragment_intensity'] = sub['fragment_intensity'].sum()
            out.append(row)
    result = pd.DataFrame(out).reset_index(drop=True)
    logger.info(f"Preprocessed MS2 data: {len(ms2)} -> {len(result)} fragments")
    return result


def matchMS2(
    ms2_extracted_data: pd.DataFrame,
    precursor_matched_data: pd.DataFrame,
    precursor_composition: str,
    fragment_mass_tol: float = 0.02,
    intensity_threshold: float = 1e2,
    ppm_tol: float = 10,
    db_path: str = None,
    ion_series: str = "ABCXYZ",
    max_cleavages: int = 2
) -> pd.DataFrame:
    """
    Vectorized MS2 matching using cKDTree on theoretical fragments generated
    from the N-glycan IUPAC database.

    For each numerical composition, all parsable IUPAC candidate structures are
    fragmented, observed MS2 fragments are matched to those theoretical
    fragments, and the best-scoring IUPAC structure is returned.
    NH4 adducts are only included if glycan starts with '2'.

    :param fragment_mass_tol: Fragment m/z matching tolerance in Da.
    :param ppm_tol: Precursor matching tolerance in ppm.
    :param db_path: Path to N-glycan database CSV/Excel file. If None, uses default.
    """
    if db_path is None:
        db_path = DEFAULT_PRECURSOR_DB

    ion_series = _normalize_ion_series(ion_series)
    max_cleavages = _normalize_max_cleavages(max_cleavages)
    precursor_composition = _normalize_glycan(precursor_composition)
    dbf = _build_fragment_db(
        db_path,
        precursor_composition,
        ion_series=ion_series,
        max_cleavages=max_cleavages
    )
    if dbf.empty:
        logger.info(f"No fragmentable IUPAC database entries for {precursor_composition}")
        return pd.DataFrame()

    adduct_df = _build_adduct_table(dbf, precursor_composition)
    if adduct_df.empty:
        logger.info(f"No theoretical fragment adducts for {precursor_composition}")
        return pd.DataFrame()

    theo_mzs = adduct_df['Theo_mz'].to_numpy()
    tree = cKDTree(theo_mzs.reshape(-1, 1))

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
        logger.info(f"No MS1 precursor matches for {precursor_composition}")
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
        logger.info(f"No MS2 data after filtering for {precursor_composition}")
        return pd.DataFrame()
    ms2f = pd.concat(filt, ignore_index=True)

    # 4) cluster fragments within each spectrum
    ms2f = preprocess_ms2_data(ms2f, fragment_mass_tol=0.1)

    # 5) perform vectorized adduct matching against all candidate IUPAC structures
    obs_mzs = ms2f['fragment_mz'].to_numpy()
    neighbors = tree.query_ball_point(obs_mzs.reshape(-1, 1), r=fragment_mass_tol)

    matched = []
    base_cols = ['scan_number', 'rt', 'precursor_mz', 'fragment_intensity', 'Glycan', 'PrecursorAdduct']
    base_data = ms2f[base_cols].to_dict('records')

    for i, nbrs in enumerate(neighbors):
        if not nbrs:
            continue
        hits = adduct_df.iloc[nbrs].copy()
        hits['mz_diff'] = np.abs(hits['Theo_mz'].to_numpy() - obs_mzs[i])
        hits = hits.sort_values('mz_diff').drop_duplicates(subset=['IUPAC'], keep='first')

        for _, hit in hits.iterrows():
            diff = float(hit['mz_diff'])
            row = dict(base_data[i])
            row.update({
                'Fragment': hit['Fragment'],
                'FragmentType': hit['FragmentType'],
                'Charge': int(hit['Charge']),
                'Fragment_mz': float(hit['Theo_mz']),
                'mz_diff': diff,
                'ppm_error': float(diff / hit['Theo_mz'] * 1e6),
                'Adduct': hit['Adduct'],
                'IUPAC': hit['IUPAC'],
                'NumericalComposition': hit['NumericalComposition'],
                'Composition': hit['Composition'],
            })
            matched.append(row)

    if not matched:
        logger.info(f"No fragments matched for {precursor_composition}")
        return pd.DataFrame()

    all_matches = pd.DataFrame(matched)
    scores = (
        all_matches.groupby('IUPAC')
        .agg(
            match_count=('Fragment', 'size'),
            unique_fragments=('Fragment', 'nunique'),
            total_intensity=('fragment_intensity', 'sum'),
            mean_abs_ppm=('ppm_error', lambda s: float(np.mean(np.abs(s)))),
        )
        .reset_index()
        .sort_values(
            ['match_count', 'unique_fragments', 'total_intensity', 'mean_abs_ppm'],
            ascending=[False, False, False, True]
        )
    )
    best = scores.iloc[0]
    out = all_matches.loc[all_matches['IUPAC'] == best['IUPAC']].copy()
    out['IUPAC_match_count'] = int(best['match_count'])
    out['IUPAC_unique_fragments'] = int(best['unique_fragments'])
    out['IUPAC_total_intensity'] = float(best['total_intensity'])

    logger.info(f"Selected IUPAC for {precursor_composition}: {best['IUPAC']}")
    logger.info(
        f"Structure score: {int(best['match_count'])} matches, "
        f"{int(best['unique_fragments'])} unique fragments"
    )
    logger.info(f"Matched {len(out)} fragments for {precursor_composition}:")
    logger.info("\n%s", out['Charge'].value_counts().rename_axis('Charge').to_string())
    return out
