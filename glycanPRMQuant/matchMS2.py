import hashlib
import logging

import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from .constants import DEFAULT_PRECURSOR_DB, METHANOL_MASS
from .fragment_structure import fragment_glycan
from .mass_tolerance import (
    conservative_query_radius_da,
    mass_error_ppm,
    normalize_tolerance_unit,
    validate_tolerance,
    within_tolerance,
)
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


def _build_adduct_table(
    dbf: pd.DataFrame,
    precursor_adducts,
) -> pd.DataFrame:
    """Build product-ion forms allowed by each observed precursor adduct."""
    if isinstance(precursor_adducts, str):
        precursor_adducts = [precursor_adducts]
    precursor_adducts = sorted(
        {str(adduct) for adduct in precursor_adducts if pd.notna(adduct)}
    )
    if not precursor_adducts:
        precursor_adducts = [""]

    fragment_iupac = (
        dbf['fragment_iupac']
        if 'fragment_iupac' in dbf
        else pd.Series("", index=dbf.index)
    )
    contains_neuac = (
        dbf['contains_neuac'].fillna(False).astype(bool)
        if 'contains_neuac' in dbf
        else pd.Series(False, index=dbf.index)
    )
    common = {
        'Fragment': dbf['Fragment'],
        'FragmentType': dbf['FragmentType'],
        'IUPAC': dbf['IUPAC'],
        'NumericalComposition': dbf['NumericalComposition'],
        'Composition': dbf['Composition'],
        'fragment_iupac': fragment_iupac,
        'contains_neuac': contains_neuac,
    }
    table = []
    for precursor_adduct in precursor_adducts:
        product_ions = [
            ('[M+H]+', 1, '+H'),
            ('[M+2H]2+', 2, '+2H'),
        ]
        if precursor_adduct == 'H+NH4':
            product_ions.extend([
                ('[M+NH4]+', 1, '+NH4'),
                ('[M+NH4+H]2+', 2, '+H+NH4'),
            ])
        elif precursor_adduct == '2NH4':
            product_ions.extend([
                ('[M+NH4]+', 1, '+NH4'),
                ('[M+2NH4]2+', 2, '+2NH4'),
            ])

        for mass_column, charge, fragment_adduct in product_ions:
            if mass_column not in dbf.columns:
                continue
            table.append(
                pd.DataFrame({
                    **common,
                    'Theo_mz': dbf[mass_column],
                    'Charge': charge,
                    'Adduct': fragment_adduct,
                    '_PrecursorAdduct': precursor_adduct,
                })
            )

    adduct_df = pd.concat(table, ignore_index=True)
    adduct_df = adduct_df.dropna(subset=['Theo_mz'])
    adduct_df['neutral_loss'] = ""
    adduct_df['fragment_annotation'] = adduct_df['Fragment'].astype(str)

    methanol_losses = adduct_df.loc[adduct_df['contains_neuac']].copy()
    if not methanol_losses.empty:
        methanol_losses['Theo_mz'] = (
            methanol_losses['Theo_mz']
            - METHANOL_MASS / methanol_losses['Charge']
        )
        methanol_losses['neutral_loss'] = "CH3OH"
        methanol_losses['fragment_annotation'] = (
            methanol_losses['Fragment'].astype(str) + "-CH3OH"
        )
        adduct_df = pd.concat([adduct_df, methanol_losses], ignore_index=True)
    return adduct_df


def _build_decoy_adduct_table(
    target_adducts: pd.DataFrame,
    fragment_mass_tolerance: float,
    decoy_seed: int,
    fragment_mass_tolerance_unit: str = "Da",
) -> pd.DataFrame:
    """Create one deterministic, charge-aware shifted ion per target ion."""
    if target_adducts.empty:
        return target_adducts.copy()
    fragment_mass_tolerance, fragment_mass_tolerance_unit = validate_tolerance(
        fragment_mass_tolerance, fragment_mass_tolerance_unit
    )
    target_mzs = np.sort(target_adducts['Theo_mz'].to_numpy(dtype=float))

    def collides_with_target(mz: float) -> bool:
        radius = float(
            conservative_query_radius_da(
                np.asarray([mz]),
                fragment_mass_tolerance,
                fragment_mass_tolerance_unit,
            )[0]
        )
        left = int(np.searchsorted(target_mzs, mz - radius, side="left"))
        right = int(np.searchsorted(target_mzs, mz + radius, side="right"))
        neighbors = target_mzs[left:right]
        return bool(
            neighbors.size
            and np.any(
                within_tolerance(
                    mz,
                    neighbors,
                    fragment_mass_tolerance,
                    fragment_mass_tolerance_unit,
                )
            )
        )

    decoys = target_adducts.copy()
    shifted_mzs = []
    neutral_shifts = []
    for _, row in decoys.iterrows():
        stable_key = "|".join(
            [
                str(decoy_seed),
                str(row.get('IUPAC', '')),
                str(row.get('Fragment', '')),
                str(row.get('FragmentType', '')),
                str(row.get('fragment_annotation', '')),
                str(row.get('Adduct', '')),
            ]
        )
        charge = max(int(row['Charge']), 1)
        shifted_mz = np.nan
        neutral_shift = np.nan
        for attempt in range(32):
            digest = hashlib.sha256(
                f"{stable_key}|{attempt}".encode("utf-8")
            ).digest()
            unit_interval = int.from_bytes(digest[:8], "big") / float(2**64)
            candidate_shift = 1.0 + 29.0 * unit_interval
            candidate_mz = float(row['Theo_mz']) + candidate_shift / charge
            if not collides_with_target(candidate_mz):
                shifted_mz = candidate_mz
                neutral_shift = candidate_shift
                break
        if not np.isfinite(shifted_mz):
            # A dense target library can occasionally consume all hashed
            # proposals. Scan the allowed window so a paired ion is never
            # silently dropped from the decoy library.
            for candidate_shift in np.linspace(1.0, 30.0, 29001):
                candidate_mz = float(row['Theo_mz']) + candidate_shift / charge
                if not collides_with_target(candidate_mz):
                    shifted_mz = candidate_mz
                    neutral_shift = float(candidate_shift)
                    break
        if not np.isfinite(shifted_mz):
            raise ValueError(
                "Could not place a shifted decoy ion outside the target "
                "fragment tolerance within the 1-30 Da shift window"
            )
        shifted_mzs.append(shifted_mz)
        neutral_shifts.append(neutral_shift)

    decoys['Theo_mz'] = shifted_mzs
    decoys['decoy_neutral_mass_shift'] = neutral_shifts
    decoys = decoys.dropna(subset=['Theo_mz']).reset_index(drop=True)
    decoys['is_decoy'] = True
    return decoys


def _match_adduct_library(
    ms2f: pd.DataFrame,
    adduct_df: pd.DataFrame,
    fragment_mass_tolerance: float,
    is_decoy: bool,
    fragment_mass_tolerance_unit: str = "Da",
) -> pd.DataFrame:
    """Match one target or decoy product-ion library to the retained spectra."""
    if adduct_df.empty:
        return pd.DataFrame()
    fragment_mass_tolerance, fragment_mass_tolerance_unit = validate_tolerance(
        fragment_mass_tolerance, fragment_mass_tolerance_unit
    )
    theo_mzs = adduct_df['Theo_mz'].to_numpy(dtype=float)
    tree = cKDTree(theo_mzs.reshape(-1, 1))
    obs_mzs = ms2f['fragment_mz'].to_numpy(dtype=float)
    query_radii = conservative_query_radius_da(
        obs_mzs, fragment_mass_tolerance, fragment_mass_tolerance_unit
    )
    neighbors = tree.query_ball_point(
        obs_mzs.reshape(-1, 1), r=query_radii
    )
    base_cols = [
        'scan_number', 'rt', 'precursor_mz', 'precursor_intensity',
        'fragment_intensity', 'Glycan', 'PrecursorAdduct', 'database_mz',
        'precursor_ppm_error'
    ]
    base_data = ms2f[base_cols].to_dict('records')
    matched = []
    for i, nbrs in enumerate(neighbors):
        if not nbrs:
            continue
        hits = adduct_df.iloc[nbrs].copy()
        hits = hits.loc[
            hits['_PrecursorAdduct'].eq(str(base_data[i]['PrecursorAdduct']))
        ]
        if hits.empty:
            continue
        hits['mz_diff'] = np.abs(hits['Theo_mz'].to_numpy() - obs_mzs[i])
        hits['ppm_error'] = mass_error_ppm(
            obs_mzs[i], hits['Theo_mz'].to_numpy(dtype=float)
        )
        hits = hits.loc[
            within_tolerance(
                obs_mzs[i],
                hits['Theo_mz'].to_numpy(dtype=float),
                fragment_mass_tolerance,
                fragment_mass_tolerance_unit,
            )
        ]
        if hits.empty:
            continue
        hits = hits.sort_values('mz_diff').drop_duplicates(
            subset=['IUPAC'], keep='first'
        )
        for _, hit in hits.iterrows():
            diff = float(hit['mz_diff'])
            row = dict(base_data[i])
            row.update({
                'Fragment': hit['fragment_annotation'],
                'BaseFragment': hit['Fragment'],
                'FragmentType': hit['FragmentType'],
                'fragment_annotation': hit['fragment_annotation'],
                'fragment_iupac': hit['fragment_iupac'],
                'contains_neuac': bool(hit['contains_neuac']),
                'neutral_loss': hit['neutral_loss'],
                'Charge': int(hit['Charge']),
                'Fragment_mz': float(hit['Theo_mz']),
                'observed_fragment_mz': float(obs_mzs[i]),
                'theoretical_fragment_mz': float(hit['Theo_mz']),
                'mz_diff': diff,
                'ppm_error': float(hit['ppm_error']),
                'fragment_mass_tolerance_value': fragment_mass_tolerance,
                'fragment_mass_tolerance_unit': fragment_mass_tolerance_unit,
                'Adduct': hit['Adduct'],
                'IUPAC': hit['IUPAC'],
                'NumericalComposition': hit['NumericalComposition'],
                'Composition': hit['Composition'],
                'is_decoy': bool(is_decoy),
                'decoy_neutral_mass_shift': hit.get(
                    'decoy_neutral_mass_shift', np.nan
                ),
            })
            matched.append(row)
    return pd.DataFrame(matched)


def _select_best_iupac(matches: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series | None]:
    """Apply the same structure-level ranking to target or decoy matches."""
    if matches.empty:
        return matches.copy(), None
    scores = (
        matches.groupby('IUPAC')
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
    selected = matches.loc[matches['IUPAC'] == best['IUPAC']].copy()
    selected['IUPAC_match_count'] = int(best['match_count'])
    selected['IUPAC_unique_fragments'] = int(best['unique_fragments'])
    selected['IUPAC_total_intensity'] = float(best['total_intensity'])
    return selected, best


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

def cluster_1d(
    mzs: np.ndarray,
    tol: float,
    tolerance_unit: str = "Da",
    maximum_cluster_gap_da: float | None = None,
) -> np.ndarray:
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
    tol, tolerance_unit = validate_tolerance(tol, tolerance_unit)
    for i in range(1, len(sorted_mzs)):
        gap = sorted_mzs[i] - sorted_mzs[i - 1]
        same_cluster = bool(
            within_tolerance(
                sorted_mzs[i], sorted_mzs[i - 1], tol, tolerance_unit
            )
        )
        if maximum_cluster_gap_da is not None:
            same_cluster = same_cluster and gap <= maximum_cluster_gap_da
        if not same_cluster:
            cl += 1
        labels[i] = cl
    # invert the sort
    inv = np.empty_like(idx)
    inv[idx] = np.arange(len(idx))
    return labels[inv]


def preprocess_ms2_data(
    ms2: pd.DataFrame,
    fragment_mass_tol: float = 0.1,
    fragment_mass_tol_unit: str = "Da",
) -> pd.DataFrame:
    """
    For each (scan_number, precursor_mz) group, cluster fragment_mz via cluster_1d
    and collapse to one row per cluster (mean m/z, summed intensity).
    """
    out = []
    grouped = ms2.groupby(['scan_number', 'precursor_mz'])
    for (scan, prec), grp in grouped:
        mzs = grp['fragment_mz'].to_numpy()
        labels = cluster_1d(
            mzs,
            fragment_mass_tol,
            fragment_mass_tol_unit,
            maximum_cluster_gap_da=0.1,
        )
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
    max_cleavages: int = 2,
    generate_decoys: bool = False,
    decoy_seed: int = 1729,
    fragment_mass_tol_unit: str = "Da",
) -> pd.DataFrame:
    """
    Vectorized MS2 matching using cKDTree on theoretical fragments generated
    from the N-glycan IUPAC database.

    For each numerical composition, all parsable IUPAC candidate structures are
    fragmented, observed MS2 fragments are matched to those theoretical
    fragments, and the best-scoring IUPAC structure is returned.
    Protonated product ions are always considered. Ammonium-bearing product
    ions are additionally considered for precursor assignments labeled H+NH4
    or 2NH4.

    :param fragment_mass_tol: Numeric fragment m/z matching tolerance.
    :param fragment_mass_tol_unit: Fragment tolerance unit, ``Da`` or ``ppm``.
    :param ppm_tol: Precursor matching tolerance in ppm.
    :param db_path: Path to N-glycan database CSV/Excel file. If None, uses default.
    :param generate_decoys: Also return matches to a paired shifted-ion library.
    :param decoy_seed: Reproducible seed used to construct fragment mass shifts.
    """
    if db_path is None:
        db_path = DEFAULT_PRECURSOR_DB

    fragment_mass_tol, fragment_mass_tol_unit = validate_tolerance(
        fragment_mass_tol, fragment_mass_tol_unit
    )
    fragment_mass_tol_unit = normalize_tolerance_unit(fragment_mass_tol_unit)
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

    # 3) select MS2 rows whose precursor_mz matches an MS1 precursor (track adduct per precursor)
    precursor_matched_data = precursor_matched_data.copy()
    precursor_matched_data['Glycan'] = (
        precursor_matched_data['Glycan']
        .astype(str)
        .apply(_normalize_glycan)
    )
    precursor_columns = ['precursor_mz', 'Adduct']
    for optional_column in ('database_mz', 'ppm_error'):
        if optional_column in precursor_matched_data.columns:
            precursor_columns.append(optional_column)
    prec_rows = precursor_matched_data[
        precursor_matched_data['Glycan'] == precursor_composition
    ][precursor_columns]
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
        ].copy()
        if not sel.empty:
            if 'precursor_intensity' not in sel.columns:
                sel['precursor_intensity'] = 0.0
            sel['Glycan'] = precursor_composition
            sel['PrecursorAdduct'] = adduct_label
            sel['database_mz'] = row.get('database_mz', np.nan)
            sel['precursor_ppm_error'] = row.get('ppm_error', np.nan)
            sel['_precursor_assignment_error'] = abs(
                float(row.get('ppm_error', np.inf))
            )
            filt.append(sel)
    if not filt:
        logger.info(f"No MS2 data after filtering for {precursor_composition}")
        return pd.DataFrame()
    ms2f = pd.concat(filt, ignore_index=True)
    # A scan can satisfy more than one precursor-database row. Preserve only
    # the closest precursor assignment before fragment clustering so repeated
    # selection does not multiply the observed fragment intensity.
    ms2f = (
        ms2f.sort_values('_precursor_assignment_error')
        .drop_duplicates(
            subset=['scan_number', 'precursor_mz', 'fragment_mz'], keep='first'
        )
        .drop(columns='_precursor_assignment_error')
        .reset_index(drop=True)
    )

    # 4) cluster fragments within each spectrum
    ms2f = preprocess_ms2_data(
        ms2f,
        fragment_mass_tol=fragment_mass_tol,
        fragment_mass_tol_unit=fragment_mass_tol_unit,
    )

    # 5) perform vectorized adduct matching against all candidate IUPAC structures
    adduct_df = _build_adduct_table(dbf, ms2f['PrecursorAdduct'].unique())
    if adduct_df.empty:
        logger.info(f"No theoretical fragment adducts for {precursor_composition}")
        return pd.DataFrame()

    target_matches = _match_adduct_library(
        ms2f,
        adduct_df,
        fragment_mass_tol,
        is_decoy=False,
        fragment_mass_tolerance_unit=fragment_mass_tol_unit,
    )
    target_out, target_best = _select_best_iupac(target_matches)

    decoy_out = pd.DataFrame()
    if generate_decoys:
        decoy_adducts = _build_decoy_adduct_table(
            adduct_df,
            fragment_mass_tolerance=fragment_mass_tol,
            decoy_seed=int(decoy_seed),
            fragment_mass_tolerance_unit=fragment_mass_tol_unit,
        )
        decoy_matches = _match_adduct_library(
            ms2f,
            decoy_adducts,
            fragment_mass_tol,
            is_decoy=True,
            fragment_mass_tolerance_unit=fragment_mass_tol_unit,
        )
        decoy_out, _ = _select_best_iupac(decoy_matches)

    if target_out.empty and decoy_out.empty:
        logger.info(f"No target or decoy fragments matched for {precursor_composition}")
        return pd.DataFrame()
    if target_best is not None:
        logger.info(
            f"Selected IUPAC for {precursor_composition}: {target_best['IUPAC']}"
        )
        logger.info(
            f"Structure score: {int(target_best['match_count'])} matches, "
            f"{int(target_best['unique_fragments'])} unique fragments"
        )
        logger.info(f"Matched {len(target_out)} target fragments for {precursor_composition}")
    if generate_decoys:
        logger.info(
            f"Matched {len(decoy_out)} paired decoy fragments for {precursor_composition}"
        )
    return pd.concat([target_out, decoy_out], ignore_index=True, sort=False)
