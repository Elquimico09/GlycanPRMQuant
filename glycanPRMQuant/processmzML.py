import os
import logging
import pandas as pd
import numpy as np

from glycanPRMQuant.spectra import extract_ms2
from glycanPRMQuant.matchMS1     import matchMS1
from glycanPRMQuant.matchMS2     import matchMS2
from glycanPRMQuant.calculateAUC import calculateAUC
from glycanPRMQuant.candidate_scoring import (
    CandidateScoringConfig,
    score_and_resolve_candidates,
)
from glycanPRMQuant.plotFragmentIntensity import plot_ms2_fragments, plot_total_chromatogram_with_window
from glycanPRMQuant.plotMS2spectrum   import plotMS2spectrum
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

logger = logging.getLogger(__name__)
ISOBARIC_RESOLUTION_PPM = 20.0

def _smooth_signal(y, method: str, window: int):
    if not window or window <= 0:
        return y
    method = (method or "gaussian").lower()
    if method in ("gaussian", "gauss"):
        return gaussian_filter1d(y, sigma=window, mode='nearest')
    if method in ("savgol", "sav-gol", "savitzky-golay", "sg"):
        n = len(y)
        if n < 3:
            return y
        win = int(window)
        if win % 2 == 0:
            win += 1
        if win < 3:
            win = 3
        if win > n:
            win = n if n % 2 == 1 else n - 1
            if win < 3:
                return y
        return savgol_filter(y, window_length=win, polyorder=2, mode='nearest')
    return y

def _resample_uniform(rt, y):
    rt = np.asarray(rt, dtype=float)
    y = np.asarray(y, dtype=float)
    if rt.size < 3:
        return rt, y
    order = np.argsort(rt)
    rt = rt[order]
    y = y[order]
    diffs = np.diff(rt)
    step = np.median(diffs[diffs > 0]) if np.any(diffs > 0) else None
    if step is None or step <= 0:
        return rt, y
    grid = np.arange(rt.min(), rt.max() + step * 0.5, step)
    y_interp = np.interp(grid, rt, y)
    return grid, y_interp

def _build_common_grid(rts):
    rts = np.asarray(rts, dtype=float)
    if rts.size < 3:
        return rts
    rts = np.sort(rts)
    diffs = np.diff(rts)
    step = np.median(diffs[diffs > 0]) if np.any(diffs > 0) else None
    if step is None or step <= 0:
        return rts
    return np.arange(rts.min(), rts.max() + step * 0.5, step)

def _build_total_trace_like_total_plot(
    sub: pd.DataFrame,
    smoothing_window: int,
    smoothing_method: str
):
    """
    Build the total trace using the same per-adduct interpolation + sum approach
    as plot_ms2_fragments(..., group_col=None) and plot_total_chromatogram_with_window.
    """
    if 'PrecursorAdduct' in sub.columns:
        agg = (
            sub.groupby(['scan_number', 'rt', 'PrecursorAdduct'])
               .agg(sum_intensity=('fragment_intensity', 'sum'))
               .reset_index()
        )
        all_rts = agg['rt'].values
        x = _build_common_grid(all_rts)
        y = np.zeros_like(x, dtype=float)
        for g in agg['PrecursorAdduct'].unique().tolist():
            ad = agg[agg['PrecursorAdduct'] == g].sort_values('rt')
            xi = ad['rt'].values
            yi = ad['sum_intensity'].values
            if xi.size == 0:
                continue
            yi = np.interp(x, xi, yi, left=0.0, right=0.0)
            if smoothing_window and smoothing_window > 0:
                yi = _smooth_signal(yi, smoothing_method, smoothing_window)
            y += yi
        return x, y

    total = (
        sub.groupby(['scan_number', 'rt'])['fragment_intensity']
           .sum()
           .reset_index()
           .sort_values('rt')
    )
    x = total['rt'].values
    y = total['fragment_intensity'].values
    if smoothing_window and smoothing_window > 0:
        xg, yg = _resample_uniform(x, y)
        y = _smooth_signal(yg, smoothing_method, smoothing_window)
        x = xg
    return x, y

def _compute_total_window_boundaries(
    all_df: pd.DataFrame,
    smoothing_window: int,
    smoothing_method: str,
    rel_height: float,
    rel_height_mode: str
) -> pd.DataFrame:
    """
    Compute AUC boundaries from the same total trace model used for the displayed
    total chromatogram (per-adduct interpolation + sum), then reuse calculateAUC
    boundary logic on that constructed trace.
    """
    rows = []
    if all_df.empty:
        return pd.DataFrame(columns=['Glycan', 'peak_rt', 'start_rt', 'end_rt', 'AUC'])

    for glycan, sub in all_df.groupby('Glycan'):
        x, y = _build_total_trace_like_total_plot(sub, smoothing_window, smoothing_method)
        if x.size == 0:
            continue
        trace_df = pd.DataFrame({
            'Glycan': [glycan] * x.size,
            'PrecursorAdduct': ['ALL'] * x.size,
            'scan_number': np.arange(x.size, dtype=int),
            'rt': x,
            'fragment_intensity': y,
        })
        per_adduct_df, _ = calculateAUC(
            trace_df,
            adduct_col='PrecursorAdduct',
            smoothing_window=0,
            rel_height=rel_height,
            rel_height_mode=rel_height_mode
        )
        if per_adduct_df.empty:
            continue
        rows.append(per_adduct_df.iloc[0].to_dict())

    return pd.DataFrame(rows)

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

def _compute_glycan_apex_rt(all_df: pd.DataFrame, smoothing_window: int, smoothing_method: str) -> dict:
    apex = {}
    if all_df.empty:
        return apex
    summed = (
        all_df.groupby(['Glycan', 'scan_number'])
              .agg(rt=('rt', 'first'),
                   summed_intensity=('fragment_intensity', 'sum'))
              .reset_index()
    )
    for glycan, sub in summed.groupby('Glycan'):
        sub = sub.sort_values('rt')
        x = sub['rt'].to_numpy()
        y = sub['summed_intensity'].to_numpy()
        if smoothing_window and smoothing_window > 0:
            xg, yg = _resample_uniform(x, y)
            y = _smooth_signal(yg, smoothing_method, smoothing_window)
            x = xg
        if y.size == 0:
            continue
        apex[glycan] = float(x[int(y.argmax())])
    return apex

def _assign_precursor_clusters(all_df: pd.DataFrame, ppm_tol: float = ISOBARIC_RESOLUTION_PPM) -> pd.Series:
    """Cluster precursor m/z values within a ppm window for conflict resolution."""
    precursors = pd.to_numeric(all_df['precursor_mz'], errors='coerce')
    unique_mzs = np.sort(precursors.dropna().unique())
    cluster_by_mz = {}
    cluster_id = -1
    anchor_mz = None

    for mz in unique_mzs:
        if anchor_mz is None or abs(mz - anchor_mz) > anchor_mz * ppm_tol / 1e6:
            cluster_id += 1
            anchor_mz = mz
        cluster_by_mz[mz] = cluster_id

    clusters = precursors.map(cluster_by_mz)
    if clusters.isna().any():
        next_id = cluster_id + 1
        for idx in clusters[clusters.isna()].index:
            clusters.loc[idx] = next_id
            next_id += 1
    return clusters.astype(int)


def _resolve_precursor_conflicts(
    all_df: pd.DataFrame,
    ppm_tol: float = ISOBARIC_RESOLUTION_PPM
) -> pd.DataFrame:
    """
    For precursor m/z clusters that map to multiple glycans, keep the glycan with
    the most fragment matches and (tie-breaker) highest total fragment intensity.
    """
    if all_df.empty:
        return all_df
    work = all_df.copy()
    work['_precursor_cluster'] = _assign_precursor_clusters(work, ppm_tol=ppm_tol)
    scores = (
        work.groupby(['_precursor_cluster', 'Glycan'])
              .agg(
                  frag_count=('fragment_mz', 'size'),
                  total_intensity=('fragment_intensity', 'sum')
              )
              .reset_index()
    )
    scores = scores.sort_values(
        ['_precursor_cluster', 'frag_count', 'total_intensity'],
        ascending=[True, False, False]
    )
    winners = scores.drop_duplicates(subset=['_precursor_cluster'], keep='first')
    winner_map = winners.set_index('_precursor_cluster')['Glycan'].to_dict()
    keep = work['Glycan'] == work['_precursor_cluster'].map(winner_map)
    return all_df.loc[keep].copy()


def _count_precursor_conflicts(
    all_df: pd.DataFrame,
    ppm_tol: float = ISOBARIC_RESOLUTION_PPM
) -> int:
    """Return the number of precursor m/z clusters assigned to multiple glycans."""
    if all_df.empty:
        return 0
    work = all_df.copy()
    work['_precursor_cluster'] = _assign_precursor_clusters(work, ppm_tol=ppm_tol)
    return int(
        work.groupby('_precursor_cluster')['Glycan']
            .nunique()
            .gt(1)
            .sum()
    )


def _filter_ms1_to_resolved_assignments(
    ms1_results: pd.DataFrame,
    resolved_ms2: pd.DataFrame,
    ppm_tol: float = ISOBARIC_RESOLUTION_PPM
) -> pd.DataFrame:
    """Keep MS1 rows whose precursor cluster and glycan survived MS2 resolution."""
    if ms1_results.empty or resolved_ms2.empty:
        return ms1_results.iloc[0:0].copy()

    ms1_work = ms1_results.copy()
    ms2_work = resolved_ms2.copy()
    combined = pd.concat(
        [
            ms1_work[['precursor_mz']].assign(_source='ms1'),
            ms2_work[['precursor_mz']].assign(_source='ms2'),
        ],
        ignore_index=True
    )
    combined['_precursor_cluster'] = _assign_precursor_clusters(combined, ppm_tol=ppm_tol)
    ms1_work['_precursor_cluster'] = combined.loc[combined['_source'] == 'ms1', '_precursor_cluster'].to_numpy()
    ms2_work['_precursor_cluster'] = combined.loc[combined['_source'] == 'ms2', '_precursor_cluster'].to_numpy()
    ms1_work['_glycan_norm'] = ms1_work['Glycan'].apply(_normalize_glycan)
    ms2_work['_glycan_norm'] = ms2_work['Glycan'].apply(_normalize_glycan)

    keep_pairs = set(
        ms2_work[['_precursor_cluster', '_glycan_norm']]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    keep = [
        (cluster, glycan) in keep_pairs
        for cluster, glycan in zip(ms1_work['_precursor_cluster'], ms1_work['_glycan_norm'])
    ]
    return ms1_work.loc[keep, ms1_results.columns].copy()

def process_mzml_pipeline(
    mzml_file: str,
    output_dir: str,
    ppm_ms1_tol: float = 10,
    mz_offset: float = 0.0,
    mass_offset: float = 0.0,
    intensity_threshold: float = 1e2,
    fragment_mass_tol: float = 0.02,
    smoothing_window: int = 11,
    smoothing_method: str = "gaussian",
    fragment_top_n: int = 10,
    spectrum_window_minutes: float = 2.0,
    enable_adduct_plots: bool = True,
    enable_total_plots: bool = True,
    rel_height: float = 0.7,
    rel_height_mode: str = "prominence",
    skyline_transition: bool = False,
    enable_smoothing: bool = True,
    fragment_ion_series: str = "ABCXYZ",
    fragment_max_cleavages: int = 2,
    precursor_db_path: str = None,
    structure_db_path: str = None,
    resolve_isobaric_conflicts: bool = True,
    candidate_min_fragments: int = 2,
    candidate_min_explained_intensity: float = 0.005,
    candidate_min_score: float = 35.0,
    candidate_min_evidence_difference: float = 4.0,
    candidate_mass_outlier_min_delta: float = 2.0,
):
    base_name = os.path.splitext(os.path.basename(mzml_file))[0]
    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, 'images')
    os.makedirs(images_dir, exist_ok=True)
    logger.info(
        "Isobaric precursor conflict resolution: %s (%g ppm)",
        "enabled" if resolve_isobaric_conflicts else "disabled",
        ISOBARIC_RESOLUTION_PPM
    )

    # 1) Extract MS2
    logger.info("Extracting MS2 data from %s", mzml_file)
    ms2_data = extract_ms2(mzml_file, min_intensity=intensity_threshold)
    logger.info(f"Extracted {len(ms2_data)} points")

    # 2) Match MS1
    logger.info("Matching MS1 precursors…")
    ms1_results = matchMS1(
        ms2_data,
        ppm_tol=ppm_ms1_tol,
        mz_offset=mz_offset,
        mass_offset=mass_offset,
        db_path=precursor_db_path
    )
    ms1_out = os.path.join(output_dir, "ms1_results.csv")
    ms1_results.to_csv(ms1_out, index=False)
    logger.info(f"Wrote {len(ms1_results)} MS1 matches")

    if ms1_results.empty:
        logger.info("No MS1 matches; done.")
        return

    # Pre-filter MS2 by matched precursor m/z values (speed-up)
    precs = ms1_results['precursor_mz'].unique()
    if precs.size > 0:
        mask = pd.Series(False, index=ms2_data.index)
        for p in precs:
            tol = p * ppm_ms1_tol / 1e6
            mask |= ms2_data['precursor_mz'].between(p - tol, p + tol)
        ms2_data = ms2_data.loc[mask].reset_index(drop=True)
        logger.info(f"Pre-filtered MS2 to {len(ms2_data)} rows for matched precursors")

    all_matched = []

    # build clean glycan list
    glycans = set()
    for entry in ms1_results['Glycan'].dropna().astype(str):
        for comp in entry.split(';'):
            c = _normalize_glycan(comp)
            if c and c.lower() != 'nan':
                glycans.add(c)

    # 3) For each glycan: match MS2 (collect only)
    adduct_charge = {'2H':2,'3H':3,'4H':4,'H+NH4':2,'2NH4':2}
    effective_window = smoothing_window if enable_smoothing else 0

    for glycan in sorted(glycans):
        logger.info(f"Processing glycan {glycan!r}...")
        matched_ms2 = matchMS2(
            ms2_data,
            ms1_results,
            precursor_composition=glycan,
            ppm_tol=ppm_ms1_tol,
            fragment_mass_tol=fragment_mass_tol,
            intensity_threshold=intensity_threshold,
            ion_series=fragment_ion_series,
            max_cleavages=fragment_max_cleavages,
            db_path=structure_db_path
        )
        if matched_ms2.empty:
            logger.info(f"No MS2 fragments for {glycan!r}")
            continue

        # **Normalize** Fragment_mz -> fragment_mz so downstream code sees it
        if 'Fragment_mz' in matched_ms2.columns:
            matched_ms2 = matched_ms2.rename(columns={'Fragment_mz': 'fragment_mz'})

        all_matched.append(matched_ms2)

    if not all_matched:
        logger.info("No MS2 matches; done.")
        return

    # Score and resolve precursor conflicts across glycans. Ambiguous and
    # possibly coisolated candidates are retained rather than forcing a winner.
    all_df = pd.concat(all_matched, ignore_index=True)
    conflict_count = _count_precursor_conflicts(all_df, ppm_tol=ISOBARIC_RESOLUTION_PPM)
    rows_before = len(all_df)
    glycans_before = all_df['Glycan'].nunique()
    scoring_config = CandidateScoringConfig(
        precursor_cluster_ppm=ISOBARIC_RESOLUTION_PPM,
        fragment_mass_tolerance=fragment_mass_tol,
        minimum_distinct_fragments=candidate_min_fragments,
        minimum_explained_intensity=candidate_min_explained_intensity,
        minimum_candidate_score=candidate_min_score,
        minimum_discriminative_evidence_difference=candidate_min_evidence_difference,
        mass_outlier_min_delta_ppm=candidate_mass_outlier_min_delta,
    )
    scoring_result = score_and_resolve_candidates(
        all_df,
        scan_data=ms2_data,
        config=scoring_config,
        resolve=resolve_isobaric_conflicts,
    )
    all_df = scoring_result.resolved_rows
    score_path = os.path.join(output_dir, "candidate_scores.csv")
    scoring_result.candidate_scores.to_csv(score_path, index=False)
    not_quantified = scoring_result.reported_rows.loc[
        ~scoring_result.reported_rows['selected']
    ].copy()
    not_quantified_path = os.path.join(output_dir, "candidate_rows_not_quantified.csv")
    not_quantified.to_csv(not_quantified_path, index=False)
    status_counts = scoring_result.candidate_scores['resolution_status'].value_counts().to_dict()
    logger.info("Wrote candidate scoring audit table to %s", score_path)
    logger.info(
        "Wrote %d matched row(s) excluded from composition-level quantification to %s",
        len(not_quantified),
        not_quantified_path,
    )
    logger.info(
        "Scored %d precursor conflict(s) within %g ppm: %d -> %d rows, %d -> %d glycans; statuses=%s",
        conflict_count,
        ISOBARIC_RESOLUTION_PPM,
        rows_before,
        len(all_df),
        glycans_before,
        all_df['Glycan'].nunique(),
        status_counts,
    )
    ms1_resolved = _filter_ms1_to_resolved_assignments(
        ms1_results,
        all_df,
        ppm_tol=ISOBARIC_RESOLUTION_PPM
    )
    ms1_resolved_out = os.path.join(output_dir, "ms1_results_resolved.csv")
    ms1_resolved.to_csv(ms1_resolved_out, index=False)
    logger.info(f"Wrote {len(ms1_resolved)} resolved MS1 matches to {ms1_resolved_out}")

    # 4) For each glycan after conflict resolution: save CSV, plot
    for glycan, sub in all_df.groupby('Glycan'):
        csv_path = os.path.join(output_dir, f"ms2_{glycan}.csv")
        sub.to_csv(csv_path, index=False)
        logger.info(f"Wrote {len(sub)} MS2 matches to {csv_path}")

        # Chromatograms per fragment (legacy view)
        chrom_frag_pdf = os.path.join(images_dir, f"ms2_{glycan}.pdf")
        try:
            plot_ms2_fragments(csv_path, window=effective_window,
                               top_n=fragment_top_n, save_path=chrom_frag_pdf,
                               figsize=(6.2, 4),
                               group_col='Fragment',
                               smoothing_method=smoothing_method)
            logger.info(f"Saved fragment-level chromatogram to {chrom_frag_pdf}")
        except Exception as e:
            logger.warning(f"Fragment chromatogram failed: {e}")

        # Chromatogram per precursor adduct
        if enable_adduct_plots:
            chrom_adduct_pdf = os.path.join(images_dir, f"ms2_{glycan}_by_precursor_adduct.pdf")
            try:
                plot_ms2_fragments(csv_path, window=effective_window,
                                   top_n=None, save_path=chrom_adduct_pdf,
                                   group_col='PrecursorAdduct',
                                   smoothing_method=smoothing_method)
                logger.info(f"Saved precursor-adduct chromatogram to {chrom_adduct_pdf}")
            except Exception as e:
                logger.warning(f"Precursor-adduct chromatogram failed: {e}")

        # Total chromatogram (all fragments/adducts summed)
        if enable_total_plots:
            chrom_total_pdf = os.path.join(images_dir, f"ms2_{glycan}_total.pdf")
            try:
                plot_ms2_fragments(csv_path, window=effective_window,
                                   top_n=None, save_path=chrom_total_pdf,
                                   group_col=None,
                                   smoothing_method=smoothing_method)
                logger.info(f"Saved total chromatogram to {chrom_total_pdf}")
            except Exception as e:
                logger.warning(f"Total chromatogram failed: {e}")

        # averaged spectrum
        spec_pdf = os.path.join(images_dir, f"ms2_{glycan}_ms2spectrum.pdf")
        try:
            plotMS2spectrum(csv_path, window_minutes=spectrum_window_minutes,
                            top_n=fragment_top_n, save_path=spec_pdf)
            logger.info(f"Saved spectrum to {spec_pdf}")
        except Exception as e:
            logger.warning(f"Spectrum plot failed: {e}")
    # 4) AUC
    if not all_df.empty:
        logger.info("Calculating AUC...")
        per_adduct_df, total_df = calculateAUC(
            all_df,
            smoothing_window=effective_window,
            smoothing_method=smoothing_method,
            adduct_col='PrecursorAdduct',
            rel_height=rel_height,
            rel_height_mode=rel_height_mode
        )

        # Total AUC (backward compatible filename)
        auc_path = os.path.join(output_dir, f"{base_name}_auc_values.csv")
        total_df.to_csv(auc_path, index=False)
        logger.info(f"Wrote total AUC (summed across adducts) to {auc_path}")

        # Per-adduct detail
        auc_adduct_path = os.path.join(output_dir, f"{base_name}_auc_values_by_adduct.csv")
        per_adduct_df.to_csv(auc_adduct_path, index=False)
        logger.info(f"Wrote per-adduct AUC values to {auc_adduct_path}")

        # Compute total-window boundaries for shaded AUC plots
        total_window_df = _compute_total_window_boundaries(
            all_df,
            smoothing_window=effective_window,
            smoothing_method=smoothing_method,
            rel_height=rel_height,
            rel_height_mode=rel_height_mode
        )
        if not total_window_df.empty:
            total_window_df = total_window_df.set_index('Glycan')
        else:
            total_window_df = pd.DataFrame(columns=['peak_rt', 'start_rt', 'end_rt', 'AUC'])

    # 5) Total chromatogram with AUC window shading
    if not all_df.empty:
        for glycan, sub in all_df.groupby('Glycan'):
            if glycan not in total_window_df.index:
                continue
            start_rt = float(total_window_df.loc[glycan, 'start_rt'])
            end_rt = float(total_window_df.loc[glycan, 'end_rt'])
            csv_path = os.path.join(output_dir, f"ms2_{glycan}.csv")
            shaded_pdf = os.path.join(images_dir, f"ms2_{glycan}_total_auc.pdf")
            try:
                plot_total_chromatogram_with_window(
                    csv_path,
                    window=effective_window,
                    save_path=shaded_pdf,
                    start_rt=start_rt,
                    end_rt=end_rt
                )
                logger.info(f"Saved total chromatogram with AUC window to {shaded_pdf}")
            except Exception as e:
                logger.warning(f"Total AUC chromatogram failed: {e}")

    # 6) Skyline transition export (unique fragments, apex RT per glycan)
    if skyline_transition and not all_df.empty:
        apex_rt = _compute_glycan_apex_rt(all_df, effective_window, smoothing_method)
        charge = all_df['PrecursorAdduct'].map(adduct_charge).fillna(0).astype(int)
        trans_df = pd.DataFrame({
            'Glycan': all_df['Glycan'],
            'Adduct': all_df['PrecursorAdduct'],
            'Precursor m/z': all_df['precursor_mz'],
            'Precursor charge': charge,
            'Fragment m/z': all_df['fragment_mz'],
            'Fragment charge': all_df['Charge'],
        })
        trans_df = trans_df.drop_duplicates(subset=[
            'Glycan', 'Adduct', 'Precursor m/z', 'Precursor charge', 'Fragment m/z', 'Fragment charge'
        ])
        trans_df['Retention time'] = trans_df['Glycan'].map(apex_rt)
        trans_path = os.path.join(output_dir, f"{base_name}_skyline_transitions.xlsx")
        trans_df.to_excel(trans_path, index=False)
        logger.info(f"Wrote Skyline transition list to {trans_path}")

    logger.info("Done.")
