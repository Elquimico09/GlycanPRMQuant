import os
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

logger = logging.getLogger(__name__)

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

def _find_width_at_height(y: np.ndarray, peak_idx: int, height: float):
    y = np.asarray(y, dtype=float)
    n = y.size
    if n == 0:
        return 0.0, float(max(n - 1, 0))
    peak_idx = int(np.clip(peak_idx, 0, n - 1))
    height = float(height)

    left_candidates = np.where(y[:peak_idx + 1] <= height)[0]
    if left_candidates.size == 0:
        left_ip = 0.0
    else:
        li = left_candidates[-1]
        if li == peak_idx:
            left_ip = float(li)
        else:
            y1, y2 = y[li], y[li + 1]
            if y2 == y1:
                left_ip = float(li)
            else:
                frac = (height - y1) / (y2 - y1)
                left_ip = li + float(frac)

    right_candidates = np.where(y[peak_idx:] <= height)[0]
    if right_candidates.size == 0:
        right_ip = float(n - 1)
    else:
        ri = peak_idx + right_candidates[0]
        if ri == peak_idx:
            right_ip = float(ri)
        else:
            y1, y2 = y[ri - 1], y[ri]
            if y2 == y1:
                right_ip = float(ri)
            else:
                frac = (height - y1) / (y2 - y1)
                right_ip = (ri - 1) + float(frac)

    return left_ip, right_ip

def calculateAUC(
    ms2_input,
    glycan_col: str = 'Glycan',
    scan_col: str = 'scan_number',
    rt_col: str = 'rt',
    intensity_col: str = 'fragment_intensity',
    adduct_col: str = 'Adduct',
    rel_height: float = 0.7,
    rel_height_mode: str = "prominence",
    prominence: float = None,
    smoothing_window: int = 30,
    smoothing_method: str = "gaussian",
    plot: bool = False,
    save_path: str = None,
    window = 0
) -> pd.DataFrame:
    """
    Calculate AUC for each glycan/adduct by optionally smoothing the summed fragment-intensity
    chromatogram, detecting the main peak, determining its boundaries at a relative height,
    and integrating the (smoothed or raw) intensity between those boundaries. Also returns a
    glycan-level total that sums AUCs across all adducts for that glycan.

    smoothing_window <= 0 will skip smoothing.

    Parameters
    ----------
    ms2_input : pd.DataFrame or str
        Matched MS2 DataFrame or path to CSV/Excel.
    glycan_col : str
        Column name for glycan composition.
    scan_col : str
        Column name for scan number.
    rt_col : str
        Column name for retention time.
    intensity_col : str
        Column for fragment intensity.
    rel_height : float
        Relative height (0–1) for width calculation (e.g. 0.5 for half-height).
    prominence : float or None
        Minimum peak prominence passed to find_peaks.
    smoothing_window : int
        Smoothing window. For Gaussian, this is sigma. For Sav-Gol, this is window length.
    smoothing_method : str
        "gaussian" (default) or "savgol".
    plot : bool
        If True, plot smoothed vs raw chromatogram and integration window.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        (per_adduct_df, total_df)
        per_adduct_df columns: [glycan_col, adduct_col, 'peak_rt', 'start_rt', 'end_rt', 'AUC']
        total_df      columns: [glycan_col, 'AUC'] (AUC summed across adducts per glycan)
    """
    # load data
    if isinstance(ms2_input, str):
        ext = os.path.splitext(ms2_input)[1].lower()
        if ext == '.csv':
            df = pd.read_csv(ms2_input)
        elif ext in ('.xlsx', '.xls'):
            df = pd.read_excel(ms2_input)
        else:
            raise ValueError("Unsupported file type.")
    else:
        df = ms2_input.copy()

    # validate
    missing = {glycan_col, scan_col, rt_col, intensity_col} - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # If adduct column is absent, treat all signal as one pseudo-adduct so grouping works
    if adduct_col not in df.columns:
        df = df.copy()
        df[adduct_col] = 'ALL'

    # sum per scan for each glycan
    summed = (
        df.groupby([glycan_col, adduct_col, scan_col])
          .agg(rt=(rt_col, 'first'),
               summed_intensity=(intensity_col, 'sum'))
          .reset_index()
    )

    results = []
    for (glycan, adduct), sub in summed.groupby([glycan_col, adduct_col]):
        sub = sub.sort_values('rt')
        x = sub['rt'].to_numpy()
        y = sub['summed_intensity'].to_numpy()

        # apply smoothing if requested
        if smoothing_window and smoothing_window > 0:
            xg, yg = _resample_uniform(x, y)
            y_smooth = _smooth_signal(yg, smoothing_method, smoothing_window)
            x = xg
        else:
            y_smooth = y

        # detect peaks on y_smooth
        peaks, props = find_peaks(y_smooth, prominence=prominence)
        if len(peaks) == 0:
            main_idx = np.argmax(y_smooth)
        else:
            main_idx = peaks[np.argmax(y_smooth[peaks])]

        # compute width at rel_height
        mode = (rel_height_mode or "prominence").lower()
        if mode in ("height", "peak", "absolute"):
            peak_y = float(y_smooth[main_idx])
            height_level = peak_y * (1.0 - float(rel_height))
            left_ip, right_ip = _find_width_at_height(y_smooth, main_idx, height_level)
        else:
            widths, h_eval, left_ips, right_ips = peak_widths(
                y_smooth, [main_idx], rel_height=rel_height
            )
            left_ip, right_ip = left_ips[0], right_ips[0]

        # map to retention time
        idxs = np.arange(len(x))
        start_rt = np.interp(left_ip, idxs, x)
        end_rt   = np.interp(right_ip, idxs, x)
        peak_rt  = x[main_idx]

        # Integrate with interpolated boundary points so narrow windows that
        # fall between scans do not collapse to a single apex sample.
        interior_mask = (x > start_rt) & (x < end_rt)
        x_auc = np.concatenate(([start_rt], x[interior_mask], [end_rt]))
        y_auc = np.interp(x_auc, x, y_smooth)
        auc = np.trapezoid(y_auc, x_auc)

        logger.info(
            f"Glycan {glycan!r}: peak RT={peak_rt:.2f}, "
            f"window=[{start_rt:.2f}, {end_rt:.2f}], AUC={auc:.2f}"
        )
        results.append({
            glycan_col: glycan,
            adduct_col: adduct,
            'peak_rt': peak_rt,
            'start_rt': start_rt,
            'end_rt': end_rt,
            'AUC': auc
        })
        if plot:
            # The plotting style is optional and should not be required for
            # non-plotting AUC calculations (including feature-level scoring).
            try:
                import scienceplots  # noqa: F401

                plt.style.use(['science', 'no-latex'])
            except (ImportError, OSError):
                logger.warning(
                    "SciencePlots style is unavailable; using Matplotlib defaults"
                )
            plt.rcParams['font.family'] = 'Arial'
            fig, ax = plt.subplots(figsize=(4.8, 4))
            ax.plot(x, y_smooth, label=(
                f'smoothed ({smoothing_method}, w={smoothing_window})'
                if smoothing_window and smoothing_window > 0 else 'raw'
            ))
            ax.axvspan(start_rt, end_rt, color='red', alpha=0.1,
                       label='integration window')
            ax.set_xlabel('RT (min)')
            ax.set_ylabel('Intensity')
            plt.xlim(x.min()-window, x.max()+window)
            plt.ylim(0, y_smooth.max() * 1.1)
            ax.set_title(f"{glycan} ({adduct}): Integration Window")
            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300)
                logger.info(f"Saved plot to {save_path}")
            else:
                plt.show()

    per_adduct_df = pd.DataFrame(results)
    total_df = per_adduct_df.groupby(glycan_col, as_index=False)['AUC'].sum()
    return per_adduct_df, total_df


def calculate_feature_auc(
    ms2_input,
    glycan_col: str = "Glycan",
    precursor_cluster_col: str = "precursor_cluster",
    feature_col: str = "rt_feature_id",
    adduct_col: str = "PrecursorAdduct",
    rel_height: float = 0.7,
    rel_height_mode: str = "prominence",
    smoothing_window: int = 30,
    smoothing_method: str = "gaussian",
) -> pd.DataFrame:
    """Calculate an independent AUC for every selected scoring feature.

    Unlike :func:`calculateAUC`, this function does not collapse separate
    chromatographic features that share a glycan composition. Candidate and
    target-decoy audit fields are copied onto each feature-level result so the
    features can subsequently be aligned and ranked across runs.
    """
    if isinstance(ms2_input, str):
        ext = os.path.splitext(ms2_input)[1].lower()
        if ext == ".csv":
            df = pd.read_csv(ms2_input)
        elif ext in (".xlsx", ".xls"):
            df = pd.read_excel(ms2_input)
        else:
            raise ValueError("Unsupported file type.")
    else:
        df = ms2_input.copy()

    grouping_columns = [glycan_col, precursor_cluster_col, feature_col]
    missing = set(grouping_columns) - set(df.columns)
    if missing:
        raise ValueError(f"Feature-level AUC requires columns: {sorted(missing)}")

    audit_columns = [
        "feature_start_rt",
        "feature_apex_rt",
        "feature_end_rt",
        "feature_scan_count",
        "feature_left_flank_scans",
        "feature_right_flank_scans",
        "chromatographic_peak_valid",
        "chromatographic_peak_rejection_reason",
        "candidate_score",
        "discriminative_score",
        "distinct_fragment_count",
        "candidate_specific_fragment_count",
        "explained_intensity_fraction",
        "resolution_status",
        "assignment_q_value",
        "target_decoy_score_margin",
        "target_decoy_pass",
        "quantification_source",
    ]
    quantification_summary_columns = [
        "accepted_transition_count",
        "quantification_scan_count",
        "quantification_trace_point_count",
        "detected_trace_point_count",
        "zero_trace_point_count",
    ]
    rows = []
    for keys, feature_data in df.groupby(grouping_columns, dropna=False):
        per_adduct, _ = calculateAUC(
            feature_data,
            glycan_col=glycan_col,
            adduct_col=adduct_col,
            rel_height=rel_height,
            rel_height_mode=rel_height_mode,
            smoothing_window=smoothing_window,
            smoothing_method=smoothing_method,
            plot=False,
        )
        if per_adduct.empty:
            continue

        metadata = dict(zip(grouping_columns, keys))
        for column in audit_columns:
            if column not in feature_data.columns:
                continue
            values = feature_data[column].dropna()
            metadata[column] = values.iloc[0] if not values.empty else np.nan

        metadata["quantification_scan_count"] = int(
            feature_data["scan_number"].nunique()
        )
        metadata["quantification_trace_point_count"] = int(len(feature_data))
        theoretical_column = next(
            (
                column
                for column in ("theoretical_fragment_mz", "fragment_mz")
                if column in feature_data.columns
            ),
            None,
        )
        if theoretical_column is not None:
            transition_columns = [theoretical_column]
            transition_columns.extend(
                column
                for column in ("Fragment", "Charge", "Adduct")
                if column in feature_data.columns
            )
            metadata["accepted_transition_count"] = int(
                feature_data[transition_columns].drop_duplicates().shape[0]
            )
        if "quantification_peak_detected" in feature_data.columns:
            detected = feature_data["quantification_peak_detected"].fillna(False).astype(bool)
            metadata["detected_trace_point_count"] = int(detected.sum())
            metadata["zero_trace_point_count"] = int((~detected).sum())

        for result in per_adduct.to_dict("records"):
            rows.append({**metadata, **result})

    if not rows:
        return pd.DataFrame(
            columns=grouping_columns
            + [adduct_col, "peak_rt", "start_rt", "end_rt", "AUC"]
            + audit_columns
            + quantification_summary_columns
        )
    return pd.DataFrame(rows)
