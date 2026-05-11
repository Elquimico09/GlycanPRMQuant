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
        plt.style.use(['science', 'no-latex'])
        plt.rcParams['font.family'] = 'Arial'

        if plot:
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
