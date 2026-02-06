import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
import scienceplots

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

def plot_ms2_fragments(ms2_csv_file, window=11, top_n=None, save_path=None,
                       figsize=(10, 5), group_col='Fragment', smoothing_method: str = "gaussian"):
    """
    Reads an MS2 results CSV file, aggregates fragment data per scan grouped by `group_col`,
    applies smoothing on a uniform RT grid if window > 0, and plots summed intensity vs. RT
    with one smoothed line per group plus the total intensity line, shading under each curve.

    Parameters
    ----------
    ms2_csv_file : str
        Path to the MS2 CSV file.
    window : int
        Smoothing window. For Gaussian, this is sigma. For Sav-Gol, this is window length.
        If <= 0, no smoothing.
    smoothing_method : str
        "gaussian" (default) or "savgol".
    top_n : int or None
        If provided and group_col is not None, plot only the top N groups by total unsmoothed intensity.
    save_path : str or None
        Path to save the figure (SVG/PNG). If None, shows the plot interactively.
    figsize : tuple
        Figure size in inches (width, height).
    group_col : str or None
        Column to group chromatograms by. Use 'Fragment' (default) for fragment-ion traces,
        'Adduct' to compare adduct-level traces, or None/'all' to plot a single total trace.
    """
    # Read data
    df = pd.read_csv(ms2_csv_file)
    plt.style.use(['science', 'no-latex'])
    plt.rcParams['font.family'] = 'Arial'
    
    # Validate required columns
    required_cols = {'scan_number', 'rt', 'fragment_intensity'}
    if group_col not in (None, 'all'):
        required_cols.add(group_col)
        if group_col == 'Fragment':
            required_cols.add('fragment_mz')  # used for mean_mz labels
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")

    # Handle total-only plot (match per-adduct aggregation when possible)
    if group_col in (None, 'all'):
        total_group = 'PrecursorAdduct' if 'PrecursorAdduct' in df.columns else None
        if total_group:
            agg = (
                df.groupby(['scan_number', 'rt', total_group])
                  .agg(sum_intensity=('fragment_intensity', 'sum'))
                  .reset_index()
            )
            all_rts = agg['rt'].values
            x = _build_common_grid(all_rts)
            y = np.zeros_like(x, dtype=float)
            groups = agg[total_group].unique().tolist()
            for g in groups:
                sub = agg[agg[total_group] == g].sort_values('rt')
                xi = sub['rt'].values
                yi = sub['sum_intensity'].values
                if xi.size == 0:
                    continue
                yi = np.interp(x, xi, yi, left=0.0, right=0.0)
                if window and window > 0:
                    yi = _smooth_signal(yi, smoothing_method, window)
                y += yi
            title_suffix = "" if (window and window > 0) else " (raw)"
        else:
            total = (
                df.groupby(['scan_number', 'rt'])['fragment_intensity']
                  .sum()
                  .reset_index()
                  .sort_values('rt')
            )
            x = total['rt']
            y = total['fragment_intensity']
            if window and window > 0:
                xg, yg = _resample_uniform(x.values, y.values)
                y = _smooth_signal(yg, smoothing_method, window)
                x = xg
                title_suffix = ""
            else:
                title_suffix = " (raw)"

        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(x, y, label='Total intensity', linewidth=1.5, color='black')
        ax.fill_between(x, y, alpha=0.15, color='gray')
        ax.set_xlabel('Retention Time (RT)', fontdict={'fontsize': 13})
        ax.set_ylabel('Intensity', fontdict={'fontsize': 13})
        ax.set_title(f'MS2 Extracted Ion Chromatogram{title_suffix}', fontdict={'fontsize': 14})
        ax.legend(loc='upper right')
        ax.set_ylim(bottom=0)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=300)
            print(f"Saved plot to {save_path}")
        else:
            plt.show()
        return

    # Aggregate per scan + group (no zero-filling across missing scans)
    agg = (
        df.groupby(['scan_number', 'rt', group_col])
          .agg(sum_intensity=('fragment_intensity', 'sum'))
          .reset_index()
    )

    # Determine top N groups by total observed intensity
    totals_by_group = agg.groupby(group_col)['sum_intensity'].sum().sort_values(ascending=False)
    if top_n and top_n > 0:
        top_groups = totals_by_group.head(top_n).index.tolist()
    else:
        top_groups = totals_by_group.index.tolist()

    # Build a common RT grid across all groups to compute total consistently
    all_rts = agg['rt'].values
    total_x = _build_common_grid(all_rts)
    if window and window > 0:
        title_suffix = ""
    else:
        title_suffix = " (raw)"

    # Compute total as sum of per-group traces on the common grid
    total_y = np.zeros_like(total_x, dtype=float)
    for g in top_groups:
        sub = agg[agg[group_col] == g].sort_values('rt')
        x = sub['rt'].values
        y = sub['sum_intensity'].values
        if x.size == 0:
            continue
        y_interp = np.interp(total_x, x, y, left=0.0, right=0.0)
        if window and window > 0:
            y_interp = _smooth_signal(y_interp, smoothing_method, window)
        total_y += y_interp

    # Plotting (each group on the common RT grid so totals align)
    fig, ax = plt.subplots(figsize=figsize)
    labels = {g: f"{g}" for g in top_groups}
    for g in sorted(top_groups):
        sub = agg[agg[group_col] == g].sort_values('rt')
        x = sub['rt'].values
        y = sub['sum_intensity'].values
        if x.size == 0:
            continue
        y = np.interp(total_x, x, y, left=0.0, right=0.0)
        if window and window > 0:
            y = _smooth_signal(y, smoothing_method, window)
        x = total_x
        line, = ax.plot(
            x,
            y,
            label=labels[g],
            linewidth=1.25
        )
        ax.fill_between(
            x,
            y,
            alpha=0.2,
            color=line.get_color()
        )

    # Plot total intensity
    ax.plot(
        total_x,
        total_y,
        linestyle='--',
        label='Total Intensity',
        linewidth=1.25,
        color='black'
    )

    # Labels and styling
    ax.set_xlabel('Retention Time (RT)', fontdict={'fontsize': 13})
    ax.set_ylabel('Intensity', fontdict={'fontsize': 13})
    title_core = 'Chromatograms by ' + group_col
    ax.set_title(f'MS2 {title_core}{title_suffix}', fontdict={'fontsize': 14})
    ax.legend(title=group_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim(bottom=0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
        print(f"Saved plot to {save_path}")
    else:
        plt.show()
        plt.close(fig)


def plot_total_chromatogram_with_window(
    ms2_csv_file,
    window=11,
    save_path=None,
    figsize=(10, 5),
    smoothing_method: str = "gaussian",
    start_rt: float = None,
    end_rt: float = None
):
    """
    Plot total intensity chromatogram and optionally shade the integration window.
    Uses the same per-adduct aggregation as the total plot when PrecursorAdduct exists.
    """
    df = pd.read_csv(ms2_csv_file)
    plt.style.use(['science', 'no-latex'])
    plt.rcParams['font.family'] = 'Arial'

    required = {'scan_number', 'rt', 'fragment_intensity'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")

    if 'PrecursorAdduct' in df.columns:
        agg = (
            df.groupby(['scan_number', 'rt', 'PrecursorAdduct'])
              .agg(sum_intensity=('fragment_intensity', 'sum'))
              .reset_index()
        )
        all_rts = agg['rt'].values
        x = _build_common_grid(all_rts)
        y = np.zeros_like(x, dtype=float)
        for g in agg['PrecursorAdduct'].unique().tolist():
            sub = agg[agg['PrecursorAdduct'] == g].sort_values('rt')
            xi = sub['rt'].values
            yi = sub['sum_intensity'].values
            if xi.size == 0:
                continue
            yi = np.interp(x, xi, yi, left=0.0, right=0.0)
            if window and window > 0:
                yi = _smooth_signal(yi, smoothing_method, window)
            y += yi
    else:
        total = (
            df.groupby(['scan_number', 'rt'])['fragment_intensity']
              .sum()
              .reset_index()
              .sort_values('rt')
        )
        x = total['rt'].values
        y = total['fragment_intensity'].values
        if window and window > 0:
            xg, yg = _resample_uniform(x, y)
            y = _smooth_signal(yg, smoothing_method, window)
            x = xg

    title_suffix = "" if (window and window > 0) else " (raw)"
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(x, y, label='Total intensity', linewidth=1.5, color='black')

    if start_rt is not None and end_rt is not None:
        mask = (x >= start_rt) & (x <= end_rt)
        ax.fill_between(x[mask], y[mask], alpha=0.25, color='red', label='AUC window')

    ax.set_xlabel('Retention Time (RT)', fontdict={'fontsize': 13})
    ax.set_ylabel('Intensity', fontdict={'fontsize': 13})
    ax.set_title(f'MS2 Total Chromatogram{title_suffix}', fontdict={'fontsize': 14})
    ax.legend(loc='upper right')
    ax.set_ylim(bottom=0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
        print(f"Saved plot to {save_path}")
    else:
        plt.show()
        plt.close(fig)

