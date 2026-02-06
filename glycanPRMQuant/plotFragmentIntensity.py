import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
import scienceplots

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
    applies smoothing if window > 0, and plots summed intensity vs. RT
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

    # Handle total-only plot
    if group_col in (None, 'all'):
        total = (
            df.groupby(['scan_number', 'rt'])['fragment_intensity']
              .sum()
              .reset_index()
              .sort_values('rt')
        )
        x = total['rt']
        y = total['fragment_intensity']
        if window and window > 0:
            y = _smooth_signal(y, smoothing_method, window)
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

    # Aggregate per scan + group
    agg = (
        df.groupby(['scan_number', 'rt', group_col])
          .agg(
              mean_mz=('fragment_mz', lambda x: round(x.mean(), 4)) if 'fragment_mz' in df.columns else ('fragment_intensity', 'sum'),
              sum_intensity=('fragment_intensity', 'sum')
          )
          .reset_index()
    )

    # Pivot to wide format
    pivot = agg.pivot(index='rt', columns=group_col, values='sum_intensity').fillna(0)

    # Determine top N groups
    if top_n and top_n > 0:
        top_groups = pivot.sum(axis=0).nlargest(top_n).index.tolist()
    else:
        top_groups = pivot.columns.tolist()
    data = pivot[top_groups]

    # Compute total intensity
    total_intensity = data.sum(axis=1)

    # Apply Gaussian smoothing if requested
    if window and window > 0:
        smoothed = data.apply(
            lambda col: _smooth_signal(col.values, smoothing_method, window),
            axis=0
        )
        total_smoothed = _smooth_signal(total_intensity.values, smoothing_method, window)
        title_suffix = ""
    else:
        smoothed = data
        total_smoothed = total_intensity.values
        title_suffix = " (raw)"

    # Prepare labels
    labels = {g: f"{g}" for g in top_groups}

    # Plotting
    fig, ax = plt.subplots(figsize=figsize)
    for g in sorted(top_groups):
        y = smoothed[g]
        line, = ax.plot(
            smoothed.index,
            y,
            label=labels[g],
            linewidth=1.25
        )
        ax.fill_between(
            smoothed.index,
            y,
            alpha=0.2,
            color=line.get_color()
        )
    # Plot total intensity
    ax.plot(
        smoothed.index,
        total_smoothed,
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

