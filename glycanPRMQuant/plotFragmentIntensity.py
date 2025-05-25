import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import scienceplots

def plot_ms2_fragments(ms2_csv_file, window=11, top_n=None, save_path=None,
                       figsize=(10, 5)):
    """
    Reads an MS2 results CSV file, aggregates fragment data per scan grouped by Fragment,
    applies Gaussian smoothing if window > 0, and plots summed intensity vs. RT
    with one smoothed line per fragment type plus the total intensity line, shading under each curve.

    Parameters
    ----------
    ms2_csv_file : str
        Path to the MS2 CSV file.
    window : int
        Sigma for Gaussian smoothing (in number of scans). If <= 0, no smoothing.
    top_n : int or None
        If provided, plot only the top N fragments by total unsmoothed intensity.
    save_path : str or None
        Path to save the figure (SVG/PNG). If None, shows the plot interactively.
    figsize : tuple
        Figure size in inches (width, height).
    """
    # Read data
    df = pd.read_csv(ms2_csv_file)
    plt.style.use(['science', 'no-latex'])
    plt.rcParams['font.family'] = 'Arial'
    
    # Validate required columns
    required_cols = {'scan_number', 'rt', 'Fragment', 'fragment_mz', 'fragment_intensity'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")
    
    # Aggregate per scan + fragment
    agg = (
        df.groupby(['scan_number', 'rt', 'Fragment'])
          .agg(
              mean_mz=('fragment_mz', lambda x: round(x.mean(), 4)),
              sum_intensity=('fragment_intensity', 'sum')
          )
          .reset_index()
    )
    
    # Pivot to wide format
    pivot = agg.pivot(index='rt', columns='Fragment', values='sum_intensity').fillna(0)
    
    # Determine top N fragments
    if top_n and top_n > 0:
        top_frags = pivot.sum(axis=0).nlargest(top_n).index.tolist()
    else:
        top_frags = pivot.columns.tolist()
    data = pivot[top_frags]
    
    # Compute total intensity
    total_intensity = data.sum(axis=1)
    
    # Apply Gaussian smoothing if requested
    if window and window > 0:
        smoothed = data.apply(
            lambda col: gaussian_filter1d(col, sigma=window, mode='nearest'),
            axis=0
        )
        total_smoothed = gaussian_filter1d(total_intensity.values, sigma=window, mode='nearest')
        title_suffix = ""
    else:
        smoothed = data
        total_smoothed = total_intensity.values
        title_suffix = " (raw)"
    
    # Prepare labels
    mean_mz_dict = agg.groupby('Fragment')['mean_mz'].first().to_dict()
    labels = {frag: f"{frag}" for frag in top_frags}
    
    # Plotting
    fig, ax = plt.subplots(figsize=figsize)
    for frag in sorted(top_frags):
        y = smoothed[frag]
        line, = ax.plot(
            smoothed.index,
            y,
            label=labels[frag],
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
    ax.set_title(f'MS2 Fragment Ion Chromatograms{title_suffix}', fontdict={'fontsize': 14})
    ax.legend(title='Fragment (mean m/z)', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim(bottom=0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300)
        print(f"Saved plot to {save_path}")
    else:
        plt.show()

