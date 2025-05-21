import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_ms2_fragments(ms2_csv_file, window=20):
    """
    Reads an MS2 results CSV file, aggregates fragment data per scan grouped by Fragment,
    applies moving average smoothing, and plots summed intensity vs. RT
    with one smoothed line per fragment type.

    Parameters
    ----------
    ms2_csv_file : str
        Path to the MS2 CSV file.
    window : int
        Window size (in number of scans) for the moving average.
    """
    # Read the CSV file
    df = pd.read_csv(ms2_csv_file)
    
    # Check required columns
    required_cols = {'scan_number', 'rt', 'Fragment', 'fragment_mz', 'fragment_intensity'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")
    
    # Aggregate: mean m/z (rounded to 4 decimals) and sum intensity per scan + fragment
    agg = (
        df
        .groupby(['scan_number', 'rt', 'Fragment'])
        .agg(
            mean_mz=('fragment_mz', lambda x: round(x.mean(), 4)),
            sum_intensity=('fragment_intensity', 'sum')
        )
        .reset_index()
    )
    
    # Pivot: index=rt, columns=Fragment, values=sum_intensity
    pivot = agg.pivot(index='rt', columns='Fragment', values='sum_intensity').fillna(0)
    
    # Apply moving average smoothing along the RT axis
    pivot_smoothed = pivot.rolling(window=window, center=True, min_periods=1).mean()
    
    # Prepare labels including mean m/z for each fragment
    mean_mz_dict = agg.groupby('Fragment')['mean_mz'].first().to_dict()
    labels = {frag: f"{frag} ({mean_mz_dict[frag]:.4f})" for frag in pivot_smoothed.columns}
    
    # Plotting
    plt.figure(figsize=(10, 6))
    for frag in sorted(pivot_smoothed.columns):
        plt.plot(
            pivot_smoothed.index,
            pivot_smoothed[frag],
            label=labels[frag],
            linewidth=2
        )
    plt.xlabel('Retention Time (RT)')
    plt.ylabel('Smoothed Summed Fragment Intensity')
    plt.title(f'MS2 Fragment Ion Chromatograms (MA window={window})')
    plt.legend(title='Fragment (mean m/z)', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
