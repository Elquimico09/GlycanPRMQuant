import os
import pandas as pd
import matplotlib.pyplot as plt

import os
import pandas as pd
import matplotlib.pyplot as plt

def plotMS2spectrum(file_path, top_n=10):
    """
    Reads an MS2 results file (CSV or Excel), accounts for multiple charges per fragment,
    computes the average intensity per fragment+charge across all scans, and plots
    a single averaged MS2 spectrum showing only the top N fragments by average intensity,
    with fragment type + m/z and charge annotations above each peak.

    Returns a DataFrame with columns ['Fragment', 'Charge', 'mz', 'avg_intensity'] for the top N.
    """
    # Load data
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.csv':
        df = pd.read_csv(file_path)
    elif ext in ('.xlsx', '.xls'):
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file type. Provide a .csv or Excel file.")

    # Validate required columns
    required = {'scan_number', 'Fragment', 'Charge', 'fragment_mz', 'fragment_intensity'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Sum intensities per scan + fragment + charge, and compute mean m/z
    scan_frag_charge = (
        df
        .groupby(['scan_number', 'Fragment', 'Charge'])
        .agg(
            sum_intensity=('fragment_intensity', 'sum'),
            mean_mz=('fragment_mz', 'mean')
        )
        .reset_index()
    )

    # Average across all scans per fragment+charge
    frag_charge_avg = (
        scan_frag_charge
        .groupby(['Fragment', 'Charge'])
        .agg(
            mz=('mean_mz', 'mean'),
            avg_intensity=('sum_intensity', 'mean')
        )
        .reset_index()
        .sort_values('avg_intensity', ascending=False)
    )

    # Select top N fragments by average intensity
    top_fragments = frag_charge_avg.head(top_n)

    # Plot averaged spectrum as a stick plot with annotations
    fig, ax = plt.subplots(figsize=(14, 6))
    markerline, stemlines, baseline = ax.stem(
        top_fragments['mz'],
        top_fragments['avg_intensity']
    )
    baseline.set_visible(False)
    
    # Annotate each peak with fragment type on first line, then m/z+charge
    for frag, x, y, charge in zip(
        top_fragments['Fragment'],
        top_fragments['mz'],
        top_fragments['avg_intensity'],
        top_fragments['Charge']
    ):
        ax.text(
            x, y * 1.02,  # slightly above the peak
            f"{frag}\n{x:.4f}+{int(charge)}",
            ha='center', va='bottom',
            rotation=0,
            fontsize=10
        )

    ax.set_xlabel('m/z')
    ax.set_ylabel('Average Intensity')
    ax.set_ylim(0, top_fragments['avg_intensity'].max() * 1.1)
    ax.set_xlim(top_fragments['mz'].min() - 50, top_fragments['mz'].max() * 1.05)
    ax.set_title(f'Averaged MS2 Spectrum (Top {top_n} Fragment+Charge)')
    plt.tight_layout()
    plt.show()

    return top_fragments
