import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots

def plotMS2spectrum(file_path, window_minutes=2, top_n=20,
                    save_path=None, figsize=(15, 5)):
    """
    Reads an MS2 results file (CSV or Excel), identifies the scan with the maximum total fragment intensity,
    averages fragment intensities within a time window around that scan, and plots a single averaged MS2 spectrum
    with fragment type and m/z+charge annotations.

    Parameters
    ----------
    file_path : str
        Path to the MS2 results file (.csv, .xlsx, .xls).
    window_minutes : float
        Total window length in minutes around the peak scan to include in the average (default 2 minutes).
    top_n : int or None
        If provided, display only the top N fragments by average intensity.

    Returns
    -------
    pandas.DataFrame
        Columns ['Fragment', 'Charge', 'mz', 'avg_intensity'] for fragments averaged within the window.
    """
    # Load data
    plt.style.use(['science', 'no-latex'])
    plt.rcParams['font.family'] = 'Arial'
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.csv':
        df = pd.read_csv(file_path)
    elif ext in ('.xlsx', '.xls'):
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file type. Provide a .csv or Excel file.")

    # Validate required columns
    required = {'scan_number', 'rt', 'Fragment', 'Charge', 'fragment_mz', 'fragment_intensity'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Compute total intensity per scan
    scan_intensity = (
        df.groupby('scan_number')
          .agg(rt=('rt', 'first'),
               total_intensity=('fragment_intensity', 'sum'))
          .reset_index()
    )
    # Find the peak scan
    peak = scan_intensity.loc[scan_intensity['total_intensity'].idxmax()]
    rt_peak = peak['rt']
    half_window = window_minutes / 2
    rt_min, rt_max = rt_peak - half_window, rt_peak + half_window

    print(f"Peak scan at RT={rt_peak:.2f} min (total intensity={peak['total_intensity']:.1f}); "
          f"averaging ±{half_window:.2f} min → RT window [{rt_min:.2f}, {rt_max:.2f}]")

    # Filter to window
    df_win = df[(df['rt'] >= rt_min) & (df['rt'] <= rt_max)]
    print(f"  → {df_win['scan_number'].nunique()} scans selected out of {df['scan_number'].nunique()}")

    # Aggregate per scan + fragment + charge
    scan_frag_charge = (
        df_win
        .groupby(['scan_number', 'Fragment', 'Charge'])
        .agg(
            sum_intensity=('fragment_intensity', 'sum'),
            mean_mz=('fragment_mz', 'mean')
        )
        .reset_index()
    )

    # Average across selected scans per fragment+charge
    frag_charge_avg = (
        scan_frag_charge
        .groupby(['Fragment', 'Charge'])
        .agg(
            mz=('mean_mz', 'mean'),
            avg_intensity=('sum_intensity', 'mean')
        )
        .reset_index()
    )

    # Select top N fragments if requested
    if top_n is not None and top_n < len(frag_charge_avg):
        total_frag = len(frag_charge_avg)
        frag_charge_avg = frag_charge_avg.nlargest(top_n, 'avg_intensity')
        print(f"Displaying top {top_n} fragments out of {total_frag} by average intensity.")
    # Sort by m/z for consistent plotting order
    frag_charge_avg = frag_charge_avg.sort_values('mz').reset_index(drop=True)
    max_intensity = frag_charge_avg['avg_intensity'].max()
    if max_intensity > 0:
        frag_charge_avg['relative_intensity'] = frag_charge_avg['avg_intensity'] / max_intensity
    else:
        frag_charge_avg['relative_intensity'] = 0.0

    # Plot averaged spectrum
    fig, ax = plt.subplots(figsize=figsize)
    markerline, stemlines, baseline = ax.stem(
        frag_charge_avg['mz'],
        frag_charge_avg['relative_intensity']
    )
    baseline.set_visible(False)

    # Annotate peaks with fragment type + m/z+charge
    for frag, x, y, charge in zip(
        frag_charge_avg['Fragment'],
        frag_charge_avg['mz'],
        frag_charge_avg['relative_intensity'],
        frag_charge_avg['Charge']
    ):
        ax.text(
            x, y * 1.05,
            f"{frag}\n{x:.4f}+{int(charge)}",
            ha='center', va='bottom',
            rotation=90, fontsize=8
        )

    ax.set_xlabel('m/z')
    ax.set_ylabel('Relative Intensity')
    # Add max intensity text to top left of the plot
    ax.text(
        0.01, 0.95,
        f"Base Peak Intensity: {max_intensity:.0f}",
        transform=ax.transAxes,
        ha='left', va='top',
        fontsize=12,
        fontdict={'weight': 'bold', 'color': 'black'}
    )
    ax.set_ylim(0, 1.4)
    ax.set_xlim(frag_charge_avg['mz'].min() - 100, frag_charge_avg['mz'].max() * 1.05)
    title = f"Averaged MS2 Spectrum"
    if top_n is not None:
        title += f" (Top {top_n})"
    ax.set_title(title)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300)
        print(f"Saved plot to {save_path}")
        plt.close(fig)
    else:
        plt.show()
        plt.close(fig)

    return frag_charge_avg
