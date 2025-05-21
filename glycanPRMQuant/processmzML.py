import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1 import matchMS1
from glycanPRMQuant.matchMS2 import matchMS2

def process_mzml_pipeline(
    mzml_file: str,
    output_dir: str,
    # MS1 matching parameters
    ppm_ms1_tol: float = 10,
    mz_min: float = 400,
    mz_max: float = 2000,
    # MS2 matching parameters
    intensity_threshold: float = 1e2,
    ppm_ms2_tol: float = 10,
    mz_tol: float = 0.02,
    smoothing_window: int = 20
):
    """
    Full MS1→MS2 pipeline for one .mzML file, using built-in glycan DBs,
    and saves per-glycan MS2 plots into output_dir/images/.

    Parameters
    ----------
    mzml_file : str
        Path to the .mzML file to process.
    output_dir : str
        Directory where results (CSVs + images/) will be written.
    ppm_ms1_tol : float
        Tolerance (ppm) for MS1 precursor matching.
    mz_min : float
        Minimum m/z to consider in MS1.
    mz_max : float
        Maximum m/z to consider in MS1.
    intensity_threshold : float
        Minimum fragment intensity for MS2 extraction.
    ppm_ms2_tol : float
        Tolerance (ppm) for filtering MS2 scans by precursor m/z.
    mz_tol : float
        Absolute Da tolerance for fragment-ion matching in MS2.
    smoothing_window : int
        Window size (in scans) for the moving‐average smoothing on the fragment traces.
    """
    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, 'images')
    os.makedirs(images_dir, exist_ok=True)

    # 1) Extract all MS2 scans
    print(f"Extracting MS2 data (min_intensity={intensity_threshold})…")
    ms2_data = extractMS2(mzml_file, min_intensity=intensity_threshold)
    print(f" → Extracted {len(ms2_data)} MS2 points.")

    # 2) Match MS1 precursors
    print(f"Matching MS1 precursors (±{ppm_ms1_tol} ppm, m/z {mz_min}-{mz_max})…")
    ms1_results = matchMS1(
        ms2_data,
        ppm_tol=ppm_ms1_tol,
        mz_min=mz_min,
        mz_max=mz_max
    )
    ms1_out = os.path.join(output_dir, "ms1_results.csv")
    ms1_results.to_csv(ms1_out, index=False)
    print(f" → Wrote {len(ms1_results)} MS1 matches to {ms1_out}")

    if ms1_results.empty:
        print("No MS1 matches found; skipping MS2 matching.")
        return

    # 3) For each glycan, match MS2 and then plot+save
    glycans = {
        comp
        for comps in ms1_results['Glycan']
        for comp in str(comps).split(';')
        if comp.strip()
    }

    for glycan in sorted(glycans):
        print(f"Processing glycan {glycan!r}…")
        # match MS2
        matched_ms2 = matchMS2(
            ms2_data,
            ms1_results,
            precursor_composition=glycan,
            ppm_tol=ppm_ms2_tol,
            mz_tol=mz_tol,
            intensity_threshold=intensity_threshold
        )
        if matched_ms2.empty:
            print(f"  → No MS2 fragments for {glycan!r}.")
            continue

        # write CSV
        csv_path = os.path.join(output_dir, f"ms2_{glycan}.csv")
        matched_ms2.to_csv(csv_path, index=False)
        print(f"  → Wrote {len(matched_ms2)} MS2 matches to {csv_path}")

        # now plot & save
        #  a) aggregate by scan_number, rt, Fragment
        agg = (
            matched_ms2
            .groupby(['scan_number', 'rt', 'Fragment'])
            .agg(
                mean_mz=('fragment_mz', lambda x: round(x.mean(), 4)),
                sum_intensity=('fragment_intensity', 'sum')
            )
            .reset_index()
        )
        #  b) pivot to RT × Fragment
        pivot = agg.pivot(index='rt', columns='Fragment', values='sum_intensity').fillna(0)
        #  c) smooth
        pivot_smoothed = pivot.rolling(window=smoothing_window, center=True, min_periods=1).mean()

        #  d) build labels
        mean_mz_map = agg.groupby('Fragment')['mean_mz'].first().to_dict()
        labels = {frag: f"{frag} ({mean_mz_map[frag]:.4f})" for frag in pivot_smoothed.columns}

        #  e) plot
        fig, ax = plt.subplots(figsize=(10, 6))
        for frag in sorted(pivot_smoothed.columns):
            ax.plot(
                pivot_smoothed.index,
                pivot_smoothed[frag],
                label=labels[frag],
                linewidth=2
            )
        ax.set_xlabel('Retention Time (RT)')
        ax.set_ylabel('Smoothed Summed Fragment Intensity')
        ax.set_title(f'Glycan {glycan}: Fragment Chromatograms (MA window={smoothing_window})')
        ax.legend(title='Fragment (mean m/z)', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(alpha=0.3)
        fig.tight_layout()

        #  f) save figure
        img_path = os.path.join(images_dir, f"ms2_{glycan}.png")
        fig.savefig(img_path)
        plt.close(fig)
        print(f"  → Saved plot to {img_path}")

    print("Pipeline complete.")
