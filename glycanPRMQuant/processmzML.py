import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1 import matchMS1
from glycanPRMQuant.matchMS2 import matchMS2
from glycanPRMQuant.calculateAUC import calculateAUC

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
    smoothing_window: int = 20,
    mz_offset: float = 0.0,
    mass_offset: float = 0.0,
):
    """
    Full MS1→MS2 pipeline for one .mzML file, using built-in glycan DBs,
    saves per-glycan MS2 plots into output_dir/images/, shows top 3 fragments
    + total, and computes AUC per glycan, writing to
    <mzML_basename>_auc_values.csv in output_dir.
    """
    base_name = os.path.splitext(os.path.basename(mzml_file))[0]

    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, 'images')
    os.makedirs(images_dir, exist_ok=True)

    # 1) Extract all MS2 scans
    print(f"Extracting MS2 data from {mzml_file} (min_intensity={intensity_threshold})…")
    ms2_data = extractMS2(mzml_file, min_intensity=intensity_threshold)
    print(f" → Extracted {len(ms2_data)} MS2 points.")

    # 2) Match MS1 precursors
    print(f"Matching MS1 precursors (±{ppm_ms1_tol} ppm, m/z {mz_min}-{mz_max})…")
    ms1_results = matchMS1(
        ms2_data,
        ppm_tol=ppm_ms1_tol,
        mz_min=mz_min,
        mz_max=mz_max,
        mz_offset=mz_offset,
        mass_offset=mass_offset
    )
    ms1_out = os.path.join(output_dir, "ms1_results.csv")
    ms1_results.to_csv(ms1_out, index=False)
    print(f" → Wrote {len(ms1_results)} MS1 matches to {ms1_out}")

    if ms1_results.empty:
        print("No MS1 matches found; skipping MS2 matching and AUC.")
        return

    # prepare to collect all matched MS2
    all_matched = []

    # 3) For each glycan, match MS2 and then plot+save
    glycans = {
        comp
        for comps in ms1_results['Glycan']
        for comp in str(comps).split(';')
        if comp.strip()
    }

    for glycan in sorted(glycans):
        print(f"Processing glycan {glycan!r}…")
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

        # collect for AUC
        all_matched.append(matched_ms2)

        # write per-glycan CSV
        csv_path = os.path.join(output_dir, f"ms2_{glycan}.csv")
        matched_ms2.to_csv(csv_path, index=False)
        print(f"  → Wrote {len(matched_ms2)} MS2 matches to {csv_path}")

        # aggregate by scan_number, rt, Fragment
        agg = (
            matched_ms2
            .groupby(['scan_number', 'rt', 'Fragment'])
            .agg(
                mean_mz=('fragment_mz', lambda x: round(x.mean(), 4)),
                sum_intensity=('fragment_intensity', 'sum')
            )
            .reset_index()
        )
        pivot = agg.pivot(index='rt', columns='Fragment', values='sum_intensity').fillna(0)
        if smoothing_window > 0:
            pivot_smoothed = pivot.apply(
                lambda x: gaussian_filter1d(x, sigma=smoothing_window, mode='nearest'),
                axis=0
            )
        else:
            pivot_smoothed = pivot

        # top 3 fragments
        top3 = pivot.sum(axis=0).nlargest(3).index.tolist()
        total_per_rt = pivot_smoothed.sum(axis=1)
        mean_mz_map = agg.groupby('Fragment')['mean_mz'].first().to_dict()
        labels = {frag: f"{frag} ({mean_mz_map[frag]:.4f})" for frag in top3}

        # plotting
        fig, ax = plt.subplots(figsize=(10, 6))
        for frag in top3:
            ax.plot(pivot_smoothed.index, pivot_smoothed[frag], label=labels[frag], linewidth=2)
        ax.plot(total_per_rt.index, total_per_rt.values, linestyle='--', label='Total Intensity', linewidth=2)
        ax.set_xlabel('Retention Time (RT)')
        ax.set_ylabel('Smoothed Summed Fragment Intensity')
        ax.set_title(f'Glycan {glycan}: Top 3 + Total (MA window={smoothing_window})')
        ax.legend(title='Fragment (mean m/z)', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(alpha=0.3)
        fig.tight_layout()

        img_path = os.path.join(images_dir, f"ms2_{glycan}.png")
        fig.savefig(img_path)
        plt.close(fig)
        print(f"  → Saved plot to {img_path}")

    # 4) Compute and save AUC values if any matched MS2 exist
    if all_matched:
        all_df = pd.concat(all_matched, ignore_index=True)
        print("Calculating AUC values for all glycans…")
        auc_df = calculateAUC(all_df, smoothing_window=smoothing_window)
        auc_path = os.path.join(output_dir, f"{base_name}_auc_values.csv")
        auc_df.to_csv(auc_path, index=False)
        print(f" → Wrote AUC values to {auc_path}")

    print("Pipeline complete.")
