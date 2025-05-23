import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1 import matchMS1
from glycanPRMQuant.matchMS2 import matchMS2
from glycanPRMQuant.calculateAUC import calculateAUC
from glycanPRMQuant.plotFragmentIntensity import plot_ms2_fragments
from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum

def process_mzml_pipeline(
    mzml_file: str,
    output_dir: str,
    # MS1 matching parameters
    ppm_ms1_tol: float = 10,
    mz_min: float = 400,
    mz_max: float = 2000,
    mz_offset: float = 0.0,
    mass_offset: float = 0.0,
    # MS2 matching parameters
    intensity_threshold: float = 1e2,
    ppm_ms2_tol: float = 10,
    mz_tol: float = 0.02,
    smoothing_window: int = 11,
    fragment_top_n: int = 10,
    spectrum_window_minutes: float = 2.0
):
    """
    Full MS1→MS2 pipeline for one .mzML file, using built-in glycan DBs.
    Saves per-glycan MS2 CSVs and plots (chromatograms + averaged spectra) into output_dir/images/,
    and computes AUC per glycan, writing to <basename>_auc_values.csv.
    """
    base_name = os.path.splitext(os.path.basename(mzml_file))[0]
    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, 'images')
    os.makedirs(images_dir, exist_ok=True)

    # 1) Extract MS2
    print(f"Extracting MS2 data from {mzml_file} (min_intensity={intensity_threshold})…")
    ms2_data = extractMS2(mzml_file, min_intensity=intensity_threshold)
    print(f" → Extracted {len(ms2_data)} MS2 points.")

    # 2) Match MS1
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

    all_matched = []

    # 3) MS2 match + plotting
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

        all_matched.append(matched_ms2)

        # save CSV
        csv_path = os.path.join(output_dir, f"ms2_{glycan}.csv")
        matched_ms2.to_csv(csv_path, index=False)
        print(f"  → Wrote {len(matched_ms2)} MS2 matches to {csv_path}")

        # plot fragment chromatograms
        chrom_svg = os.path.join(images_dir, f"ms2_{glycan}.svg")
        plot_ms2_fragments(
            ms2_csv_file=csv_path,
            window=smoothing_window,
            top_n=fragment_top_n,
            save_path=chrom_svg,
            figsize=(9, 4)
        )
        print(f"  → Saved chromatogram to {chrom_svg}")

        # plot averaged spectrum
        spec_svg = os.path.join(images_dir, f"ms2_{glycan}_ms2spectrum.svg")
        _ = plotMS2spectrum(
            file_path=csv_path,
            window_minutes=spectrum_window_minutes,
            top_n=fragment_top_n,
            save_path=spec_svg,
            figsize=(8, 4)
        )
        print(f"  → Saved averaged spectrum to {spec_svg}")

    # 4) AUC calculation
    if all_matched:
        all_df = pd.concat(all_matched, ignore_index=True)
        print("Calculating AUC values for all glycans…")
        auc_df = calculateAUC(all_df, smoothing_window=smoothing_window)
        auc_path = os.path.join(output_dir, f"{base_name}_auc_values.csv")
        auc_df.to_csv(auc_path, index=False)
        print(f" → Wrote AUC values to {auc_path}")

    print("Pipeline complete.")
