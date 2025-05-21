import os
import pandas as pd
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
    ppm_ms2_tol: float = 50,
    mz_tol: float = 0.05
):
    """
    Full MS1→MS2 pipeline for one .mzML file, using built-in glycan DBs.

    Parameters
    ----------
    mzml_file : str
        Path to the .mzML file to process.
    output_dir : str
        Directory where results will be written.

    # MS1 settings
    ppm_ms1_tol : float
        Tolerance (ppm) for MS1 precursor matching.
    mz_min : float
        Minimum m/z to consider in MS1.
    mz_max : float
        Maximum m/z to consider in MS1.

    # MS2 settings
    intensity_threshold : float
        Minimum fragment intensity for MS2 extraction.
    ppm_ms2_tol : float
        Tolerance (ppm) for filtering MS2 scans by precursor m/z.
    mz_tol : float
        Absolute Da tolerance for fragment‐ion matching in MS2.
    """
    os.makedirs(output_dir, exist_ok=True)

    # 1) Extract all MS2 scans
    print(f"Extracting MS2 data (min_intensity={intensity_threshold})…")
    ms2_data = extractMS2(mzml_file, min_intensity=intensity_threshold)
    print(f" → Extracted {len(ms2_data)} MS2 points.")

    # 2) Match MS1 precursors (uses ppm_ms1_tol)
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

    # 3) For each glycan, match MS2 (uses ppm_ms2_tol & mz_tol)
    glycans = {
        comp
        for comps in ms1_results['Glycan']
        for comp in str(comps).split(';')
        if comp.strip()
    }

    for glycan in sorted(glycans):
        print(f"Matching MS2 for glycan {glycan!r} (±{ppm_ms2_tol} ppm precursor, ±{mz_tol} Da fragment)…")
        matched_ms2 = matchMS2(
            ms2_data,
            ms1_results,
            precursor_composition=glycan,
            ppm_tol=ppm_ms2_tol,           # now clearly MS2 ppm tol
            mz_tol=mz_tol,
            intensity_threshold=intensity_threshold
        )
        if matched_ms2.empty:
            print(f"  → No MS2 fragments for {glycan!r}.")
            continue

        out_file = os.path.join(output_dir, f"ms2_{glycan}.csv")
        matched_ms2.to_csv(out_file, index=False)
        print(f"  → Wrote {len(matched_ms2)} MS2 matches to {out_file}")

    print("Pipeline complete.")
