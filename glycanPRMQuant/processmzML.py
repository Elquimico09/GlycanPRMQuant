import os
import pandas as pd

from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1     import matchMS1
from glycanPRMQuant.matchMS2     import matchMS2
from glycanPRMQuant.calculateAUC import calculateAUC
from glycanPRMQuant.plotFragmentIntensity import plot_ms2_fragments
from glycanPRMQuant.plotMS2spectrum   import plotMS2spectrum

def _normalize_glycan(val):
    if pd.isna(val):
        return ""
    s = str(val).strip()
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except ValueError:
        pass
    if s.endswith(".0"):
        return s[:-2]
    return s

def process_mzml_pipeline(
    mzml_file: str,
    output_dir: str,
    ppm_ms1_tol: float = 10,
    mz_min: float = 400,
    mz_max: float = 2000,
    mz_offset: float = 0.0,
    mass_offset: float = 0.0,
    intensity_threshold: float = 1e2,
    ppm_ms2_tol: float = 10,
    mz_tol: float = 0.02,
    smoothing_window: int = 11,
    fragment_top_n: int = 10,
    spectrum_window_minutes: float = 2.0,
    enable_adduct_plots: bool = True,
    enable_total_plots: bool = True,
    rel_height: float = 0.7,
    skyline_transition: bool = False
):
    base_name = os.path.splitext(os.path.basename(mzml_file))[0]
    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, 'images')
    os.makedirs(images_dir, exist_ok=True)

    # 1) Extract MS2
    print(f"Extracting MS2 data from {mzml_file}…")
    ms2_data = extractMS2(mzml_file, min_intensity=intensity_threshold)
    print(f" → Extracted {len(ms2_data)} points")

    # 2) Match MS1
    print("Matching MS1 precursors…")
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
    print(f" → Wrote {len(ms1_results)} MS1 matches")

    if ms1_results.empty:
        print("No MS1 matches; done.")
        return

    all_matched = []

    # build clean glycan list
    glycans = set()
    for entry in ms1_results['Glycan'].dropna().astype(str):
        for comp in entry.split(';'):
            c = _normalize_glycan(comp)
            if c and c.lower() != 'nan':
                glycans.add(c)

    # 3) For each glycan: match MS2, save CSV, plot
    transitions = []
    adduct_charge = {'2H':2,'3H':3,'4H':4,'H+NH4':2,'2NH4':2}

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
            print(f"  → No MS2 fragments for {glycan!r}")
            continue

        # **Normalize** Fragment_mz → fragment_mz so downstream code sees it
        if 'Fragment_mz' in matched_ms2.columns:
            matched_ms2 = matched_ms2.rename(columns={'Fragment_mz': 'fragment_mz'})

        all_matched.append(matched_ms2)

        # save CSV
        csv_path = os.path.join(output_dir, f"ms2_{glycan}.csv")
        matched_ms2.to_csv(csv_path, index=False)
        print(f"  → Wrote {len(matched_ms2)} MS2 matches to {csv_path}")

        # Collect transitions for Skyline (one per matched fragment row)
        if skyline_transition:
            charge = matched_ms2['PrecursorAdduct'].map(adduct_charge).fillna(0)
            trans = pd.DataFrame({
                'Glycan': matched_ms2['Glycan'],
                'Adduct': matched_ms2['PrecursorAdduct'],
                'Precursor m/z': matched_ms2['precursor_mz'],
                'Precursor charge': charge.astype(int),
                'Fragment m/z': matched_ms2['fragment_mz'],
                'Fragment charge': matched_ms2['Charge'],
                'Retention time': matched_ms2['rt']
            })
            transitions.append(trans)

        # chromatogram
        # Chromatograms per fragment (legacy view)
        chrom_frag_svg = os.path.join(images_dir, f"ms2_{glycan}.svg")
        try:
            plot_ms2_fragments(csv_path, window=smoothing_window,
                               top_n=fragment_top_n, save_path=chrom_frag_svg,
                               group_col='Fragment')
            print(f"  → Saved fragment-level chromatogram to {chrom_frag_svg}")
        except Exception as e:
            print(f"  [warn] Fragment chromatogram failed: {e}")

        # Chromatogram per precursor adduct
        if enable_adduct_plots:
            chrom_adduct_svg = os.path.join(images_dir, f"ms2_{glycan}_by_precursor_adduct.svg")
            try:
                plot_ms2_fragments(csv_path, window=smoothing_window,
                                   top_n=None, save_path=chrom_adduct_svg,
                                   group_col='PrecursorAdduct')
                print(f"  → Saved precursor-adduct chromatogram to {chrom_adduct_svg}")
            except Exception as e:
                print(f"  [warn] Precursor-adduct chromatogram failed: {e}")

        # Total chromatogram (all fragments/adducts summed)
        if enable_total_plots:
            chrom_total_svg = os.path.join(images_dir, f"ms2_{glycan}_total.svg")
            try:
                plot_ms2_fragments(csv_path, window=smoothing_window,
                                   top_n=None, save_path=chrom_total_svg,
                                   group_col=None)
                print(f"  → Saved total chromatogram to {chrom_total_svg}")
            except Exception as e:
                print(f"  [warn] Total chromatogram failed: {e}")

        # averaged spectrum
        spec_svg = os.path.join(images_dir, f"ms2_{glycan}_ms2spectrum.svg")
        try:
            plotMS2spectrum(csv_path, window_minutes=spectrum_window_minutes,
                            top_n=fragment_top_n, save_path=spec_svg)
            print(f"  → Saved spectrum to {spec_svg}")
        except Exception as e:
            print(f"  [warn] Spectrum plot failed: {e}")

    # 4) AUC
    if all_matched:
        all_df = pd.concat(all_matched, ignore_index=True)
        print("Calculating AUC…")
        per_adduct_df, total_df = calculateAUC(
            all_df,
            smoothing_window=smoothing_window,
            adduct_col='PrecursorAdduct',
            rel_height=rel_height
        )

        # Total AUC (backward compatible filename)
        auc_path = os.path.join(output_dir, f"{base_name}_auc_values.csv")
        total_df.to_csv(auc_path, index=False)
        print(f" → Wrote total AUC (summed across adducts) to {auc_path}")

        # Per-adduct detail
        auc_adduct_path = os.path.join(output_dir, f"{base_name}_auc_values_by_adduct.csv")
        per_adduct_df.to_csv(auc_adduct_path, index=False)
        print(f" → Wrote per-adduct AUC values to {auc_adduct_path}")

    # 5) Skyline transition export
    if skyline_transition and transitions:
        trans_df = pd.concat(transitions, ignore_index=True)
        trans_path = os.path.join(output_dir, f"{base_name}_skyline_transitions.xlsx")
        trans_df.to_excel(trans_path, index=False)
        print(f" → Wrote Skyline transition list to {trans_path}")

    print("Done.")
