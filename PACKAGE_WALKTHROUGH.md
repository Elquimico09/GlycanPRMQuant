# GlycanPRMQuant package walkthrough

This document explains how each part of the GlycanPRMQuant Python package works, how the pieces fit together, and where you can safely modify behavior.

## High level flow

There are two major feature areas:

1) **PRM quantification pipeline**
   - Read mzML data
   - Match MS1 precursors against a glycan mass database
   - Match MS2 fragments against a fragment database
   - Plot chromatograms and spectra
   - Compute AUC values for each glycan
   - Optionally batch process multiple Thermo RAW or mzML files in parallel

2) **Glycan structure building and visualization**
   - Build glycan graphs from compositions under biosynthetic rules
   - Generate isomers for a given composition
   - Visualize structures via PyQt5 GUI or matplotlib-only fallback

## Repo layout

- `glycanPRMQuant/`  Package source
- `database/`        Precursor and fragment reference tables
- `fragments/`       Example fragment images (used in docs or demos)
- `tests/`           Unit tests for core functions
- `demo_glycan_visualization.py`  Standalone demo for glycan visualization
- `launch_glycan_builder.py`      Standalone launcher for the PyQt GUI
- `SUMMARY.md` and `GUI_*.md`      High level docs and GUI notes

## Core pipeline modules (mzML -> AUC)

### `glycanPRMQuant/msfileReader.py`

Purpose: read mzML files and extract MS1, MS2, or XIC data into DataFrames.
Thermo RAW MS2 extraction is handled separately by `thermo_raw.py` and selected
by the extension-based dispatcher in `spectra.py`.

Key functions:
- `extractMS1(mzml_file, min_intensity=1)`
  - Iterates over spectra with MS level 1.
  - Emits one row per peak above `min_intensity` with columns like `scan_number`, `rt`, `mz`, `intensity`, `tic`, `base_peak`.
- `extractXIC(mzml_file, mz, ppm_tol=10)`
  - Builds an extracted ion chromatogram by summing intensities within ppm tolerance.
  - Returns a DataFrame with `scan_number`, `rt`, `intensity`.
- `extractMS2(mzml_file, min_intensity=1)`
  - Reads MS2 spectra, captures precursor m/z and intensity, and all fragment peaks above threshold.
  - Returns a DataFrame with `scan_number`, `rt`, `precursor_mz`, `fragment_mz`, `precursor_intensity`, `fragment_intensity`.

If you want to change the mzML parsing behavior, do it here. For example, to include additional metadata, add columns to the emitted DataFrame in these functions.

### `glycanPRMQuant/matchMS1.py`

Purpose: match MS1 precursor m/z values against a theoretical glycan database with common adducts.

Key pieces:
- Constants in `glycanPRMQuant/constants.py` define proton and ammonium masses.
- `_load_precursor_db` reads `database/glycan_precursor_mz_list.xlsx` and computes adduct m/z columns.
- `matchMS1`:
  - Filters precursors to `mz_min`/`mz_max`.
  - Computes ppm tolerance per precursor.
  - Finds database rows where adduct m/z is within tolerance.
  - Returns a DataFrame with `precursor_mz`, `Glycan`, `Adduct`, `database_mz`, `ppm_error`.

If you want new adduct types, add computed columns in `_load_precursor_db` and update `adduct_cols` and `adduct_labels`.

### `glycanPRMQuant/matchMS2.py`

Purpose: match observed MS2 fragments to theoretical fragment masses for a specific glycan.

Main flow:
1) Load fragment database CSV (`database/fragment_database.csv`).
2) Filter to the target glycan composition.
3) Build a theoretical fragment list with adducts and charges.
4) Filter MS2 rows to those whose precursor matches the MS1 hits for the glycan.
5) Cluster fragments per scan using `preprocess_ms2_data` (reduces near-duplicate peaks).
6) Use a KD-tree to match observed fragment m/z to theoretical m/z within `fragment_mass_tol`.

Outputs a DataFrame with fragment assignments, including `Fragment`, `Charge`, `Fragment_mz`, `mz_diff`, and `ppm_error`.

If fragment matching is too strict or too loose, adjust `fragment_mass_tol`, the clustering tolerance in `preprocess_ms2_data`, or the adduct table.

### `glycanPRMQuant/plotFragmentIntensity.py`

Purpose: plot chromatograms from matched MS2 fragment data.

Key functions:
- `plot_ms2_fragments(...)`
  - Reads `ms2_*.csv` outputs.
  - Builds intensity traces per fragment or per adduct (or a single total trace).
  - Supports smoothing (Gaussian or Savitzky-Golay) and common RT grids.
- `plot_total_chromatogram_with_window(...)`
  - Same total chromatogram logic but can shade a chosen AUC window.

If plot layout or style changes are needed, edit here. For example, to change fonts or line widths, update the plotting code.

### `glycanPRMQuant/plotMS2spectrum.py`

Purpose: create an averaged MS2 spectrum around the most intense scan.

Algorithm:
- Find the scan with max total fragment intensity.
- Average fragments in a time window around that scan.
- Optionally keep only top N fragments.
- Plot stem spectrum and annotate fragment labels and m/z.

### `glycanPRMQuant/calculateAUC.py`

Purpose: calculate AUC per glycan and adduct from matched MS2 data.

Steps:
- Aggregate fragment intensity per scan.
- Optionally smooth the chromatogram.
- Detect the main peak using `find_peaks`.
- Compute an integration window from `rel_height`.
- Integrate using trapezoidal rule.

Returns two DataFrames:
- Per-adduct AUC details
- Per-glycan total AUC (sum across adducts)

If you need a different peak selection strategy, modify the peak detection and window logic here.

### `glycanPRMQuant/processmzML.py`

Purpose: orchestrate a full run for a single mzML file.

Pipeline steps:
1) Extract MS2 data through `spectra.extract_ms2`, which dispatches to
   Pyteomics for mzML or AlphaRaw for Thermo RAW.
2) Match precursors in MS1 via `matchMS1`.
3) For each glycan, match MS2 fragments via `matchMS2`.
4) Resolve conflicts where a precursor m/z maps to multiple glycans.
5) Write per-glycan MS2 CSV outputs.
6) Generate chromatogram plots and averaged MS2 spectra.
7) Compute AUC (per adduct and total).
8) Optionally export Skyline transition lists.

If you want to change pipeline ordering, add new outputs, or skip steps, this is the file to edit.

### `glycanPRMQuant/parallelProcess.py`

Purpose: run the pipeline on many Thermo RAW or mzML files in parallel.

Main functions:
- `_process_one_file(...)` runs `process_mzml_pipeline` on a single file and returns status.
- `run_parallel_pipeline(...)` manages a process pool, logging, and optional final AUC consolidation.

The GUI calls `run_parallel_pipeline`, and you can also call it from scripts. If you need different worker behavior (for example, chunking or retries), change it here.

### `glycanPRMQuant/consolidateAUC.py`

Purpose: merge multiple `<sample>_auc_values.csv` files into a single table with one row per glycan and one column per sample.

This is used after batch processing.

### `glycanPRMQuant/constants.py`

Purpose: shared mass constants and default database paths used by MS1/MS2 matching.

## Plotting and analysis helpers

### `glycanPRMQuant/msPlotter.py`

Utility for plotting MS data from a DataFrame with columns `mz` and `intensity`.

### `glycanPRMQuant/intensityBarplot.py`

Plots boxplots of glycan AUC distributions across samples from a consolidated AUC table.

### `glycanPRMQuant/glycantypeBarplot.py`

Computes relative abundances grouped by glycan class or type and plots mean plus SEM.

### `glycanPRMQuant/performPCA.py`

Performs PCA on consolidated AUC data (samples as rows) and plots the first two PCs.

### `glycanPRMQuant/skylineTransition.py`

Clusters fragment m/z values per precursor and writes a deduplicated transition list, useful for Skyline.

## Glycan building and visualization

### `glycanPRMQuant/glycanBuilder.py`

This file contains:
- `build_glycan(...)` for building a NetworkX graph with custom N-glycan rules.
- Visualization helpers (SNFG-like layout).
- Fragmentation functions (`fragment_glycans`, `fragment_glycans_yy`, etc.) that enumerate b/y fragments and compute masses.

This is the original, rule-based builder. It is feature rich but large. If you are changing fragmentation logic or the rule set for N-glycans, this is where most of that code lives.

### `glycanPRMQuant/glycanBuilder_improved.py`

This is the newer strict rule builder. It focuses on:
- Deterministic layering rules
- Validation of structures
- Enumerating all possible isomers for a composition

Primary functions:
- `build_nglycan_strict(...)`
- `validate_nglycan_structure(...)`
- `generate_all_isomers(...)` (with helpers for high mannose vs complex/hybrid)

The GUI uses this file, so changes here affect the interactive builder.

### `glycanPRMQuant/glycanBuilderGUI.py`

PyQt5 GUI for building and visualizing glycans. It:
- Collects counts of HexNAc, Hex, Fuc, Sia
- Builds all isomers via `generate_all_isomers`
- Displays them in a custom canvas with SNFG-style shapes
- Provides navigation between isomers and summary info

The GUI is launched via the console script `glycan-builder-gui` or `launch_glycan_builder.py`.

### `demo_glycan_visualization.py`

Standalone demo script that:
- Runs the PyQt5 GUI if available
- Falls back to a matplotlib-only visualization demo

If you want a quick smoke test of the glycan builder without the full GUI, run this script.

## Classification and class-based plots

### `glycanPRMQuant/glycanClassification.py`

Takes a consolidated AUC table and classifies each glycan based on 5-digit composition strings.
Adds `Class` (high mannose, fucosylated, etc.) and `Type` (high mannose, complex, hybrid) columns.

### `glycanPRMQuant/glycanClassificationUI.py`

Tkinter GUI with two tabs:
- Classification: runs `classifyGlycan` and writes output CSV
- Barplot: plots glycan class/type abundances using `glycantypeBarplot.py`

## Databases and static inputs

- `database/glycan_precursor_mz_list.xlsx` is the MS1 precursor database used by `matchMS1`.
- `database/fragment_database.csv` is the fragment database used by `matchMS2`.

If you add new glycans or fragments, update these files. You may also need to update adduct logic to keep matching consistent.

## Tests

`tests/` contains unit tests for:
- MS1/MS2 matching
- AUC computation and consolidation
- PCA and plotting utilities
- Fragmentation rules

These are a good reference for expected inputs/outputs and are useful for regression checks after changes.

## Where to modify specific behavior

- Change MS1/precursor matching logic: `glycanPRMQuant/matchMS1.py`
- Change MS2/fragment matching logic: `glycanPRMQuant/matchMS2.py`
- Change smoothing, peak selection, or AUC integration: `glycanPRMQuant/calculateAUC.py`
- Add new outputs to the pipeline: `glycanPRMQuant/processmzML.py`
- Change batch processing: `glycanPRMQuant/parallelProcess.py`
- Change glycan building rules: `glycanPRMQuant/glycanBuilder.py` or `glycanPRMQuant/glycanBuilder_improved.py`
- Change GUI behavior or look: `glycanPRMQuant/glycanBuilderGUI.py` and `glycanPRMQuant/pipelineGUI.py`
- Update plots: `glycanPRMQuant/plotFragmentIntensity.py`, `glycanPRMQuant/plotMS2spectrum.py`, `glycanPRMQuant/glycantypeBarplot.py`, `glycanPRMQuant/intensityBarplot.py`

## Notes on data shapes and naming

Most pipeline steps pass pandas DataFrames. The key columns are:
- MS2 extracted data: `scan_number`, `rt`, `precursor_mz`, `fragment_mz`, `fragment_intensity`
- MS1 matches: `precursor_mz`, `Glycan`, `Adduct`
- MS2 matches: `Glycan`, `PrecursorAdduct`, `Fragment`, `Charge`, `fragment_mz`, `fragment_intensity`
- AUC tables: `Glycan` and `AUC` (plus per-adduct details)

Keeping these column names consistent is important for downstream functions to work without modification.

## Setup and entry points

- `setup.py` defines dependencies and installs a console entry point:
  - `glycan-builder-gui` -> `glycanPRMQuant.glycanBuilderGUI:launch_gui`

If you add new command-line tools, you can register them in `setup.py`.
