<p align="center">
  <img src="GPQ%20Logo.png" alt="glycanPRMQuant logo" width="560">
</p>

# glycanPRMQuant

<p align="center">
  <a href="https://pypi.org/project/glycanPRMQuant/"><img src="https://img.shields.io/pypi/v/glycanPRMQuant.svg" alt="PyPI version"></a>
  <a href="https://pypi.org/project/glycanPRMQuant/"><img src="https://img.shields.io/pypi/pyversions/glycanPRMQuant.svg" alt="Python versions"></a>
  <a href="https://github.com/Elquimico09/GlycanPRMQuant"><img src="https://img.shields.io/github/stars/Elquimico09/GlycanPRMQuant?style=social" alt="GitHub stars"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="MIT license"></a>
</p>

`glycanPRMQuant` is a Python package for targeted PRM glycomics analysis from
`.mzML` data. It extracts MS2 spectra, matches precursor ions to N-glycan
compositions, generates theoretical fragments from IUPAC structures, resolves
likely structures, plots chromatograms/spectra, and quantifies glycan signal by
AUC.

The package can be run from a Tkinter GUI for batch processing or called
programmatically from Python.

## What It Does

- Reads vendor-converted `.mzML` files with `pyteomics`.
- Matches MS1 precursor m/z values against glycan compositions.
- Calculates precursor neutral masses from the bundled `N_glycan_db.csv` using
  `glypy`, grouped once per `Composition`.
- Generates theoretical MS2 fragments from each candidate `Condensed IUPAC`
  structure for a matched numerical composition.
- Scores candidate IUPAC structures and returns the most likely structure with
  the numerical composition.
- Supports configurable fragment ion series, maximum cleavage count, m/z
  tolerances, intensity thresholds, smoothing, and AUC boundary logic.
- Produces per-glycan MS2 CSV files, chromatograms, spectra, AUC tables, and
  optional Skyline transition lists.
- Runs one file or many files in parallel.

## Repository Layout

- `glycanPRMQuant/processmzML.py`  
  Single-file end-to-end pipeline: extraction, MS1 matching, MS2 matching,
  plotting, AUC, and optional Skyline export.
- `glycanPRMQuant/parallelProcess.py`  
  Parallel multi-file runner used by the GUI and programmatic batch workflows.
- `glycanPRMQuant/pipelineGUI.py`  
  Tkinter GUI for selecting input files, output folder, matching parameters,
  plotting options, DB overrides, and batch execution.
- `glycanPRMQuant/matchMS1.py`  
  Precursor matching. Uses the N-glycan database by default and calculates
  neutral masses from grouped IUPAC compositions.
- `glycanPRMQuant/matchMS2.py`  
  Fragment matching. Generates fragments from IUPAC candidates, matches
  observed fragments, and selects the best IUPAC structure.
- `glycanPRMQuant/fragment_structure.py`  
  `glypy`-based theoretical glycan fragmentation.
- `glycanPRMQuant/calculateAUC.py`  
  Peak picking, integration windows, smoothing, and AUC summarization.
- `glycanPRMQuant/plotFragmentIntensity.py` and `plotMS2spectrum.py`  
  Chromatogram and spectrum plotting utilities.
- `glycanPRMQuant/database/N_glycan_db.csv`  
  Default structure database with `Condensed IUPAC`, `Composition`, and
  `Numerical Composition` columns.

## Installation

(Optional) Create a new python environment:

```bash
python -m venv .venv
```

(Optional) Activate the environment:

```bash
# Windows
.venv\Scripts\activate

# macOS/Linux
source .venv/bin/activate
```

Install from PyPi:

```bash
python -m pip install --upgrade pip
pip install glycanprmquant
```

Check the command-line entry point and bundled database:

```bash
glycan-prmquant --help
python -c "from glycanPRMQuant.constants import DEFAULT_PRECURSOR_DB; import os; print(os.path.exists(DEFAULT_PRECURSOR_DB), DEFAULT_PRECURSOR_DB)"
```

The package expects Python `>=3.12`.

### Development Install

For local development, clone the repository and install it in editable mode:

```bash
git clone https://github.com/Elquimico09/GlycanPRMQuant.git
cd GlycanPRMQuant
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

On Windows, activate the environment with:

```bash
.venv\Scripts\activate
```

## Dependencies

Installed from `pyproject.toml`:

- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `seaborn`
- `statsmodels`
- `scikit-learn`
- `openpyxl`
- `scienceplots`
- `pyteomics`
- `psims`
- `glypy`
- `lxml`

External requirement:

- Input data must be in `.mzML` format. Convert vendor files with ProteoWizard
  `msconvert` before running the pipeline.

## Development Checks

Install the development extra and run the tests:

```bash
pip install -e ".[dev]"
python -m pytest
python -m build
python -m twine check dist/*
```

## Quick Start: GUI

Run:

```bash
glycan-prmquant gui
```

In the GUI:

1. Select one or more `.mzML` files.
2. Select an output folder.
3. Optionally provide custom precursor/structure DB files. Leave blank to use
   the bundled `N_glycan_db.csv`.
4. Set MS1/MS2 tolerances and intensity thresholds.
5. Set fragment options:
   - `Fragment ion series`: any combination of `A`, `B`, `C`, `X`, `Y`, `Z`.
     Default: `ABCXYZ`.
   - `Max cleavages`: maximum number of cleavages used during theoretical
     fragmentation. Default: `2`.
6. Choose output options and run.

You can also launch the GUI as a module:

```bash
python -m glycanPRMQuant.pipelineGUI
```

## Quick Start: Command Line

Process one file:

```bash
glycan-prmquant run path/to/sample.mzML path/to/output_dir \
  --ppm-ms1-tol 10 \
  --ppm-ms2-tol 10 \
  --mz-tol 0.02 \
  --fragment-ion-series BY \
  --fragment-max-cleavages 2
```

Process a folder of `.mzML` files:

```bash
glycan-prmquant batch \
  --input-dir path/to/mzml_folder \
  --output-root path/to/results \
  --workers 4
```

Process specific files:

```bash
glycan-prmquant batch \
  --input-files path/to/file1.mzML path/to/file2.mzML \
  --output-root path/to/results \
  --workers 2
```

Useful CLI flags:

- `--precursor-db-path` and `--structure-db-path` override the bundled
  `N_glycan_db.csv`.
- `--skyline-transition` writes Skyline transition lists.
- `--disable-smoothing` disables chromatogram/AUC smoothing.
- `--disable-isobaric-resolution` keeps all glycan assignments whose precursor
  m/z values fall within the 20 ppm isobaric-resolution window instead of
  selecting the assignment with stronger fragment evidence.
- `--quiet` shows warnings/errors only.
- `-v` and `-vv` increase logging verbosity.

## Quick Start: Single File

```python
from glycanPRMQuant.processmzML import process_mzml_pipeline

process_mzml_pipeline(
    mzml_file="path/to/sample.mzML",
    output_dir="path/to/output_dir",
    ppm_ms1_tol=10,
    mz_min=400,
    mz_max=2000,
    intensity_threshold=1e2,
    ppm_ms2_tol=10,
    mz_tol=0.02,
    fragment_ion_series="BY",
    fragment_max_cleavages=2,
)
```

## Quick Start: Multiple Files

On Windows, keep the `if __name__ == "__main__"` guard for multiprocessing.

```python
import multiprocessing
from glycanPRMQuant.parallelProcess import run_parallel_pipeline

if __name__ == "__main__":
    multiprocessing.freeze_support()
    run_parallel_pipeline(
        input_files=[
            r"path\to\file1.mzML",
            r"path\to\file2.mzML",
        ],
        output_root=r"path\to\results",
        n_workers=4,
        ppm_ms1_tol=10,
        ppm_ms2_tol=10,
        mz_tol=0.02,
        fragment_ion_series="ABCXYZ",
        fragment_max_cleavages=2,
    )
```

## Custom Databases

By default, both MS1 and MS2 use the bundled `N_glycan_db.csv`.

You can override the database paths:

```python
process_mzml_pipeline(
    mzml_file="path/to/sample.mzML",
    output_dir="path/to/output_dir",
    precursor_db_path="path/to/N_glycan_db.csv",
    structure_db_path="path/to/N_glycan_db.csv",
)
```

The N-glycan structure database should include:

- `Condensed IUPAC`
- `Composition`
- `Numerical Composition`

`matchMS1` groups by `Composition` and calculates mass once per composition.
`matchMS2` groups by `Numerical Composition` and fragments each candidate IUPAC
structure for that composition.

## Matching Details

### MS1

`matchMS1` calculates neutral masses from the first parsable IUPAC structure for
each unique `Composition`, then generates precursor adduct m/z values:

- `2H`
- `3H`
- `4H`
- `H+NH4`
- `2NH4`

The output includes:

- `precursor_mz`
- `Glycan` using the numerical composition ID when available
- `Adduct`
- `database_mz`
- `ppm_error`

### MS2

`matchMS2` uses the matched numerical composition to find all candidate IUPAC
structures, generates theoretical fragments, and matches observed fragments by
m/z tolerance. It scores candidate structures by:

1. Total matched fragment count
2. Unique matched fragment count
3. Total matched fragment intensity
4. Mean absolute ppm error

The returned rows are restricted to the selected best-scoring IUPAC and include:

- `Glycan`
- `NumericalComposition`
- `Composition`
- `IUPAC`
- `Fragment`
- `FragmentType`
- `fragment_mz`
- `fragment_intensity`
- `Charge`
- `Adduct`
- `IUPAC_match_count`
- `IUPAC_unique_fragments`
- `IUPAC_total_intensity`

## Important Parameters

- `ppm_ms1_tol`: precursor matching tolerance in ppm.
- `mz_min`, `mz_max`: precursor m/z search range.
- `mz_offset`: offset applied to calculated precursor adduct m/z values.
- `mass_offset`: offset applied to neutral masses before precursor adduct
  calculation.
- `intensity_threshold`: minimum MS2 fragment intensity used during extraction
  and matching.
- `ppm_ms2_tol`: tolerance used to associate MS2 scans with matched precursors.
- `mz_tol`: fragment m/z tolerance in Da.
- `fragment_ion_series`: allowed theoretical fragment ion series. Use any
  combination of `A`, `B`, `C`, `X`, `Y`, `Z`.
- `fragment_max_cleavages`: maximum number of cleavages during theoretical
  fragmentation.
- `smoothing_window`: smoothing strength/window for chromatograms and AUC.
- `smoothing_method`: `gaussian` or `savgol`.
- `rel_height`: AUC boundary relative height.
- `rel_height_mode`: `prominence` or `height`.
- `skyline_transition`: write a Skyline transition list when `True`.

## Outputs

Each sample output directory can include:

- `ms1_results.csv`  
  Matched precursor assignments.
- `ms1_results_resolved.csv`  
  MS1 precursor assignments that survived MS2-based isobaric resolution.
- `ms2_<glycan>.csv`  
  Matched MS2 rows for a numerical glycan composition, including selected IUPAC
  structure information.
- `<sample>_auc_values.csv`  
  Glycan-level total AUC.
- `<sample>_auc_values_by_adduct.csv`  
  Per-adduct AUC values.
- `<sample>_skyline_transitions.xlsx`  
  Optional Skyline transition export.
- `images/*.pdf`  
  Fragment chromatograms, precursor-adduct chromatograms, total chromatograms,
  shaded AUC plots, and averaged MS2 spectra.

For multi-file runs:

- `combined_auc_values.csv` is written at the output root when more than one
  file is processed.

## Notes For Packaging

Default database paths are resolved through `glycanPRMQuant.resources`, which
supports both source-tree execution and PyInstaller-style bundled resources.
When building an executable, include `glycanPRMQuant/database/` as bundled data.

## Data Availability

Development and benchmarking data are available through MassIVE: `MSV000101208`.

The package is archived on Zenodo:
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19189798-blue)](https://doi.org/10.5281/zenodo.19189798)

## To-Do List / Work in Progress

- [ ] Finalize the visualization module allowing for gui-based visualization of quality-control data.
- [ ] Test the module on permethylated O-glycans, database already supports usage with O-glycans.
- [ ] Test the module with native glycans in negative mode, database already supports usage but needs to be validated.
- [ ] Support for Isomeric Analysis using PGC/MGC.

## Contact

For issues regarding this package, bug reports, feature requests, or questions
about usage, please contact:

- Vishal Sandilya
- Email: [vis.sandilya@gmail.com](mailto:vis.sandilya@gmail.com)
- ORCID: [0009-0006-5834-7404](https://orcid.org/0009-0006-5834-7404)

## License

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-green.svg)](https://www.python.org)
