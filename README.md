# glycanPRMQuant

`glycanPRMQuant` is a Python package for targeted PRM glycomics analysis from mzML data.  
It supports glycan precursor matching, fragment matching, chromatogram generation, peak-boundary/AUC quantification, and batch processing with optional GUI workflows.

## Package Information

- Extracts MS2 data from `.mzML` files
- Matches precursor ions to glycan/adduct hypotheses (MS1)
- Matches fragment ions to glycan-specific fragments (MS2)
- Resolves ambiguous precursor-to-glycan assignments
- Generates chromatogram and spectrum figures
- Quantifies glycan signals with configurable AUC boundary logic (`rel_height`, mode)
- Processes files in parallel and consolidates AUC outputs across samples

## Main Components

- `glycanPRMQuant/processmzML.py`  
  Single-file end-to-end processing (`process_mzml_pipeline`).
- `glycanPRMQuant/parallelProcess.py`  
  Parallel multi-file processing (`run_parallel_pipeline`).
- `glycanPRMQuant/calculateAUC.py`  
  Peak boundary detection and AUC integration.
- `glycanPRMQuant/pipelineGUI.py`  
  Tkinter GUI for batch pipeline runs.
- `glycanPRMQuant/glycanBuilderGUI.py`  
  GUI for glycan structure building/visualization.

## Dependencies

### Core Python Dependencies

From `setup.py`:

- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `seaborn`
- `statsmodels`
- `scikit-learn`
- `openpyxl`
- `pyteomics`
- `networkx`

### Optional / Workflow-Specific

- `PyQt5` (for glycan builder GUI): install via `.[gui]`
- `scienceplots` (used by multiple plotting modules for figure style)

### External Requirement

- Input files should be in `.mzML` format.  
  Vendor specific files should be converted with ProteoWizard `msconvert` before running this package.

## Installation

### 1) Clone Repository

```bash
git clone https://github.com/Elquimico09/GlycanPRMQuant.git
cd GlycanPRMQuant
```

### 2) Create Environment (Recommended)

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate
```

### 3) Install Package

```bash
pip install -e .
```

### 4) Install Optional Dependencies

```bash
# Plot styling used in several figure scripts
pip install scienceplots

# Glycan builder GUI support
pip install -e .[gui]
```

## Quick Start

### Run Batch Pipeline GUI

```bash
python -m glycanPRMQuant.pipelineGUI
```

### Run Programmatically (Single File)

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
)
```

### Run Programmatically (Parallel Files)

```python
import multiprocessing
from glycanPRMQuant.parallelProcess import run_parallel_pipeline

if __name__ == "__main__":
    multiprocessing.freeze_support()  # important on Windows
    run_parallel_pipeline(
        input_files=[
            r"path\to\file1.mzML",
            r"path\to\file2.mzML",
            r"path\to\file3.mzML",
        ],
        output_root=r"path\to\results",
        n_workers=4,
    )
```

## Typical Outputs

Per sample:

- `ms1_results.csv`
- `ms2_<glycan>.csv`
- `<sample>_auc_values.csv` (glycan-level total AUC)
- `<sample>_auc_values_by_adduct.csv` (per-adduct AUC)
- `images/*.pdf` chromatograms and spectra

Multi-sample:

- `combined_auc_values.csv` (written after parallel runs with multiple files)

## Data Availability

The sample data used for development and benchmarking is available via MassIVE: [MSV000101208]
The pacakge is also archived on Zenodo: [![DOI](https://zenodo.org/badge/945763571.svg)](https://doi.org/10.5281/zenodo.19189798)

## License

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-green.svg)](https://www.python.org)
