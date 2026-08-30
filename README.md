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
Thermo `.raw` or `.mzML` data. It extracts MS2 spectra, matches precursor ions to N-glycan
compositions, generates theoretical fragments from IUPAC structures, resolves
likely structures, plots chromatograms/spectra, and quantifies glycan signal by
AUC.

The package can be run from a Tkinter GUI for batch processing or called
programmatically from Python.

## What It Does

- Reads Thermo `.raw` files directly with `AlphaRaw`, or `.mzML` files with
  `pyteomics`.
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
- `glycanPRMQuant/spectra.py` and `glycanPRMQuant/thermo_raw.py`
  Input-format dispatch and direct Thermo RAW extraction through AlphaRaw.
- `glycanPRMQuant/parallelProcess.py`  
  Parallel multi-file runner used by the GUI and programmatic batch workflows.
- `glycanPRMQuant/pipelineGUI.py`  
  Tkinter GUI for selecting input files, output folder, a shared glycan
  database, matching parameters, plotting options, and batch execution.
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

- `alpharaw`
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

Thermo RAW files are read directly; no mzML conversion is performed. AlphaRaw
initializes its .NET runtime only when a `.raw` file is read, so mzML-only runs
continue to use Pyteomics without loading .NET. Convert unsupported vendor
formats to `.mzML` with a tool such as ProteoWizard `msconvert`.

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

1. Select one or more Thermo `.raw` files or `.mzML` files. Do not mix formats
   within one batch.
2. Select an output folder.
3. Optionally select a custom glycan database. Leave the field blank to use the
   bundled `N_glycan_db.csv`. One custom database is used for both precursor
   and fragment matching and must contain nonblank `Condensed IUPAC`,
   `Composition`, and `Numerical Composition` columns.
4. Set MS1/MS2 tolerances and intensity thresholds.
5. Set fragment options:
   - `Fragment ion series`: any combination of `A`, `B`, `C`, `X`, `Y`, `Z`.
     Default: `ABCXYZ`.
   - `Max cleavages`: maximum number of cleavages used during theoretical
     fragmentation. Default: `2`.
6. For a replicate batch, set the cross-run peak-consensus controls. The
   default aligned RT tolerance is `±0.3` minutes and a peak group must occur in
   at least `0.8` (80%) of the input runs. All files selected together are
   treated as one comparison set.
7. Choose output options and run.

You can also launch the GUI as a module:

```bash
python -m glycanPRMQuant.pipelineGUI
```

## Quick Start: Command Line

Process one file:

```bash
glycan-prmquant run path/to/sample.raw path/to/output_dir \
  --ppm-ms1-tol 10 \
  --fragment-mass-tol 0.02 \
  --figure-filetype pdf \
  --fragment-ion-series BY \
  --fragment-max-cleavages 2
```

Process a folder of Thermo `.raw` files:

```bash
glycan-prmquant batch \
  --input-dir path/to/raw_folder \
  --output-root path/to/results \
  --workers 4
```

Process specific files:

```bash
glycan-prmquant batch \
  --input-files path/to/file1.raw path/to/file2.raw \
  --output-root path/to/results \
  --workers 2
```

Useful CLI flags:

- `--precursor-db-path` and `--structure-db-path` override the bundled
  `N_glycan_db.csv`.
- `--skyline-transition` writes Skyline transition lists.
- `--disable-smoothing` disables chromatogram/AUC smoothing.
- `--disable-isobaric-resolution` reports all glycan assignments whose precursor
  m/z values fall within the 20 ppm isobaric-resolution window. Candidate
  scores are still calculated, but contested assignments are not included in
  composition-level AUC because no winner was selected.
- `--candidate-min-fragments`, `--candidate-min-explained-intensity`,
  `--candidate-min-score`, and `--candidate-min-evidence-difference` control
  assignment acceptance and isobaric resolution.
  `--candidate-mass-outlier-min-delta` sets the minimum precursor-error
  separation used by audited pruning.
- `--fragment-mass-tol` supplies the numeric product-ion tolerance and
  `--fragment-mass-tol-unit {da,ppm}` selects its unit. The defaults remain
  `0.02 Da`; `20 ppm` is a typical high-resolution starting point, while a
  wider Da tolerance should be validated for low-resolution data.
- `--figure-filetype {png,pdf,svg}` selects the format for every chromatogram
  and MS2 spectrum generated by the pipeline. The default is `pdf`.
- Target-decoy validation is enabled by default. `--candidate-max-q-value`
  sets the acceptance threshold, `--target-decoy-seed` makes the paired decoy
  library reproducible, and `--disable-target-decoy` disables this validation.
- Multi-file runs perform cross-run RT alignment and consensus peak selection
  by default. `--consensus-rt-tolerance` controls the allowed aligned apex ΔRT
  in either direction (for example, `0.5` means `±0.5 min`),
  `--consensus-min-replicate-fraction` controls the required
  run coverage, and `--disable-consensus-peak-selection` restores legacy
  glycan-level AUC consolidation.
- `--quiet` shows warnings/errors only.
- `-v` and `-vv` increase logging verbosity.

## Quick Start: Single File

```python
from glycanPRMQuant.processmzML import process_mzml_pipeline

process_mzml_pipeline(
    mzml_file="path/to/sample.raw",
    output_dir="path/to/output_dir",
    ppm_ms1_tol=10,
    intensity_threshold=1e2,
    fragment_mass_tol=20,
    fragment_mass_tol_unit="ppm",
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
            r"path\to\file1.raw",
            r"path\to\file2.raw",
        ],
        output_root=r"path\to\results",
        n_workers=4,
        ppm_ms1_tol=10,
        fragment_mass_tol=0.02,
        fragment_mass_tol_unit="Da",
        fragment_ion_series="ABCXYZ",
        fragment_max_cleavages=2,
    )
```

## Custom Databases

By default, both MS1 and MS2 use the bundled `N_glycan_db.csv`.

You can override the database paths:

```python
process_mzml_pipeline(
    mzml_file="path/to/sample.raw",
    output_dir="path/to/output_dir",
    precursor_db_path="path/to/N_glycan_db.csv",
    structure_db_path="path/to/N_glycan_db.csv",
)
```

The N-glycan structure database should include:

- `Condensed IUPAC`
- `Composition`
- `Numerical Composition`

`Numerical Composition` is validated against every condensed IUPAC structure
as the concatenated `HexNAc`, `Hex`, `Fuc`, `Neu5Ac`, and `Neu5Gc` residue
counts, in that order.

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
m/z tolerance.

Protonated `+H` and `+2H` product ions are considered for every precursor.
For an `H+NH4` precursor assignment, `+NH4` and `+H+NH4` product ions are also
considered. For a `2NH4` precursor assignment, `+NH4` and `+2NH4` product ions
are also considered. Product-ion matching is restricted to the adduct forms
allowed by that scan's assigned precursor adduct.

For each theoretical fragment whose substructure contains NeuAc, matching also
includes a methanol neutral-loss variant (`-CH3OH`). The neutral mass loss is
32.026215 Da, so the fragment m/z shift is divided by its charge.

It scores candidate structures by:

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
- `BaseFragment`
- `FragmentType`
- `fragment_annotation`
- `fragment_iupac`
- `contains_neuac`
- `neutral_loss`
- `fragment_mz`
- `fragment_intensity`
- `Charge`
- `Adduct`
- `IUPAC_match_count`
- `IUPAC_unique_fragments`
- `IUPAC_total_intensity`

### Isobaric composition scoring

After structure-level matching, numerical-composition candidates within 20 ppm
are scored independently inside each detected retention-time feature. The
composition scorer:

1. Deduplicates repeated matches to the same observed peak within a scan.
2. Counts distinct fragment transitions instead of matched table rows.
3. Down-weights observed peaks shared by multiple composition candidates.
4. Uses the specificity-weighted explained MS2 intensity fraction.
5. Estimates a robust signed precursor-error center and sigma from the run,
   then converts each candidate's error into a Gaussian likelihood relative to
   the best precursor match in the same feature.
6. Measures chromatographic coherence against the precursor trace, a consensus
   of fragments shared by the conflict set, or the MS2 TIC fallback using
   peak-shape correlation, apex agreement, and trace overlap.
7. Uses chromatographic coherence only when it is evaluable for every candidate
   in a contested feature. If any candidate is missing it, every candidate gets
   the same neutral coelution component.
8. Requires the chromatographic apex to have at least two acquired scans on
   both its rising and falling flanks. Boundary or sparsely sampled features
   remain auditable but are not eligible for assignment or quantification.
9. Searches a paired charge/adduct-matched shifted-fragment decoy library,
   performs feature-level target-decoy competition, and estimates assignment
   FDRs and q-values.

Fragment mass error is measured in the selected tolerance unit. Its score is
still Gaussian in the median error divided by the tolerance, but its influence
decreases as the effective tolerance widens. At the feature's median observed
fragment m/z, the selected tolerance is converted to equivalent ppm and the
reliability is:

```text
R = min(1, sqrt(20 / equivalent_tolerance_ppm))
```

Thus 20 ppm or tighter receives full reliability. For example, 0.5 Da at m/z
500 is 1000 ppm and receives `R = 0.141`. The raw mass-accuracy weight is
`0.15 * R`, and all bounded-score weights are renormalized so the candidate
score still spans 0-100:

```text
W = 0.15 * R

100 / (0.85 + W) * (
    0.35 * explained-intensity component
  + 0.20 * distinct-fragment support
  + W    * fragment mass accuracy
  + 0.10 * within-feature precursor relative likelihood
  + 0.20 * feature-symmetric chromatographic coherence
)
```

The bounded score remains an absolute evidence-quality check. Candidate ranking
and runner-up separation use a separate discriminative score:

```text
S = 90 / (75 + 15 * R)

S * (35 * explained-intensity component
   + 20 * distinct-fragment support
   + 15 * R * fragment mass accuracy
   + 20 * feature-symmetric chromatographic coherence)
+  2 * within-feature precursor log-likelihood
```

This score is centered within each feature, so evidence shared equally by all
candidates cancels exactly. Resolution uses the top-minus-runner-up
discriminative evidence difference and no longer uses a ratio of total scores.

A candidate with zero candidate-specific fragments is pruned before ranking only
when its calibrated precursor error is also farther from the best candidate by
more than both 2 ppm and four calibrated sigmas. The decision is recorded in
`mass_outlier_pruned`, `mass_outlier_threshold_ppm`, and
`candidate_rejection_reason`; the row is never deleted from the audit output. If
no candidate in a contested feature has a distinguishing fragment, the feature
receives `no_discriminating_fragment_evidence` and no winner is selected.

Every assignment, including an uncontested feature, must contain an interior,
adequately sampled chromatographic apex and pass the minimum fragment,
explained-intensity, and bounded-score checks. Contested assignments must
additionally pass the discriminative-difference rule. Otherwise they receive
`no_chromatographic_peak`, `insufficient_evidence`, `ambiguous`,
`possible_coisolation`, or `no_discriminating_fragment_evidence` and remain
unselected. A sole surviving candidate after audited pruning is reported as
`resolved_after_mass_pruning`.

For target-decoy validation, each target product ion receives a reproducible
random neutral-mass shift between 1 and 30 Da; the m/z shift is divided by ion
charge. Shifts overlapping any target theoretical ion within the selected Da
or ppm fragment tolerance are rejected and regenerated. Target and decoy
libraries preserve
the same precursor assignments, fragment counts, structures, ion series,
charges, adducts, and neutral-loss counts, and are searched and structure-ranked
separately so decoys do not change target specificity weights. Decoys reuse
the target run's precursor-error calibration because each paired decoy has the
same precursor hypothesis and differs only in its shifted fragment evidence.
For a shared feature, decoys also reuse the target feature's tolerance
reliability so both sides of the competition give mass accuracy equal weight.

The selected target and selected decoy compete within each RT feature using the
bounded candidate score. At score threshold `s`, the estimated FDR is:

```text
(number of decoy-winning features at or above s + 1)
/ max(number of target-winning features at or above s, 1)
```

Monotone q-values are calculated across score thresholds. A preliminarily
selected target is removed when the decoy out-scores it or its q-value exceeds
the configured maximum. The `+1` correction is deliberately conservative; at
the default q-value of `0.05`, at least 20 target-winning features are needed
before any assignment can pass. Non-quantified fragment rows are kept in
`candidate_rows_not_quantified.csv`.

### Cross-run chromatographic peak consensus

Each selected composition/RT feature is now integrated independently before
results are combined across files. Nearby adduct-level features for the same
glycan are first joined into run-local chromatographic peaks. Reproducible,
high-scoring glycans provide RT landmarks, and each run is mapped to the batch
median retention-time scale with a robust Theil-Sen linear fit. When too few
landmarks are available, the method uses a median RT shift; with no landmarks,
it records an identity alignment rather than inventing a correction.

After alignment, peaks for the same glycan are grouped within the configured RT
tolerance, with at most one peak from each run in a group. Each group receives:

```text
replicate coverage = number of runs containing the peak / number of batch runs
consensus score    = median candidate score - candidate-score MAD
```

The median absolute deviation (MAD) penalizes a peak whose candidate score is
unstable across runs. A group is eligible only when both its run coverage and,
when target-decoy validation is available, its target-decoy pass fraction meet
the minimum replicate fraction. Eligible groups are ranked first by coverage,
then by consensus score, then by aligned-RT consistency. Only the highest
ranked group is used for that glycan in `combined_auc_values.csv`; abundance is
not part of the ranking, so a large inconsistent peak cannot win merely because
it has the largest AUC. Alternative peak groups are retained in audit tables.

## Important Parameters

- `ppm_ms1_tol`: tolerance in ppm for precursor database matching and for
  associating MS2 scans with matched precursors.
- `mz_offset`: offset applied to calculated precursor adduct m/z values.
- `mass_offset`: offset applied to neutral masses before precursor adduct
  calculation.
- `intensity_threshold`: minimum MS2 fragment intensity used during extraction
  and matching.
- `fragment_mass_tol`: numeric fragment m/z tolerance value; default `0.02`.
- `fragment_mass_tol_unit`: `Da` or `ppm`; default `Da`. The same unit/value is
  used for target matching, decoy matching, overlap prevention, near-duplicate
  peak clustering, and fragment mass-accuracy scoring.
- `fragment_ion_series`: allowed theoretical fragment ion series. Use any
  combination of `A`, `B`, `C`, `X`, `Y`, `Z`.
- `fragment_max_cleavages`: maximum number of cleavages during theoretical
  fragmentation.
- `smoothing_window`: smoothing strength/window for chromatograms and AUC.
- `smoothing_method`: `gaussian` or `savgol`.
- `rel_height`: AUC boundary relative height.
- `rel_height_mode`: `prominence` or `height`.
- `skyline_transition`: write a Skyline transition list when `True`.
- `candidate_min_fragments`: minimum distinct fragments needed to accept any
  assignment; default `2`.
- `candidate_min_explained_intensity`: minimum specificity-weighted explained
  intensity fraction; default `0.01` (1%).
- `candidate_min_score`: minimum bounded candidate score; default `35`.
- `candidate_min_evidence_difference`: minimum top-minus-runner-up difference
  in discriminative evidence; default `4`. Common-mode evidence does not affect
  this difference.
- `candidate_mass_outlier_min_delta`: absolute floor, in ppm, for pruning a
  zero-specific-evidence candidate. The effective threshold is the larger of
  this value and four run-calibrated precursor sigmas; default `2`.
- `enable_target_decoy`: generate and search the paired shifted-fragment decoy
  library; default `True`.
- `candidate_max_q_value`: maximum assignment q-value; default `0.05`.
- `target_decoy_seed`: reproducible decoy-generation seed; default `1729`.
- `minimum_peak_flank_scans`: minimum acquired scans required before and after
  a feature apex for assignment and quantification; programmatic default `2`.
- `enable_consensus_peak_selection`: align and choose one reproducible peak
  group per glycan across a multi-file batch; default `True`.
- `consensus_rt_tolerance`: allowed absolute difference between an aligned peak
  apex and its consensus-group center, in minutes. For example, `0.5` accepts
  peaks within `±0.5 min`; default `0.3`.
- `consensus_min_replicate_fraction`: minimum fraction of all batch runs that
  must contain a peak group (and pass target-decoy validation when available);
  default `0.8`.

## Outputs

Each sample output directory can include:

- `ms1_results.csv`  
  Matched precursor assignments.
- `ms1_results_resolved.csv`  
  MS1 precursor assignments that survived MS2-based isobaric resolution.
- `candidate_scores.csv`
  Candidate-level component scores, RT-feature boundaries, runner-up metrics,
  precursor calibration and pruning audit fields, coelution comparability,
  discriminative evidence, selected tolerance and equivalent ppm, mass-accuracy
  reliability/effective weights, feature scan/flank counts,
  `chromatographic_peak_valid`, `reported`/`selected` flags, quantification
  weight, and resolution status for every composition considered. Here,
  `selected=True` specifically means that the candidate is eligible for
  composition-level quantification.
- `decoy_fragment_matches.csv`
  Raw observed-fragment matches to the shifted-ion decoy library.
- `decoy_candidate_scores.csv`
  Candidate-level scores calculated independently from decoy matches.
- `target_decoy_competitions.csv`
  Complete feature-level target and decoy scores, competition winner, estimated
  FDR, q-value, and pass/fail result, including decoy-only winning features.
- `candidate_rows_not_quantified.csv`
  Matched fragment rows for losing or unresolved candidate hypotheses. These
  rows remain available for review but are excluded from glycan AUC outputs.
- `ms2_<glycan>.csv`  
  Matched MS2 rows for a numerical glycan composition, including selected IUPAC
  structure information.
- `<sample>_auc_values.csv`  
  Glycan-level total AUC.
- `<sample>_auc_values_by_adduct.csv`  
  Per-adduct AUC values.
- `<sample>_feature_auc_values.csv`
  Independent AUC and scoring/target-decoy audit fields for every selected
  chromatographic feature; this is the input to cross-run peak consensus.
- `<sample>_skyline_transitions.xlsx`  
  Optional Skyline transition export.
- `images/*.{png,pdf,svg}`
  Fragment chromatograms, precursor-adduct chromatograms, total chromatograms,
  shaded AUC plots, and averaged MS2 spectra in the selected figure format.

For multi-file runs:

- `combined_auc_values.csv`
  One consensus-selected peak group per glycan, with a separate AUC column for
  every run.
- `consensus_peak_groups.csv`
  Coverage, median score, score MAD, consensus score, aligned RT MAD, rank,
  eligibility, and selection decision for every possible peak group.
- `aligned_feature_auc_values.csv`
  Run-local peaks with raw/aligned RTs, alignment parameters, group assignment,
  source features/adducts, and all consensus metrics.
- `combined_all_feature_auc_values.csv`
  Wide AUC table containing selected and alternative peak groups.
- `retention_time_alignment.csv`
  Per-run slope, intercept, landmark count, residual MAD, and alignment method.

Cross-run consensus should be run on files that form a meaningful comparison
set (for example technical replicates or samples expected to share the same LC
method). If unrelated batches should not constrain one another, process them
separately or disable consensus selection.

## Notes For Packaging

Default database paths are resolved through `glycanPRMQuant.resources`, which
supports both source-tree execution and PyInstaller-style bundled resources.
The included `GlycanPRMQuant.spec` bundles the glycan database, AlphaRaw's
Thermo assemblies, and the Python/.NET bridge needed for direct RAW access.

On Windows, build the complete GUI application from the `defaultenv` Conda
environment with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\build_windows.ps1
```

The script produces a versioned ZIP and SHA-256 checksum under
`.packaging/release/`. Distribute the ZIP as the GitHub Release asset because
the executable depends on the adjacent files in its `onedir` bundle.

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
