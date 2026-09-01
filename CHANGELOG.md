# Changelog

## v1.3.0 - 2026-09-01

- Added target-decoy candidate scoring and reporting.
- Added MS2 noise detection and filtering.
- Added retention-time alignment and consensus feature quantification across runs.
- Added the interior-apex rule for excluding decaying signal windows.
- Improved GUI usability and automated analysis-directory generation.

## v1.2.1 - 2026-07-08

- Added CLI and GUI controls to enable or disable MS2-based isobaric precursor conflict resolution.
- Updated isobaric precursor conflict resolution to group precursor m/z values.
- Added `ms1_results_resolved.csv`, which reports the MS1 precursor assignments that survive MS2-based isobaric resolution.
- Added logging that reports whether isobaric resolution is enabled and how many precursor conflicts are resolved or preserved.
- Documented the new `--disable-isobaric-resolution` CLI flag and resolved MS1 output.
- Added focused tests for the CLI flag, 20 ppm precursor clustering, and resolved MS1 filtering.

## v1.2.0 - 2026-05-12

- Published `GlycanPRMQuant` on PyPI.
- Added modern `pyproject.toml` packaging metadata and package data handling.
- Bundled `N_glycan_db.csv` inside the Python package for installed builds.
- Removed legacy precursor and fragment database files from the root `database/` directory.
- Added the `glycan-prmquant` command-line interface with `run`, `batch`, and `gui` commands.
- Updated MS1 matching to calculate masses from grouped N-glycan compositions.
- Updated MS2 matching to dynamically fragment candidate IUPAC structures and report the best-scoring structure.
- Added GUI controls for fragment ion series and maximum cleavages.
- Improved GUI logging, including worker-process log forwarding and final runtime summary.
- Added `lxml` as a direct dependency for XML-backed mzML workflows.
- Added pytest coverage for packaging resources, MS1/MS2 matching, AUC, plotting helpers, consolidation, Skyline export, and classification utilities.
- Added README logo, badges, PyPI installation instructions, and development checks.
