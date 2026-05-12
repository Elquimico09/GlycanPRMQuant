# Changelog

## v1.2.0 - 2026-05-12

- Published `glycanPRMQuant` on PyPI.
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
