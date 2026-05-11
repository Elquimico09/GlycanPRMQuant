"""
Shared constants for glycanPRMQuant package.

This module contains all physical and chemical constants used throughout
the package to ensure consistency across modules.
"""

from .resources import resource_path

# Mass constants (in Daltons)
PROTON_MASS = 1.007276
NH4_MASS = 18.033826  # Ammonium adduct mass

# Default file paths
DEFAULT_PRECURSOR_DB = resource_path("database/N_glycan_db.csv")
