"""
glycanPRMQuant: A package for glycan PRM (Parallel Reaction Monitoring) quantification.

This package provides tools for processing mass spectrometry data of glycans,
including MS1/MS2 matching, fragmentation analysis, and quantification.
"""

__version__ = "1.2.1"
__author__ = "Vishal Sandilya"
__email__ = "vishal.sandilya@ttu.edu"

__all__ = [
    "matchMS1",
    "matchMS2",
    "process_mzml_pipeline",
    "calculateAUC",
]
