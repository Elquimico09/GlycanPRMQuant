#!/usr/bin/env python
"""
Launch script for N-Glycan Builder GUI

This script launches the N-Glycan Builder graphical interface.
"""

import sys
import os

# Add the package directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from glycanPRMQuant.glycanBuilderGUI import launch_gui
    launch_gui()
except ImportError as e:
    print(f"Error importing GUI: {e}")
    print("\nMissing dependencies. Please install:")
    print("  pip install PyQt5 matplotlib networkx")
    sys.exit(1)
