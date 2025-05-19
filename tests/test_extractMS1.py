from glycanPRMQuant.msfileReader import extractMS2
import pandas as pd
import numpy as np
from pyteomics import mzml
import os

ms2_data = extractMS2("sample_data/mzML/AZ_0_5ug_R1.mzML", min_intensity=1e2)
print(ms2_data.head())
# print the total number of unique precursors
print("Total number of unique precursors:", len(ms2_data['precursor_mz'].unique()))