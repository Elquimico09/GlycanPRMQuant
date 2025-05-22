import pandas as pd
import numpy as np
from pyteomics import mzml
from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1 import matchMS1

test_mzml_file = "c:\\Users\\Vishal\\Documents\\Sherif_MCI\\PRM_Glycomics_AD_001.mzML"

extracted_ms2_data = extractMS2(test_mzml_file, min_intensity=1e2)

matched_results = matchMS1(extracted_ms2_data, ppm_tol=10, mz_min=400, mz_max=2000, mz_offset=0.02)
print("Matched Results:", matched_results) 
matched_results.to_csv("c:\\Users\\Vishal\\Documents\\Sherif_MCI\\matched_results.csv", index=False)