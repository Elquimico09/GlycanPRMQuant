import pandas as pd
import numpy as np
from pyteomics import mzml
from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1 import matchMS1

test_mzml_file = "sample_data/mzML/AZ_0_5ug_R1.mzML"

extracted_ms2_data = extractMS2(test_mzml_file, min_intensity=1e2)

matched_results = matchMS1(extracted_ms2_data, ppm_tol=10, mz_min=400, mz_max=2000)
print("Matched Results:", matched_results)