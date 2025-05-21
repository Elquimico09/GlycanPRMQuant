import pandas as pd
import numpy as np
from glycanPRMQuant.msfileReader import extractMS1, extractMS2
from glycanPRMQuant.matchMS1 import matchMS1
from glycanPRMQuant.matchMS2 import matchMS2
from glycanPRMQuant.processmzML import process_mzml_pipeline

mzml_path = "sample_data/mzML/AZ_0_5ug_R1.mzML"
output_dir = "sample_data/mzML/AZ_0_5ug_R1_matched"
process_mzml_pipeline(
    mzml_file=mzml_path,
    output_dir=output_dir,
    ppm_ms1_tol=10,
    mz_min=400,
    mz_max=2000,
    intensity_threshold=1e2,
    ppm_ms2_tol=50,
    mz_tol=0.05
)