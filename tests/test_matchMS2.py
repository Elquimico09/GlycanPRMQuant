import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pyteomics import mzml
from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1 import matchMS1
from glycanPRMQuant.matchMS2 import matchMS2
from sklearn.cluster import DBSCAN

test_mzml_file = "sample_data/mzML/AZ_0_5ug_R1.mzML"
write_path = "sample_data/mzML/AZ_0_5ug_R1_matched.csv"
    
    # Extract MS2 data
print("Extracting MS2 data from mzML file...")
extracted_ms2_data = extractMS2(test_mzml_file, min_intensity=1e2)
print(f"Extracted {len(extracted_ms2_data)} MS2 data points")
    
    # Match MS1 data
print("Matching MS1 data...")
matched_results = matchMS1(extracted_ms2_data, ppm_tol=10, mz_min=400, mz_max=2000)
print(f"Found {len(matched_results)} matched MS1 precursors")
    
    # Match MS2 data with the updated function
print("Matching MS2 data with Gaussian fitting...")
matched_ms2 = matchMS2(
    extracted_ms2_data, 
    matched_results, 
    precursor_composition="53304", 
    mz_tol=0.05, 
    intensity_threshold=1e2, 
    ppm_tol=10
)
    
    # Save the results
matched_ms2.to_csv(write_path, index=False)
print(f"Matched MS2 data saved to {write_path}")

ms2_df = pd.read_csv(write_path)
agg = (
    ms2_df
    .groupby('scan_number')
    .agg({
        'rt': 'first',                    # or 'mean' if rt varies
        'fragment_intensity': 'sum'       # sum all fragments per scan
    })
    .reset_index()
)

window = 30
agg['smoothed_intensity'] = (
    agg['fragment_intensity']
       .rolling(window=window, center=True)
       .mean()
)

# plot raw vs smoothed
plt.figure(figsize=(8,5))
plt.plot(agg['rt'], agg['fragment_intensity'], linestyle='-', alpha=0.4, label='Raw')
plt.plot(agg['rt'], agg['smoothed_intensity'], linestyle='-', label=f'{window}-Scan MA')
plt.xlabel('Retention Time (min)')
plt.ylabel('Sum of Fragment Intensities')
plt.title('Total Fragment Ion Signal per Scan vs. RT')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
