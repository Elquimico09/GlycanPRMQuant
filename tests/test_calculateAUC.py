from glycanPRMQuant.calculateAUC import calculateAUC

## Write a loop to iterate over all CSV files in the directory
## and calculate AUC for each file
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from glycanPRMQuant.calculateAUC import calculateAUC
# # Define the directory containing the CSV files
# directory = "sample_data/results_parallel/AZ_2ug_R3/"
# # Get a list of all CSV files in the directory
# csv_files = glob.glob(os.path.join(directory, "*.csv"))
# # remove any files that are not MS2 files
# csv_files = [f for f in csv_files if "ms2" in f]
# # Loop through each CSV file and calculate AUC
# for csv_file in csv_files:
#     # Define the save path for the AUC plot
#     save_path = csv_file.replace(".csv", "_AUC.svg")
#     # Calculate AUC and save the plot
#     auc = calculateAUC(csv_file,
#         glycan_col='Glycan',
#         rt_col='rt',
#         intensity_col='fragment_intensity',
#         scan_col='scan_number',
#         rel_height=0.9,
#         smoothing_window=11,
#         plot=True,
#         save_path=save_path
#     )
#     print(f"AUC for {csv_file}: {auc}")



csv_file = "sample_data/results_parallel_re/BZ_2ug_R1/ms2_29000.csv"
save_path = csv_file.replace(".csv", "_AUC.svg")
auc = calculateAUC(csv_file,
    glycan_col='Glycan',
    rt_col='rt',
    intensity_col='fragment_intensity',
    scan_col='scan_number',
    rel_height=0.9,
    smoothing_window=20,
    plot=True,
    save_path=save_path
)