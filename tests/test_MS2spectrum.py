
from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum

# Example usage:
file_path = "sample_data/results_parallel/AZ_2ug_R1/ms2_63312.csv"
svg_path = file_path.replace(".csv", "_ms2spectrum.svg")
avg_df = plotMS2spectrum(file_path=file_path, window_minutes=2,
                         save_path=svg_path,
                         figsize=(15, 4))
print(avg_df)
