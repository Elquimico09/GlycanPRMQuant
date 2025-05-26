
from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum

# Example usage:
file_path = "C:\\Users\\Vishal\\Documents\\Pompe_PRM\\Processed\\G_10_2_10252024\\ms2_27000.csv"
svg_path = file_path.replace(".csv", "_ms2spectrum.svg")
avg_df = plotMS2spectrum(file_path=file_path, window_minutes=2,
                         save_path=svg_path,
                         top_n=8,
                         figsize=(6.3, 4))
print(avg_df)
