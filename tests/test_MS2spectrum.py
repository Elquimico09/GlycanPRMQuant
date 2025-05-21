
from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum

# Example usage:
avg_df = plotMS2spectrum("sample_data/results_parallel/BZ_2ug_R2/ms2_25000.csv", top_n=5)
print(avg_df)
