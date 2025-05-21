
from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum

# Example usage:
avg_df = plotMS2spectrum("sample_data/results_parallel/BZ_2ug_R2/ms2_26000.csv", top_n=10)
print(avg_df)
