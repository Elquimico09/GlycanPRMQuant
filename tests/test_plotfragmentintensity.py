import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glycanPRMQuant.plotFragmentIntensity import plot_ms2_fragments
from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum
file_path = "C:\\Users\\Vishal\\Documents\\Pompe_PRM\\Pompe\\Processed\\Pompe_1_1_10252024\\ms2_44211.csv"
svg_path = file_path.replace(".csv", ".svg")
svg_path_2 = file_path.replace(".csv", "_ms2spectrum.svg")
plot_ms2_fragments(file_path, window = 20, top_n=5, save_path=svg_path,
                   figsize=(7, 3))
avg_df = plotMS2spectrum(file_path=file_path, window_minutes=2,
                         top_n=15,
                         save_path=svg_path_2,
                         figsize=(10, 5))