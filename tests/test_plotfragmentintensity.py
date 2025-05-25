import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glycanPRMQuant.plotFragmentIntensity import plot_ms2_fragments
from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum
file_path = "c:\\Users\\Vishal\\Documents\\Pompe_PRM\\Processed\\Control_11_1_10252024\\ms2_27000.csv"
svg_path = file_path.replace(".csv", ".svg")
svg_path_2 = file_path.replace(".csv", "_ms2spectrum.svg")
plot_ms2_fragments(file_path, window = 20, top_n=5, save_path=svg_path,
                   figsize=(7, 3))
avg_df = plotMS2spectrum(file_path=file_path, window_minutes=2,
                         top_n=7,
                         save_path=svg_path_2,
                         figsize=(7, 3.7))