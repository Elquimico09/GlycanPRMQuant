import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glycanPRMQuant.plotFragmentIntensity import plot_ms2_fragments

path_to_csv = "sample_data/mzML/AZ_0_5ug_R1_matched/ms2_43201.csv"
plot_ms2_fragments(path_to_csv, window = 10)