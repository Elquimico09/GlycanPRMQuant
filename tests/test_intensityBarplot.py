from glycanPRMQuant.intensityBarplot import plot_glycan_intensity_boxplot

# Test the function with a sample consolidated CSV file
plot_glycan_intensity_boxplot(consolidated_csv="sample_data/Fetuin/AUCs.csv",
                              save_path="sample_data/Fetuin/AUCs_boxplot.svg")