from glycanPRMQuant.performPCA import plot_pca_2d
file_path = "sample_data/Fetuin/AUCs.csv"
save_path = file_path.replace(".csv", "_PCA.svg")
pca_plot = plot_pca_2d(consolidated_csv=file_path,
                       save_path = save_path)