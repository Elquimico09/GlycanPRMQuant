from glycanPRMQuant.calculateAUC import calculateAUC

csv_file = "sample_data/results_parallel/AZ_1ug_R3/ms2_53301.csv"
auc = calculateAUC(csv_file,
    glycan_col='Glycan',
    rt_col='rt',
    intensity_col='fragment_intensity',
    scan_col='scan_number',
    rel_height=0.9,
    smoothing_window=40,
    plot=True
)