from glycanPRMQuant.glycanClassification import classifyGlycan
from glycanPRMQuant.glycantypeBarplot import plot_barplot

auc_file = "sample_data/Fetuin/AUCs.csv"
save_path = auc_file.replace(".csv", "_glycan_type_barplot.svg")
# Classify glycans and plot the barplot
plot_barplot(auc_file, figsize=(5.5,4), save_path=save_path)