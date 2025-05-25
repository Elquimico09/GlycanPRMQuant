import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import scienceplots

def plot_pca_2d(consolidated_csv: str,
                save_path: str = None):
    """
    Reads a consolidated glycan AUC CSV, transposes it so that samples (files) are observations,
    performs PCA (2 components) on the samples, and displays a 2D PCA scatter plot with sample labels.

    Parameters
    ----------
    consolidated_csv : str
        Path to the consolidated CSV file with 'Glycan' and sample columns.

    Returns
    -------
    pca : sklearn.decomposition.PCA
        Fitted PCA model.
    scores_df : pandas.DataFrame
        PCA scores DataFrame with sample names as index and ['PC1','PC2'] columns.
    """
    # Set plot style
    plt.style.use(['science', 'no-latex'])
    plt.rcParams['font.family'] = 'Arial'
    # Load and prepare data
    df = pd.read_csv(consolidated_csv).fillna(0)
    # set glycan names as columns, samples become rows after transpose
    df = df.set_index('Glycan')
    df_t = df.T  # rows = samples, cols = glycans

    # Scale features (glycan AUCs)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df_t)

    # PCA with 2 components
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X_scaled)

    # Build scores DataFrame
    scores_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'], index=df_t.index)

    # Plot
    plt.figure(figsize=(5, 4))
    plt.scatter(scores_df['PC1'], scores_df['PC2'])
    for sample, (x, y) in scores_df.iterrows():
        plt.text(x, y, sample, fontsize=8, ha='right', va='bottom')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('Unsupervised Principal Component Analysis')
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.xlim(scores_df['PC1'].min() - 2, scores_df['PC1'].max() + 2)
    plt.ylim(scores_df['PC2'].min() - 0.5, scores_df['PC2'].max() + 0.5)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi = 600)
    else:
        plt.show()

    return pca, scores_df

# Example usage:
# pca_model, scores = plot_pca_2d_samples("all_glycan_auc_summary.csv")
# print(scores.head())
