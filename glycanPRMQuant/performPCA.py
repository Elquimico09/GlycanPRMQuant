import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def plot_pca_2d(consolidated_csv: str):
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
    plt.figure(figsize=(8, 6))
    plt.scatter(scores_df['PC1'], scores_df['PC2'])
    for sample, (x, y) in scores_df.iterrows():
        plt.text(x, y, sample, fontsize=8, ha='right', va='bottom')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('2D PCA of Samples (using glycan AUCs)')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return pca, scores_df

# Example usage:
# pca_model, scores = plot_pca_2d_samples("all_glycan_auc_summary.csv")
# print(scores.head())
