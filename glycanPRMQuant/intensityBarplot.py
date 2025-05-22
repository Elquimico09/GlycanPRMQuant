import pandas as pd
import matplotlib.pyplot as plt

def plot_glycan_intensity_boxplot(consolidated_csv: str):
    """
    Reads the consolidated AUC CSV and plots a boxplot of glycan intensities per sample.

    Parameters
    ----------
    consolidated_csv : str
        Path to the consolidated CSV with 'Glycan' and one column per sample.
    """
    # Load data
    df = pd.read_csv(consolidated_csv)
    
    # Melt to long format: columns = ['Glycan', 'Sample', 'AUC']
    long_df = df.melt(id_vars='Glycan', var_name='Sample', value_name='AUC').fillna(0)
    
    # Prepare data: list of lists of AUCs per sample in sorted sample order
    samples = sorted(long_df['Sample'].unique())
    data = [long_df.loc[long_df['Sample'] == s, 'AUC'].values for s in samples]
    
    # Plot boxplot
    plt.figure(figsize=(10, 6))
    plt.boxplot(data, labels=samples, showfliers=False)
    plt.xlabel('Sample')
    plt.ylabel('Glycan AUC')
    plt.title('Distribution of Glycan Intensities per Sample')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
