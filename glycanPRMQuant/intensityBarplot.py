import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

def plot_glycan_intensity_boxplot(consolidated_csv: str,
                                  save_path: str = None) -> None:
    """
    Reads the consolidated AUC CSV and plots a boxplot of glycan intensities per sample,
    preserving the original sample column order.
    """
    plt.style.use(['science', 'no-latex'])
    plt.rcParams['font.family'] = 'Arial'
    
    # Load data
    df = pd.read_csv(consolidated_csv)
    
    # Identify samples in original order (all columns except 'Glycan')
    samples = [col for col in df.columns if col != 'Glycan']
    
    # Melt to long format
    long_df = df.melt(
        id_vars='Glycan',
        value_vars=samples,
        var_name='Sample',
        value_name='AUC'
    ).fillna(0)
    
    # Prepare data: one list of AUCs per sample, in original order
    data = [long_df.loc[long_df['Sample'] == s, 'AUC'].values for s in samples]
    
    # Plot
    plt.figure(figsize=(4.8, 4))
    plt.boxplot(data, labels=samples, showfliers=False)
    plt.ylabel('Glycan Abundance', fontsize=14)
    plt.title('Distribution of Glycan Abundance per Sample')
    plt.xticks(rotation=45, ha='right', fontsize = 12)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Saved boxplot to {save_path}")
    else:
        plt.show()
