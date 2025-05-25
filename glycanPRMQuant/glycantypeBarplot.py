import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

from glycanPRMQuant.glycanClassification import classifyGlycan  # adjust as needed

def plot_barplot(consolidated_csv: str,
                 figsize = (6, 4),
                                 save_path: str = None) -> pd.Series:
    # 1) classify
    df = classifyGlycan(consolidated_csv)

    # 2) get sample columns
    meta = {'Glycan', 'pos1','pos2','pos3','pos4','pos5','Class'}
    sample_cols = [c for c in df.columns if c not in meta]

    # 3) compute relative abundances per sample
    rel = df.copy()
    rel[sample_cols] = rel[sample_cols].div(rel[sample_cols].sum(axis=0), axis=1)

    # 4) sum by class and then compute mean + SEM across samples
    summed = rel.groupby('Class')[sample_cols].sum()
    means = summed.mean(axis=1)
    sems  = summed.std(axis=1) / np.sqrt(summed.shape[1])

    # 5) keep only our four classes in the desired order
    class_order = ['high mannose','sialylated','fucosylated','sialofucosylated']
    means = means.reindex(class_order).fillna(0)
    sems  = sems.reindex(class_order).fillna(0)

    # 6) plot with colored fill, black outlines, and error bars
    colors = {
      'high mannose':     'green',
      'sialofucosylated': 'blue',
      'sialylated':       'magenta',
      'fucosylated':      'red'
    }

    plt.style.use(['science','no-latex'])
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=figsize)

    bars = ax.bar(
        class_order,
        means.values,
        yerr=sems.values,
        capsize=5,                           # little caps on the error bars
        edgecolor='black',                   # bar borders
        linewidth=1.2,                       # bar border width
        color=[colors[c] for c in class_order],
        error_kw={
          'elinewidth':1.5,                   # error‐bar line width
          'ecolor':'black'                    # error‐bar color
        }
    )

    ax.set_xlabel('Glycan Class', fontsize=12)
    ax.set_ylabel('Mean Relative Abundance', fontsize=12)
    ax.set_title('Average Relative Abundance by Glycan Class', fontsize=14)
    ax.set_ylim(0, (means+sems).max() * 1.1)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()

    return means
