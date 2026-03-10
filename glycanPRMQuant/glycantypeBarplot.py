import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from pathlib import Path                 # Path for filesystem paths
from typing import Tuple, Optional  

from glycanPRMQuant.glycanClassification import classifyGlycan  # adjust as needed

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

def plot_barplot(consolidated_csv: str,
                 figsize=(4.8, 4),
                 save_path: str | None = None) -> pd.Series:
    # 1) classify
    df = classifyGlycan(consolidated_csv)

    # 2) identify *only* the numeric sample columns
    meta = {'Glycan',
            'pos1', 'pos2', 'pos3', 'pos4', 'pos5',
            'Class', 'Type'}                       #  ← ensure Type is excluded
    sample_cols = [c for c in df.columns if c not in meta]

    # force numeric dtype (strings, blanks → NaN → float)
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors='coerce')

    # 3) compute relative abundances per sample
    rel = df.copy()
    rel[sample_cols] = rel[sample_cols].div(rel[sample_cols].sum(axis=0), axis=1)

    # 4) sum by class, then mean ± SEM
    summed = rel.groupby('Class')[sample_cols].sum()
    means = summed.mean(axis=1)
    sems  = summed.std(axis=1) / np.sqrt(summed.shape[1])

    # 5) keep only the four classes in order
    class_order = ['high mannose', 'sialylated',
                   'fucosylated', 'sialofucosylated']
    means = means.reindex(class_order).fillna(0)
    sems  = sems.reindex(class_order).fillna(0)

    # 6) plot
    colors = {'high mannose':     'green',
              'sialofucosylated': 'blue',
              'sialylated':       'magenta',
              'fucosylated':      'red'}

    plt.style.use(['science', 'no-latex'])
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=figsize)

    ax.bar(class_order,
           means.values,
           yerr=sems.values,
           capsize=5,
           edgecolor='black',
           linewidth=1.2,
           color=[colors[c] for c in class_order],
           error_kw={'elinewidth': 1.5, 'ecolor': 'black'})

    ax.set_xlabel('Glycan Class', fontsize=12)
    ax.set_ylabel('Mean Relative Abundance', fontsize=12)
    ax.set_title('Average Relative Abundance by Glycan Class', fontsize=14)
    ax.set_ylim(0, (means + sems).max() * 1.1)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()

    return means


# --------------------------------------------------------------------------
#  helper – ensures the 'Type' values are capital-ised as the user expects
# --------------------------------------------------------------------------
_TYPE_ORDER = ["High Mannose", "Complex", "Hybrid"]
_TYPE_COLORS = {
    "High Mannose": "green",
    "Complex":      "steelblue",
    "Hybrid":       "orange",
}


def _clean_type_series(s: pd.Series) -> pd.Series:
    """
    Normalise the 'Type' strings coming out of classifyGlycan to an
    exact ('High Mannose' | 'Complex' | 'Hybrid') spelling, else NaN.
    """
    mapping = {
        "high mannose": "High Mannose",
        "high_mannose": "High Mannose",
        "highmannose":  "High Mannose",
        "complex":      "Complex",
        "hybrid":       "Hybrid",
    }
    return (
        s.str.strip()
         .str.lower()
         .map(mapping, na_action="ignore")  # unknowns → NaN
    )


# --------------------------------------------------------------------------
#  main plotting routine
# --------------------------------------------------------------------------
def plot_type_barplot(
    consolidated_csv: str | Path,
    figsize: Tuple[int, int] = (4.8, 4),
    save_path: Optional[str | Path] = None,
) -> pd.Series:
    """
    Bar-plot of *relative* glycan abundance grouped by **Type**
    (High Mannose / Complex / Hybrid).

    Parameters
    ----------
    consolidated_csv : str or Path
        Consolidated AUC table.
    figsize : tuple, default (6, 4)
        Matplotlib figure size in inches.
    save_path : str or Path, optional
        If provided, write PNG/PDF; otherwise show interactively.

    Returns
    -------
    pandas.Series
        Mean relative abundance for each Type (index ordered
        High Mannose → Complex → Hybrid).
    """
    # 1) classify & clean
    df = classifyGlycan(str(consolidated_csv))
    df["Type"] = _clean_type_series(df["Type"])
    df = df.dropna(subset=["Type"])

    meta = {
        "Glycan", "pos1", "pos2", "pos3", "pos4", "pos5",
        "Class", "Type",
    }
    sample_cols = [c for c in df.columns if c not in meta]

    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce")

    # 2) relative abundances
    rel = df.copy()
    rel[sample_cols] = rel[sample_cols].div(rel[sample_cols].sum(axis=0), axis=1)

    # 3) aggregate → mean ± SEM
    summed = rel.groupby("Type")[sample_cols].sum()
    means = summed.mean(axis=1).reindex(_TYPE_ORDER).fillna(0)
    sems  = (summed.std(axis=1) / np.sqrt(summed.shape[1])
             ).reindex(_TYPE_ORDER).fillna(0)

    # 4) plot
    plt.style.use(["science", "no-latex"])
    plt.rcParams["font.family"] = "Arial"

    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(
        _TYPE_ORDER,
        means.values,
        yerr=sems.values,
        capsize=5,
        edgecolor="black",
        linewidth=1.2,
        color=[_TYPE_COLORS[t] for t in _TYPE_ORDER],
        error_kw={"elinewidth": 1.5, "ecolor": "black"},
    )
    ax.set_xlabel("Glycan Type")
    ax.set_ylabel("Mean Relative Abundance")
    ax.set_title("Average Relative Abundance by Glycan Type")
    ax.set_ylim(0, (means + sems).max() * 1.1)
    ax.tick_params(axis="x", rotation=30)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()

    return means
