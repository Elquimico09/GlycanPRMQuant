import os
import pandas as pd

def consolidate_auc_results(results_root: str, output_csv: str):
    """
    Consolidate all <mzML_basename>_auc_values.csv files under `results_root`
    into one CSV file.

    The output will have one row per glycan, a 'Glycan' column plus one column
    per mzML file (named by the folder/mzML basename) containing its AUC values.

    Parameters
    ----------
    results_root : str
        Path to the directory containing per-file subfolders (each with *_auc_values.csv).
    output_csv : str
        Path to write the consolidated CSV file (e.g. "all_auc_summary.csv").
    """
    auc_dfs = []
    for sub in sorted(os.listdir(results_root)):
        subdir = os.path.join(results_root, sub)
        if not os.path.isdir(subdir):
            continue
        auc_file = os.path.join(subdir, f"{sub}_auc_values.csv")
        if not os.path.isfile(auc_file):
            print(f"Warning: no AUC file found for {sub} (looking for {auc_file})")
            continue

        df = pd.read_csv(auc_file)
        # Identify AUC column(s) (case-insensitive contains 'auc')
        auc_cols = [c for c in df.columns if 'auc' in c.lower()]
        if not auc_cols:
            print(f"Warning: no AUC column found in {auc_file}")
            continue
        # Prefer exact 'AUC' match if present
        auc_col = 'AUC' if 'AUC' in auc_cols else auc_cols[0]
        # Rename that column to the sub (folder) name
        df = df.rename(columns={auc_col: sub})
        # Keep only Glycan and the renamed AUC column
        df = df[['Glycan', sub]]
        auc_dfs.append(df)

    if not auc_dfs:
        raise RuntimeError(f"No _auc_values.csv files found in {results_root}")

    # Merge all on 'Glycan' using outer join
    merged = auc_dfs[0]
    for df in auc_dfs[1:]:
        merged = pd.merge(merged, df, on='Glycan', how='outer')

    merged = merged.sort_values('Glycan').reset_index(drop=True)

    # Write to CSV
    merged.to_csv(output_csv, index=False)
    print(f"Wrote consolidated AUC summary to {output_csv}")
