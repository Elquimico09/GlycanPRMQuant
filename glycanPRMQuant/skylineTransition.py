import pandas as pd
from sklearn.cluster import DBSCAN

def dedupe_ms2_fragments(
    ms2_csv_file: str,
    output_csv: str,
    mz_tol: float = 0.02
) -> pd.DataFrame:
    """
    Read an MS2 CSV, DBSCAN‐cluster the fragment m/z values for each precursor,
    and write out one representative fragment m/z and its charge per cluster.

    Parameters
    ----------
    ms2_csv_file : str
        Path to the input MS2 CSV (must contain 'precursor_mz', 'fragment_mz', and 'Charge').
    output_csv : str
        Path where the deduplicated CSV will be written.
    mz_tol : float
        Mass tolerance (Da) for clustering fragment m/z values.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['precursor_mz', 'fragment_mz', 'Charge'] listing one
        fragment_mz+Charge per cluster (cluster mean m/z, representative charge).
    """
    # Load
    df = pd.read_csv(ms2_csv_file)
    required = {'precursor_mz', 'fragment_mz', 'Charge'}
    if not required.issubset(df.columns):
        raise ValueError(f"Input must contain columns: {required}")

    out_rows = []
    # Loop over each precursor
    for prec, grp in df.groupby('precursor_mz'):
        mzs = grp['fragment_mz'].values.reshape(-1, 1)
        labels = DBSCAN(eps=mz_tol, min_samples=1).fit_predict(mzs)
        grp = grp.assign(_cluster=labels)

        # For each cluster, compute mean m/z and take the first Charge
        for cid, sub in grp.groupby('_cluster'):
            mean_mz = sub['fragment_mz'].mean()
            charge = sub['Charge'].iloc[0]
            out_rows.append({
                'precursor_mz': prec,
                'fragment_mz': float(mean_mz),
                'Charge': int(charge)
            })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(output_csv, index=False)
    print(f"Wrote {len(out_df)} unique fragments to {output_csv}")
    return out_df