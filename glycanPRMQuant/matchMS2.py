import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyteomics import mzml
from sklearn.cluster import DBSCAN

def preprocess_ms2_data(ms2_extracted_data, mz_tol=0.1):
    """
    Preprocess MS2 data by clustering fragment m/z values within mz_tol,
    taking the average m/z per cluster, and summing intensities.

    Parameters
    ----------
    ms2_extracted_data : pandas.DataFrame
        Must contain columns:
            ['scan_number', 'precursor_mz', 
             'fragment_mz', 'fragment_intensity']
    mz_tol : float, optional
        Mass tolerance in Da for clustering fragment m/z values.
        Default is 0.02 Da.

    Returns
    -------
    pandas.DataFrame
        Consolidated MS2 data with:
        - fragment_mz = mean cluster m/z
        - fragment_intensity = summed intensities per cluster
    """
    df = ms2_extracted_data.copy()
    result_rows = []

    # Process each spectrum separately
    for (scan, precursor), group in df.groupby(['scan_number', 'precursor_mz']):
        # cluster on the fragment_mz column
        mz_array = group['fragment_mz'].values.reshape(-1, 1)
        clustering = DBSCAN(eps=mz_tol, min_samples=1).fit(mz_array)
        group = group.assign(cluster=clustering.labels_)

        # for each cluster, compute mean m/z and summed intensity
        for cid, sub in group.groupby('cluster'):
            mean_mz = sub['fragment_mz'].mean()
            total_intensity = sub['fragment_intensity'].sum()
            # take one representative row and overwrite m/z & intensity
            row = sub.iloc[0].copy()
            row['fragment_mz'] = mean_mz
            row['fragment_intensity'] = total_intensity
            result_rows.append(row)

    result = pd.DataFrame(result_rows).reset_index(drop=True)
    print(f"Preprocessed MS2 data: reduced from "
          f"{len(ms2_extracted_data)} to {len(result)} fragment ions.")
    return result


# Update the MS2 matching function
def matchMS2(ms2_extracted_data, precursor_matched_data, precursor_composition,
              mz_tol = 0.02, intensity_threshold = 1e2, ppm_tol = 10):
    """
    Match MS2 data with glycan database.
    :param ms2_extracted_data: DataFrame containing extracted MS2 data with 'scan_number', 'rt', 'precursor_mz', 
                              'fragment_mz', 'precursor_intensity', and 'fragment_intensity' columns.
    :param precursor_matched_data: DataFrame containing matched MS1 data with 'precursor_mz' and 'Glycan' columns.
                                  The Glycan column may contain multiple compositions separated by ';'.
    :param precursor_composition: String specifying the glycan composition to filter for (e.g., "2500").
    :param ms2_database: DataFrame containing MS2 database with fragment information.
    :param mz_tol: Mass tolerance in Da for fragment matching (default: 0.02).
    :param ppm_tol: Mass tolerance in ppm for precursor matching (default: 10).
    :param intensity_threshold: Minimum intensity threshold for MS2 data (default: 1e2).
    :return: DataFrame with matched MS2 data and corresponding glycan information.
    """
    # Constants for mass calculations
    PROTON_MASS = 1.007276  # Mass of a proton in Da
    
    print(f"Starting MS2 matching with m/z tolerance: {mz_tol} Da")
    print("Loading MS2 database...")
    # Load the MS2 database
    ms2_database_path = "database/fragment_database.csv"
    ms2_database = pd.read_csv(ms2_database_path)
    ms2_database['Glycan'] = ms2_database['Glycan'].astype(str)

    # Filter MS2 database for the specific glycan composition
    ms2_database_filtered = ms2_database[ms2_database['Glycan'] == precursor_composition].copy()
    print(f"Loaded MS2 database with {len(ms2_database_filtered)} entries for {precursor_composition}")
    
    # Calculate [M+H]+ and [M+2H]2+ m/z values for each fragment
    ms2_database_filtered['Fragment_mz_1plus'] = ms2_database_filtered['Fragment Mass'] + PROTON_MASS
    ms2_database_filtered['Fragment_mz_2plus'] = (ms2_database_filtered['Fragment Mass'] + 2 * PROTON_MASS) / 2
    
    # Filter precursor_matched_data for the specific glycan composition
    # Handle semicolon-separated compositions in the Glycan column
    glycan_matches = precursor_matched_data[
        precursor_matched_data['Glycan'].apply(
            lambda x: precursor_composition in [comp.strip() for comp in str(x).split(';')]
        )
    ]
    
    if glycan_matches.empty:
        print(f"No precursor matches found for glycan composition: {precursor_composition}")
        return pd.DataFrame()
    
    # Get the precursor m/z values for the specified glycan composition
    target_precursor_mzs = glycan_matches['precursor_mz'].values
    
    # Create an empty list to store filtered MS2 data
    filtered_ms2_data = []
    
    # Filter MS2 data for each matching precursor m/z (within ppm tolerance)
    for precursor_mz in target_precursor_mzs:
        # Calculate absolute tolerance based on ppm
        abs_tol = precursor_mz * ppm_tol / 1e6
        
        # Filter MS2 data for this precursor m/z (within tolerance)
        precursor_filtered = ms2_extracted_data[
            (ms2_extracted_data['precursor_mz'] >= precursor_mz - abs_tol) & 
            (ms2_extracted_data['precursor_mz'] <= precursor_mz + abs_tol) &
            (ms2_extracted_data['fragment_intensity'] >= intensity_threshold)
        ]
        
        if not precursor_filtered.empty:
            # Add glycan information to the filtered data
            precursor_filtered = precursor_filtered.assign(Glycan=precursor_composition)
            filtered_ms2_data.append(precursor_filtered)
    
    if not filtered_ms2_data:
        print(f"No MS2 data found for glycan composition: {precursor_composition}")
        return pd.DataFrame()
    
    # Combine all filtered MS2 data
    ms2_filtered_df = pd.concat(filtered_ms2_data, ignore_index=True)
    print(f"Found {len(ms2_filtered_df)} MS2 spectra matching precursor m/z for {precursor_composition}")
    
    # Preprocess MS2 data to group similar m/z values
    ms2_filtered_df = preprocess_ms2_data(ms2_filtered_df, mz_tol=0.1)
    
    # Match fragment ions against the MS2 database
    matched_fragments = []
    
    for _, fragment_row in ms2_filtered_df.iterrows():
        fragment_mz = fragment_row['fragment_mz']
        
        # Try matching against both [M+H]+ and [M+2H]2+ fragments
        matching_db_fragments_1plus = ms2_database_filtered[
            (ms2_database_filtered['Fragment_mz_1plus'] >= fragment_mz - mz_tol) &
            (ms2_database_filtered['Fragment_mz_1plus'] <= fragment_mz + mz_tol)
        ].copy()
        
        matching_db_fragments_2plus = ms2_database_filtered[
            (ms2_database_filtered['Fragment_mz_2plus'] >= fragment_mz - mz_tol) &
            (ms2_database_filtered['Fragment_mz_2plus'] <= fragment_mz + mz_tol)
        ].copy()
        
        # Add charge state information
        if not matching_db_fragments_1plus.empty:
            matching_db_fragments_1plus['Charge'] = 1
            matching_db_fragments_1plus['Fragment_mz'] = matching_db_fragments_1plus['Fragment_mz_1plus']
            matching_db_fragments_1plus['mz_diff'] = abs(matching_db_fragments_1plus['Fragment_mz_1plus'] - fragment_mz)
        
        if not matching_db_fragments_2plus.empty:
            matching_db_fragments_2plus['Charge'] = 2
            matching_db_fragments_2plus['Fragment_mz'] = matching_db_fragments_2plus['Fragment_mz_2plus']
            matching_db_fragments_2plus['mz_diff'] = abs(matching_db_fragments_2plus['Fragment_mz_2plus'] - fragment_mz)
        
        # Combine both matches
        all_matching_fragments = pd.concat([matching_db_fragments_1plus, matching_db_fragments_2plus], ignore_index=True)
        
        if not all_matching_fragments.empty:
            # Get the closest match
            best_match = all_matching_fragments.loc[all_matching_fragments['mz_diff'].idxmin()]
            
            # Create a new row with combined information
            matched_row = fragment_row.copy()
            for col in best_match.index:
                if col not in matched_row.index:
                    matched_row[col] = best_match[col]
            
            # Calculate and add ppm error
            matched_row['ppm_error'] = (matched_row['fragment_mz'] - matched_row['Fragment_mz']) / matched_row['Fragment_mz'] * 1e6
            
            matched_fragments.append(matched_row)
    
    if not matched_fragments:
        print(f"No matching fragments found in the database for glycan composition: {precursor_composition}")
        return pd.DataFrame()
    
    # Create the final dataframe with matched MS2 data
    matched_ms2_df = pd.DataFrame(matched_fragments)
    print(f"Successfully matched {len(matched_ms2_df)} fragment ions for glycan composition: {precursor_composition}")
    
    # Summarize results by charge state
    charge_counts = matched_ms2_df['Charge'].value_counts()
    for charge, count in charge_counts.items():
        print(f"  - Matched {count} fragments with charge {charge}+")
    
    return matched_ms2_df