import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pyteomics import mzml
from glycanPRMQuant.msfileReader import extractMS2
from glycanPRMQuant.matchMS1 import matchMS1

# First, let's define the Gaussian fit function that we created
def gaussian_fit_mz_values(mz_values, intensities=None, plot=False):
    """
    Apply a Gaussian fit to a set of m/z values that should represent the same ion.
    
    Parameters:
    -----------
    mz_values : array-like
        A list or array of m/z values for the same ion.
    intensities : array-like, optional
        Corresponding intensity values. If None, all points are weighted equally.
    plot : bool, optional
        If True, create a plot of the data and the Gaussian fit.
        
    Returns:
    --------
    center : float
        The center of the Gaussian fit (mean).
    std_dev : float
        The standard deviation of the Gaussian fit.
    """
    # Convert to numpy arrays if needed
    mz_values = np.array(mz_values)
    
    # If no intensities are provided, assume equal weights
    if intensities is None:
        weights = None
        intensities = np.ones_like(mz_values)
    else:
        intensities = np.array(intensities)
        # Normalize weights to sum to 1
        weights = intensities / np.sum(intensities)
    
    # Define Gaussian function
    def gaussian(x, mean, sigma, amplitude):
        return amplitude * np.exp(-(x - mean)**2 / (2 * sigma**2))
    
    # Initial guesses for parameters
    # Use weighted statistics if weights are provided
    if weights is not None:
        mean_guess = np.average(mz_values, weights=weights)
        # Calculate weighted variance
        variance = np.average((mz_values - mean_guess)**2, weights=weights)
        sigma_guess = np.sqrt(variance)
    else:
        mean_guess = np.mean(mz_values)
        sigma_guess = np.std(mz_values)
    
    amplitude_guess = np.max(intensities)
    
    # Perform the curve fit
    try:
        popt, pcov = curve_fit(gaussian, mz_values, intensities, 
                               p0=[mean_guess, sigma_guess, amplitude_guess],
                               sigma=None if weights is None else 1/np.sqrt(weights+1e-10),
                               absolute_sigma=True)
        
        center, std_dev, amplitude = popt
        
        # Create a plot if requested
        if plot:
            plt.figure(figsize=(10, 6))
            
            # Plot the original data points
            plt.scatter(mz_values, intensities, color='blue', label='Observed m/z values')
            
            # Plot the Gaussian fit
            x_fit = np.linspace(min(mz_values), max(mz_values), 1000)
            y_fit = gaussian(x_fit, center, std_dev, amplitude)
            plt.plot(x_fit, y_fit, 'r-', label=f'Gaussian fit (μ={center:.6f}, σ={std_dev:.6f})')
            
            # Add vertical line at the center
            plt.axvline(x=center, color='green', linestyle='--', 
                       label=f'Center: {center:.6f}')
            
            plt.xlabel('m/z')
            plt.ylabel('Intensity')
            plt.title('Gaussian Fit of m/z Values')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()
        
        return center, std_dev
    
    except RuntimeError:
        print("Error: Curve fitting failed. Using median and MAD instead.")
        center = np.median(mz_values)
        # Calculate Median Absolute Deviation (MAD)
        mad = np.median(np.abs(mz_values - center))
        # Convert MAD to standard deviation (assuming normal distribution)
        std_dev = mad * 1.4826
        
        # Create a plot if requested
        if plot:
            plt.figure(figsize=(10, 6))
            plt.scatter(mz_values, intensities, color='blue', label='Observed m/z values')
            plt.axvline(x=center, color='red', linestyle='--', 
                       label=f'Median: {center:.6f}')
            plt.xlabel('m/z')
            plt.ylabel('Intensity')
            plt.title('Median of m/z Values (Gaussian Fit Failed)')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()
        
        return center, std_dev

# Define the preprocessing function
def preprocess_ms2_data(ms2_extracted_data, mz_tol=0.02):
    """
    Preprocess MS2 data by grouping similar m/z values and applying Gaussian fits.
    
    Parameters:
    -----------
    ms2_extracted_data : DataFrame
        DataFrame containing MS2 data with fragment_mz and fragment_intensity columns.
    mz_tol : float, optional
        Tolerance for grouping m/z values (in Da).
        
    Returns:
    --------
    processed_data : DataFrame
        Processed MS2 data with consolidated m/z values.
    """
    # Make a copy to avoid modifying the original data
    processed_data = ms2_extracted_data.copy()
    
    # Group the data by scan_number and precursor_mz
    grouped = processed_data.groupby(['scan_number', 'precursor_mz'])
    
    result_dfs = []
    
    for (scan, precursor), group in grouped:
        # Sort by fragment_mz
        group = group.sort_values('fragment_mz')
        
        # Initialize variables for grouping
        current_group_mzs = []
        current_group_intensities = []
        current_group_indices = []
        
        # Process each row
        for i, row in group.iterrows():
            mz = row['fragment_mz']
            intensity = row['fragment_intensity']
            
            # If this is the first point or it's close to the previous group
            if not current_group_mzs or abs(mz - current_group_mzs[-1]) <= mz_tol:
                current_group_mzs.append(mz)
                current_group_intensities.append(intensity)
                current_group_indices.append(i)
            else:
                # Process the current group if it has at least 1 point
                if len(current_group_mzs) >= 1:
                    # If there's only one point, keep it as is
                    if len(current_group_mzs) == 1:
                        # No need to modify this single point
                        pass
                    else:
                        # Apply Gaussian fit to multiple points
                        center, _ = gaussian_fit_mz_values(current_group_mzs, current_group_intensities)
                        
                        # Update the m/z values for all points in this group
                        for idx in current_group_indices:
                            processed_data.at[idx, 'fragment_mz'] = center
                
                # Start a new group
                current_group_mzs = [mz]
                current_group_intensities = [intensity]
                current_group_indices = [i]
        
        # Process the last group
        if len(current_group_mzs) > 1:
            center, _ = gaussian_fit_mz_values(current_group_mzs, current_group_intensities)
            for idx in current_group_indices:
                processed_data.at[idx, 'fragment_mz'] = center
    
    # Now consolidate rows with identical fragment_mz values (after Gaussian fitting)
    consolidated = []
    
    for (scan, precursor), group in processed_data.groupby(['scan_number', 'precursor_mz']):
        # Find unique fragment_mz values after fitting
        unique_mzs = group['fragment_mz'].unique()
        
        for unique_mz in unique_mzs:
            # Get all rows with this m/z
            same_mz_rows = group[group['fragment_mz'] == unique_mz]
            
            # Sum the intensities
            total_intensity = same_mz_rows['fragment_intensity'].sum()
            
            # Take the first row and update the intensity
            consolidated_row = same_mz_rows.iloc[0].copy()
            consolidated_row['fragment_intensity'] = total_intensity
            
            consolidated.append(consolidated_row)
    
    # Create the consolidated DataFrame
    result = pd.DataFrame(consolidated)
    
    print(f"Preprocessed MS2 data: reduced from {len(ms2_extracted_data)} to {len(result)} fragment ions after grouping similar m/z values")
    
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
    ms2_filtered_df = preprocess_ms2_data(ms2_filtered_df, mz_tol=mz_tol/2)
    
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

# Now let's create the main test script
if __name__ == "__main__":
    # Define the paths
    test_mzml_file = "sample_data/mzML/AZ_0_5ug_R1.mzML"
    write_path = "sample_data/mzML/AZ_0_5ug_R1_matched.csv"
    
    # Extract MS2 data
    print("Extracting MS2 data from mzML file...")
    extracted_ms2_data = extractMS2(test_mzml_file, min_intensity=1e2)
    print(f"Extracted {len(extracted_ms2_data)} MS2 data points")
    
    # Match MS1 data
    print("Matching MS1 data...")
    matched_results = matchMS1(extracted_ms2_data, ppm_tol=10, mz_min=400, mz_max=2000)
    print(f"Found {len(matched_results)} matched MS1 precursors")
    
    # Match MS2 data with the updated function
    print("Matching MS2 data with Gaussian fitting...")
    matched_ms2 = matchMS2(
        extracted_ms2_data, 
        matched_results, 
        precursor_composition="25000", 
        mz_tol=0.02, 
        intensity_threshold=1e2, 
        ppm_tol=10
    )
    
    # Save the results
    matched_ms2.to_csv(write_path, index=False)
    print(f"Matched MS2 data saved to {write_path}")