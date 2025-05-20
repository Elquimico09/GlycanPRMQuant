import pandas as pd
import numpy as np


def matchMS1(ms1_data, ppm_tol=10, mz_min=400, mz_max=2000):
  """
  Match MS1 data with glycan database.

  :param glycan_database: DataFrame containing glycan database with various adduct columns.
  :param ms1_data: DataFrame containing MS1 data with 'precursor_mz' column.
  :param ppm_tol: Mass tolerance in parts per million.
  :param mz_min: Minimum m/z value to consider (default: 400).
  :param mz_max: Maximum m/z value to consider (default: 2000).
  :return: DataFrame with matched MS1 data and corresponding glycan information.
  """
  print(f"Starting MS1 matching with ppm tolerance: {ppm_tol}")
  
  print("Loading glycan database...")
  glycan_database = pd.read_excel("database/glycan_precursor_mz_list.xlsx")
  print(f"Loaded glycan database with {len(glycan_database)} entries")

  # Filter MS1 data to only include precursors within the specified m/z range
  print(f"Filtering MS1 data with m/z range {mz_min}-{mz_max}...")
  filtered_ms1_data = ms1_data[(ms1_data['precursor_mz'] >= mz_min) & 
                              (ms1_data['precursor_mz'] <= mz_max)]
  print(f"Filtered MS1 data: {len(filtered_ms1_data)} entries (from {len(ms1_data)} total)")
  
  # Extract unique precursor m/z values to avoid redundant matching
  unique_precursors = filtered_ms1_data['precursor_mz'].unique()
  print(f"Found {len(unique_precursors)} unique precursor m/z values")
  
  # Identify adduct columns in the glycan database
  adduct_columns = [col for col in glycan_database.columns 
            if col.startswith('[M+') and col.endswith('+')]
  print(f"Found {len(adduct_columns)} adduct columns in database: {adduct_columns}")
  
  # Create a dictionary to store matches for each unique precursor m/z
  precursor_matches = {}
  match_count = 0
  
  # For each unique precursor m/z
  print("Starting matching process...")
  for i, precursor_mz in enumerate(unique_precursors):
    if i % 100 == 0 and i > 0:
      print(f"  Processed {i}/{len(unique_precursors)} precursors, found {match_count} matches so far")
    
    # Calculate tolerance in m/z units
    tolerance = (ppm_tol * precursor_mz) / 1e6
    
    # Check each adduct column for matches
    for adduct_col in adduct_columns:
      # Find matches within tolerance and m/z range
      matches = glycan_database[
        (glycan_database[adduct_col] >= precursor_mz - tolerance) & 
        (glycan_database[adduct_col] <= precursor_mz + tolerance) &
        (glycan_database[adduct_col] >= mz_min) &
        (glycan_database[adduct_col] <= mz_max)
      ]
      
      # Add matches to precursor_matches dictionary
      for _, match_row in matches.iterrows():
        if precursor_mz not in precursor_matches:
          precursor_matches[precursor_mz] = []
        
        precursor_matches[precursor_mz].append({
          'matched_glycan': str(match_row.get('Composition', '')),
          'matched_adduct': adduct_col,
          'database_mz': match_row[adduct_col],
          'ppm_error': ((precursor_mz - match_row[adduct_col]) / 
              match_row[adduct_col]) * 1e6
        })
        match_count += 1
  print(f"Matching completed. Found matches for {len(precursor_matches)} unique precursor m/z values")
  print(f"Total number of matches: {match_count}")
  
  # Create a list to store matched results
  matched_results = []
  
  # Create a dataframe with unique precursor m/z values and their matches
  print("Creating final results table with one row per unique precursor...")
  unique_results = []
  
  for precursor_mz, matches in precursor_matches.items():
    # Collect all matched glycans for this precursor
    matched_glycans = [match['matched_glycan'] for match in matches]
    # Join them with a separator to create a single string of all matches
    matched_compositions = '; '.join(matched_glycans)
    
    # Create a result row with just the precursor m/z and matched compositions
    result_row = {
      'precursor_mz': precursor_mz,
      'Glycan': matched_compositions
    }
    unique_results.append(result_row)
  
  # Return results as DataFrame with one row per unique precursor
  result_df = pd.DataFrame(unique_results) if unique_results else pd.DataFrame()
  print(f"Final matched results: {len(result_df)} entries")
  return result_df

