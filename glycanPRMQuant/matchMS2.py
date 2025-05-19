import pandas as pd
import numpy as np

def matchMS2(ms2_extracted_data, precursor_matched_data, precursor_composition,
              mz_tol = 0.02, intensity_threshold = 1e2):
        """
        Match MS2 data with glycan database.
        :param ms2_extracted_data: DataFrame containing extracted MS2 data with 'scan_number', 'rt', 'intensity' columns.
        :param precursor_matched_data: DataFrame containing matched MS1 data with 'precursor_mz' column.
        :param
        :param ms2_database: DataFrame containing MS2 database with glycan information.
        :param mz_tol: Mass tolerance in Da (default: 0.02).
        :param
        :param intensity_threshold: Minimum intensity threshold for MS2 data (default: 1e2).
        :return: DataFrame with matched MS2 data and corresponding glycan information.
        """
        print(f"Starting MS2 matching with m/z tolerance: {mz_tol} Da")
        print("Loading MS2 database...")
        # Load the MS2 database
        ms2_database_path = "database/fragment_database.csv"
        ms2_database = pd.read_csv(ms2_database_path)
        print(f"Loaded MS2 database with {len(ms2_database)} entries")
