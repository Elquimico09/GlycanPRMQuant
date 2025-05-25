from glycanPRMQuant.skylineTransition import dedupe_ms2_fragments
import pandas as pd

deduped = dedupe_ms2_fragments(
    ms2_csv_file="sample_data/results_parallel/AZ_2ug_R1/ms2_43202.csv",
    output_csv="sample_data/results_parallel/AZ_2ug_R1/ms2_43202_unique_fragments.csv",
    mz_tol=0.08
)