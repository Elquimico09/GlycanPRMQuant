import multiprocessing
from glycanPRMQuant.parallelProcess import run_parallel_pipeline

if __name__ == '__main__':
    # On Windows, this is required for any multiprocessing spawn
    multiprocessing.freeze_support()

    run_parallel_pipeline(
        input_dir="sample_data/mzML",
        output_root="sample_data/results_parallel",
        n_workers=8,
        ppm_ms1_tol=40,
        mz_min=400,
        mz_max=2000,
        intensity_threshold=1e2,
        ppm_ms2_tol=50,
        mz_tol=0.05,
        smoothing_window=20
    )