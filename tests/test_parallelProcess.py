import multiprocessing
from glycanPRMQuant.parallelProcess import run_parallel_pipeline

if __name__ == '__main__':
    # On Windows, this is required for any multiprocessing spawn
    multiprocessing.freeze_support()

    run_parallel_pipeline(
        input_dir="sample_data/RiboB",
        output_root="sample_data/RiboB/processed/",
        n_workers=4,
        ppm_ms1_tol=20,
        mz_min=400,
        mz_max=2000,
        intensity_threshold=1e2,
        ppm_ms2_tol=20,
        mz_tol=0.05,        
        smoothing_window=11,
        mass_offset=0,
        mz_offset=0.0
    )