# glycanPRMQuant/parallelProcess.py

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from glycanPRMQuant.processmzML import process_mzml_pipeline

def _process_one_file(
    mzml_path: str,
    output_root: str,
    ppm_ms1_tol: float,
    mz_min: float,
    mz_max: float,
    intensity_threshold: float,
    ppm_ms2_tol: float,
    mz_tol: float,
    smoothing_window: int
) -> str:
    """
    Worker function to process a single .mzML file.
    Returns the basename (without extension) on success.
    """
    base = os.path.splitext(os.path.basename(mzml_path))[0]
    out_dir = os.path.join(output_root, base)
    os.makedirs(out_dir, exist_ok=True)

    process_mzml_pipeline(
        mzml_file=mzml_path,
        output_dir=out_dir,
        ppm_ms1_tol=ppm_ms1_tol,
        mz_min=mz_min,
        mz_max=mz_max,
        intensity_threshold=intensity_threshold,
        ppm_ms2_tol=ppm_ms2_tol,
        mz_tol=mz_tol,
        smoothing_window=smoothing_window
    )

    return base


def run_parallel_pipeline(
    input_dir: str,
    output_root: str,
    n_workers: int = None,
    # MS1 settings
    ppm_ms1_tol: float = 10,
    mz_min: float = 400,
    mz_max: float = 2000,
    # MS2 settings
    intensity_threshold: float = 1e2,
    ppm_ms2_tol: float = 50,
    mz_tol: float = 0.05,
    smoothing_window: int = 20
):
    """
    Discover all .mzML files in `input_dir` and process them in parallel.
    Each file gets its own subfolder under `output_root` named after the file (no extension).

    Parameters
    ----------
    input_dir : str
        Folder containing .mzML files.
    output_root : str
        Root folder where per‐file output directories will be created.
    n_workers : int, optional
        Number of parallel workers (defaults to os.cpu_count()).
    ... (other parameters are passed to process_mzml_pipeline)
    """
    os.makedirs(output_root, exist_ok=True)

    # find all .mzML files
    mzml_files = [
        os.path.join(input_dir, fn)
        for fn in os.listdir(input_dir)
        if fn.lower().endswith('.mzml')
    ]
    if not mzml_files:
        print("No .mzML files found in", input_dir)
        return

    # submit each to the pool
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                _process_one_file,
                mzml_path,
                output_root,
                ppm_ms1_tol,
                mz_min,
                mz_max,
                intensity_threshold,
                ppm_ms2_tol,
                mz_tol,
                smoothing_window
            ): mzml_path
            for mzml_path in mzml_files
        }

        for fut in as_completed(futures):
            mzml_path = futures[fut]
            try:
                base = fut.result()
                print(f"[✓] Finished processing {base}")
            except Exception as e:
                print(f"[✗] Error processing {mzml_path}: {e}")
