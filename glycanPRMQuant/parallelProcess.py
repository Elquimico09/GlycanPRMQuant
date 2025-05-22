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
):
    """
    Worker wrapper: skips processing if AUC file already exists.
    Returns (basename, status, message) where status is 'done', 'skipped', or 'error'.
    """
    base = os.path.splitext(os.path.basename(mzml_path))[0]
    out_dir = os.path.join(output_root, base)
    auc_file = os.path.join(out_dir, f"{base}_auc_values.csv")

    # Skip if AUC file exists
    if os.path.isfile(auc_file):
        return base, 'skipped', 'AUC file already exists'

    os.makedirs(out_dir, exist_ok=True)

    try:
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
        return base, 'done', None
    except Exception as e:
        return base, 'error', str(e)

def run_parallel_pipeline(
    input_dir: str,
    output_root: str,
    n_workers: int = None,
    ppm_ms1_tol: float = 10,
    mz_min: float = 400,
    mz_max: float = 2000,
    intensity_threshold: float = 1e2,
    ppm_ms2_tol: float = 50,
    mz_tol: float = 0.05,
    smoothing_window: int = 20
):
    """
    Discover all .mzML files in `input_dir` and process them in parallel.
    Skips any sample whose AUC file already exists.
    """
    os.makedirs(output_root, exist_ok=True)
    mzml_files = [
        os.path.join(input_dir, fn)
        for fn in os.listdir(input_dir)
        if fn.lower().endswith('.mzml')
    ]
    if not mzml_files:
        print("No .mzML files found in", input_dir)
        return

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                _process_one_file,
                path,
                output_root,
                ppm_ms1_tol,
                mz_min,
                mz_max,
                intensity_threshold,
                ppm_ms2_tol,
                mz_tol,
                smoothing_window
            ): path for path in mzml_files
        }
        for fut in as_completed(futures):
            base, status, msg = fut.result()
            if status == 'done':
                print(f"[✓] Finished processing {base}")
            elif status == 'skipped':
                print(f"[→] Skipped {base}: {msg}")
            else:
                print(f"[✗] Error processing {base}: {msg}")
