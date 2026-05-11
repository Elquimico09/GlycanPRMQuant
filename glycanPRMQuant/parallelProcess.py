import os
import logging
from contextlib import redirect_stderr, redirect_stdout
from concurrent.futures import ProcessPoolExecutor, as_completed
from glycanPRMQuant.consolidateAUC import consolidate_auc_results
from glycanPRMQuant.logging_utils import configure_logging

logger = logging.getLogger(__name__)


class QueueWriter:
    """File-like writer that forwards stdout/stderr text to the GUI log queue."""

    def __init__(self, queue):
        self.queue = queue

    def write(self, data):
        if data:
            self.queue.put(data)

    def flush(self):
        pass

def _process_one_file(
    mzml_path: str,
    output_root: str,
    ppm_ms1_tol: float,
    mz_min: float,
    mz_max: float,
    intensity_threshold: float,
    ppm_ms2_tol: float,
    mz_tol: float,
    fragment_ion_series: str,
    fragment_max_cleavages: int,
    smoothing_window: int,
    smoothing_method: str = "gaussian",
    mz_offset: float = 0.0,
    mass_offset: float = 0.0,
    overwrite: bool = False,
    enable_adduct_plots: bool = True,
    enable_total_plots: bool = True,
    dry_run: bool = False,
    rel_height: float = 0.7,
    rel_height_mode: str = "prominence",
    skyline_transition: bool = False,
    enable_smoothing: bool = True,
    precursor_db_path: str = None,
    structure_db_path: str = None,
    log_queue=None
):
    """
    Worker wrapper: skips processing if AUC file already exists.
    Returns (basename, status, message) where status is 'done', 'skipped', or 'error'.
    """
    base = os.path.splitext(os.path.basename(mzml_path))[0]
    out_dir = os.path.join(output_root, base)
    auc_file = os.path.join(out_dir, f"{base}_auc_values.csv")

    # Skip if AUC file exists
    if not overwrite and os.path.isfile(auc_file):
        return base, 'skipped', 'AUC file already exists'

    os.makedirs(out_dir, exist_ok=True)

    if dry_run:
        if log_queue is not None:
            configure_logging(log_queue=log_queue, force=True)
        return base, 'dry-run', None

    try:
        configure_logging(log_queue=log_queue, force=log_queue is not None)
        from glycanPRMQuant.processmzML import process_mzml_pipeline

        def _run_pipeline():
            process_mzml_pipeline(
                mzml_file=mzml_path,
                output_dir=out_dir,
                ppm_ms1_tol=ppm_ms1_tol,
                mz_min=mz_min,
                mz_max=mz_max,
                intensity_threshold=intensity_threshold,
                ppm_ms2_tol=ppm_ms2_tol,
                mz_tol=mz_tol,
                fragment_ion_series=fragment_ion_series,
                fragment_max_cleavages=fragment_max_cleavages,
                smoothing_window=smoothing_window,
                smoothing_method=smoothing_method,
                mz_offset=mz_offset,
                mass_offset=mass_offset,
                enable_adduct_plots=enable_adduct_plots,
                enable_total_plots=enable_total_plots,
                rel_height=rel_height,
                rel_height_mode=rel_height_mode,
                skyline_transition=skyline_transition,
                enable_smoothing=enable_smoothing,
                precursor_db_path=precursor_db_path,
                structure_db_path=structure_db_path
            )

        if log_queue is None:
            _run_pipeline()
        else:
            writer = QueueWriter(log_queue)
            with redirect_stdout(writer), redirect_stderr(writer):
                _run_pipeline()
        return base, 'done', None
    except Exception as e:
        return base, 'error', str(e)

def run_parallel_pipeline(
    input_dir: str = None,
    input_files: list = None,
    output_root: str = None,
    n_workers: int = None,
    ppm_ms1_tol: float = 10,
    mz_min: float = 400,
    mz_max: float = 2000,
    intensity_threshold: float = 1e2,
    ppm_ms2_tol: float = 50,
    mz_tol: float = 0.05,
    fragment_ion_series: str = "ABCXYZ",
    fragment_max_cleavages: int = 2,
    smoothing_window: int = 20,
    smoothing_method: str = "gaussian",
    mz_offset: float = 0.0,
    mass_offset: float = 0.0,
    overwrite: bool = False,
    enable_adduct_plots: bool = True,
    enable_total_plots: bool = True,
    dry_run: bool = False,
    rel_height: float = 0.7,
    rel_height_mode: str = "prominence",
    skyline_transition: bool = False,
    enable_smoothing: bool = True,
    precursor_db_path: str = None,
    structure_db_path: str = None,
    log_queue=None,
    progress_queue=None
):
    """
    Discover all .mzML files in `input_dir` (or use explicit list) and process them in parallel.
    Skips any sample whose AUC file already exists.
    """
    if output_root is None:
        raise ValueError("output_root must be provided")

    os.makedirs(output_root, exist_ok=True)

    if input_files:
        mzml_files = list(input_files)
    else:
        if not input_dir:
            raise ValueError("Either input_files or input_dir must be provided")
        mzml_files = [
            os.path.join(input_dir, fn)
            for fn in os.listdir(input_dir)
            if fn.lower().endswith('.mzml')
        ]
    if not mzml_files:
        logger.warning("No .mzML files found")
        return

    # Clamp worker count to a safe upper bound for Windows (ProcessPool has a hard cap)
    if n_workers is None or n_workers <= 0:
        max_workers = os.cpu_count() or 1
    else:
        max_workers = n_workers
    max_workers = min(max_workers, 61)

    def _run():
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
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
                    fragment_ion_series,
                    fragment_max_cleavages,
                    smoothing_window,
                    smoothing_method,
                    mz_offset,
                    mass_offset,
                    overwrite,
                    enable_adduct_plots,
                    enable_total_plots,
                    dry_run,
                    rel_height,
                    rel_height_mode,
                    skyline_transition,
                    enable_smoothing,
                    precursor_db_path,
                    structure_db_path,
                    log_queue
                ): path for path in mzml_files
            }
            for fut in as_completed(futures):
                base, status, msg = fut.result()
                if progress_queue:
                    progress_queue.put((base, status, msg))
                if status == 'done':
                    logger.info("[✓] Finished processing %s", base)
                elif status == 'skipped':
                    logger.info("[→] Skipped %s: %s", base, msg)
                elif status == 'dry-run':
                    logger.info("[i] Planned (dry-run) %s", base)
                else:
                    logger.error("[✗] Error processing %s: %s", base, msg)
        if progress_queue:
            progress_queue.put(None)

    if log_queue is None:
        configure_logging()
        _run()
    else:
        configure_logging(log_queue=log_queue, force=True)
        _run()

    # After processing, write combined AUC summary if applicable
    if (not dry_run) and len(mzml_files) > 1:
        try:
            combined_path = os.path.join(output_root, "combined_auc_values.csv")
            consolidate_auc_results(output_root, combined_path)
            logger.info("[✓] Wrote combined AUC table to %s", combined_path)
        except Exception as e:
            logger.error("[✗] Failed to write combined AUC table: %s", e)

    if log_queue is not None:
        log_queue.put(None)  # sentinel
