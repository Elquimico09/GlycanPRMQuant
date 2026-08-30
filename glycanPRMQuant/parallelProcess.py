import os
import logging
from contextlib import redirect_stderr, redirect_stdout
from concurrent.futures import ProcessPoolExecutor, as_completed
from glycanPRMQuant.consolidateAUC import consolidate_auc_results
from glycanPRMQuant.logging_utils import configure_logging
from glycanPRMQuant.retention_alignment import consolidate_consensus_peak_results
from glycanPRMQuant.spectra import SUPPORTED_SUFFIXES, validate_input_file_types

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
    intensity_threshold: float,
    fragment_mass_tol: float,
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
    log_queue=None,
    resolve_isobaric_conflicts: bool = True,
    candidate_min_fragments: int = 2,
    candidate_min_explained_intensity: float = 0.01,
    candidate_min_score: float = 35.0,
    candidate_min_evidence_difference: float = 4.0,
    candidate_mass_outlier_min_delta: float = 2.0,
    enable_target_decoy: bool = True,
    candidate_max_q_value: float = 0.05,
    target_decoy_seed: int = 1729,
    fragment_mass_tol_unit: str = "Da",
    figure_filetype: str = "pdf",
):
    """
    Worker wrapper: skips processing if AUC file already exists.
    Returns (basename, status, message) where status is 'done', 'skipped', or 'error'.
    """
    base = os.path.splitext(os.path.basename(mzml_path))[0]
    out_dir = os.path.join(output_root, base)
    auc_file = os.path.join(out_dir, f"{base}_auc_values.csv")
    feature_auc_file = os.path.join(out_dir, f"{base}_feature_auc_values.csv")

    # Newer batch consolidation also requires feature-level AUCs. Older output
    # folders that contain only the composition-level table are reprocessed.
    if (
        not overwrite
        and os.path.isfile(auc_file)
        and os.path.isfile(feature_auc_file)
    ):
        return base, 'skipped', 'AUC and feature-level AUC files already exist'

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
                intensity_threshold=intensity_threshold,
                fragment_mass_tol=fragment_mass_tol,
                fragment_mass_tol_unit=fragment_mass_tol_unit,
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
                resolve_isobaric_conflicts=resolve_isobaric_conflicts,
                candidate_min_fragments=candidate_min_fragments,
                candidate_min_explained_intensity=candidate_min_explained_intensity,
                candidate_min_score=candidate_min_score,
                candidate_min_evidence_difference=candidate_min_evidence_difference,
                candidate_mass_outlier_min_delta=candidate_mass_outlier_min_delta,
                enable_target_decoy=enable_target_decoy,
                candidate_max_q_value=candidate_max_q_value,
                target_decoy_seed=target_decoy_seed,
                precursor_db_path=precursor_db_path,
                structure_db_path=structure_db_path,
                figure_filetype=figure_filetype,
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
    intensity_threshold: float = 1e2,
    fragment_mass_tol: float = 0.05,
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
    progress_queue=None,
    resolve_isobaric_conflicts: bool = True,
    candidate_min_fragments: int = 2,
    candidate_min_explained_intensity: float = 0.01,
    candidate_min_score: float = 35.0,
    candidate_min_evidence_difference: float = 4.0,
    candidate_mass_outlier_min_delta: float = 2.0,
    enable_target_decoy: bool = True,
    candidate_max_q_value: float = 0.05,
    target_decoy_seed: int = 1729,
    fragment_mass_tol_unit: str = "Da",
    figure_filetype: str = "pdf",
    enable_consensus_peak_selection: bool = True,
    consensus_rt_tolerance: float = 0.3,
    consensus_min_replicate_fraction: float = 0.8,
):
    """
    Discover all Thermo .raw or .mzML files in `input_dir` (or use an explicit list) and process them.
    Skips a sample only when both its legacy glycan-level AUC and feature-level
    AUC files already exist. Multi-file runs align feature retention times and
    select the most reproducibly scored peak group by default.
    """
    if output_root is None:
        raise ValueError("output_root must be provided")
    if consensus_rt_tolerance <= 0:
        raise ValueError("Consensus RT tolerance must be positive")
    if not 0 < consensus_min_replicate_fraction <= 1:
        raise ValueError("Consensus minimum replicate fraction must be in (0, 1]")

    os.makedirs(output_root, exist_ok=True)

    if input_files:
        ms_files = list(input_files)
    else:
        if not input_dir:
            raise ValueError("Either input_files or input_dir must be provided")
        ms_files = [
            os.path.join(input_dir, fn)
            for fn in os.listdir(input_dir)
            if os.path.splitext(fn)[1].lower() in SUPPORTED_SUFFIXES
        ]
    if not ms_files:
        logger.warning("No Thermo .raw or .mzML files found")
        return

    validate_input_file_types(ms_files)

    # Clamp worker count to a safe upper bound for Windows (ProcessPool has a hard cap)
    if n_workers is None or n_workers <= 0:
        max_workers = os.cpu_count() or 1
    else:
        max_workers = n_workers
    max_workers = min(max_workers, 61)
    run_statuses = {}
    expected_samples = [
        os.path.splitext(os.path.basename(path))[0] for path in ms_files
    ]

    def _run():
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    _process_one_file,
                    path,
                    output_root,
                    ppm_ms1_tol,
                    intensity_threshold,
                    fragment_mass_tol,
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
                    log_queue,
                    resolve_isobaric_conflicts,
                    candidate_min_fragments,
                    candidate_min_explained_intensity,
                    candidate_min_score,
                    candidate_min_evidence_difference,
                    candidate_mass_outlier_min_delta,
                    enable_target_decoy,
                    candidate_max_q_value,
                    target_decoy_seed,
                    fragment_mass_tol_unit,
                    figure_filetype,
                ): path for path in ms_files
            }
            for fut in as_completed(futures):
                base, status, msg = fut.result()
                run_statuses[base] = status
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
    if (not dry_run) and len(ms_files) > 1:
        combined_path = os.path.join(output_root, "combined_auc_values.csv")
        failed_samples = sorted(
            sample for sample, status in run_statuses.items() if status == "error"
        )
        if enable_consensus_peak_selection and failed_samples:
            logger.error(
                "[✗] Cross-run consensus was not calculated because processing "
                "failed for: %s",
                ", ".join(failed_samples),
            )
        else:
            try:
                if enable_consensus_peak_selection:
                    consensus = consolidate_consensus_peak_results(
                        output_root,
                        combined_path,
                        rt_tolerance_minutes=consensus_rt_tolerance,
                        minimum_replicate_fraction=consensus_min_replicate_fraction,
                        expected_samples=expected_samples,
                    )
                    logger.info(
                        "[✓] Wrote consensus AUC table with %d selected peak group(s) to %s",
                        int(consensus.peak_groups["consensus_selected"].sum()),
                        combined_path,
                    )
                    logger.info(
                        "[✓] Wrote RT alignment and alternative-peak audits to %s",
                        output_root,
                    )
                else:
                    consolidate_auc_results(output_root, combined_path)
                    logger.info("[✓] Wrote combined AUC table to %s", combined_path)
            except Exception as e:
                stage = (
                    "Consensus peak selection"
                    if enable_consensus_peak_selection
                    else "AUC consolidation"
                )
                logger.error("[✗] %s failed: %s", stage, e)

    if log_queue is not None:
        log_queue.put(None)  # sentinel
