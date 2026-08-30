"""Command-line interface for glycanPRMQuant."""

import argparse
import logging
import multiprocessing

from glycanPRMQuant.logging_utils import configure_logging


def _add_common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--ppm-ms1-tol", type=float, default=10)
    parser.add_argument("--mz-offset", type=float, default=0.0)
    parser.add_argument("--mass-offset", type=float, default=0.0)
    parser.add_argument("--intensity-threshold", type=float, default=1e2)
    parser.add_argument(
        "--fragment-mass-tol",
        type=float,
        default=0.02,
        help="Numeric fragment matching tolerance (default: 0.02)",
    )
    parser.add_argument(
        "--fragment-mass-tol-unit",
        type=str.lower,
        choices=["da", "ppm"],
        default="da",
        help="Unit for --fragment-mass-tol (default: da)",
    )
    parser.add_argument("--fragment-ion-series", default="BY")
    parser.add_argument("--fragment-max-cleavages", type=int, default=2)
    parser.add_argument("--smoothing-window", type=int, default=11)
    parser.add_argument("--smoothing-method", choices=["gaussian", "savgol"], default="gaussian")
    parser.add_argument("--disable-smoothing", action="store_true")
    parser.add_argument("--rel-height", type=float, default=0.7)
    parser.add_argument("--rel-height-mode", choices=["prominence", "height"], default="prominence")
    parser.add_argument("--precursor-db-path")
    parser.add_argument("--structure-db-path")
    parser.add_argument("--skyline-transition", action="store_true")
    parser.add_argument("--disable-isobaric-resolution", action="store_true")
    parser.add_argument("--candidate-min-fragments", type=int, default=2)
    parser.add_argument("--candidate-min-explained-intensity", type=float, default=0.005)
    parser.add_argument("--candidate-min-score", type=float, default=35.0)
    parser.add_argument("--candidate-min-evidence-difference", type=float, default=4.0)
    parser.add_argument("--candidate-mass-outlier-min-delta", type=float, default=2.0)
    parser.add_argument("--candidate-max-q-value", type=float, default=0.05)
    parser.add_argument("--disable-target-decoy", action="store_true")
    parser.add_argument("--target-decoy-seed", type=int, default=1729)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("-v", "--verbose", action="count", default=0)


def _log_level(args: argparse.Namespace) -> int:
    if args.quiet:
        return logging.WARNING
    if args.verbose >= 2:
        return logging.DEBUG
    return logging.INFO


def _run_one(args: argparse.Namespace) -> int:
    configure_logging(_log_level(args), force=True)
    from glycanPRMQuant.processmzML import process_mzml_pipeline

    process_mzml_pipeline(
        mzml_file=args.input_file,
        output_dir=args.output_dir,
        ppm_ms1_tol=args.ppm_ms1_tol,
        mz_offset=args.mz_offset,
        mass_offset=args.mass_offset,
        intensity_threshold=args.intensity_threshold,
        fragment_mass_tol=args.fragment_mass_tol,
        fragment_mass_tol_unit=args.fragment_mass_tol_unit,
        smoothing_window=args.smoothing_window,
        smoothing_method=args.smoothing_method,
        enable_smoothing=not args.disable_smoothing,
        rel_height=args.rel_height,
        rel_height_mode=args.rel_height_mode,
        skyline_transition=args.skyline_transition,
        fragment_ion_series=args.fragment_ion_series,
        fragment_max_cleavages=args.fragment_max_cleavages,
        resolve_isobaric_conflicts=not args.disable_isobaric_resolution,
        candidate_min_fragments=args.candidate_min_fragments,
        candidate_min_explained_intensity=args.candidate_min_explained_intensity,
        candidate_min_score=args.candidate_min_score,
        candidate_min_evidence_difference=args.candidate_min_evidence_difference,
        candidate_mass_outlier_min_delta=args.candidate_mass_outlier_min_delta,
        enable_target_decoy=not args.disable_target_decoy,
        candidate_max_q_value=args.candidate_max_q_value,
        target_decoy_seed=args.target_decoy_seed,
        precursor_db_path=args.precursor_db_path,
        structure_db_path=args.structure_db_path,
    )
    return 0


def _run_batch(args: argparse.Namespace) -> int:
    configure_logging(_log_level(args), force=True)
    multiprocessing.freeze_support()
    from glycanPRMQuant.parallelProcess import run_parallel_pipeline

    run_parallel_pipeline(
        input_dir=args.input_dir,
        input_files=args.input_files,
        output_root=args.output_root,
        n_workers=args.workers,
        ppm_ms1_tol=args.ppm_ms1_tol,
        mz_offset=args.mz_offset,
        mass_offset=args.mass_offset,
        intensity_threshold=args.intensity_threshold,
        fragment_mass_tol=args.fragment_mass_tol,
        fragment_mass_tol_unit=args.fragment_mass_tol_unit,
        fragment_ion_series=args.fragment_ion_series,
        fragment_max_cleavages=args.fragment_max_cleavages,
        smoothing_window=args.smoothing_window,
        smoothing_method=args.smoothing_method,
        enable_smoothing=not args.disable_smoothing,
        rel_height=args.rel_height,
        rel_height_mode=args.rel_height_mode,
        skyline_transition=args.skyline_transition,
        precursor_db_path=args.precursor_db_path,
        structure_db_path=args.structure_db_path,
        resolve_isobaric_conflicts=not args.disable_isobaric_resolution,
        candidate_min_fragments=args.candidate_min_fragments,
        candidate_min_explained_intensity=args.candidate_min_explained_intensity,
        candidate_min_score=args.candidate_min_score,
        candidate_min_evidence_difference=args.candidate_min_evidence_difference,
        candidate_mass_outlier_min_delta=args.candidate_mass_outlier_min_delta,
        enable_target_decoy=not args.disable_target_decoy,
        candidate_max_q_value=args.candidate_max_q_value,
        target_decoy_seed=args.target_decoy_seed,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
    )
    return 0


def _run_gui(args: argparse.Namespace) -> int:
    configure_logging(_log_level(args), force=True)
    from glycanPRMQuant.pipelineGUI import PipelineGUI

    app = PipelineGUI()
    app.mainloop()
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="glycan-prmquant")
    sub = parser.add_subparsers(dest="command", required=True)

    run_parser = sub.add_parser("run", help="Process one Thermo RAW or mzML file")
    run_parser.add_argument("input_file")
    run_parser.add_argument("output_dir")
    _add_common_options(run_parser)
    run_parser.set_defaults(func=_run_one)

    batch_parser = sub.add_parser("batch", help="Process multiple Thermo RAW or mzML files")
    source = batch_parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input-dir")
    source.add_argument("--input-files", nargs="+")
    batch_parser.add_argument("--output-root", required=True)
    batch_parser.add_argument("--workers", type=int)
    batch_parser.add_argument("--overwrite", action="store_true")
    batch_parser.add_argument("--dry-run", action="store_true")
    _add_common_options(batch_parser)
    batch_parser.set_defaults(func=_run_batch)

    gui_parser = sub.add_parser("gui", help="Launch the Tkinter GUI")
    gui_parser.add_argument("--quiet", action="store_true")
    gui_parser.add_argument("-v", "--verbose", action="count", default=0)
    gui_parser.set_defaults(func=_run_gui)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
