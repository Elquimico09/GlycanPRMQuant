import argparse
import json
import multiprocessing
import os
import time
from pathlib import Path

from glycanPRMQuant.parallelProcess import run_parallel_pipeline


def _build_smoothing_configs():
    configs = [
        {
            "name": "no_smoothing",
            "enable_smoothing": False,
            "smoothing_method": "gaussian",
            "smoothing_window": 0,
        }
    ]
    for w in (2, 4, 6, 8):
        configs.append(
            {
                "name": f"gaussian_{w}",
                "enable_smoothing": True,
                "smoothing_method": "gaussian",
                "smoothing_window": w,
            }
        )
    for w in (2, 4, 6, 8):
        configs.append(
            {
                "name": f"savgol_{w}",
                "enable_smoothing": True,
                "smoothing_method": "savgol",
                "smoothing_window": w,
            }
        )
    return configs


def _validate_inputs(files):
    missing = [f for f in files if not Path(f).is_file()]
    if missing:
        joined = "\n".join(missing)
        raise FileNotFoundError(f"These input files were not found:\n{joined}")


def _evaluate_run(config_output_root: Path, input_files, dry_run: bool):
    if dry_run:
        return "dry-run", []

    missing_auc = []
    for f in input_files:
        stem = Path(f).stem
        auc_file = config_output_root / stem / f"{stem}_auc_values.csv"
        if not auc_file.is_file():
            missing_auc.append(str(auc_file))

    if missing_auc:
        return "incomplete", missing_auc
    return "ok", []


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark glycanPRMQuant across smoothing modes/windows on the same input files."
        )
    )
    parser.add_argument(
        "--input-files",
        nargs="+",
        required=True,
        help="Paths to mzML files (e.g., 3 files you want to benchmark).",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        help="Root output directory for benchmark runs.",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=3,
        help="Parallel workers per benchmark run.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Plan only; do not process files.")

    # Forwarded pipeline parameters (keep defaults aligned with GUI/pipeline defaults)
    parser.add_argument("--ppm-ms1-tol", type=float, default=10.0)
    parser.add_argument("--mz-min", type=float, default=400.0)
    parser.add_argument("--mz-max", type=float, default=2000.0)
    parser.add_argument("--intensity-threshold", type=float, default=1e2)
    parser.add_argument("--ppm-ms2-tol", type=float, default=10.0)
    parser.add_argument("--mz-tol", type=float, default=0.02)
    parser.add_argument("--mz-offset", type=float, default=0.0)
    parser.add_argument("--mass-offset", type=float, default=0.0)
    parser.add_argument("--rel-height", type=float, default=0.7)
    parser.add_argument(
        "--rel-height-mode",
        choices=["prominence", "height"],
        default="prominence",
    )
    parser.add_argument("--enable-adduct-plots", action="store_true", default=True)
    parser.add_argument("--disable-adduct-plots", action="store_true")
    parser.add_argument("--enable-total-plots", action="store_true", default=True)
    parser.add_argument("--disable-total-plots", action="store_true")
    parser.add_argument("--skyline-transition", action="store_true", default=False)

    args = parser.parse_args()

    input_files = [str(Path(f).resolve()) for f in args.input_files]
    _validate_inputs(input_files)

    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    enable_adduct_plots = args.enable_adduct_plots and not args.disable_adduct_plots
    enable_total_plots = args.enable_total_plots and not args.disable_total_plots

    configs = _build_smoothing_configs()
    results = []

    print(f"Running {len(configs)} smoothing configurations on {len(input_files)} files...")
    for i, cfg in enumerate(configs, start=1):
        cfg_root = output_root / cfg["name"]
        cfg_root.mkdir(parents=True, exist_ok=True)
        print(
            f"[{i}/{len(configs)}] {cfg['name']} | method={cfg['smoothing_method']} "
            f"window={cfg['smoothing_window']} enable_smoothing={cfg['enable_smoothing']}"
        )

        t0 = time.perf_counter()
        run_parallel_pipeline(
            input_files=input_files,
            output_root=str(cfg_root),
            n_workers=args.n_workers,
            ppm_ms1_tol=args.ppm_ms1_tol,
            mz_min=args.mz_min,
            mz_max=args.mz_max,
            intensity_threshold=args.intensity_threshold,
            ppm_ms2_tol=args.ppm_ms2_tol,
            mz_tol=args.mz_tol,
            smoothing_window=cfg["smoothing_window"],
            smoothing_method=cfg["smoothing_method"],
            mz_offset=args.mz_offset,
            mass_offset=args.mass_offset,
            overwrite=args.overwrite,
            enable_adduct_plots=enable_adduct_plots,
            enable_total_plots=enable_total_plots,
            dry_run=args.dry_run,
            rel_height=args.rel_height,
            rel_height_mode=args.rel_height_mode,
            skyline_transition=args.skyline_transition,
            enable_smoothing=cfg["enable_smoothing"],
        )
        elapsed = time.perf_counter() - t0

        status, missing_auc = _evaluate_run(cfg_root, input_files, args.dry_run)
        results.append(
            {
                "config": cfg["name"],
                "enable_smoothing": cfg["enable_smoothing"],
                "smoothing_method": cfg["smoothing_method"],
                "smoothing_window": cfg["smoothing_window"],
                "elapsed_seconds": round(elapsed, 3),
                "status": status,
                "output_root": str(cfg_root),
                "missing_auc_files": missing_auc,
            }
        )

    summary_json = output_root / "benchmark_summary.json"
    with summary_json.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    summary_csv = output_root / "benchmark_summary.csv"
    with summary_csv.open("w", encoding="utf-8") as f:
        f.write(
            "config,enable_smoothing,smoothing_method,smoothing_window,elapsed_seconds,status,output_root\n"
        )
        for r in results:
            f.write(
                f"{r['config']},{r['enable_smoothing']},{r['smoothing_method']},"
                f"{r['smoothing_window']},{r['elapsed_seconds']},{r['status']},"
                f"\"{r['output_root']}\"\n"
            )

    print(f"\nBenchmark complete. Summary files:")
    print(f"  {summary_json}")
    print(f"  {summary_csv}")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
