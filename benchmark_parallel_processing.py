import multiprocessing
import random
import time
from pathlib import Path

import pandas as pd

from glycanPRMQuant.parallelProcess import run_parallel_pipeline


INPUT_FILES = [
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_5_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_5_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_6_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_6_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_7_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_7_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_8_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_8_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_9_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\G_9_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_1_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_1_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_2_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_2_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_3_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_3_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_4_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_4_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_5_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_5_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_6_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_6_2_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_7_1_10252024.mzML",
    r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\Pompe_7_2_10252024.mzML"
]

OUTPUT_ROOT = Path(r"C:\Users\Vishal\Documents\PRM Pompe Glycomics\parallel_bench_output_22_files")
WORKER_COUNTS = list(range(1, 9))
REPEATS = 3


def _validate_inputs(paths):
    missing = [p for p in paths if not Path(p).is_file()]
    if missing:
        missing_str = "\n".join(missing)
        raise FileNotFoundError(f"Missing input mzML file(s):\n{missing_str}")


def main():
    _validate_inputs(INPUT_FILES)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    jobs = [(workers, rep) for workers in WORKER_COUNTS for rep in range(1, REPEATS + 1)]
    random.shuffle(jobs)

    rows = []
    for worker_count, repeat in jobs:
        print(f"Running with {worker_count} workers, repeat {repeat}...")
        out_dir = OUTPUT_ROOT / f"workers_{worker_count}" / f"repeat_{repeat}"
        out_dir.mkdir(parents=True, exist_ok=True)

        t0 = time.perf_counter()
        run_parallel_pipeline(
            input_files=INPUT_FILES,
            output_root=str(out_dir),
            n_workers=worker_count,
            ppm_ms1_tol=10,
            mz_min=400,
            mz_max=2000,
            intensity_threshold=1000,
            ppm_ms2_tol=20,
            mz_tol=0.02,
            mass_offset= 1,
            smoothing_window=4,
            smoothing_method="savgol",
            overwrite=True,
            enable_adduct_plots=True,
            enable_total_plots=True,
            dry_run=False,
            rel_height=0.9,
            rel_height_mode="height",
            skyline_transition=False,
            enable_smoothing=True,
        )
        elapsed = time.perf_counter() - t0
        print(f"Finished with {worker_count} workers, repeat {repeat}: {elapsed:.2f} s")
        rows.append(
            {
                "worker_count": worker_count,
                "repeat": repeat,
                "elapsed_seconds": round(elapsed, 3),
            }
        )

    raw_df = pd.DataFrame(rows).sort_values(["worker_count", "repeat"]).reset_index(drop=True)
    raw_path = OUTPUT_ROOT / "parallel_benchmark_results.csv"
    raw_df.to_csv(raw_path, index=False)

    summary_df = (
        raw_df.groupby("worker_count", as_index=False)["elapsed_seconds"]
        .agg(["mean", "median", "min", "max"])
        .reset_index()
    )
    baseline = summary_df.loc[summary_df["worker_count"] == 1, "median"].iloc[0]
    summary_df["speedup_vs_1worker"] = baseline / summary_df["median"]
    summary_path = OUTPUT_ROOT / "parallel_benchmark_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    print(f"\nSaved raw timings: {raw_path}")
    print(f"Saved summary: {summary_path}")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
