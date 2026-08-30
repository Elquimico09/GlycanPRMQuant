"""Cross-run retention-time alignment and consensus peak selection."""

from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import theilslopes


logger = logging.getLogger(__name__)


@dataclass
class ConsensusPeakResult:
    """Audit tables produced by cross-run consensus peak selection."""

    aligned_features: pd.DataFrame
    peak_groups: pd.DataFrame
    selected_auc: pd.DataFrame
    all_feature_auc: pd.DataFrame
    alignment_models: pd.DataFrame


def _median_absolute_deviation(values) -> float:
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if numeric.empty:
        return float("nan")
    center = float(numeric.median())
    return float((numeric - center).abs().median())


def _boolean_series(values, default: bool = True) -> pd.Series:
    series = pd.Series(values)
    if series.empty:
        return pd.Series(dtype=bool)
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(default).astype(bool)
    normalized = series.astype(str).str.strip().str.lower()
    mapped = normalized.map(
        {
            "true": True,
            "1": True,
            "yes": True,
            "false": False,
            "0": False,
            "no": False,
            "nan": default,
            "none": default,
            "": default,
        }
    )
    return mapped.fillna(default).astype(bool)


def _numeric_boundary(
    frame: pd.DataFrame, column: str, operation: str, fallback: float
) -> float:
    """Return a finite feature boundary, falling back to the peak apex."""
    if column not in frame.columns:
        return fallback
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    if values.empty:
        return fallback
    return float(values.min() if operation == "min" else values.max())


def _cluster_rows_by_rt(group: pd.DataFrame, rt_column: str, tolerance: float):
    """Assign rows to compact RT groups without single-linkage chaining."""
    clusters: list[list[int]] = []
    for index in group.sort_values(rt_column).index:
        rt = float(group.at[index, rt_column])
        candidates = []
        for cluster_index, members in enumerate(clusters):
            center = float(group.loc[members, rt_column].median())
            distance = abs(rt - center)
            if distance <= tolerance:
                candidates.append((distance, cluster_index))
        if candidates:
            _, cluster_index = min(candidates)
            clusters[cluster_index].append(index)
        else:
            clusters.append([index])
    return clusters


def build_within_run_peaks(
    feature_auc: pd.DataFrame,
    rt_tolerance_minutes: float,
) -> pd.DataFrame:
    """Combine adduct-level feature AUCs that represent the same run-local peak."""
    required = {"sample", "Glycan", "peak_rt", "AUC", "candidate_score"}
    missing = required - set(feature_auc.columns)
    if missing:
        raise ValueError(
            f"Consensus peak selection requires columns: {sorted(missing)}"
        )
    if rt_tolerance_minutes <= 0:
        raise ValueError("Consensus RT tolerance must be positive")

    work = feature_auc.copy()
    work["sample"] = work["sample"].astype(str)
    work["Glycan"] = work["Glycan"].astype(str).str.replace(r"\.0$", "", regex=True)
    work["peak_rt"] = pd.to_numeric(work["peak_rt"], errors="coerce")
    work["AUC"] = pd.to_numeric(work["AUC"], errors="coerce").fillna(0.0)
    work["candidate_score"] = pd.to_numeric(
        work["candidate_score"], errors="coerce"
    )
    work = work.dropna(subset=["peak_rt", "candidate_score"])

    records = []
    for (sample, glycan), group in work.groupby(["sample", "Glycan"], sort=True):
        clusters = _cluster_rows_by_rt(group, "peak_rt", rt_tolerance_minutes)
        ordered_clusters = sorted(
            clusters,
            key=lambda indices: float(group.loc[indices, "peak_rt"].median()),
        )
        for peak_number, indices in enumerate(ordered_clusters, start=1):
            peak = group.loc[indices]
            q_values = pd.to_numeric(
                peak.get("assignment_q_value", pd.Series(np.nan, index=peak.index)),
                errors="coerce",
            )
            pass_values = _boolean_series(
                peak.get("target_decoy_pass", pd.Series(True, index=peak.index))
            )
            adduct_column = (
                "PrecursorAdduct"
                if "PrecursorAdduct" in peak.columns
                else "Adduct"
                if "Adduct" in peak.columns
                else None
            )
            adducts = (
                ";".join(sorted(peak[adduct_column].dropna().astype(str).unique()))
                if adduct_column
                else ""
            )
            source_features = []
            if {"precursor_cluster", "rt_feature_id"}.issubset(peak.columns):
                source_features = sorted(
                    {
                        f"{cluster}:{feature}"
                        for cluster, feature in peak[
                            ["precursor_cluster", "rt_feature_id"]
                        ].itertuples(index=False, name=None)
                    }
                )
            records.append(
                {
                    "sample": sample,
                    "Glycan": glycan,
                    "run_peak_id": f"{glycan}_run_peak_{peak_number}",
                    "raw_peak_rt": float(peak["peak_rt"].median()),
                    "raw_start_rt": _numeric_boundary(
                        peak,
                        "start_rt",
                        "min",
                        float(peak["peak_rt"].median()),
                    ),
                    "raw_end_rt": _numeric_boundary(
                        peak,
                        "end_rt",
                        "max",
                        float(peak["peak_rt"].median()),
                    ),
                    "AUC": float(peak["AUC"].sum()),
                    "peak_candidate_score": float(peak["candidate_score"].max()),
                    "within_peak_score_mad": _median_absolute_deviation(
                        peak["candidate_score"]
                    ),
                    "assignment_q_value": (
                        float(q_values.max()) if q_values.notna().any() else np.nan
                    ),
                    "target_decoy_pass": bool(pass_values.all()),
                    "source_feature_count": int(len(peak)),
                    "source_features": ";".join(source_features),
                    "adduct_count": int(
                        peak[adduct_column].nunique() if adduct_column else 0
                    ),
                    "adducts": adducts,
                }
            )
    return pd.DataFrame(records)


def _build_alignment_anchors(
    run_peaks: pd.DataFrame,
    sample_count: int,
    minimum_replicate_fraction: float,
    rt_tolerance_minutes: float,
) -> pd.DataFrame:
    """Choose reproducible high-scoring peaks as preliminary RT landmarks."""
    minimum_samples = max(
        2, int(math.ceil(minimum_replicate_fraction * sample_count))
    )
    best_per_run = (
        run_peaks.sort_values(
            ["sample", "Glycan", "peak_candidate_score", "raw_peak_rt"],
            ascending=[True, True, False, True],
        )
        .drop_duplicates(["sample", "Glycan"], keep="first")
        .copy()
    )
    anchor_limit = max(1.0, 3.0 * rt_tolerance_minutes)
    anchors = []
    for glycan, group in best_per_run.groupby("Glycan"):
        if group["sample"].nunique() < minimum_samples:
            continue
        reference_rt = float(group["raw_peak_rt"].median())
        if float((group["raw_peak_rt"] - reference_rt).abs().max()) > anchor_limit:
            continue
        for row in group.itertuples(index=False):
            anchors.append(
                {
                    "sample": row.sample,
                    "Glycan": glycan,
                    "observed_rt": float(row.raw_peak_rt),
                    "reference_rt": reference_rt,
                }
            )
    return pd.DataFrame(anchors)


def _fit_alignment_models(
    samples: list[str], anchors: pd.DataFrame
) -> pd.DataFrame:
    models = []
    for sample in samples:
        sample_anchors = anchors.loc[anchors.get("sample", pd.Series(dtype=str)) == sample]
        source = "identity_no_anchors"
        slope = 1.0
        intercept = 0.0
        if not sample_anchors.empty:
            observed = sample_anchors["observed_rt"].to_numpy(dtype=float)
            reference = sample_anchors["reference_rt"].to_numpy(dtype=float)
            if len(np.unique(observed)) >= 3:
                fitted = theilslopes(reference, observed)
                proposed_slope = float(fitted.slope)
                proposed_intercept = float(fitted.intercept)
                if 0.8 <= proposed_slope <= 1.2:
                    slope = proposed_slope
                    intercept = proposed_intercept
                    source = "robust_linear_anchors"
                else:
                    intercept = float(np.median(reference - observed))
                    source = "median_shift_slope_guard"
            else:
                intercept = float(np.median(reference - observed))
                source = "median_shift_anchors"
        if sample_anchors.empty:
            residual_mad = np.nan
        else:
            predicted = slope * sample_anchors["observed_rt"] + intercept
            residual_mad = _median_absolute_deviation(
                predicted - sample_anchors["reference_rt"]
            )
        models.append(
            {
                "sample": sample,
                "alignment_slope": slope,
                "alignment_intercept": intercept,
                "alignment_anchor_count": int(len(sample_anchors)),
                "alignment_residual_mad": residual_mad,
                "alignment_source": source,
            }
        )
    return pd.DataFrame(models)


def _assign_consensus_groups(
    aligned: pd.DataFrame,
    rt_tolerance_minutes: float,
) -> pd.DataFrame:
    assignments = pd.Series(index=aligned.index, dtype=object)
    for glycan, glycan_peaks in aligned.groupby("Glycan", sort=True):
        groups: list[list[int]] = []
        ordered_indices = glycan_peaks.sort_values(
            ["peak_candidate_score", "aligned_peak_rt"],
            ascending=[False, True],
        ).index
        for index in ordered_indices:
            row = aligned.loc[index]
            candidates = []
            for group_index, members in enumerate(groups):
                member_samples = set(aligned.loc[members, "sample"])
                if row["sample"] in member_samples:
                    continue
                center = float(aligned.loc[members, "aligned_peak_rt"].median())
                distance = abs(float(row["aligned_peak_rt"]) - center)
                if distance <= rt_tolerance_minutes:
                    candidates.append((distance, group_index))
            if candidates:
                _, group_index = min(candidates)
                groups[group_index].append(index)
            else:
                groups.append([index])

        groups = sorted(
            groups,
            key=lambda members: float(
                aligned.loc[members, "aligned_peak_rt"].median()
            ),
        )
        for number, members in enumerate(groups, start=1):
            assignments.loc[members] = f"{glycan}_peak_{number}"
    result = aligned.copy()
    result["consensus_peak_id"] = assignments
    return result


def score_consensus_peak_groups(
    aligned_features: pd.DataFrame,
    sample_count: int,
    minimum_replicate_fraction: float,
) -> pd.DataFrame:
    """Score RT-matched peak groups by coverage and cross-run score stability."""
    records = []
    for (glycan, peak_id), group in aligned_features.groupby(
        ["Glycan", "consensus_peak_id"], sort=True
    ):
        median_score = float(group["peak_candidate_score"].median())
        score_mad = _median_absolute_deviation(group["peak_candidate_score"])
        rt_mad = _median_absolute_deviation(group["aligned_peak_rt"])
        coverage = float(group["sample"].nunique() / max(sample_count, 1))
        pass_fraction = float(_boolean_series(group["target_decoy_pass"]).mean())
        records.append(
            {
                "Glycan": glycan,
                "consensus_peak_id": peak_id,
                "replicate_count": int(group["sample"].nunique()),
                "replicate_coverage": coverage,
                "median_candidate_score": median_score,
                "candidate_score_mad": score_mad,
                "consensus_score": median_score - score_mad,
                "median_aligned_rt": float(group["aligned_peak_rt"].median()),
                "aligned_rt_mad": rt_mad,
                "target_decoy_pass_fraction": pass_fraction,
                "median_assignment_q_value": float(
                    pd.to_numeric(group["assignment_q_value"], errors="coerce").median()
                ),
                "total_auc": float(group["AUC"].sum()),
                "consensus_eligible": bool(
                    coverage >= minimum_replicate_fraction
                    and pass_fraction >= minimum_replicate_fraction
                ),
            }
        )
    groups = pd.DataFrame(records)
    if groups.empty:
        return groups

    groups["consensus_rank"] = np.nan
    groups["consensus_selected"] = False
    for _, indices in groups.groupby("Glycan").groups.items():
        ordered = groups.loc[indices].sort_values(
            [
                "consensus_eligible",
                "replicate_coverage",
                "consensus_score",
                "aligned_rt_mad",
                "median_aligned_rt",
            ],
            ascending=[False, False, False, True, True],
        )
        groups.loc[ordered.index, "consensus_rank"] = np.arange(
            1, len(ordered) + 1
        )
        if bool(ordered.iloc[0]["consensus_eligible"]):
            groups.loc[ordered.index[0], "consensus_selected"] = True
    return groups.sort_values(["Glycan", "consensus_rank"]).reset_index(drop=True)


def build_consensus_peak_result(
    feature_auc: pd.DataFrame,
    rt_tolerance_minutes: float = 0.3,
    minimum_replicate_fraction: float = 0.8,
    expected_samples: list[str] | None = None,
) -> ConsensusPeakResult:
    """Align feature AUCs and choose one reproducibly supported peak per glycan."""
    if not 0 < minimum_replicate_fraction <= 1:
        raise ValueError("Minimum replicate fraction must be in (0, 1]")
    observed_samples = set(feature_auc["sample"].dropna().astype(str).unique())
    if expected_samples is None:
        samples = sorted(observed_samples)
    else:
        samples = sorted({str(sample) for sample in expected_samples})
        unexpected = observed_samples - set(samples)
        if unexpected:
            raise ValueError(
                "Feature-level AUC table contains samples outside the expected "
                f"batch: {sorted(unexpected)}"
            )
    if len(samples) < 2:
        raise ValueError("Consensus peak selection requires at least two runs")
    if feature_auc.empty:
        raise ValueError("Consensus peak selection received no quantified features")

    run_peaks = build_within_run_peaks(feature_auc, rt_tolerance_minutes)
    anchors = _build_alignment_anchors(
        run_peaks,
        len(samples),
        minimum_replicate_fraction,
        rt_tolerance_minutes,
    )
    models = _fit_alignment_models(samples, anchors)
    aligned = run_peaks.merge(models, on="sample", how="left")
    aligned["aligned_peak_rt"] = (
        aligned["alignment_slope"] * aligned["raw_peak_rt"]
        + aligned["alignment_intercept"]
    )
    aligned["aligned_start_rt"] = (
        aligned["alignment_slope"] * aligned["raw_start_rt"]
        + aligned["alignment_intercept"]
    )
    aligned["aligned_end_rt"] = (
        aligned["alignment_slope"] * aligned["raw_end_rt"]
        + aligned["alignment_intercept"]
    )
    aligned = _assign_consensus_groups(aligned, rt_tolerance_minutes)
    groups = score_consensus_peak_groups(
        aligned, len(samples), minimum_replicate_fraction
    )
    aligned = aligned.merge(
        groups,
        on=["Glycan", "consensus_peak_id"],
        how="left",
        validate="many_to_one",
    )

    all_feature_auc = (
        aligned.pivot_table(
            index=["Glycan", "consensus_peak_id"],
            columns="sample",
            values="AUC",
            aggfunc="sum",
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    selected = aligned.loc[aligned["consensus_selected"]].copy()
    if selected.empty:
        selected_auc = pd.DataFrame(columns=["Glycan", *samples])
    else:
        selected_auc = (
            selected.pivot_table(
                index="Glycan",
                columns="sample",
                values="AUC",
                aggfunc="sum",
            )
            .reindex(columns=samples)
            .reset_index()
            .rename_axis(columns=None)
        )
    return ConsensusPeakResult(
        aligned_features=aligned.sort_values(
            ["Glycan", "consensus_rank", "sample"]
        ).reset_index(drop=True),
        peak_groups=groups,
        selected_auc=selected_auc,
        all_feature_auc=all_feature_auc,
        alignment_models=models,
    )


def consolidate_consensus_peak_results(
    results_root: str,
    output_csv: str,
    rt_tolerance_minutes: float = 0.3,
    minimum_replicate_fraction: float = 0.8,
    expected_samples: list[str] | None = None,
) -> ConsensusPeakResult:
    """Load per-run feature AUCs, align them, and write consensus audit files."""
    frames = []
    samples_with_feature_output = []
    if expected_samples is None:
        candidate_samples = [
            sample
            for sample in sorted(os.listdir(results_root))
            if os.path.isdir(os.path.join(results_root, sample))
        ]
    else:
        candidate_samples = sorted({str(sample) for sample in expected_samples})

    for sample in candidate_samples:
        sample_dir = os.path.join(results_root, sample)
        if not os.path.isdir(sample_dir):
            logger.warning("No output directory found for expected sample %s", sample)
            continue
        path = os.path.join(sample_dir, f"{sample}_feature_auc_values.csv")
        if not os.path.isfile(path):
            logger.warning("No feature-level AUC file found for %s", sample)
            continue
        samples_with_feature_output.append(sample)
        frame = pd.read_csv(path, low_memory=False)
        if frame.empty:
            continue
        frame["sample"] = sample
        frames.append(frame)
    coverage_samples = (
        candidate_samples
        if expected_samples is not None
        else samples_with_feature_output
    )
    if len(coverage_samples) < 2:
        raise RuntimeError(
            "At least two samples are required for consensus peak selection"
        )
    if not frames:
        raise RuntimeError("No quantified chromatographic features were found")

    result = build_consensus_peak_result(
        pd.concat(frames, ignore_index=True),
        rt_tolerance_minutes=rt_tolerance_minutes,
        minimum_replicate_fraction=minimum_replicate_fraction,
        expected_samples=coverage_samples,
    )
    result.selected_auc.to_csv(output_csv, index=False)
    result.aligned_features.to_csv(
        os.path.join(results_root, "aligned_feature_auc_values.csv"), index=False
    )
    result.peak_groups.to_csv(
        os.path.join(results_root, "consensus_peak_groups.csv"), index=False
    )
    result.all_feature_auc.to_csv(
        os.path.join(results_root, "combined_all_feature_auc_values.csv"),
        index=False,
    )
    result.alignment_models.to_csv(
        os.path.join(results_root, "retention_time_alignment.csv"), index=False
    )
    return result
