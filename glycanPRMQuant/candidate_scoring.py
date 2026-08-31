"""Evidence-based scoring of isobaric glycan-composition candidates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from glycanPRMQuant.mass_tolerance import (
    mass_accuracy_reliability,
    mass_error_ppm,
    validate_tolerance,
)


@dataclass(frozen=True)
class CandidateScoringConfig:
    """Thresholds and weights used by composition-level candidate scoring."""

    precursor_cluster_ppm: float = 20.0
    fragment_mass_tolerance: float = 0.02
    fragment_mass_tolerance_unit: str = "Da"
    fragment_mass_accuracy_reference_ppm: float = 20.0
    minimum_distinct_fragments: int = 2
    minimum_explained_intensity: float = 0.01
    minimum_candidate_score: float = 35.0
    minimum_discriminative_evidence_difference: float = 4.0
    maximum_assignment_q_value: float = 0.05
    precursor_likelihood_default_sigma_ppm: float = 1.0
    precursor_likelihood_min_sigma_ppm: float = 0.25
    precursor_calibration_min_points: int = 10
    precursor_log_likelihood_weight: float = 2.0
    mass_outlier_min_delta_ppm: float = 2.0
    mass_outlier_sigma_multiplier: float = 4.0
    minimum_coelution_transitions: int = 2
    minimum_transition_scans: int = 2
    coelution_pass_score: float = 0.70
    coisolation_score: float = 0.75
    feature_smoothing_sigma: float = 1.0
    feature_prominence_fraction: float = 0.10
    minimum_peak_flank_scans: int = 2

    def __post_init__(self):
        tolerance_value, tolerance_unit = validate_tolerance(
            self.fragment_mass_tolerance,
            self.fragment_mass_tolerance_unit,
        )
        object.__setattr__(
            self,
            "fragment_mass_tolerance",
            tolerance_value,
        )
        object.__setattr__(
            self,
            "fragment_mass_tolerance_unit",
            tolerance_unit,
        )
        if self.precursor_cluster_ppm <= 0:
            raise ValueError("Precursor cluster tolerance must be positive")
        if self.fragment_mass_accuracy_reference_ppm <= 0:
            raise ValueError("Fragment mass-accuracy reference ppm must be positive")
        if self.minimum_distinct_fragments < 1:
            raise ValueError("Minimum distinct fragments must be at least 1")
        if not 0 <= self.minimum_explained_intensity <= 1:
            raise ValueError("Minimum explained intensity must be between 0 and 1")
        if (
            self.minimum_candidate_score < 0
            or self.minimum_discriminative_evidence_difference < 0
        ):
            raise ValueError("Candidate score thresholds cannot be negative")
        if not 0 < self.maximum_assignment_q_value <= 1:
            raise ValueError("Maximum assignment q-value must be in (0, 1]")
        if (
            self.precursor_likelihood_default_sigma_ppm <= 0
            or self.precursor_likelihood_min_sigma_ppm <= 0
        ):
            raise ValueError("Precursor likelihood sigma values must be positive")
        if self.precursor_calibration_min_points < 1:
            raise ValueError("Precursor calibration requires at least one point")
        if self.precursor_log_likelihood_weight <= 0:
            raise ValueError("Precursor log-likelihood weight must be positive")
        if self.mass_outlier_min_delta_ppm < 0 or self.mass_outlier_sigma_multiplier <= 0:
            raise ValueError("Mass-outlier thresholds must be non-negative")
        if self.minimum_peak_flank_scans < 1:
            raise ValueError("Minimum peak flank scans must be at least 1")


@dataclass
class CandidateScoringResult:
    """Quantifiable rows, all reported rows, and the candidate audit table."""

    resolved_rows: pd.DataFrame
    reported_rows: pd.DataFrame
    candidate_scores: pd.DataFrame
    decoy_scores: pd.DataFrame | None = None
    target_decoy_competitions: pd.DataFrame | None = None


def _normalize_glycan(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    try:
        number = float(text)
        if number.is_integer():
            return str(int(number))
    except ValueError:
        pass
    return text[:-2] if text.endswith(".0") else text


def assign_precursor_clusters(values: pd.Series, ppm_tolerance: float) -> pd.Series:
    """Assign anchor-based precursor-m/z clusters while preserving the input index."""
    numeric = pd.to_numeric(values, errors="coerce")
    unique_mz = np.sort(numeric.dropna().unique())
    cluster_by_mz: dict[float, int] = {}
    cluster_id = -1
    anchor = None

    for mz in unique_mz:
        if anchor is None or abs(mz - anchor) > anchor * ppm_tolerance / 1e6:
            cluster_id += 1
            anchor = mz
        cluster_by_mz[float(mz)] = cluster_id

    clusters = numeric.map(cluster_by_mz)
    next_id = cluster_id + 1
    for index in clusters[clusters.isna()].index:
        clusters.loc[index] = next_id
        next_id += 1
    return clusters.astype(int)


def _build_scan_table(
    matched: pd.DataFrame,
    scan_data: pd.DataFrame | None,
    config: CandidateScoringConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build one row per MS2 scan and assign clusters consistently to both tables."""
    work = matched.copy()
    work["Glycan"] = work["Glycan"].map(_normalize_glycan)

    if scan_data is not None and {
        "scan_number", "rt", "precursor_mz", "fragment_intensity"
    }.issubset(scan_data.columns):
        raw = scan_data.copy()
        raw["fragment_intensity"] = pd.to_numeric(
            raw["fragment_intensity"], errors="coerce"
        ).fillna(0.0)
        if "precursor_intensity" not in raw.columns:
            raw["precursor_intensity"] = 0.0
        aggregation = {
            "precursor_intensity": ("precursor_intensity", "max"),
            "ms2_tic": (
                "denoised_ms2_tic"
                if "denoised_ms2_tic" in raw.columns
                else "fragment_intensity",
                "max" if "denoised_ms2_tic" in raw.columns else "sum",
            ),
        }
        if "raw_ms2_tic" in raw.columns:
            aggregation["raw_ms2_tic"] = ("raw_ms2_tic", "max")
        scans = raw.groupby(
            ["scan_number", "rt", "precursor_mz"], dropna=False
        ).agg(**aggregation).reset_index()
        if "raw_ms2_tic" not in scans.columns:
            scans["raw_ms2_tic"] = scans["ms2_tic"]
    else:
        fallback = work.copy()
        if "scan_number" not in fallback.columns:
            fallback["scan_number"] = np.arange(len(fallback), dtype=int)
        if "rt" not in fallback.columns:
            fallback["rt"] = 0.0
        if "precursor_intensity" not in fallback.columns:
            fallback["precursor_intensity"] = 0.0
        scans = (
            fallback.groupby(["scan_number", "rt", "precursor_mz"], dropna=False)
            .agg(
                precursor_intensity=("precursor_intensity", "max"),
                ms2_tic=("fragment_intensity", "max"),
            )
            .reset_index()
        )
        scans["raw_ms2_tic"] = scans["ms2_tic"]

    combined = pd.concat(
        [
            work[["precursor_mz"]].assign(_source="matched"),
            scans[["precursor_mz"]].assign(_source="scan"),
        ],
        ignore_index=True,
    )
    combined["_precursor_cluster"] = assign_precursor_clusters(
        combined["precursor_mz"], config.precursor_cluster_ppm
    )
    work["_precursor_cluster"] = combined.loc[
        combined["_source"] == "matched", "_precursor_cluster"
    ].to_numpy()
    scans["_precursor_cluster"] = combined.loc[
        combined["_source"] == "scan", "_precursor_cluster"
    ].to_numpy()
    return work, scans


def _feature_labels(signal: np.ndarray, config: CandidateScoringConfig) -> np.ndarray:
    """Split a trace at valleys between reproducible chromatographic peaks."""
    n = signal.size
    if n < 5 or not np.isfinite(signal).any():
        return np.zeros(n, dtype=int)
    y = np.nan_to_num(signal.astype(float), nan=0.0, posinf=0.0, neginf=0.0)
    if np.max(y) <= 0:
        return np.zeros(n, dtype=int)
    smooth = gaussian_filter1d(
        y, sigma=max(config.feature_smoothing_sigma, 0.0), mode="nearest"
    )
    prominence = max(float(np.max(smooth)) * config.feature_prominence_fraction, 1e-12)
    peaks, _ = find_peaks(smooth, prominence=prominence, distance=2)

    endpoint_peaks = []
    if smooth[0] > smooth[1] and smooth[0] >= prominence:
        endpoint_peaks.append(0)
    if smooth[-1] > smooth[-2] and smooth[-1] >= prominence:
        endpoint_peaks.append(n - 1)
    peaks = np.unique(np.concatenate([peaks, np.asarray(endpoint_peaks, dtype=int)]))
    peaks = peaks[smooth[peaks] >= 0.10 * float(np.max(smooth))]
    if peaks.size <= 1:
        return np.zeros(n, dtype=int)

    labels = np.zeros(n, dtype=int)
    start = 0
    for feature_id, (left_peak, right_peak) in enumerate(zip(peaks[:-1], peaks[1:])):
        valley = int(left_peak + np.argmin(smooth[left_peak:right_peak + 1]))
        labels[start:valley + 1] = feature_id
        start = valley + 1
    labels[start:] = peaks.size - 1
    return labels


def _assign_rt_features(
    work: pd.DataFrame,
    scans: pd.DataFrame,
    config: CandidateScoringConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scans = scans.copy()
    summaries = []
    feature_offset = 0

    for cluster, indices in scans.groupby("_precursor_cluster").groups.items():
        sub = scans.loc[indices].sort_values(["rt", "scan_number"])
        precursor = pd.to_numeric(sub["precursor_intensity"], errors="coerce").fillna(0.0)
        use_precursor = precursor.gt(0).sum() >= 3 and precursor.std() > 0
        signal = precursor.to_numpy() if use_precursor else sub["ms2_tic"].to_numpy(dtype=float)
        local_labels = _feature_labels(signal, config)
        global_labels = local_labels + feature_offset
        scans.loc[sub.index, "_rt_feature"] = global_labels

        for local_id in np.unique(local_labels):
            mask = local_labels == local_id
            feature_signal = np.asarray(signal)[mask]
            feature_rt = sub.loc[mask, "rt"].to_numpy(dtype=float)
            apex_index = int(np.argmax(feature_signal)) if feature_signal.size else 0
            left_flank_scans = apex_index
            right_flank_scans = max(int(feature_signal.size) - apex_index - 1, 0)
            chromatographic_peak_valid = bool(
                left_flank_scans >= config.minimum_peak_flank_scans
                and right_flank_scans >= config.minimum_peak_flank_scans
            )
            summaries.append(
                {
                    "_precursor_cluster": int(cluster),
                    "_rt_feature": int(local_id + feature_offset),
                    "feature_start_rt": float(np.min(feature_rt)),
                    "feature_apex_rt": float(feature_rt[apex_index]),
                    "feature_end_rt": float(np.max(feature_rt)),
                    "feature_reference": "precursor" if use_precursor else "ms2_tic",
                    "feature_scan_count": int(feature_signal.size),
                    "feature_apex_scan_index": apex_index,
                    "feature_left_flank_scans": left_flank_scans,
                    "feature_right_flank_scans": right_flank_scans,
                    "minimum_peak_flank_scans": config.minimum_peak_flank_scans,
                    "chromatographic_peak_valid": chromatographic_peak_valid,
                    "chromatographic_peak_rejection_reason": (
                        ""
                        if chromatographic_peak_valid
                        else "apex_lacks_required_scans_on_both_flanks"
                    ),
                }
            )
        feature_offset += int(local_labels.max()) + 1 if local_labels.size else 1

    scans["_rt_feature"] = scans["_rt_feature"].astype(int)
    scan_feature = (
        scans.sort_values("rt")
        .drop_duplicates(["_precursor_cluster", "scan_number"], keep="first")
        .set_index(["_precursor_cluster", "scan_number"])["_rt_feature"]
        .to_dict()
    )

    features = []
    for _, row in work.iterrows():
        key = (int(row["_precursor_cluster"]), row.get("scan_number"))
        feature = scan_feature.get(key)
        if feature is None:
            cluster_scans = scans.loc[scans["_precursor_cluster"] == key[0]]
            if cluster_scans.empty or "rt" not in row:
                feature = 0
            else:
                nearest = (cluster_scans["rt"] - float(row["rt"])).abs().idxmin()
                feature = int(cluster_scans.loc[nearest, "_rt_feature"])
        features.append(feature)
    work["_rt_feature"] = np.asarray(features, dtype=int)
    return work, scans, pd.DataFrame(summaries)


def _prepare_evidence(work: pd.DataFrame) -> pd.DataFrame:
    evidence = work.copy()
    if "scan_number" not in evidence.columns:
        evidence["scan_number"] = np.arange(len(evidence), dtype=int)
    if "rt" not in evidence.columns:
        evidence["rt"] = 0.0

    observed_col = (
        "observed_fragment_mz"
        if "observed_fragment_mz" in evidence.columns
        else "fragment_mz"
    )
    annotation_col = (
        "fragment_annotation"
        if "fragment_annotation" in evidence.columns
        else "Fragment"
    )
    evidence["_observed_fragment_mz"] = pd.to_numeric(
        evidence[observed_col], errors="coerce"
    )
    theoretical_col = (
        "theoretical_fragment_mz"
        if "theoretical_fragment_mz" in evidence.columns
        else "fragment_mz"
    )
    evidence["_theoretical_fragment_mz"] = pd.to_numeric(
        evidence.get(theoretical_col, evidence["_observed_fragment_mz"]),
        errors="coerce",
    )
    evidence["_event_mz"] = evidence["_observed_fragment_mz"].round(5)
    evidence["_annotation"] = evidence[annotation_col].astype(str)
    charge = evidence.get("Charge", pd.Series(0, index=evidence.index)).astype(str)
    adduct = evidence.get("Adduct", pd.Series("", index=evidence.index)).astype(str)
    evidence["_transition"] = evidence["_annotation"] + "|z" + charge + "|" + adduct

    if "mz_diff" in evidence.columns:
        evidence["_abs_mass_error_da"] = pd.to_numeric(
            evidence["mz_diff"], errors="coerce"
        ).abs()
    else:
        ppm = pd.to_numeric(
            evidence.get("ppm_error", pd.Series(np.nan, index=evidence.index)),
            errors="coerce",
        ).abs()
        evidence["_abs_mass_error_da"] = (
            ppm * evidence["_theoretical_fragment_mz"] / 1e6
        )
    evidence["_abs_mass_error_ppm"] = mass_error_ppm(
        evidence["_observed_fragment_mz"].to_numpy(dtype=float),
        evidence["_theoretical_fragment_mz"].to_numpy(dtype=float),
    )
    fallback_ppm = pd.to_numeric(
        evidence.get("ppm_error", pd.Series(np.nan, index=evidence.index)),
        errors="coerce",
    ).abs()
    evidence["_abs_mass_error_ppm"] = pd.Series(
        evidence["_abs_mass_error_ppm"], index=evidence.index
    ).fillna(fallback_ppm)
    evidence["_abs_mass_error"] = evidence["_abs_mass_error_da"]

    evidence["fragment_intensity"] = pd.to_numeric(
        evidence["fragment_intensity"], errors="coerce"
    ).fillna(0.0)
    evidence = evidence.sort_values("_abs_mass_error", na_position="last")
    evidence = evidence.drop_duplicates(
        [
            "_precursor_cluster",
            "_rt_feature",
            "Glycan",
            "scan_number",
            "_event_mz",
        ],
        keep="first",
    )

    event_group = ["_precursor_cluster", "_rt_feature", "scan_number", "_event_mz"]
    event_frequency = evidence.groupby(event_group)["Glycan"].nunique().rename("_event_frequency")
    evidence = evidence.join(event_frequency, on=event_group)
    candidate_count = (
        evidence.groupby(["_precursor_cluster", "_rt_feature"])["Glycan"]
        .nunique()
        .rename("_candidate_count")
    )
    evidence = evidence.join(candidate_count, on=["_precursor_cluster", "_rt_feature"])

    max_idf = np.log((evidence["_candidate_count"] + 1.0) / 2.0)
    idf = np.log(
        (evidence["_candidate_count"] + 1.0) / (evidence["_event_frequency"] + 1.0)
    )
    normalized = np.divide(
        idf,
        max_idf,
        out=np.ones(len(evidence), dtype=float),
        where=max_idf.to_numpy() > 0,
    )
    evidence["fragment_specificity"] = np.clip(0.20 + 0.80 * normalized, 0.20, 1.0)
    return evidence


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        return float("nan")
    values = values[valid]
    weights = weights[valid]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cutoff = 0.5 * weights.sum()
    return float(values[np.searchsorted(np.cumsum(weights), cutoff, side="left")])


def _coelution_components(
    group: pd.DataFrame,
    feature_scans: pd.DataFrame,
    config: CandidateScoringConfig,
    common_fragment_reference: np.ndarray | None = None,
) -> dict:
    scan_numbers = feature_scans.sort_values(["rt", "scan_number"])["scan_number"].tolist()
    if len(scan_numbers) < 3:
        return {
            "coelution_score": np.nan,
            "coelution_evaluable": False,
            "coeluting_transition_count": 0,
            "evaluable_transition_count": 0,
            "median_shape_correlation": np.nan,
            "median_apex_delta": np.nan,
            "coelution_reference": "unavailable",
        }

    ordered_scans = feature_scans.sort_values(["rt", "scan_number"])
    rts = ordered_scans["rt"].to_numpy(dtype=float)
    precursor = pd.to_numeric(
        ordered_scans["precursor_intensity"], errors="coerce"
    ).fillna(0.0).to_numpy()
    if np.count_nonzero(precursor) >= 3 and np.std(precursor) > 0:
        reference = precursor
        reference_name = "precursor"
    elif (
        common_fragment_reference is not None
        and np.count_nonzero(common_fragment_reference) >= 3
        and np.std(common_fragment_reference) > 0
    ):
        reference = common_fragment_reference
        reference_name = "common_fragments"
    else:
        reference = ordered_scans["ms2_tic"].to_numpy(dtype=float)
        reference_name = "ms2_tic"
    reference = np.nan_to_num(reference, nan=0.0)
    if reference.max(initial=0.0) <= 0 or np.std(reference) == 0:
        return {
            "coelution_score": np.nan,
            "coelution_evaluable": False,
            "coeluting_transition_count": 0,
            "evaluable_transition_count": 0,
            "median_shape_correlation": np.nan,
            "median_apex_delta": np.nan,
            "coelution_reference": reference_name,
        }

    reference = gaussian_filter1d(
        reference, sigma=max(config.feature_smoothing_sigma, 0.0), mode="nearest"
    )
    ref_norm = reference / reference.max()
    ref_apex = int(np.argmax(reference))
    positive_steps = np.diff(rts)
    positive_steps = positive_steps[positive_steps > 0]
    scan_spacing = float(np.median(positive_steps)) if positive_steps.size else 0.0
    half_height = np.flatnonzero(ref_norm >= 0.5)
    if half_height.size >= 2:
        peak_width = float(rts[half_height[-1]] - rts[half_height[0]])
    else:
        peak_width = float(rts[-1] - rts[0]) if len(rts) > 1 else 0.0
    allowed_apex_delta = max(2.0 * scan_spacing, 0.20 * peak_width, 1e-9)

    transition_rows = []
    for transition, transition_group in group.groupby("_transition"):
        trace = (
            transition_group.groupby("scan_number")["fragment_intensity"]
            .max()
            .reindex(scan_numbers, fill_value=0.0)
            .to_numpy(dtype=float)
        )
        if np.count_nonzero(trace) < config.minimum_transition_scans or trace.max(initial=0.0) <= 0:
            continue
        trace = gaussian_filter1d(
            trace, sigma=max(config.feature_smoothing_sigma, 0.0), mode="nearest"
        )
        if np.std(trace) == 0:
            continue
        trace_norm = trace / trace.max()
        shape = float(np.corrcoef(np.log1p(trace), np.log1p(reference))[0, 1])
        shape = float(np.clip(np.nan_to_num(shape, nan=0.0), 0.0, 1.0))
        apex_delta = abs(float(rts[int(np.argmax(trace))] - rts[ref_apex]))
        apex_score = float(np.exp(-0.5 * (apex_delta / allowed_apex_delta) ** 2))
        union = np.maximum(trace_norm, ref_norm).sum()
        overlap = float(np.minimum(trace_norm, ref_norm).sum() / union) if union > 0 else 0.0
        score = 0.50 * shape + 0.30 * apex_score + 0.20 * overlap
        transition_rows.append(
            {
                "transition": transition,
                "score": score,
                "shape": shape,
                "apex_delta": apex_delta,
                "weight": float(
                    transition_group["fragment_specificity"].mean()
                    * np.sqrt(transition_group["fragment_intensity"].sum())
                ),
            }
        )

    if len(transition_rows) < config.minimum_coelution_transitions:
        return {
            "coelution_score": np.nan,
            "coelution_evaluable": False,
            "coeluting_transition_count": sum(
                row["score"] >= config.coelution_pass_score for row in transition_rows
            ),
            "evaluable_transition_count": len(transition_rows),
            "median_shape_correlation": np.nan,
            "median_apex_delta": np.nan,
            "coelution_reference": reference_name,
        }

    scores = np.asarray([row["score"] for row in transition_rows], dtype=float)
    shapes = np.asarray([row["shape"] for row in transition_rows], dtype=float)
    apex_deltas = np.asarray([row["apex_delta"] for row in transition_rows], dtype=float)
    weights = np.asarray([row["weight"] for row in transition_rows], dtype=float)
    return {
        "coelution_score": _weighted_median(scores, weights),
        "coelution_evaluable": True,
        "coeluting_transition_count": int((scores >= config.coelution_pass_score).sum()),
        "evaluable_transition_count": len(transition_rows),
        "median_shape_correlation": _weighted_median(shapes, weights),
        "median_apex_delta": _weighted_median(apex_deltas, weights),
        "coelution_reference": reference_name,
    }


def _score_candidates(
    evidence: pd.DataFrame,
    scans: pd.DataFrame,
    feature_summary: pd.DataFrame,
    config: CandidateScoringConfig,
    precursor_calibration: dict | None = None,
    fragment_reliability_by_feature: dict[tuple[int, int], dict] | None = None,
) -> pd.DataFrame:
    rows = []
    keys = ["_precursor_cluster", "_rt_feature"]
    feature_lookup = feature_summary.set_index(keys).to_dict("index") if not feature_summary.empty else {}

    for (cluster, feature), conflict in evidence.groupby(keys):
        feature_scans = scans.loc[
            (scans["_precursor_cluster"] == cluster)
            & (scans["_rt_feature"] == feature)
        ]
        denominator = float(feature_scans["ms2_tic"].sum())
        if denominator <= 0:
            denominator = float(
                conflict.drop_duplicates(["scan_number", "_event_mz"])["fragment_intensity"].sum()
            )

        ordered_scan_numbers = feature_scans.sort_values(
            ["rt", "scan_number"]
        )["scan_number"].tolist()
        number_of_candidates = int(conflict["Glycan"].nunique())
        supplied_reliability = (fragment_reliability_by_feature or {}).get(
            (int(cluster), int(feature))
        )
        if supplied_reliability is None:
            feature_events = conflict.drop_duplicates(["scan_number", "_event_mz"])
            representative_fragment_mz = float(
                feature_events["_observed_fragment_mz"].median()
            )
            equivalent_ppm, accuracy_reliability = mass_accuracy_reliability(
                config.fragment_mass_tolerance,
                config.fragment_mass_tolerance_unit,
                representative_fragment_mz,
                config.fragment_mass_accuracy_reference_ppm,
            )
        else:
            representative_fragment_mz = float(
                supplied_reliability["fragment_tolerance_representative_mz"]
            )
            equivalent_ppm = float(
                supplied_reliability["fragment_tolerance_equivalent_ppm"]
            )
            accuracy_reliability = float(
                supplied_reliability["fragment_mass_accuracy_reliability"]
            )
        common_events = conflict.loc[
            conflict["_event_frequency"] >= number_of_candidates
        ]
        if common_events.empty:
            common_reference = None
        else:
            common_reference = (
                common_events.groupby(["scan_number", "_event_mz"])["fragment_intensity"]
                .max()
                .groupby("scan_number")
                .sum()
                .reindex(ordered_scan_numbers, fill_value=0.0)
                .to_numpy(dtype=float)
            )

        for glycan, group in conflict.groupby("Glycan"):
            transition_specificity = group.groupby("_transition")["fragment_specificity"].mean()
            distinct_fragments = int(transition_specificity.size)
            candidate_specific = int((transition_specificity >= 0.95).sum())
            effective_fragments = float(transition_specificity.sum())
            weighted_intensity = float(
                (group["fragment_intensity"] * group["fragment_specificity"]).sum()
            )
            explained_fraction = (
                float(np.clip(weighted_intensity / denominator, 0.0, 1.0))
                if denominator > 0
                else 0.0
            )
            explained_component = float(np.clip(explained_fraction / 0.25, 0.0, 1.0))
            fragment_support = float(1.0 - np.exp(-effective_fragments / 2.0))

            mass_error_column = (
                "_abs_mass_error_ppm"
                if config.fragment_mass_tolerance_unit == "ppm"
                else "_abs_mass_error_da"
            )
            normalized_error = (
                pd.to_numeric(group[mass_error_column], errors="coerce")
                / max(config.fragment_mass_tolerance, 1e-12)
            ).to_numpy(dtype=float)
            mass_weights = np.sqrt(group["fragment_intensity"].to_numpy(dtype=float))
            median_normalized_error = _weighted_median(normalized_error, mass_weights)
            median_mass_error_da = _weighted_median(
                pd.to_numeric(
                    group["_abs_mass_error_da"], errors="coerce"
                ).to_numpy(dtype=float),
                mass_weights,
            )
            median_mass_error_ppm = _weighted_median(
                pd.to_numeric(
                    group["_abs_mass_error_ppm"], errors="coerce"
                ).to_numpy(dtype=float),
                mass_weights,
            )
            fragment_accuracy = (
                float(np.exp(-0.5 * median_normalized_error ** 2))
                if np.isfinite(median_normalized_error)
                else 0.5
            )

            signed_precursor_error = pd.to_numeric(
                group.get("precursor_ppm_error", pd.Series(np.nan, index=group.index)),
                errors="coerce",
            ).to_numpy(dtype=float)
            median_signed_precursor_error = _weighted_median(
                signed_precursor_error, mass_weights
            )
            median_absolute_precursor_error = _weighted_median(
                np.abs(signed_precursor_error), mass_weights
            )

            specific_transitions = group.loc[group["fragment_specificity"] >= 0.50]
            if (
                specific_transitions["_transition"].nunique()
                >= config.minimum_coelution_transitions
            ):
                coelution_group = specific_transitions
            else:
                coelution_group = group
            coelution = _coelution_components(
                coelution_group,
                feature_scans,
                config,
                common_fragment_reference=common_reference,
            )
            feature_info = feature_lookup.get((cluster, feature), {})
            rows.append(
                {
                    "_precursor_cluster": int(cluster),
                    "_rt_feature": int(feature),
                    "Glycan": glycan,
                    "explained_intensity_fraction": explained_fraction,
                    "explained_intensity_score": explained_component,
                    "specificity_weighted_intensity": weighted_intensity,
                    "distinct_fragment_count": distinct_fragments,
                    "candidate_specific_fragment_count": candidate_specific,
                    "effective_fragment_count": effective_fragments,
                    "fragment_support_score": fragment_support,
                    "median_fragment_mass_error": median_mass_error_da,
                    "median_fragment_mass_error_ppm": median_mass_error_ppm,
                    "fragment_mass_accuracy_score": fragment_accuracy,
                    "fragment_mass_tolerance_value": config.fragment_mass_tolerance,
                    "fragment_mass_tolerance_unit": config.fragment_mass_tolerance_unit,
                    "fragment_tolerance_representative_mz": representative_fragment_mz,
                    "fragment_tolerance_equivalent_ppm": equivalent_ppm,
                    "fragment_mass_accuracy_reliability": accuracy_reliability,
                    "median_signed_precursor_ppm_error": median_signed_precursor_error,
                    "median_precursor_ppm_error": median_absolute_precursor_error,
                    **coelution,
                    **feature_info,
                }
            )
    return _finalize_candidate_scores(
        pd.DataFrame(rows), config, precursor_calibration=precursor_calibration
    )


def _finalize_candidate_scores(
    scores: pd.DataFrame,
    config: CandidateScoringConfig,
    precursor_calibration: dict | None = None,
) -> pd.DataFrame:
    """Calibrate precursor evidence and make optional terms feature-comparable."""
    if scores.empty:
        return scores
    finalized = scores.copy()
    keys = ["_precursor_cluster", "_rt_feature"]

    # Hand-built score tables used by programmatic callers may predate the
    # chromatographic gate. Production scores always receive these fields from
    # _assign_rt_features; missing fields remain backward-compatible.
    if "chromatographic_peak_valid" not in finalized.columns:
        finalized["chromatographic_peak_valid"] = True
    if "chromatographic_peak_rejection_reason" not in finalized.columns:
        finalized["chromatographic_peak_rejection_reason"] = ""

    if "fragment_mass_accuracy_reliability" not in finalized.columns:
        representative_mz = 500.0
        equivalent_ppm, reliability = mass_accuracy_reliability(
            config.fragment_mass_tolerance,
            config.fragment_mass_tolerance_unit,
            representative_mz,
            config.fragment_mass_accuracy_reference_ppm,
        )
        finalized["fragment_mass_tolerance_value"] = config.fragment_mass_tolerance
        finalized["fragment_mass_tolerance_unit"] = config.fragment_mass_tolerance_unit
        finalized["fragment_tolerance_representative_mz"] = representative_mz
        finalized["fragment_tolerance_equivalent_ppm"] = equivalent_ppm
        finalized["fragment_mass_accuracy_reliability"] = reliability

    # Use one value per precursor-cluster/candidate pair so long RT features do
    # not receive more influence on calibration than short ones. Iteratively
    # choose the candidate closest to the run-level error center in each
    # cluster, then estimate that center and spread robustly.
    if precursor_calibration is not None:
        calibration_center = float(precursor_calibration["center_ppm"])
        sigma = float(precursor_calibration["sigma_ppm"])
        calibration_point_count = int(precursor_calibration["points"])
        calibration_source = str(precursor_calibration["source"])
    else:
        calibration_candidates = (
            finalized.groupby(["_precursor_cluster", "Glycan"], as_index=False)
            ["median_signed_precursor_ppm_error"]
            .median()
            .dropna(subset=["median_signed_precursor_ppm_error"])
        )
        calibration_center = 0.0
        calibration_values = np.asarray([], dtype=float)
        if not calibration_candidates.empty:
            for _ in range(2):
                distance = (
                    calibration_candidates["median_signed_precursor_ppm_error"]
                    - calibration_center
                ).abs()
                closest_indices = distance.groupby(
                    calibration_candidates["_precursor_cluster"]
                ).idxmin()
                calibration_values = calibration_candidates.loc[
                    closest_indices, "median_signed_precursor_ppm_error"
                ].to_numpy(dtype=float)
                calibration_center = float(np.median(calibration_values))

        calibration_point_count = int(calibration_values.size)
        if calibration_values.size >= config.precursor_calibration_min_points:
            q1, q3 = np.quantile(calibration_values, [0.25, 0.75])
            empirical_sigma = float((q3 - q1) / 1.349)
            if not np.isfinite(empirical_sigma) or empirical_sigma <= 0:
                empirical_sigma = float(
                    1.4826
                    * np.median(np.abs(calibration_values - calibration_center))
                )
            sigma = max(config.precursor_likelihood_min_sigma_ppm, empirical_sigma)
            calibration_source = "empirical_iqr"
        else:
            calibration_center = 0.0
            sigma = config.precursor_likelihood_default_sigma_ppm
            calibration_source = "configured_fallback"

    calibrated_error = (
        finalized["median_signed_precursor_ppm_error"] - calibration_center
    ).abs()
    # Retain compatibility with inputs that contain only absolute precursor
    # errors. Production matches contain the signed error used above.
    calibrated_error = calibrated_error.fillna(
        finalized["median_precursor_ppm_error"]
    )
    best_error = calibrated_error.groupby(
        [finalized[key] for key in keys]
    ).transform("min")
    delta = calibrated_error - best_error
    valid_delta = calibrated_error.notna() & best_error.notna()
    absolute_log_likelihood = pd.Series(0.0, index=finalized.index)
    absolute_log_likelihood.loc[valid_delta] = -0.5 * (
        calibrated_error.loc[valid_delta] / max(sigma, 1e-12)
    ) ** 2
    best_log_likelihood = absolute_log_likelihood.groupby(
        [finalized[key] for key in keys]
    ).transform("max")
    log_likelihood = absolute_log_likelihood - best_log_likelihood
    relative_likelihood = pd.Series(0.5, index=finalized.index)
    relative_likelihood.loc[valid_delta] = np.exp(
        np.clip(log_likelihood.loc[valid_delta], -50.0, 0.0)
    )

    finalized["best_precursor_ppm_error"] = finalized.groupby(keys)[
        "median_precursor_ppm_error"
    ].transform("min")
    finalized["best_calibrated_precursor_error_ppm"] = best_error
    finalized["calibrated_precursor_error_ppm"] = calibrated_error
    finalized["precursor_error_delta_ppm"] = delta
    finalized["precursor_absolute_log_likelihood"] = absolute_log_likelihood
    finalized["precursor_log_likelihood"] = log_likelihood
    finalized["precursor_relative_likelihood"] = relative_likelihood
    finalized["precursor_mass_accuracy_score"] = relative_likelihood
    finalized["precursor_likelihood_center_ppm"] = calibration_center
    finalized["precursor_likelihood_sigma_ppm"] = sigma
    finalized["precursor_calibration_points"] = calibration_point_count
    finalized["precursor_calibration_source"] = calibration_source

    candidate_count = finalized.groupby(keys)["Glycan"].transform("size")
    all_coelution_evaluable = finalized.groupby(keys)["coelution_evaluable"].transform("all")
    coelution_comparable = (candidate_count == 1) | all_coelution_evaluable
    finalized["coelution_comparable"] = coelution_comparable
    finalized["coelution_component_used"] = np.where(
        coelution_comparable & finalized["coelution_evaluable"],
        finalized["coelution_score"],
        0.5,
    )

    raw_mass_accuracy_weight = (
        0.15 * finalized["fragment_mass_accuracy_reliability"]
    )
    bounded_weight_total = 0.85 + raw_mass_accuracy_weight
    finalized["effective_fragment_mass_accuracy_weight"] = (
        raw_mass_accuracy_weight / bounded_weight_total
    )
    finalized["candidate_score"] = 100.0 * (
        0.35 * finalized["explained_intensity_score"]
        + 0.20 * finalized["fragment_support_score"]
        + raw_mass_accuracy_weight * finalized["fragment_mass_accuracy_score"]
        + 0.10 * finalized["precursor_relative_likelihood"]
        + 0.20 * finalized["coelution_component_used"]
    ) / bounded_weight_total

    # The discriminative score is centered within each feature. Equal/common
    # components therefore cancel exactly. Precursor evidence enters as a
    # calibrated Gaussian log-likelihood instead of a flat tolerance score.
    discriminative_mass_weight = (
        15.0 * finalized["fragment_mass_accuracy_reliability"]
    )
    discriminative_scale = 90.0 / (75.0 + discriminative_mass_weight)
    finalized["effective_fragment_mass_accuracy_discriminative_weight"] = (
        discriminative_mass_weight * discriminative_scale
    )
    finalized["raw_discriminative_evidence"] = (
        discriminative_scale
        * (
            35.0 * finalized["explained_intensity_score"]
            + 20.0 * finalized["fragment_support_score"]
            + discriminative_mass_weight
            * finalized["fragment_mass_accuracy_score"]
            + 20.0 * finalized["coelution_component_used"]
        )
        + config.precursor_log_likelihood_weight
        * finalized["precursor_log_likelihood"]
    )
    feature_mean = finalized.groupby(keys)["raw_discriminative_evidence"].transform("mean")
    finalized["discriminative_score"] = (
        finalized["raw_discriminative_evidence"] - feature_mean
    )

    outlier_threshold = max(
        config.mass_outlier_min_delta_ppm,
        config.mass_outlier_sigma_multiplier * sigma,
    )
    finalized["mass_outlier_threshold_ppm"] = outlier_threshold
    finalized["mass_outlier_pruned"] = (
        (candidate_count > 1)
        & (finalized["candidate_specific_fragment_count"] == 0)
        & valid_delta
        & (delta > outlier_threshold)
    )
    invalid_chromatographic_peak = ~finalized[
        "chromatographic_peak_valid"
    ].fillna(False).astype(bool)
    finalized["candidate_rejection_reason"] = np.select(
        [
            invalid_chromatographic_peak,
            finalized["mass_outlier_pruned"],
        ],
        [
            "no_interior_chromatographic_apex",
            "no_specific_evidence_and_precursor_mass_outlier",
        ],
        default="",
    )
    finalized["eligible_for_resolution"] = ~(
        finalized["mass_outlier_pruned"] | invalid_chromatographic_peak
    )
    return finalized


def _apply_resolution(
    scores: pd.DataFrame,
    config: CandidateScoringConfig,
    resolve: bool,
) -> pd.DataFrame:
    scored = scores.copy()
    scored["candidate_rank"] = np.nan
    scored["top_candidate"] = scored["Glycan"]
    scored["top_score"] = scored["candidate_score"]
    scored["top_discriminative_score"] = scored["discriminative_score"]
    scored["runner_up_candidate"] = None
    scored["runner_up_score"] = np.nan
    scored["runner_up_discriminative_score"] = np.nan
    scored["score_margin"] = np.nan
    scored["discriminative_evidence_difference"] = np.nan
    scored["resolution_status"] = "uncontested"
    scored["reported"] = True
    scored["selected"] = True

    keys = ["_precursor_cluster", "_rt_feature"]
    for _, indices in scored.groupby(keys).groups.items():
        feature = scored.loc[indices]
        active = feature.loc[feature["eligible_for_resolution"]]
        pruned = feature.loc[~feature["eligible_for_resolution"]]
        active_ordered = active.sort_values(
            ["discriminative_score", "Glycan"], ascending=[False, True]
        )
        pruned_ordered = pruned.sort_values(
            ["discriminative_score", "Glycan"], ascending=[False, True]
        )
        ordered = pd.concat([active_ordered, pruned_ordered])
        scored.loc[ordered.index, "candidate_rank"] = np.arange(1, len(ordered) + 1)
        if not bool(feature["chromatographic_peak_valid"].fillna(False).all()):
            scored.loc[ordered.index, "resolution_status"] = (
                "no_chromatographic_peak"
            )
            scored.loc[ordered.index, "selected"] = False
            continue
        if len(ordered) == 1:
            only_index = ordered.index[0]
            only = scored.loc[only_index]
            enough_evidence = (
                only["distinct_fragment_count"] >= config.minimum_distinct_fragments
                and only["explained_intensity_fraction"]
                >= config.minimum_explained_intensity
                and only["candidate_score"] >= config.minimum_candidate_score
            )
            if not enough_evidence:
                scored.loc[only_index, "resolution_status"] = "insufficient_evidence"
                scored.loc[only_index, "selected"] = False
            continue
        # Contested candidates are reportable, but none is eligible for
        # composition-level quantification unless the conflict is resolved.
        scored.loc[ordered.index, "selected"] = False

        if not resolve:
            scored.loc[ordered.index, "resolution_status"] = "resolution_disabled"
            continue

        if active_ordered.empty:
            scored.loc[ordered.index, "resolution_status"] = "insufficient_evidence"
            continue

        top_index = active_ordered.index[0]
        top_score = float(active_ordered.iloc[0]["candidate_score"])
        top_discriminative = float(active_ordered.iloc[0]["discriminative_score"])
        scored.loc[ordered.index, "top_candidate"] = active_ordered.iloc[0]["Glycan"]
        scored.loc[ordered.index, "top_score"] = top_score
        scored.loc[ordered.index, "top_discriminative_score"] = top_discriminative

        runner = None
        if len(active_ordered) > 1:
            runner_index = active_ordered.index[1]
            runner = scored.loc[runner_index]
            runner_score = float(active_ordered.iloc[1]["candidate_score"])
            runner_discriminative = float(
                active_ordered.iloc[1]["discriminative_score"]
            )
            evidence_difference = top_discriminative - runner_discriminative
            scored.loc[ordered.index, "runner_up_candidate"] = active_ordered.iloc[1][
                "Glycan"
            ]
            scored.loc[ordered.index, "runner_up_score"] = runner_score
            scored.loc[ordered.index, "runner_up_discriminative_score"] = (
                runner_discriminative
            )
            scored.loc[ordered.index, "score_margin"] = top_score - runner_score
            scored.loc[
                ordered.index, "discriminative_evidence_difference"
            ] = evidence_difference

        if feature["candidate_specific_fragment_count"].max() == 0:
            scored.loc[
                ordered.index, "resolution_status"
            ] = "no_discriminating_fragment_evidence"
            continue

        top = scored.loc[top_index]
        enough_evidence = (
            top["distinct_fragment_count"] >= config.minimum_distinct_fragments
            and top["explained_intensity_fraction"] >= config.minimum_explained_intensity
            and top_score >= config.minimum_candidate_score
        )
        possible_coisolation = (
            runner is not None
            and bool(top["coelution_comparable"])
            and bool(top["coelution_evaluable"])
            and bool(runner["coelution_evaluable"])
            and top["candidate_specific_fragment_count"] > 0
            and runner["candidate_specific_fragment_count"] > 0
            and top["coelution_score"] >= config.coisolation_score
            and runner["coelution_score"] >= config.coisolation_score
        )

        if not enough_evidence:
            status = "insufficient_evidence"
        elif runner is None:
            status = "resolved_after_mass_pruning"
            scored.loc[top_index, "selected"] = True
        elif possible_coisolation:
            status = "possible_coisolation"
        elif (
            evidence_difference
            < config.minimum_discriminative_evidence_difference
        ):
            status = "ambiguous"
        else:
            status = "resolved"
            scored.loc[top_index, "selected"] = True
        scored.loc[ordered.index, "resolution_status"] = status
    scored["quantification_weight"] = scored["selected"].astype(float)
    return scored


def _apply_target_decoy_statistics(
    target_scores: pd.DataFrame,
    decoy_scores: pd.DataFrame | None,
    config: CandidateScoringConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach feature-level target-decoy competition and monotone q-values."""
    scored = target_scores.copy()
    scored["target_decoy_evaluated"] = decoy_scores is not None
    scored["target_assignment_score"] = np.nan
    scored["best_decoy_score"] = np.nan
    scored["target_decoy_score_margin"] = np.nan
    scored["target_decoy_winner"] = True
    scored["assignment_fdr"] = np.nan
    scored["assignment_q_value"] = np.nan
    scored["target_decoy_pass"] = True
    if decoy_scores is None or scored.empty:
        return scored, pd.DataFrame()

    keys = ["_precursor_cluster", "_rt_feature"]
    target_best = (
        scored.sort_values(keys + ["candidate_rank"])
        .drop_duplicates(keys, keep="first")
        .set_index(keys)["candidate_score"]
    )
    if decoy_scores.empty:
        decoy_best = pd.Series(dtype=float, name="candidate_score")
        decoy_best.index = pd.MultiIndex.from_arrays([[], []], names=keys)
    else:
        decoy_best = (
            decoy_scores.sort_values(keys + ["candidate_rank"])
            .drop_duplicates(keys, keep="first")
            .set_index(keys)["candidate_score"]
        )

    feature_index = target_best.index.union(decoy_best.index)
    competition = pd.DataFrame(index=feature_index)
    competition["target_score"] = target_best.reindex(feature_index).fillna(0.0)
    competition["decoy_score"] = decoy_best.reindex(feature_index).fillna(0.0)
    competition["target_winner"] = (
        competition["target_score"] >= competition["decoy_score"]
    )
    competition["winner_score"] = competition[
        ["target_score", "decoy_score"]
    ].max(axis=1)

    target_thresholds = np.sort(target_best.dropna().unique())[::-1]
    fdr_by_threshold: dict[float, float] = {}
    for threshold in target_thresholds:
        passing = competition["winner_score"] >= threshold
        target_count = int(
            (passing & competition["target_winner"]).sum()
        )
        decoy_count = int(
            (passing & ~competition["target_winner"]).sum()
        )
        fdr_by_threshold[float(threshold)] = min(
            1.0, (decoy_count + 1.0) / max(target_count, 1)
        )

    q_by_threshold: dict[float, float] = {}
    running_min = 1.0
    for threshold in target_thresholds[::-1]:
        running_min = min(running_min, fdr_by_threshold[float(threshold)])
        q_by_threshold[float(threshold)] = running_min

    target_feature_stats = pd.DataFrame(index=target_best.index)
    target_feature_stats["target_assignment_score"] = target_best
    target_feature_stats["best_decoy_score"] = decoy_best.reindex(
        target_best.index
    ).fillna(0.0)
    target_feature_stats["target_decoy_score_margin"] = (
        target_feature_stats["target_assignment_score"]
        - target_feature_stats["best_decoy_score"]
    )
    target_feature_stats["target_decoy_winner"] = (
        target_feature_stats["target_assignment_score"]
        >= target_feature_stats["best_decoy_score"]
    )
    target_feature_stats["assignment_fdr"] = target_feature_stats[
        "target_assignment_score"
    ].map(fdr_by_threshold)
    target_feature_stats["assignment_q_value"] = target_feature_stats[
        "target_assignment_score"
    ].map(q_by_threshold)
    target_feature_stats["target_decoy_pass"] = (
        target_feature_stats["target_decoy_winner"]
        & (
            target_feature_stats["assignment_q_value"]
            <= config.maximum_assignment_q_value
        )
    )

    competition = competition.join(
        target_feature_stats[
            ["assignment_fdr", "assignment_q_value", "target_decoy_pass"]
        ],
        how="left",
    )
    competition["competition_winner"] = np.where(
        competition["target_winner"], "target", "decoy"
    )

    for feature_key, stats in target_feature_stats.iterrows():
        feature_mask = (
            (scored["_precursor_cluster"] == feature_key[0])
            & (scored["_rt_feature"] == feature_key[1])
        )
        for column, value in stats.items():
            scored.loc[feature_mask, column] = value

        selected_mask = feature_mask & scored["selected"]
        if selected_mask.any() and not bool(stats["target_decoy_pass"]):
            scored.loc[selected_mask, "selected"] = False
            status = (
                "target_decoy_q_failed"
                if bool(stats["target_decoy_winner"])
                else "decoy_outcompeted"
            )
            scored.loc[feature_mask, "resolution_status"] = status

    scored["target_decoy_evaluated"] = True
    scored["quantification_weight"] = scored["selected"].astype(float)
    return scored, competition.reset_index()


def score_and_resolve_candidates(
    matched_data: pd.DataFrame,
    scan_data: pd.DataFrame | None = None,
    config: CandidateScoringConfig | None = None,
    resolve: bool = True,
    decoy_matched_data: pd.DataFrame | None = None,
) -> CandidateScoringResult:
    """Score composition candidates and resolve only well-separated conflicts."""
    config = config or CandidateScoringConfig()
    if matched_data.empty:
        return CandidateScoringResult(
            matched_data.copy(), matched_data.copy(), pd.DataFrame()
        )

    required = {"precursor_mz", "Glycan", "fragment_intensity"}
    missing = required - set(matched_data.columns)
    if missing:
        raise ValueError(f"Candidate scoring requires columns: {sorted(missing)}")

    work, scans = _build_scan_table(matched_data, scan_data, config)
    work, scans, feature_summary = _assign_rt_features(work, scans, config)
    evidence = _prepare_evidence(work)
    scores = _score_candidates(evidence, scans, feature_summary, config)
    tolerance_audit_columns = [
        "fragment_tolerance_representative_mz",
        "fragment_tolerance_equivalent_ppm",
        "fragment_mass_accuracy_reliability",
    ]
    target_fragment_reliability = (
        scores.drop_duplicates(["_precursor_cluster", "_rt_feature"])
        .set_index(["_precursor_cluster", "_rt_feature"])[tolerance_audit_columns]
        .to_dict("index")
    )
    target_precursor_calibration = {
        "center_ppm": scores["precursor_likelihood_center_ppm"].iloc[0],
        "sigma_ppm": scores["precursor_likelihood_sigma_ppm"].iloc[0],
        "points": scores["precursor_calibration_points"].iloc[0],
        "source": scores["precursor_calibration_source"].iloc[0],
    }
    scores = _apply_resolution(scores, config, resolve=resolve)
    decoy_scores = None
    if decoy_matched_data is not None:
        if decoy_matched_data.empty:
            decoy_scores = pd.DataFrame()
        else:
            decoy_work, decoy_scans = _build_scan_table(
                decoy_matched_data, scan_data, config
            )
            decoy_work, decoy_scans, decoy_feature_summary = _assign_rt_features(
                decoy_work, decoy_scans, config
            )
            decoy_evidence = _prepare_evidence(decoy_work)
            decoy_scores = _score_candidates(
                decoy_evidence,
                decoy_scans,
                decoy_feature_summary,
                config,
                precursor_calibration=target_precursor_calibration,
                fragment_reliability_by_feature=target_fragment_reliability,
            )
            decoy_scores = _apply_resolution(
                decoy_scores, config, resolve=False
            )
    scores, target_decoy_competitions = _apply_target_decoy_statistics(
        scores, decoy_scores, config
    )

    merge_columns = [
        "_precursor_cluster",
        "_rt_feature",
        "Glycan",
        "candidate_score",
        "candidate_rank",
        "top_candidate",
        "top_score",
        "top_discriminative_score",
        "runner_up_candidate",
        "runner_up_score",
        "runner_up_discriminative_score",
        "score_margin",
        "discriminative_score",
        "discriminative_evidence_difference",
        "resolution_status",
        "reported",
        "selected",
        "quantification_weight",
        "target_decoy_evaluated",
        "target_assignment_score",
        "best_decoy_score",
        "target_decoy_score_margin",
        "target_decoy_winner",
        "assignment_fdr",
        "assignment_q_value",
        "target_decoy_pass",
        "eligible_for_resolution",
        "mass_outlier_pruned",
        "candidate_rejection_reason",
        "mass_outlier_threshold_ppm",
        "explained_intensity_fraction",
        "distinct_fragment_count",
        "candidate_specific_fragment_count",
        "fragment_mass_accuracy_score",
        "median_fragment_mass_error",
        "median_fragment_mass_error_ppm",
        "fragment_mass_tolerance_value",
        "fragment_mass_tolerance_unit",
        "fragment_tolerance_representative_mz",
        "fragment_tolerance_equivalent_ppm",
        "fragment_mass_accuracy_reliability",
        "effective_fragment_mass_accuracy_weight",
        "effective_fragment_mass_accuracy_discriminative_weight",
        "median_signed_precursor_ppm_error",
        "median_precursor_ppm_error",
        "best_precursor_ppm_error",
        "best_calibrated_precursor_error_ppm",
        "calibrated_precursor_error_ppm",
        "precursor_error_delta_ppm",
        "precursor_absolute_log_likelihood",
        "precursor_log_likelihood",
        "precursor_relative_likelihood",
        "precursor_mass_accuracy_score",
        "precursor_likelihood_center_ppm",
        "precursor_likelihood_sigma_ppm",
        "precursor_calibration_points",
        "precursor_calibration_source",
        "coelution_score",
        "coelution_evaluable",
        "coelution_comparable",
        "coelution_component_used",
        "coeluting_transition_count",
        "evaluable_transition_count",
        "coelution_reference",
        "feature_start_rt",
        "feature_apex_rt",
        "feature_end_rt",
        "feature_reference",
        "feature_scan_count",
        "feature_apex_scan_index",
        "feature_left_flank_scans",
        "feature_right_flank_scans",
        "minimum_peak_flank_scans",
        "chromatographic_peak_valid",
        "chromatographic_peak_rejection_reason",
    ]
    annotated = work.merge(
        scores[merge_columns],
        on=["_precursor_cluster", "_rt_feature", "Glycan"],
        how="left",
    )
    reported = annotated.copy()
    resolved = annotated.loc[annotated["selected"]].copy()
    output_names = {
        "_precursor_cluster": "precursor_cluster",
        "_rt_feature": "rt_feature_id",
    }
    reported = reported.rename(columns=output_names)
    resolved = resolved.rename(columns=output_names)
    scores = scores.rename(
        columns={
            "_precursor_cluster": "precursor_cluster",
            "_rt_feature": "rt_feature_id",
        }
    ).sort_values(["precursor_cluster", "rt_feature_id", "candidate_rank"])
    if decoy_scores is None or decoy_scores.empty:
        decoy_scores_output = pd.DataFrame()
    else:
        decoy_scores_output = decoy_scores.rename(
            columns={
                "_precursor_cluster": "precursor_cluster",
                "_rt_feature": "rt_feature_id",
            }
        ).sort_values(["precursor_cluster", "rt_feature_id", "candidate_rank"])
    if not target_decoy_competitions.empty:
        target_decoy_competitions = target_decoy_competitions.rename(
            columns={
                "_precursor_cluster": "precursor_cluster",
                "_rt_feature": "rt_feature_id",
            }
        ).sort_values(["precursor_cluster", "rt_feature_id"])
    return CandidateScoringResult(
        resolved.reset_index(drop=True),
        reported.reset_index(drop=True),
        scores.reset_index(drop=True),
        decoy_scores_output.reset_index(drop=True),
        target_decoy_competitions.reset_index(drop=True),
    )
