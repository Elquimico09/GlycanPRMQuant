"""Evidence-based scoring of isobaric glycan-composition candidates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


@dataclass(frozen=True)
class CandidateScoringConfig:
    """Thresholds and weights used by composition-level candidate scoring."""

    precursor_cluster_ppm: float = 20.0
    fragment_mass_tolerance: float = 0.02
    minimum_distinct_fragments: int = 2
    minimum_explained_intensity: float = 0.005
    minimum_candidate_score: float = 35.0
    minimum_discriminative_evidence_difference: float = 4.0
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

    def __post_init__(self):
        if self.precursor_cluster_ppm <= 0:
            raise ValueError("Precursor cluster tolerance must be positive")
        if self.fragment_mass_tolerance <= 0:
            raise ValueError("Fragment mass tolerance must be positive")
        if self.minimum_distinct_fragments < 1:
            raise ValueError("Minimum distinct fragments must be at least 1")
        if not 0 <= self.minimum_explained_intensity <= 1:
            raise ValueError("Minimum explained intensity must be between 0 and 1")
        if (
            self.minimum_candidate_score < 0
            or self.minimum_discriminative_evidence_difference < 0
        ):
            raise ValueError("Candidate score thresholds cannot be negative")
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


@dataclass
class CandidateScoringResult:
    """Quantifiable rows, all reported rows, and the candidate audit table."""

    resolved_rows: pd.DataFrame
    reported_rows: pd.DataFrame
    candidate_scores: pd.DataFrame


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
        scans = (
            raw.groupby(["scan_number", "rt", "precursor_mz"], dropna=False)
            .agg(
                precursor_intensity=("precursor_intensity", "max"),
                ms2_tic=("fragment_intensity", "sum"),
            )
            .reset_index()
        )
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
            summaries.append(
                {
                    "_precursor_cluster": int(cluster),
                    "_rt_feature": int(local_id + feature_offset),
                    "feature_start_rt": float(np.min(feature_rt)),
                    "feature_apex_rt": float(feature_rt[apex_index]),
                    "feature_end_rt": float(np.max(feature_rt)),
                    "feature_reference": "precursor" if use_precursor else "ms2_tic",
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
    evidence["_event_mz"] = evidence["_observed_fragment_mz"].round(5)
    evidence["_annotation"] = evidence[annotation_col].astype(str)
    charge = evidence.get("Charge", pd.Series(0, index=evidence.index)).astype(str)
    adduct = evidence.get("Adduct", pd.Series("", index=evidence.index)).astype(str)
    evidence["_transition"] = evidence["_annotation"] + "|z" + charge + "|" + adduct

    if "mz_diff" in evidence.columns:
        evidence["_abs_mass_error"] = pd.to_numeric(
            evidence["mz_diff"], errors="coerce"
        ).abs()
    else:
        ppm = pd.to_numeric(
            evidence.get("ppm_error", pd.Series(np.nan, index=evidence.index)),
            errors="coerce",
        ).abs()
        evidence["_abs_mass_error"] = ppm * evidence["_observed_fragment_mz"] / 1e6

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

            normalized_error = (
                pd.to_numeric(group["_abs_mass_error"], errors="coerce")
                / max(config.fragment_mass_tolerance, 1e-12)
            ).to_numpy(dtype=float)
            mass_weights = np.sqrt(group["fragment_intensity"].to_numpy(dtype=float))
            median_normalized_error = _weighted_median(normalized_error, mass_weights)
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
                    "median_fragment_mass_error": (
                        median_normalized_error * config.fragment_mass_tolerance
                        if np.isfinite(median_normalized_error)
                        else np.nan
                    ),
                    "fragment_mass_accuracy_score": fragment_accuracy,
                    "median_signed_precursor_ppm_error": median_signed_precursor_error,
                    "median_precursor_ppm_error": median_absolute_precursor_error,
                    **coelution,
                    **feature_info,
                }
            )
    return _finalize_candidate_scores(pd.DataFrame(rows), config)


def _finalize_candidate_scores(
    scores: pd.DataFrame,
    config: CandidateScoringConfig,
) -> pd.DataFrame:
    """Calibrate precursor evidence and make optional terms feature-comparable."""
    if scores.empty:
        return scores
    finalized = scores.copy()
    keys = ["_precursor_cluster", "_rt_feature"]

    # Use one value per precursor-cluster/candidate pair so long RT features do
    # not receive more influence on calibration than short ones. Iteratively
    # choose the candidate closest to the run-level error center in each
    # cluster, then estimate that center and spread robustly.
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
    finalized["precursor_calibration_points"] = int(calibration_values.size)
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

    finalized["candidate_score"] = 100.0 * (
        0.35 * finalized["explained_intensity_score"]
        + 0.20 * finalized["fragment_support_score"]
        + 0.15 * finalized["fragment_mass_accuracy_score"]
        + 0.10 * finalized["precursor_relative_likelihood"]
        + 0.20 * finalized["coelution_component_used"]
    )

    # The discriminative score is centered within each feature. Equal/common
    # components therefore cancel exactly. Precursor evidence enters as a
    # calibrated Gaussian log-likelihood instead of a flat tolerance score.
    finalized["raw_discriminative_evidence"] = (
        35.0 * finalized["explained_intensity_score"]
        + 20.0 * finalized["fragment_support_score"]
        + 15.0 * finalized["fragment_mass_accuracy_score"]
        + config.precursor_log_likelihood_weight
        * finalized["precursor_log_likelihood"]
        + 20.0 * finalized["coelution_component_used"]
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
    finalized["candidate_rejection_reason"] = np.where(
        finalized["mass_outlier_pruned"],
        "no_specific_evidence_and_precursor_mass_outlier",
        "",
    )
    finalized["eligible_for_resolution"] = ~finalized["mass_outlier_pruned"]
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
        if len(ordered) == 1:
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


def score_and_resolve_candidates(
    matched_data: pd.DataFrame,
    scan_data: pd.DataFrame | None = None,
    config: CandidateScoringConfig | None = None,
    resolve: bool = True,
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
    scores = _apply_resolution(scores, config, resolve=resolve)

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
        "eligible_for_resolution",
        "mass_outlier_pruned",
        "candidate_rejection_reason",
        "mass_outlier_threshold_ppm",
        "explained_intensity_fraction",
        "distinct_fragment_count",
        "candidate_specific_fragment_count",
        "fragment_mass_accuracy_score",
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
    return CandidateScoringResult(
        resolved.reset_index(drop=True),
        reported.reset_index(drop=True),
        scores.reset_index(drop=True),
    )
