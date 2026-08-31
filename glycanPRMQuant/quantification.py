"""Re-extract accepted fragment transitions from unfiltered MS2 centroids."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from glycanPRMQuant.mass_tolerance import validate_tolerance


REEXTRACTION_AUDIT_COLUMNS = [
    "Glycan",
    "precursor_cluster",
    "rt_feature_id",
    "PrecursorAdduct",
    "feature_start_rt",
    "feature_apex_rt",
    "feature_end_rt",
    "accepted_transition_count",
    "quantification_scan_count",
    "quantification_trace_point_count",
    "detected_trace_point_count",
    "zero_trace_point_count",
    "raw_centroid_count_in_feature",
    "quantification_source",
]

FEATURE_TRACE_METADATA_COLUMNS = [
    "Glycan",
    "precursor_cluster",
    "rt_feature_id",
    "candidate_score",
    "discriminative_score",
    "distinct_fragment_count",
    "candidate_specific_fragment_count",
    "explained_intensity_fraction",
    "resolution_status",
    "assignment_q_value",
    "target_decoy_score_margin",
    "target_decoy_pass",
    "feature_scan_count",
    "feature_left_flank_scans",
    "feature_right_flank_scans",
    "chromatographic_peak_valid",
    "chromatographic_peak_rejection_reason",
]

TRANSITION_TRACE_METADATA_COLUMNS = [
    "Fragment",
    "BaseFragment",
    "FragmentType",
    "fragment_annotation",
    "contains_neuac",
    "neutral_loss",
    "Charge",
    "Adduct",
]


@dataclass(frozen=True)
class QuantificationReextractionResult:
    """Unfiltered accepted-transition traces and feature-level audit records."""

    peaks: pd.DataFrame
    audit: pd.DataFrame
    transitions: pd.DataFrame


def _theoretical_mz_column(data: pd.DataFrame) -> str:
    for column in ("theoretical_fragment_mz", "fragment_mz", "Fragment_mz"):
        if column in data.columns:
            return column
    raise ValueError(
        "Accepted transition re-extraction requires a theoretical fragment m/z column"
    )


def _greedy_scan_matches(
    observed_mz: np.ndarray,
    observed_intensity: np.ndarray,
    theoretical_mz: np.ndarray,
    tolerance_value: float,
    tolerance_unit: str,
) -> dict[int, int]:
    """Assign each centroid and accepted transition at most once within a scan."""
    if observed_mz.size == 0 or theoretical_mz.size == 0:
        return {}

    delta_da = np.abs(observed_mz[:, None] - theoretical_mz[None, :])
    if tolerance_unit == "ppm":
        with np.errstate(divide="ignore", invalid="ignore"):
            normalized_error = delta_da / theoretical_mz[None, :] * 1e6
    else:
        normalized_error = delta_da
    valid = np.isfinite(normalized_error) & (normalized_error <= tolerance_value)
    observed_indices, transition_indices = np.nonzero(valid)
    if observed_indices.size == 0:
        return {}

    errors = normalized_error[observed_indices, transition_indices]
    # Primary order is mass error; for exact ties prefer the more intense raw
    # centroid, then deterministic observed/transition indices.
    order = np.lexsort(
        (
            transition_indices,
            observed_indices,
            -observed_intensity[observed_indices],
            errors,
        )
    )
    used_observed: set[int] = set()
    used_transitions: set[int] = set()
    assignment: dict[int, int] = {}
    for position in order:
        observed_index = int(observed_indices[position])
        transition_index = int(transition_indices[position])
        if observed_index in used_observed or transition_index in used_transitions:
            continue
        used_observed.add(observed_index)
        used_transitions.add(transition_index)
        assignment[transition_index] = observed_index
    return assignment


def reextract_accepted_transitions(
    accepted_matches: pd.DataFrame,
    raw_ms2_data: pd.DataFrame,
    fragment_mass_tolerance: float,
    fragment_mass_tolerance_unit: str,
    precursor_ppm_tolerance: float,
) -> QuantificationReextractionResult:
    """Build raw-intensity traces for only the accepted identification transitions.

    Denoised identification rows define the candidate, transition list, precursor
    adduct, and chromatographic feature. Positive unfiltered centroid peaks are
    then assigned to those transitions within the accepted feature bounds. One
    centroid can quantify at most one transition and one transition receives at
    most one centroid per scan. Explicit zero rows retain missing scan points.
    """
    tolerance_value, tolerance_unit = validate_tolerance(
        fragment_mass_tolerance, fragment_mass_tolerance_unit
    )
    if precursor_ppm_tolerance <= 0:
        raise ValueError("Precursor ppm tolerance must be positive")

    group_columns = ["Glycan", "precursor_cluster", "rt_feature_id"]
    accepted_required = set(group_columns) | {
        "PrecursorAdduct",
        "precursor_mz",
        "feature_start_rt",
        "feature_end_rt",
    }
    raw_required = {
        "scan_number",
        "rt",
        "precursor_mz",
        "fragment_mz",
        "fragment_intensity",
    }
    missing_accepted = accepted_required.difference(accepted_matches.columns)
    missing_raw = raw_required.difference(raw_ms2_data.columns)
    if missing_accepted:
        raise ValueError(
            "Accepted transition re-extraction requires columns: "
            f"{sorted(missing_accepted)}"
        )
    if missing_raw:
        raise ValueError(
            f"Raw MS2 re-extraction requires columns: {sorted(missing_raw)}"
        )
    if accepted_matches.empty or raw_ms2_data.empty:
        return QuantificationReextractionResult(
            accepted_matches.iloc[0:0].copy(),
            pd.DataFrame(columns=REEXTRACTION_AUDIT_COLUMNS),
            pd.DataFrame(),
        )

    theoretical_column = _theoretical_mz_column(accepted_matches)
    raw = raw_ms2_data.copy()
    for column in ("rt", "precursor_mz", "fragment_mz", "fragment_intensity"):
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    raw = raw.loc[
        np.isfinite(raw["rt"])
        & np.isfinite(raw["precursor_mz"])
        & np.isfinite(raw["fragment_mz"])
        & np.isfinite(raw["fragment_intensity"])
        & raw["fragment_intensity"].gt(0)
    ].copy()
    if "precursor_intensity" not in raw.columns:
        raw["precursor_intensity"] = 0.0

    output_rows: list[dict] = []
    audit_rows: list[dict] = []
    transition_library_rows: list[dict] = []
    source_name = "unfiltered_accepted_transition_reextraction"

    for keys, feature in accepted_matches.groupby(group_columns, dropna=False):
        feature = feature.copy()
        start_rt = float(pd.to_numeric(feature["feature_start_rt"], errors="coerce").median())
        end_rt = float(pd.to_numeric(feature["feature_end_rt"], errors="coerce").median())
        apex_values = pd.to_numeric(
            feature.get("feature_apex_rt", pd.Series(dtype=float)), errors="coerce"
        ).dropna()
        apex_rt = float(apex_values.median()) if not apex_values.empty else np.nan
        if not np.isfinite(start_rt) or not np.isfinite(end_rt):
            continue
        if end_rt < start_rt:
            start_rt, end_rt = end_rt, start_rt

        feature_raw = raw.loc[raw["rt"].between(start_rt, end_rt)].copy()
        adduct_references = (
            feature.assign(
                precursor_mz=pd.to_numeric(feature["precursor_mz"], errors="coerce")
            )
            .dropna(subset=["precursor_mz", "PrecursorAdduct"])
            .groupby("PrecursorAdduct", sort=False)["precursor_mz"]
            .median()
        )
        if feature_raw.empty or adduct_references.empty:
            continue

        reference_values = adduct_references.to_numpy(dtype=float)
        raw_precursor = feature_raw["precursor_mz"].to_numpy(dtype=float)
        precursor_errors = (
            np.abs(raw_precursor[:, None] - reference_values[None, :])
            / reference_values[None, :]
            * 1e6
        )
        closest_adduct = np.argmin(precursor_errors, axis=1)
        closest_error = precursor_errors[np.arange(len(feature_raw)), closest_adduct]
        feature_raw = feature_raw.loc[closest_error <= precursor_ppm_tolerance].copy()
        if feature_raw.empty:
            continue
        feature_raw["_quant_adduct"] = adduct_references.index.to_numpy()[
            closest_adduct[closest_error <= precursor_ppm_tolerance]
        ]

        feature_record = feature.iloc[0]
        feature_metadata = {
            column: feature_record[column]
            for column in FEATURE_TRACE_METADATA_COLUMNS
            if column in feature.columns
        }
        raw_centroid_count = int(len(feature_raw))
        for precursor_adduct, adduct_feature in feature.groupby(
            "PrecursorAdduct", dropna=False, sort=False
        ):
            adduct_raw = feature_raw.loc[
                feature_raw["_quant_adduct"].eq(precursor_adduct)
            ]
            if adduct_raw.empty:
                continue

            transitions = adduct_feature.copy()
            transitions["_quant_theoretical_mz"] = pd.to_numeric(
                transitions[theoretical_column], errors="coerce"
            )
            transitions = transitions.dropna(subset=["_quant_theoretical_mz"])
            transition_identity = ["_quant_theoretical_mz"]
            for optional in ("Fragment", "Charge", "Adduct"):
                if optional in transitions.columns:
                    transition_identity.append(optional)
            transitions = transitions.drop_duplicates(transition_identity).reset_index(
                drop=True
            )
            if transitions.empty:
                continue

            glycan, precursor_cluster, feature_id = keys
            transition_ids = [
                f"{glycan}|{precursor_cluster}|{feature_id}|{precursor_adduct}|{index + 1}"
                for index in range(len(transitions))
            ]
            transitions["quantification_transition_id"] = transition_ids
            for _, transition in transitions.iterrows():
                library_row = transition.drop(
                    labels=["_quant_theoretical_mz"]
                ).to_dict()
                library_row.update(
                    {
                        "Glycan": glycan,
                        "precursor_cluster": precursor_cluster,
                        "rt_feature_id": feature_id,
                        "PrecursorAdduct": precursor_adduct,
                        "feature_start_rt": start_rt,
                        "feature_apex_rt": apex_rt,
                        "feature_end_rt": end_rt,
                        "quantification_source": source_name,
                    }
                )
                transition_library_rows.append(library_row)

            theoretical_mz = transitions["_quant_theoretical_mz"].to_numpy(dtype=float)
            detected_points = 0
            trace_points = 0
            scan_count = 0
            for scan_number, scan_peaks in adduct_raw.groupby(
                "scan_number", sort=True
            ):
                scan_count += 1
                scan_peaks = scan_peaks.sort_values("fragment_mz").reset_index(drop=True)
                observed_mz = scan_peaks["fragment_mz"].to_numpy(dtype=float)
                observed_intensity = scan_peaks["fragment_intensity"].to_numpy(
                    dtype=float
                )
                assignments = _greedy_scan_matches(
                    observed_mz,
                    observed_intensity,
                    theoretical_mz,
                    tolerance_value,
                    tolerance_unit,
                )
                scan_metadata = scan_peaks.iloc[0]
                raw_tic = float(observed_intensity.sum())
                for transition_index, transition in transitions.iterrows():
                    row = dict(feature_metadata)
                    row.update(
                        {
                            column: transition[column]
                            for column in TRANSITION_TRACE_METADATA_COLUMNS
                            if column in transition.index
                        }
                    )
                    row.update(
                        {
                            "quantification_transition_id": transition[
                                "quantification_transition_id"
                            ],
                            "scan_number": scan_number,
                            "rt": float(scan_metadata["rt"]),
                            "precursor_mz": float(scan_metadata["precursor_mz"]),
                            "precursor_intensity": float(
                                pd.to_numeric(
                                    scan_peaks["precursor_intensity"], errors="coerce"
                                ).fillna(0.0).max()
                            ),
                            "PrecursorAdduct": precursor_adduct,
                            "fragment_mz": float(theoretical_mz[transition_index]),
                            "theoretical_fragment_mz": float(
                                theoretical_mz[transition_index]
                            ),
                            "raw_ms2_tic": raw_tic,
                            "feature_start_rt": start_rt,
                            "feature_apex_rt": apex_rt,
                            "feature_end_rt": end_rt,
                            "quantification_source": source_name,
                            "quantification_peak_detected": bool(
                                transition_index in assignments
                            ),
                        }
                    )
                    observed_index = assignments.get(int(transition_index))
                    if observed_index is None:
                        row.update(
                            {
                                "fragment_intensity": 0.0,
                                "observed_fragment_mz": np.nan,
                                "mz_diff": np.nan,
                                "ppm_error": np.nan,
                            }
                        )
                    else:
                        detected_points += 1
                        observed = float(observed_mz[observed_index])
                        theoretical = float(theoretical_mz[transition_index])
                        row.update(
                            {
                                "fragment_intensity": float(
                                    observed_intensity[observed_index]
                                ),
                                "observed_fragment_mz": observed,
                                "mz_diff": abs(observed - theoretical),
                                "ppm_error": abs(observed - theoretical)
                                / theoretical
                                * 1e6,
                            }
                        )
                    output_rows.append(row)
                    trace_points += 1

            audit_rows.append(
                {
                    "Glycan": glycan,
                    "precursor_cluster": precursor_cluster,
                    "rt_feature_id": feature_id,
                    "PrecursorAdduct": precursor_adduct,
                    "feature_start_rt": start_rt,
                    "feature_apex_rt": apex_rt,
                    "feature_end_rt": end_rt,
                    "accepted_transition_count": int(len(transitions)),
                    "quantification_scan_count": int(scan_count),
                    "quantification_trace_point_count": int(trace_points),
                    "detected_trace_point_count": int(detected_points),
                    "zero_trace_point_count": int(trace_points - detected_points),
                    "raw_centroid_count_in_feature": raw_centroid_count,
                    "quantification_source": source_name,
                }
            )

    if not output_rows:
        peaks = accepted_matches.iloc[0:0].copy()
    else:
        peaks = pd.DataFrame(output_rows)
    audit = pd.DataFrame(audit_rows, columns=REEXTRACTION_AUDIT_COLUMNS)
    transitions = pd.DataFrame(transition_library_rows)
    return QuantificationReextractionResult(
        peaks.reset_index(drop=True), audit, transitions.reset_index(drop=True)
    )
