import pytest

from glycanPRMQuant.cli import build_parser


def test_run_parser_can_disable_isobaric_resolution():
    parser = build_parser()

    args = parser.parse_args(
        ["run", "sample.mzML", "out", "--disable-isobaric-resolution"]
    )

    assert args.disable_isobaric_resolution is True


def test_run_parser_uses_one_precursor_tolerance_and_named_fragment_tolerance():
    parser = build_parser()

    args = parser.parse_args(
        [
            "run",
            "sample.mzML",
            "out",
            "--ppm-ms1-tol",
            "12",
            "--fragment-mass-tol",
            "20",
            "--fragment-mass-tol-unit",
            "ppm",
        ]
    )

    assert args.ppm_ms1_tol == 12
    assert args.fragment_mass_tol == 20
    assert args.fragment_mass_tol_unit == "ppm"
    assert not hasattr(args, "ppm_ms2_tol")
    assert not hasattr(args, "mz_tol")
    assert not hasattr(args, "mz_min")
    assert not hasattr(args, "mz_max")


def test_run_parser_accepts_candidate_resolution_thresholds():
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "sample.mzML",
            "out",
            "--candidate-min-fragments",
            "3",
            "--candidate-min-explained-intensity",
            "0.02",
            "--candidate-min-score",
            "45",
            "--candidate-min-evidence-difference",
            "6",
            "--candidate-mass-outlier-min-delta",
            "2.5",
            "--candidate-max-q-value",
            "0.01",
            "--disable-target-decoy",
            "--target-decoy-seed",
            "42",
        ]
    )

    assert args.candidate_min_fragments == 3
    assert args.candidate_min_explained_intensity == 0.02
    assert args.candidate_min_score == 45
    assert args.candidate_min_evidence_difference == 6
    assert args.candidate_mass_outlier_min_delta == 2.5
    assert args.candidate_max_q_value == 0.01
    assert args.disable_target_decoy is True
    assert args.target_decoy_seed == 42


@pytest.mark.parametrize("filetype", ["png", "pdf", "svg"])
def test_run_parser_accepts_figure_filetypes(filetype):
    parser = build_parser()

    args = parser.parse_args(
        ["run", "sample.mzML", "out", "--figure-filetype", filetype]
    )

    assert args.figure_filetype == filetype


def test_run_parser_defaults_to_pdf_figures():
    parser = build_parser()

    args = parser.parse_args(["run", "sample.mzML", "out"])

    assert args.figure_filetype == "pdf"


def test_run_parser_rejects_unsupported_figure_filetype():
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            ["run", "sample.mzML", "out", "--figure-filetype", "jpg"]
        )


def test_batch_parser_defaults_to_cross_run_consensus_selection():
    parser = build_parser()

    args = parser.parse_args(
        ["batch", "--input-dir", "inputs", "--output-root", "out"]
    )

    assert args.disable_consensus_peak_selection is False
    assert args.consensus_rt_tolerance == 0.3
    assert args.consensus_min_replicate_fraction == 0.8


def test_batch_parser_accepts_consensus_peak_options():
    parser = build_parser()

    args = parser.parse_args(
        [
            "batch",
            "--input-dir",
            "inputs",
            "--output-root",
            "out",
            "--disable-consensus-peak-selection",
            "--consensus-rt-tolerance",
            "0.5",
            "--consensus-min-replicate-fraction",
            "0.6",
        ]
    )

    assert args.disable_consensus_peak_selection is True
    assert args.consensus_rt_tolerance == 0.5
    assert args.consensus_min_replicate_fraction == 0.6


@pytest.mark.parametrize(
    "removed_option",
    [
        "--ppm-ms2-tol",
        "--mz-tol",
        "--mz-min",
        "--mz-max",
        "--candidate-min-score-margin",
        "--candidate-min-score-ratio",
    ],
)
def test_run_parser_rejects_removed_tolerance_options(removed_option):
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["run", "sample.mzML", "out", removed_option, "10"])
