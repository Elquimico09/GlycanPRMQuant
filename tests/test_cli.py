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
            "0.03",
        ]
    )

    assert args.ppm_ms1_tol == 12
    assert args.fragment_mass_tol == 0.03
    assert not hasattr(args, "ppm_ms2_tol")
    assert not hasattr(args, "mz_tol")


@pytest.mark.parametrize("removed_option", ["--ppm-ms2-tol", "--mz-tol"])
def test_run_parser_rejects_removed_tolerance_options(removed_option):
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["run", "sample.mzML", "out", removed_option, "10"])
