from glycanPRMQuant.cli import build_parser


def test_run_parser_can_disable_isobaric_resolution():
    parser = build_parser()

    args = parser.parse_args(
        ["run", "sample.mzML", "out", "--disable-isobaric-resolution"]
    )

    assert args.disable_isobaric_resolution is True
