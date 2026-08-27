import logging

import pytest

from glycanPRMQuant.parallelProcess import run_parallel_pipeline


def test_run_parallel_pipeline_warns_when_no_supported_files(tmp_path, caplog):
    caplog.set_level(logging.WARNING)

    run_parallel_pipeline(input_dir=str(tmp_path), output_root=str(tmp_path / "out"))

    assert "No Thermo .raw or .mzML files found" in caplog.text


def test_run_parallel_pipeline_requires_input_source(tmp_path):
    with pytest.raises(ValueError, match="Either input_files or input_dir"):
        run_parallel_pipeline(output_root=str(tmp_path / "out"))


def test_run_parallel_pipeline_rejects_mixed_input_types(tmp_path):
    with pytest.raises(ValueError, match="cannot mix"):
        run_parallel_pipeline(
            input_files=[str(tmp_path / "one.raw"), str(tmp_path / "two.mzML")],
            output_root=str(tmp_path / "out"),
        )
