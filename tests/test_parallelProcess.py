import logging

import pytest

import glycanPRMQuant.parallelProcess as parallel_process
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


def test_run_parallel_pipeline_forwards_figure_filetype(tmp_path, monkeypatch):
    input_file = tmp_path / "sample.mzML"
    input_file.write_text("", encoding="utf-8")
    submitted_args = []

    class FinishedFuture:
        def result(self):
            return "sample", "done", None

    class RecordingExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        def submit(self, function, *args):
            submitted_args.append(args)
            return FinishedFuture()

    monkeypatch.setattr(parallel_process, "ProcessPoolExecutor", RecordingExecutor)
    monkeypatch.setattr(parallel_process, "as_completed", lambda futures: futures)

    run_parallel_pipeline(
        input_files=[str(input_file)],
        output_root=str(tmp_path / "out"),
        figure_filetype="svg",
    )

    assert submitted_args[0][-1] == "svg"
