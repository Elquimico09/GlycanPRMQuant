import logging
from types import SimpleNamespace

import pandas as pd
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


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"consensus_rt_tolerance": 0}, "RT tolerance"),
        ({"consensus_min_replicate_fraction": 1.1}, "replicate fraction"),
    ],
)
def test_run_parallel_pipeline_validates_consensus_options(
    tmp_path, kwargs, message
):
    with pytest.raises(ValueError, match=message):
        run_parallel_pipeline(
            input_dir=str(tmp_path), output_root=str(tmp_path / "out"), **kwargs
        )


def test_run_parallel_pipeline_runs_consensus_consolidation(tmp_path, monkeypatch):
    input_files = []
    for name in ["sample1.mzML", "sample2.mzML"]:
        path = tmp_path / name
        path.write_text("", encoding="utf-8")
        input_files.append(str(path))

    class FinishedFuture:
        def __init__(self, sample):
            self.sample = sample

        def result(self):
            return self.sample, "done", None

    class RecordingExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        def submit(self, function, *args):
            return FinishedFuture(args[0].rsplit("/", 1)[-1].split(".")[0])

    consensus_calls = []

    def fake_consolidate(root, output, **kwargs):
        consensus_calls.append((root, output, kwargs))
        return SimpleNamespace(
            peak_groups=pd.DataFrame({"consensus_selected": [True, False]})
        )

    monkeypatch.setattr(parallel_process, "ProcessPoolExecutor", RecordingExecutor)
    monkeypatch.setattr(parallel_process, "as_completed", lambda futures: futures)
    monkeypatch.setattr(
        parallel_process, "consolidate_consensus_peak_results", fake_consolidate
    )

    output_root = tmp_path / "out"
    run_parallel_pipeline(
        input_files=input_files,
        output_root=str(output_root),
        consensus_rt_tolerance=0.4,
        consensus_min_replicate_fraction=0.75,
    )

    assert consensus_calls == [
        (
            str(output_root),
            str(output_root / "combined_auc_values.csv"),
            {
                "rt_tolerance_minutes": 0.4,
                "minimum_replicate_fraction": 0.75,
                "expected_samples": ["sample1", "sample2"],
            },
        )
    ]
