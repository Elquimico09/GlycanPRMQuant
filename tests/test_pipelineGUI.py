import pytest

from glycanPRMQuant.pipelineGUI import _ensure_output_directory


def test_ensure_output_directory_creates_nested_path(tmp_path):
    output = tmp_path / "new" / "analysis"

    normalized = _ensure_output_directory(str(output))

    assert output.is_dir()
    assert normalized == str(output)


def test_ensure_output_directory_accepts_existing_directory(tmp_path):
    assert _ensure_output_directory(str(tmp_path)) == str(tmp_path)


def test_ensure_output_directory_rejects_blank_path():
    with pytest.raises(ValueError, match="Enter or select"):
        _ensure_output_directory("   ")


def test_ensure_output_directory_rejects_existing_file(tmp_path):
    output_file = tmp_path / "results.txt"
    output_file.write_text("not a directory", encoding="utf-8")

    with pytest.raises(NotADirectoryError, match="not a directory"):
        _ensure_output_directory(str(output_file))
