from glycanPRMQuant.processmzML import process_mzml_pipeline


def test_process_mzml_pipeline_is_importable():
    assert callable(process_mzml_pipeline)
