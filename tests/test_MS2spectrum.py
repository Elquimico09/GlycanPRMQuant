from glycanPRMQuant.plotMS2spectrum import plotMS2spectrum

import pandas as pd


def test_plot_ms2_spectrum_returns_averaged_fragments(tmp_path):
    csv_path = tmp_path / "ms2.csv"
    output_svg = tmp_path / "spectrum.svg"
    pd.DataFrame(
        {
            "scan_number": [1, 1, 2, 2],
            "rt": [1.0, 1.0, 1.5, 1.5],
            "Fragment": ["B1", "Y1", "B1", "Y1"],
            "Charge": [1, 1, 1, 1],
            "fragment_mz": [100.0, 200.0, 100.1, 200.1],
            "fragment_intensity": [50.0, 100.0, 60.0, 120.0],
        }
    ).to_csv(csv_path, index=False)

    out = plotMS2spectrum(
        file_path=str(csv_path),
        window_minutes=2,
        save_path=str(output_svg),
        top_n=2,
        figsize=(4, 3),
    )

    assert output_svg.exists()
    assert {"Fragment", "Charge", "mz", "avg_intensity"}.issubset(out.columns)
