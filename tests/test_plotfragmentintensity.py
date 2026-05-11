from glycanPRMQuant.plotFragmentIntensity import plot_ms2_fragments

import pandas as pd


def test_plot_ms2_fragments_writes_plot(tmp_path):
    csv_path = tmp_path / "ms2.csv"
    output_svg = tmp_path / "fragments.svg"
    pd.DataFrame(
        {
            "scan_number": [1, 2, 1, 2],
            "rt": [1.0, 2.0, 1.0, 2.0],
            "Fragment": ["B1", "B1", "Y1", "Y1"],
            "fragment_mz": [100.0, 100.1, 200.0, 200.1],
            "fragment_intensity": [50.0, 60.0, 100.0, 120.0],
        }
    ).to_csv(csv_path, index=False)

    plot_ms2_fragments(
        str(csv_path),
        window=0,
        top_n=2,
        save_path=str(output_svg),
        figsize=(4, 3),
    )

    assert output_svg.exists()
