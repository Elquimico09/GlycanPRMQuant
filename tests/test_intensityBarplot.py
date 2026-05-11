import pandas as pd

from glycanPRMQuant.intensityBarplot import plot_glycan_intensity_boxplot


def test_plot_glycan_intensity_boxplot_writes_plot(tmp_path):
    csv_path = tmp_path / "auc.csv"
    output_svg = tmp_path / "boxplot.svg"
    pd.DataFrame(
        {
            "Glycan": [25000, 26000],
            "sample_a": [10.0, 20.0],
            "sample_b": [15.0, 25.0],
        }
    ).to_csv(csv_path, index=False)

    plot_glycan_intensity_boxplot(str(csv_path), save_path=str(output_svg))

    assert output_svg.exists()
