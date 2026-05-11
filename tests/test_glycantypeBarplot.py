import pandas as pd

from glycanPRMQuant.glycantypeBarplot import plot_barplot


def test_plot_barplot_returns_class_means(tmp_path):
    csv_path = tmp_path / "auc.csv"
    output_svg = tmp_path / "barplot.svg"
    pd.DataFrame(
        {
            "Glycan": [25000, 43110, 43001],
            "sample_a": [10.0, 20.0, 30.0],
            "sample_b": [5.0, 15.0, 25.0],
        }
    ).to_csv(csv_path, index=False)

    means = plot_barplot(str(csv_path), figsize=(4, 3), save_path=str(output_svg))

    assert output_svg.exists()
    assert "high mannose" in means.index
