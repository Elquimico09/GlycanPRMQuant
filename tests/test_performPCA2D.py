import pandas as pd

from glycanPRMQuant.performPCA import plot_pca_2d


def test_plot_pca_2d_returns_scores(tmp_path):
    csv_path = tmp_path / "auc.csv"
    output_svg = tmp_path / "pca.svg"
    pd.DataFrame(
        {
            "Glycan": [25000, 26000, 43001],
            "sample_a": [10.0, 20.0, 5.0],
            "sample_b": [11.0, 18.0, 6.0],
            "sample_c": [30.0, 4.0, 25.0],
        }
    ).to_csv(csv_path, index=False)

    pca, scores = plot_pca_2d(str(csv_path), save_path=str(output_svg))

    assert output_svg.exists()
    assert scores.shape == (3, 2)
    assert pca.n_components == 2
