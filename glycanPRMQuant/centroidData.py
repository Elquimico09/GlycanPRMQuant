import numpy as np
import pandas as pd

def gaussian_fit(data, mz_col='mz', intensity_col='intensity', resolution=120000):
    """
    Reconstruct a mass spectrum by fitting Gaussian peaks to each (m/z, intensity) pair.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain columns:
        - scan_number
        - <mz_col>           # name of the m/z column
        - <intensity_col>    # name of the intensity column
    mz_col : str, optional
        Name of the column in `data` holding the m/z values. Default is 'mz'.
    intensity_col : str, optional
        Name of the column in `data` holding the peak intensities. Default is 'intensity'.
    resolution : float, optional
        Instrument resolution used to compute FWHM of peaks
        (FWHM = mz / resolution). Default is 120000.

    Returns
    -------
    pandas.DataFrame
        Two columns:
        - 'mz': uniformly spaced m/z grid
        - 'intensity': reconstructed spectrum
    """
    # Check for required columns
    required = {'scan_number', mz_col, intensity_col}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Input DataFrame is missing columns: {missing}")

    # Build the m/z grid
    mz_min = data[mz_col].min() - 1
    mz_max = data[mz_col].max() + 1
    num_points = 20000
    mz_grid = np.linspace(mz_min, mz_max, num_points)

    # Initialize spectrum
    reconstructed = np.zeros_like(mz_grid)

    # Sum Gaussian peaks
    for _, row in data.iterrows():
        center = row[mz_col]
        height = row[intensity_col]
        fwhm = center / resolution
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        peak = height * np.exp(-0.5 * ((mz_grid - center) / sigma) ** 2)
        reconstructed += peak

    return pd.DataFrame({'mz': mz_grid, 'intensity': reconstructed})
