import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import warnings
  
def gaussian_fit(data, resolution = 120000):
    """
    The function `gaussian_fit` takes input data with scan number, m/z, and intensity columns, fits
    Gaussian peaks to the data, and returns a reconstructed spectrum DataFrame.
    
    :param data: The `data` parameter is expected to be a DataFrame containing information about peaks
    in a mass spectrometry experiment. The DataFrame should have columns named 'scan_number', 'mz', and
    'intensity' which represent the scan number, mass-to-charge ratio, and intensity of each peak
    respectively
    :param resolution: The `resolution` parameter in the `gaussian_fit` function determines the
    resolution of the Gaussian peaks that will be generated based on the input data. It is used to
    calculate the full width at half maximum (FWHM) of the Gaussian peaks. A higher resolution value
    will result in narrower Gaussian peaks, defaults to 120000 (optional)
    :return: The function `gaussian_fit` returns a DataFrame containing the reconstructed spectrum with
    columns 'mz' and 'intensity'.
    """
    
    columns_required = ['scan_number', 'mz', 'intensity']
    if not all(col in data.columns for col in columns_required):
        raise ValueError('Input DataFrame must contain columns: scan_number, mz, intensity')
    
    # Decide on an mz grid
    mz_min = data['mz'].min() - 1
    mz_max = data['mz'].max() + 1
    num_points = 20000
    mz_grid = np.linspace(mz_min, mz_max, num_points)

    # Prepare an array to hold the summed spectrum
    reconstructed = np.zeros_like(mz_grid)

    # Build each gaussian peak and accumulate
    for idx, row in data.iterrows():
        mz_center = row['mz']
        intensity = row['intensity']

        fwhm = mz_center / resolution

        # convert FWHM to sigma for Gaussian
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))

        # Build the gaussian peak
        peak = intensity * np.exp(-0.5 * ((mz_grid - mz_center) / sigma) ** 2)

        # Add to the reconstructed spectrum
        reconstructed += peak

    # Create a DataFrame for the reconstructed spectrum
    reconstructed_df = pd.DataFrame({
        'mz': mz_grid,
        'intensity': reconstructed
    })

    return reconstructed_df