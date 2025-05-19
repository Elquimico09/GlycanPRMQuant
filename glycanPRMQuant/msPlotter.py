import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots

def plotMSData(ms_data, mz_range = None):
    """
    The function `plotMSData` plots mass spectrometry data with optional m/z range filtering.
    
    :param ms_data: The `ms_data` parameter is expected to be a DataFrame containing mass spectrometry
    data with columns 'mz' and 'intensity'. The 'mz' column represents the mass-to-charge ratio values,
    and the 'intensity' column represents the corresponding intensities of the peaks
    :param mz_range: The `mz_range` parameter in the `plotMSData` function allows you to specify a range
    of m/z values to focus on when plotting the MS data. If provided, the function will only plot the
    data points within the specified m/z range and adjust the x-axis limits accordingly. This
    :return: The function `plotMSData` returns a matplotlib figure object containing the plot of m/z
    values against intensities from the input MS data. If a specific m/z range is provided, the plot is
    limited to that range with adjusted y-axis limits based on the maximum intensity within that range.
    """

    plt.style.use(['science', 'no-latex'])
    # Check if the input DataFrame contains the required columns
    required_cols = ['mz', 'intensity']
    if not all(col in ms_data.columns for col in required_cols):
        raise ValueError('Input DataFrame must contain columns: mz, intensity')
    
    # Create a new figure
    fig, ax = plt.subplots()
    
    # Plot the m/z values against intensities
    ax.plot(ms_data['mz'], ms_data['intensity'], linestyle='-', color='b')
    
    # Set the axis labels
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    
    # Set the title
    ax.set_title('MS Data Plot')

    # Set the x-axis limits
    if mz_range:
        min_mz, max_mz = mz_range
        ax.set_xlim(min_mz, max_mz)

        data_slice = ms_data[(ms_data['mz'] >= min_mz) & (ms_data['mz'] <= max_mz)]
        
        if not data_slice.empty:
            max_intensity = data_slice['intensity'].max()
            ax.set_ylim(0, max_intensity * 1.1)
        
        else:
            pass
    
    return fig