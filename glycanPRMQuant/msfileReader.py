from pyteomics import mzml
import numpy as np
import pandas as pd

def extractMS1(mzml_file, min_intensity = 1):
    """
    This Python function extracts MS1 data from an mzML file based on a minimum intensity threshold.
    
    :param mzml_file: The `mzml_file` parameter in the `extractMS1` function is the path to the mzML
    file from which you want to extract MS1 data. This function reads the mzML file, extracts MS1 data
    (mass spectra with MS level 1), and returns a DataFrame containing
    :param min_intensity: The `min_intensity` parameter in the `extractMS1` function is used to filter
    out peaks with intensities below a certain threshold. Peaks with intensities lower than the
    specified `min_intensity` value will not be included in the extracted MS1 data, defaults to 1
    (optional)
    :return: The function `extractMS1` returns a pandas DataFrame containing extracted MS1 data from the
    provided mzML file. The DataFrame includes columns for scan number, m/z values, retention times,
    intensities, total ion current (TIC), and base peak intensities.
    """
    scan_numbers = []
    rts = []
    mzs = []
    intensities = []
    tics = []
    base_peaks = []

    with mzml.read(mzml_file) as reader:
        for i, spectrum in enumerate(reader):
            if spectrum.get('ms level', 0) == 1:
                scan_number = spectrum.get('index') + 1
                rt = spectrum.get('scanList', {}).get('scan', [{}])[0].get('scan start time')

                mz_array = spectrum.get('m/z array', [])
                intensity_array = spectrum.get('intensity array', [])

                tic = spectrum.get('total ion current', 0)
                base_peak = spectrum.get('base peak intensity', 0)
                print(f'Processing scan {scan_number} with {len(mz_array)} peaks')

                for mz, intensity in zip(mz_array, intensity_array):
                    if intensity > min_intensity:
                        scan_numbers.append(scan_number)
                        rts.append(rt)
                        mzs.append(mz)
                        intensities.append(intensity)
                        tics.append(tic)
                        base_peaks.append(base_peak)

    data = {
        'scan_number': scan_numbers,
        'mz': mzs,
        'rt': rts,
        'intensity': intensities,
        'tic': tics,
        'base_peak': base_peaks
    }

    data = pd.DataFrame(data)

    return data

def extractXIC(mzml_file, mz, ppm_tol = 10):
    """
    The function `extractXIC` reads an mzML file, extracts data points within a specified m/z tolerance
    around a given m/z value, and returns the scan numbers, retention times, and intensities of the
    extracted data points in a DataFrame.
    
    :param mzml_file: The `mzml_file` parameter in the `extractXIC` function is the path to the mzML
    file from which you want to extract extracted ion chromatogram (XIC) data. This function reads the
    mzML file, extracts data points that match the specified m/z value within the
    :param mz: The `mz` parameter in the `extractXIC` function stands for the target mass-to-charge
    ratio that you want to extract from the provided mzML file. This function reads the mzML file,
    iterates through the spectra, and extracts data points that are within a certain parts-per
    :param ppm_tol: The `ppm_tol` parameter in the `extractXIC` function stands for parts per million
    tolerance. It is used to define the acceptable deviation in parts per million (ppm) when matching
    the specified m/z value with the m/z values in the spectrum. This parameter allows for a flexible,
    defaults to 10 (optional)
    :return: The function `extractXIC` returns a pandas DataFrame containing three columns:
    'scan_number', 'rt', and 'intensity'. These columns represent the scan numbers, retention times, and
    intensities of the peaks that match the specified m/z value within the given ppm tolerance in the
    provided mzML file.
    """
    scan_numbers = []
    rts = []
    intensities = []

    with mzml.read(mzml_file) as reader:
        for i, spectrum in enumerate(reader):
            if spectrum.get('ms level', 0) == 1:
                scan_number = spectrum.get('index') + 1
                rt = spectrum.get('scanList', {}).get('scan', [{}])[0].get('scan start time')

                mz_array = spectrum.get('m/z array', [])
                intensity_array = spectrum.get('intensity array', [])

                if len(mz_array) > 0:
                    mask = np.abs(np.array(mz_array) - mz) / mz * 1e6 < ppm_tol
                    if np.any(mask):
                        matched_intensities = np.array(intensity_array)[mask]
                        total_intensity = np.sum(matched_intensities)
                        print(f'Scan = {scan_number}, RT = {rt}, Intensity = {total_intensity}')
                        scan_numbers.append(scan_number)
                        rts.append(rt)
                        intensities.append(total_intensity)

    
    data = {
        'scan_number': scan_numbers,
        'rt': rts,
        'intensity': intensities
    }

    data = pd.DataFrame(data)

    return data

def extractMS2(mzml_file, min_intensity = 1):
    """
    This Python function extracts MS2 data from an mzML file based on a minimum intensity threshold.
    
    :param mzml_file: The `mzml_file` parameter in the `extractMS2` function is the path to the mzML
    file from which you want to extract MS2 data. This function reads the mzML file, extracts MS2 data
    (mass spectra with MS level 2), and returns a DataFrame containing
    :param min_intensity: The `min_intensity` parameter in the `extractMS2` function is used to filter
    out peaks with intensities below a certain threshold. Peaks with intensities lower than the
    specified `min_intensity` value will not be included in the extracted MS2 data, defaults to 1
    (optional)
    :return: The function `extractMS2` returns a pandas DataFrame containing extracted MS2 data from the
    provided mzML file. The DataFrame includes columns for scan number, m/z values, retention times,
    intensities, total ion current (TIC), and base peak intensities.
    """
    scan_numbers = []
    rts = []
    precursor_mzs = []
    fragment_mzs = []
    precursor_intensities = []
    fragment_intensities = []

    with mzml.read(mzml_file) as reader:
        for _, spectrum in enumerate(reader):
            if spectrum.get('ms level', 0) == 2:
                scan_number = spectrum.get('index') + 1
                rt = spectrum.get('scanList', {}).get('scan', [{}])[0].get('scan start time')

                precursor_mz = spectrum.get('precursorList', {}).get('precursor', [{}])[0].get('selectedIonList', {}).get('selectedIon', [{}])[0].get('selected ion m/z')
                
                # Extract precursor intensity with fallbacks
                precursor_intensity = 0
                try:
                    precursor = spectrum.get('precursorList', {}).get('precursor', [{}])[0]
                    selected_ion = precursor.get('selectedIonList', {}).get('selectedIon', [{}])[0]
                    # Try different possible keys for intensity
                    for key in ['intensity', 'selected ion intensity', 'peak intensity']:
                        if key in selected_ion:
                            precursor_intensity = selected_ion[key]
                            break
                except (IndexError, KeyError, TypeError):
                    pass

                mz_array = spectrum.get('m/z array', [])
                intensity_array = spectrum.get('intensity array', [])

                if len(mz_array) > 0:
                    mask = np.array(intensity_array) > min_intensity
                    filtered_mzs = np.array(mz_array)[mask]
                    filtered_intensities = np.array(intensity_array)[mask]

                    print(f'Processing scan {scan_number} with {len(filtered_mzs)} fragments')

                    for frag_mz, frag_intensity in zip(filtered_mzs, filtered_intensities):
                        scan_numbers.append(scan_number)
                        rts.append(rt)
                        precursor_mzs.append(precursor_mz)
                        precursor_intensities.append(precursor_intensity)
                        fragment_mzs.append(frag_mz)
                        fragment_intensities.append(frag_intensity)
    data = {
        'scan_number': scan_numbers,
        'rt': rts,
        'precursor_mz': precursor_mzs,
        'fragment_mz': fragment_mzs,
        'precursor_intensity': precursor_intensities,
        'fragment_intensity': fragment_intensities
    }
    data = pd.DataFrame(data)
    return data
