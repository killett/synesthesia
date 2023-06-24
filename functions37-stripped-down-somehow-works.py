import csv
import os
import glob
import xarray as xr
import numpy as np
import functools
import nfft # pip install nfft Reference: https://github.com/jakevdp/nfft
import matplotlib.pyplot as plt
import datetime

from typing import Dict

from colour.colorimetry import MSDS_CMFS_STANDARD_OBSERVER
from colour import SpectralDistribution, sd_to_XYZ
from colour import XYZ_to_sRGB

outputfolder = os.path.join('.', 'output')
# Create output folder with current date
date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
outputfolder = os.path.join(outputfolder,f"{date_str}")
os.mkdir(outputfolder)
if not os.path.isdir(outputfolder):
    raise ValueError(f"!!! Problem creating {outputfolder}")

if not os.path.exists(outputfolder):
    os.mkdir(outputfolder)
elif not os.path.isdir(outputfolder):
    raise ValueError(f"{outputfolder} exists but is not a directory.")

sshafolder = os.path.join('.', 'sealevel_spectra','newest_full_grids','netCDF4')

def load_cie_functions() -> xr.Dataset:
    file = os.path.join('.', 'sealevel_spectra', 'ciexyz31_1_trimmed_400nm_700nm.csv')

    data = {
        'x': [],
        'y': [],
        'z': []
    }
    wavelengths = []

    # Open input file for reading
    try:
        with open(file, 'r') as in_fp:
            reader = csv.reader(in_fp)
            for row in reader:
                temp_long = float(row[0])  # Read wavelengths, record for all 3 functions
                wavelengths.append(temp_long)

                # Read x, y, z
                for i, val in enumerate(row[1:], start=0):
                    data[list(data.keys())[i]].append(float(val))

    except IOError:
        print(f"The CIE file, {file}, failed to open.")

    da = xr.Dataset(
        {var: ('wavelength', data[var]) for var in data},
        coords={'wavelength': wavelengths}
    )

    da.attrs['title'] = 'CIE 1931 color matching functions'
    da.coords['wavelength'].attrs['units'] = 'nm'

    return da

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def synthetic_spectrum(cie,mu,sig) -> xr.Dataset:

    wavelengths = cie.coords['wavelength'].values
    power = gaussian(wavelengths, mu, sig)
    
    # convert DataArray to Dataset
    spectrum = xr.DataArray(power, coords=[('wavelength', wavelengths)], name='power').to_dataset()

    return spectrum

def spectrum2xyz(spectrum, cie) -> xr.Dataset:
    # Create a SpectralDistribution object from the input xarray DataArray
    spd = SpectralDistribution(spectrum['power'].values, spectrum.coords['wavelength'].values)

    # Define the standard observer color matching functions.
    cmfs = MSDS_CMFS_STANDARD_OBSERVER['CIE 1931 2 Degree Standard Observer']

    # Compute the XYZ tristimulus values.
    XYZ = sd_to_XYZ(spd, cmfs)

    # Return a dataset of the XYZ values
    return xr.Dataset({'x': XYZ[0],
                       'y': XYZ[1],
                       'z': XYZ[2]})

def xyz2rgb(xyz) -> xr.Dataset:
    # Extract the XYZ tristimulus values from the xyz dataset
    XYZ = [xyz['x'].values, xyz['y'].values, xyz['z'].values]

    # Convert the XYZ tristimulus values to sRGB values
    RGB = XYZ_to_sRGB(XYZ)

    # Create the xarray Dataset
    rgb = xr.Dataset({'red': RGB[0],
                      'green': RGB[1],
                      'blue': RGB[2]})

    # Merge the original and new datasets
    result = xr.merge([xyz, rgb])
    print(f"xyz2rgb_new: {result = }")
    return result

def fix_gamut(rgb) -> xr.Dataset:
    # Extract the RGB values
    R = rgb['red'].values
    G = rgb['green'].values
    B = rgb['blue'].values
    # Combine RGB into a numpy array
    rgb_values = np.array([R, G, B])
    
    # If any values are < 0, add enough white light to make them all positive
    min_val = rgb_values.min()
    if min_val < 0:
        #Make "min" positive as in paper's appendix, and add 1/255 so the rescale value isn't 0.
        min_val = -min_val + 1.0/255.0
        #This factor rescales the luminance back to its original value.
        factor = rgb['y'].values / (rgb['y'].values + min_val)
        rgb_values = (rgb_values + min_val) * factor
    
    # Normalize to [0, 1]
    max_value = rgb_values.max()
    rgb_values /= max_value

    # Create a new xarray Dataset with the corrected RGB values
    result =  xr.Dataset({'x': rgb['x'],
                          'y': rgb['y'],
                          'z': rgb['z'],
                          'red': rgb_values[0],
                          'green': rgb_values[1],
                          'blue': rgb_values[2]})
    print(f"fix_gamut: {result = }")
    crashnow
    return result

def nfft_power(timeseries) -> xr.Dataset:
    # Convert datetime index to numeric (we use 'day' as the unit)
    x = (timeseries[x_key] - timeseries[x_key][0]).values.astype(float) / (24*3600*1e9)
    y = timeseries[y_key]

    N = len(x)
    x_min = np.min(x)
    x_range = np.max(x) - np.min(x)
    x_norm = (x - x_min) / x_range - 0.5

    # Define Fourier modes
    k = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range
    
    #Perform NFFT.
    f_k = nfft.nfft(x_norm,y)

    #Compute power spectrum, which is the square of the absolute value of the Fourier Transform
    power_spectrum = np.abs(f_k)**2

    #Only take the positive frequencies. Since the output is symmetric, this will not lose any information.
    power_spectrum = power_spectrum[N//2+1:]
    xf_half = xf[N//2+1:]

    # Create xarray DataArray with coordinates
    power_spectrum_da = xr.DataArray(power_spectrum, coords=[('frequency', xf_half)], name='power')

    # Convert this DataArray to a Dataset
    spectrum = power_spectrum_da.to_dataset()
    spectrum['frequency'].attrs['units'] = '1/days'
    
    return spectrum

def convert_spectrum_from_frequency_to_period(spectrum) -> xr.Dataset:
    freq_units = spectrum['frequency'].attrs.get('units', None)

    # Map of frequency units to period units
    units_map = {'1/days': 'days', 'Hz': 'seconds', '1/years': 'years'}

    # Calculate period as reciprocal of frequency.
    period = 1.0 / spectrum['frequency']

    # Handle units attribute
    if freq_units in units_map:
        period.attrs['units'] = units_map[freq_units]
    else:
        print(f"!!!WARNING!!! DID NOT RECOGNIZE {freq_units = }")

    # Create a new Dataset with the same variables but with an additional 'period' data variable
    new_spectrum = spectrum.assign_coords(period=('frequency', period.data))  # use .data to get the underlying numpy array

    # Drop the 'frequency' dimension and coordinate
    new_spectrum = new_spectrum.swap_dims({'frequency': 'period'}).drop('frequency')

    # Reorder the dataset so that period is increasing
    new_spectrum = new_spectrum.sortby('period')

    return new_spectrum

def map_power_spectrum(cie, power_spectrum, min_period = -1, max_period = -1) -> xr.Dataset:
    # Get existing periods
    existing_periods = power_spectrum['period'].values

    # Create a flag to track if a new period was added
    new_period_added = False

    # Check if min_period and max_period are already in existing_periods
    if min_period not in existing_periods:
        existing_periods = np.append(existing_periods, min_period)
        new_period_added = True
    if max_period not in existing_periods:
        existing_periods = np.append(existing_periods, max_period)
        new_period_added = True

    # If a new period was added, sort and reindex
    if new_period_added:
        new_periods = np.sort(existing_periods)
        power_spectrum = power_spectrum.reindex(period=new_periods, method='nearest')

    # Interpolate power values for the new periods
    power_spectrum['power'] = power_spectrum['power'].interpolate_na(dim='period')
    
    # Remove all points outside the range [min_period, max_period]
    power_spectrum = power_spectrum.where((power_spectrum['period'] >= min_period) & (power_spectrum['period'] <= max_period), drop=True)

    # Map 'power' from power_spectrum onto a new wavelength scale and interpolate
    mapped_power_values = np.interp(np.linspace(min_period, max_period, len(cie['wavelength'])), power_spectrum['period'], power_spectrum['power'])

    # Create a new Dataset that shares the 'wavelength' coordinate with cie, with 'power' as its data variable
    new_power_spectrum = xr.Dataset(
        {'power': (('wavelength',), mapped_power_values)},
        coords={'wavelength': cie['wavelength']}
    )

    return new_power_spectrum

def load_ssha_files(tskip=1) -> xr.Dataset:
    # Get the list of files
    sshafiles = sorted(glob.glob(os.path.join(sshafolder,'*.nc')))

    # Get every tskip file
    tskip_files = sshafiles[::tskip]

    # Load all files into the same dataset
    print(f"Loading {len(tskip_files)} SSHA files...")
    input_data = xr.open_mfdataset(tskip_files, combine='by_coords')
    return input_data

def timeseries_to_xyz(timeseries: xr.Dataset, cie: xr.Dataset) -> xr.Dataset:
    spectrum = synthetic_spectrum(cie, 550, 5)
    xyz = spectrum2xyz(spectrum, cie)
    return xyz

if __name__ == "__main__":
    plt.style.use('dark_background')
    plt.rcParams['font.size'] = 14  # Change the global font size
    plt.rcParams['axes.linewidth'] = 2  # Change the global linewidth
    myfigsize=(10,5)
        
    x_key = 'Time'
    y_key = 'SLA'
    lat_key = 'Latitude'
    lon_key = 'Longitude'
    min_period = 10.5
    max_period = 200
    
    debug_on = 1

    if debug_on:
        cie = load_cie_functions()
        spectrum = synthetic_spectrum(cie, 550, 5)
        xyz = spectrum2xyz(spectrum, cie)
        print(f"{xyz = }")
        rgb = xyz2rgb(xyz)
        rgb = fix_gamut(rgb)
        crashnow

    input_data = load_ssha_files(tskip=1)
    xskip = 192
    print(f"Grabbing one lat/lon point in every {xskip**2} points...")
    input_data = input_data.isel({lat_key: slice(None, None, xskip), lon_key: slice(None, None, xskip)})
    print(" done.")
    # Check if the size of x_key dimension is odd
    if input_data.dims[x_key] % 2 == 1:
        print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {input_data.dims[x_key] = }")
        print(f"!!! DELETING LAST DATA POINT!")
        # If it is, select all elements up to the second last one
        input_data = input_data.isel({x_key: slice(None, -1)})
    
    cie = load_cie_functions()

    # Stack latitude and longitude into a new single dimension latlon
    stacked = input_data.stack(latlon=[lat_key, lon_key])

    # Create a new function with min_period and max_period filled
    timeseries_to_xyz_partial = functools.partial(timeseries_to_xyz, cie=cie)

    # Apply the function to the 'y_key' variable of the stacked dataset
    result = stacked.groupby('latlon').map(timeseries_to_xyz_partial)

    #Convert to RGB, but keep XYZ values around.
    result = result.groupby('latlon').map(xyz2rgb)
    result = result.groupby('latlon').map(fix_gamut)