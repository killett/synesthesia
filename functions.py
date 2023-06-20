import csv
import os
import subprocess
import shutil
import glob
import xarray as xr
import numpy as np
import pandas as pd
import zipfile
import functools
import nfft # pip install nfft Reference: https://github.com/jakevdp/nfft
import statsmodels.api as sm #pip install statsmodels
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import datetime
import time
import timeit

from typing import Dict

outputfolder = os.path.join('.', 'output')
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
elif not os.path.isdir(outputfolder):
    raise ValueError(f"{outputfolder} exists but is not a directory.")

# Create output folder with current date
date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
outputfolder = os.path.join(outputfolder,f"{date_str}")
os.makedirs(outputfolder)
if not os.path.isdir(outputfolder):
    raise ValueError(f"!!! Problem creating {outputfolder}")

sshafolder = os.path.join('.', 'sealevel_spectra','newest_full_grids','netCDF4')

dpi_choice = 300

plot_options = {
    'outputfolder': outputfolder,
    'just_the_filenames': ['r.nc','g.nc','b.nc'],
    'output_base': 'map_parameter',
    'projection': 1, # Since the function uses plot_options.projection, a reasonable default might be 1 for the Robinson projection.
    'plot_mascons': 0,
    'coastlines': 1, #1:coast, 2:coast+InSAR.
    'blurb_disabled': 1, 
    'montage': 0, #1=left-justify titles, add (a),(b), run montage.sh.
    'region': 'global',    # Define a global region for all plots
    'frame': 'a',    # Add a frame to the plot with automatic tick intervals
    'color_scheme': 2, #1/2=white/black background
    'land_color': 'white',   # Fill continents with color white
    'sea_color': 'blue',   # Fill oceans with color blue
    'scale_digits': 2,    # Set the number of digits for the D_FORMAT
    'show_fig': False,
    'save_fig': True,
    'figsize': (10, 10),
    'dpi': dpi_choice,
    'linewidth': 2.0,
    'linestyle': '-',
    'color': 'black',
    'marker': 'o',
    'markersize': 5,
    'markerfacecolor': 'blue',
    'markeredgewidth': 1.5,
    'markeredgecolor': 'black',
    'font_size': 14,
    'font_weight': 'bold',
    'x_label': 'X-axis',
    'y_label': 'Y-axis',
    'title': 'My Plot',
    'grid': True,
    'grid_linestyle': '--',
    'grid_linewidth': 0.5,
    'grid_alpha': 0.7
}

grid = {
    'add_on': {
        'write_gmt_defs': None,
        'just_the_filenames': [],
        'color_scheme': 1,  # 1/2=white/black background
        'montage': 0,  # 1=left-justify titles, add (a),(b), run montage.sh.
        'prefixes': ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)', '(q)', '(r)', '(s)', '(t)', '(u)', '(v)', '(w)', '(x)', '(y)', '(z)'],
        'index': -1,  # Increments on each map, accesses prefixes above for montage.
        'png_options': ' -P -Tg ',  # PDF default: -E720, else 300 dpi.
        'digits': 2,  # assumed reasonable default for scale_digits.
        'output_base': '',  # empty string as placeholder, will be updated in each iteration.
        'data_name': '',  # empty string as placeholder, will be updated in each iteration.
        'write_gmt_map_data': None,
        'convert': None,
        'finish_flip_backgrounds': None,
        'finish_trim': None,
        'outputfolder': '',  # empty string as placeholder, need to be updated with real path
        'animate.sh': None,
        'montage.sh': None,
        'write_gmt_colorscale': None,
        'write_rgb_colorscale': None
    }
}

results = {
    'max_widths':2,
    'options': {'output_choice': 5},
    'latlon': {
        'lat': [10.0, 20.0, 30.0, 40.0, 50.0],
        'lon': [10.0, 20.0, 30.0, 40.0, 50.0],
        'outputs': ['output1', 'output2', 'output3'],
        'mascon_lats': [10.0, 20.0, 30.0]
    },
    'marker_lats': ['output1', 'output2', 'output3'],
    'marker_lons': ['output1', 'output2', 'output3'],
    'rgb': [1, 2, 3],
    'rgb_choice': 2,
    'units': ['unit1', 'unit2', 'unit3'],
    'maxlat': 90.0,
    'minlat': -90.0,
    'maxlon': 180.0,
    'minlon': -180.0,
    'error_bars': [0.1, 0.2, 0.3],

    # Assuming the `just_the_filenames` variable is related to the data we are dealing with, 
    # we will set it as a list with a single default value as an example. The actual values 
    # should be replaced according to your specific use case.
    'just_the_filenames': ['example_filename'],

    # Based on the code, `color_scheme` and `montage` seems to be some configuration parameters.
    # I'll assume they are integers and default them to 1.
    'color_scheme': 1,
    'montage': 1,

    # `scale_digits` appears to be controlling the formatting of some GMT parameter. 
    # I'll assume it's an integer and default it to 2 (for 2 decimal places).
    'scale_digits': 2,

    # `output_base` is being used to format filenames, I'll default it to a string 'output'.
    'output_base': 'output',

    # `titles` appears to be a list related to the `just_the_filenames`, hence it should have 
    # the same length. We will initialize it with a single default value as an example.
    'titles': ['example_title'],

    # The actual values of the keys `grid`, `plot_options`, and `outputfolder` are not clear from the provided code,
    # however, their existence can be inferred. For `grid` and `plot_options`, I'll use simple string placeholders,
    # and for `outputfolder`, I'll use a generic path.
    'grid': grid,
    'plot_options': plot_options,
    'outputfolder': '/path/to/outputfolder',
    'date': datetime.datetime.now(),
    'config': {
        'x': 0.0,
        'y': 0.0,
        'width': 0.0,
        'height': 0.0,
    },
    'xy': {
        'x_values': [[0.0]],
        'y_values': [[0.0]],
    },
    'xyz': {
        'x_values': [],
        'y_values': [],
        'z_values': [],
    },
    'min_max': {
        'x': [0.0, 0.0],
        'y': [0.0, 0.0],
        'z': [0.0, 0.0],
    },
    'misc': {
        'scale_format': '',
        'overflow': '',
        'scale_pos': '',
        'units_format': '',
        'units_pos': '',
        'misc_range': '',
        'start': '',
        'middle': '',
        'end': '',
        'scale_width': 0.0,
        'scale_x': 0.0,
        'scale_y': 0.0,
        'scale_length': 0.0,
        'plot_base': '',
    },
}

files_to_copy = ['projections.sh', 'overflow.sh', 'notation.sh']

for file in files_to_copy:
    shutil.copy(file, outputfolder)

# Determine the path of the script that is currently running
current_script_path = os.path.abspath(__file__)
# Make sure the directory exists, create if necessary
os.makedirs(outputfolder, exist_ok=True)
# Create the name for the zip file
zip_file_name = os.path.join(outputfolder, os.path.basename(current_script_path) + '.zip')
# Zip up the python script
with zipfile.ZipFile(zip_file_name, 'w') as zipf:
    zipf.write(current_script_path, arcname=os.path.basename(current_script_path))
print(f'Successfully zipped {current_script_path} to {zip_file_name}')

def load_cie_functions():
    file = os.path.join('.', 'sealevel_spectra', 'ciexyz31_1_trimmed_400nm_700nm.csv')
    #file = os.path.join('.', 'sealevel_spectra', 'ciexyz31_1_trimmed_380nm_760nm.csv') #Matches Hughes and Williams 2010

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

def synthetic_spectrum(cie,mu,sig):

    wavelengths = cie.coords['wavelength'].values
    power = gaussian(wavelengths, mu, sig)
    
    # convert DataArray to Dataset
    spectrum = xr.DataArray(power, coords=[('wavelength', wavelengths)], name='power').to_dataset()

    return spectrum

def spectrum2xyz(spectrum, cie, normalization_factor):
    xyz = {}
    wavelength_step_size = np.diff(cie['wavelength'].values).mean()  # Average wavelength step size
    for l in ['x', 'y', 'z']:
        # Multiply each spectrum value by the corresponding cie value
        temp_values = spectrum['power'].values * cie[l].values
        # Integrate over wavelengths to get a single tristimulus value
        xyz[l] = (temp_values.sum() * wavelength_step_size) / normalization_factor
    return xyz

def raise_y_to_power(xyz, power):
    power = 1 - power
    factor = pow(xyz['y'], power)
    for key in xyz:
        xyz[key] /= factor
    return xyz

def xyz2rgb(xyz):
    A = np.array([[3.2409699, -1.5373832, -0.49861079],
                  [-0.96924375, 1.8759676, 0.041555082],
                  [0.055630032, -0.20397685, 1.0569714]])

    xyz_vector = np.array([xyz['x'], xyz['y'], xyz['z']])

    # Multiply A*xyz to obtain rgb
    rgb = np.dot(A, xyz_vector)

    # If any values are < 0, add enough white light to make them all positive
    min_val = rgb.min()
    if min_val < 0:
        min_val = -min_val + 1.0/255.0
        factor = xyz['y'] / (xyz['y'] + min_val)
        rgb = (rgb + min_val) * factor

    # Normalize to [0, 1]
    rgb = {'red': rgb[0], 'green': rgb[1], 'blue': rgb[2]}
    max_value = max(rgb.values())
    for key in rgb:
        rgb[key] /= max_value

    return rgb

def gamma_correct_rgb(rgb):
    gamma_inv = 0.45
    crit = 0.018  # RGB values are gamma corrected differently below and above crit.
    h = 4.506813168
    g = -0.09914989
    f = 1.09914989
    for key in rgb:
        # Typo in Hughes and Williams 2010 equation A7, compared to Charles Poynton's GammaFAQ:
        # http://www.poynton.com/GammaFAQ.html
        if rgb[key] <= crit:
            rgb[key] *= h
        else:
            rgb[key] = f * pow(rgb[key], gamma_inv) + g
    return rgb

def synthetic_timeseries(signal='annual', signal_amplitude=1, noise='white', noise_level=0.1, 
                         temporal_resolution='monthly', time_start=datetime.datetime(2001, 1, 1), 
                         time_stop=datetime.datetime(2005, 1, 1)):
    # Generate time series
    if temporal_resolution == 'monthly':
        dates = pd.date_range(start=time_start, end=time_stop, freq='M') + pd.Timedelta(days=15)
    elif temporal_resolution == 'daily':
        dates = pd.date_range(start=time_start, end=time_stop, freq='D')
    else:
        print(f"!!!WARNING!!! {temporal_resolution = }")

    # Generate signal
    if signal == 'annual':
        t = (dates - time_start).days / 365.25
        signal_values = signal_amplitude * np.sin(2 * np.pi * t)
    else:
        print(f"!!!WARNING!!! {signal = }")

    # Generate noise
    if noise == 'white':
        noise_values = noise_level * np.random.randn(len(dates))
    else:
        print(f"!!!WARNING!!! {noise = }")

    # Combine signal and noise
    measurements = signal_values + noise_values

    # Create xarray dataset
    ds = xr.Dataset(
        {y_key: (x_key, measurements)},
        coords={x_key: dates}
    )

    return ds

def fancy_detrend(timeseries, x_key, y_key, terms=['constant', 'trend']):
    #print("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
    x = (timeseries[x_key] - timeseries[x_key][0]).values.astype(float) / (24*3600*1e9)
    y = timeseries[y_key].values

    design_matrix = []

    if 'constant' in terms:
        design_matrix.append(np.ones_like(x))

    if 'trend' in terms:
        design_matrix.append(x)

    if 'accel' in terms:
        design_matrix.append(x ** 2)

    design_matrix = np.column_stack(design_matrix)

    model = sm.OLS(y, design_matrix)
    result = model.fit()

    detrended_y = y - result.fittedvalues

    detrended_timeseries = timeseries.copy()
    detrended_timeseries[y_key].values = detrended_y

    fits = {}
    for i, term in enumerate(terms):
        fits[term] = result.params[i]

    return detrended_timeseries, fits

def turn_fits_into_timeseries(timeseries, x_key, y_key, fits):
    #print("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
    x = (timeseries[x_key] - timeseries[x_key][0]).values.astype(float) / (24*3600*1e9)
    fitted_values = np.zeros_like(x)

    for term, fit_value in fits.items():
        if term == 'constant':
            fitted_values += fit_value
        elif term == 'trend':
            fitted_values += fit_value * x
        elif term == 'accel':
            fitted_values += fit_value * x**2

    fitted_timeseries = timeseries.copy()
    fitted_timeseries[y_key].values = fitted_values

    return fitted_timeseries

def nfft_power(ds):
    # Convert datetime index to numeric (we use 'day' as the unit)
    #print("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
    x = (ds[x_key] - ds[x_key][0]).values.astype(float) / (24*3600*1e9)
    y = ds[y_key]

    #If timeseries has an odd number of points,
    #remove the last data point, then calculate min, range.
    N = 1
    while N % 2:
        N = len(x)
        x_min = np.min(x)
        x_range = np.max(x) - np.min(x)
        x_norm = (x - x_min) / x_range - 0.5
        #print(f"{N = }, {x_min = }, {x_range = }")
        if N % 2:
            print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(x) = }")
            print(f"!!! DELETING LAST DATA POINT!")
            x = np.delete(x, -1)
            y = np.delete(y, -1)

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
    ds = power_spectrum_da.to_dataset()
    ds['frequency'].attrs['units'] = '1/days'
    
    return ds

def convert_spectrum_from_frequency_to_period(ds):
    freq_units = ds['frequency'].attrs.get('units', None)

    # Map of frequency units to period units
    units_map = {'1/days': 'days', 'Hz': 'seconds', '1/years': 'years'}

    # Calculate period as reciprocal of frequency.
    period = 1.0 / ds['frequency']

    # Handle units attribute
    if freq_units in units_map:
        period.attrs['units'] = units_map[freq_units]
    else:
        print(f"!!!WARNING!!! DID NOT RECOGNIZE {freq_units = }")

    # Create a new Dataset with the same variables but with an additional 'period' data variable
    new_ds = ds.assign_coords(period=('frequency', period.data))  # use .data to get the underlying numpy array

    # Drop the 'frequency' dimension and coordinate
    new_ds = new_ds.swap_dims({'frequency': 'period'}).drop('frequency')

    # Reorder the dataset so that period is increasing
    new_ds = new_ds.sortby('period')

    return new_ds

def map_power_spectrum(cie, power_spectrum, min_period = -1, max_period = -1):
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

def plot_timeseries(ds, title):
    plt.figure(figsize=myfigsize)
    plt.plot(ds[x_key], ds[y_key], color='lime') # using a bright color for visibility
    plt.scatter(ds[x_key], ds[y_key], marker='s', color='cyan', s=10)
    plt.title(title, color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_{title.replace(' ','_')}.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

def plot_fft_spectrum(power_spectrum, title):
    plt.figure(figsize=myfigsize)
    plt.plot(power_spectrum.period, power_spectrum.power, color = 'lime', linewidth=2.0)
    plt.scatter(power_spectrum.period, power_spectrum.power, marker='s', color='cyan', s=10)
    plt.title(title)
    plt.xlabel('Period (days)')
    plt.ylabel('Power')
    #plt.grid(True, color='gray')  # set grid color to gray for visibility
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_{title.replace(' ','_')}.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

def plot_light_spectrum(power_spectrum, title):
    plt.figure(figsize=myfigsize)
    plt.plot(power_spectrum.wavelength, power_spectrum.power, color = 'lime', linewidth=2.0)
    plt.scatter(power_spectrum.wavelength, power_spectrum.power, marker='s', color='cyan', s=10)
    plt.title(title)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Power')
    #plt.grid(True, color='gray')  # set grid color to gray for visibility
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_{title.replace(' ','_')}.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

def plot_color(rgb, filename):
    fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=dpi_choice)

    # Set the facecolor using the normalized RGB values
    ax.set_facecolor(tuple(rgb[key] for key in ['red', 'green', 'blue']))

    # Remove all axes and labels
    ax.axis('off')

    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

def load_ssha_files(tskip=1):
    # Get the list of files
    sshafiles = sorted(glob.glob(os.path.join(sshafolder,'*.nc')))

    # Get every tskip file
    tskip_files = sshafiles[::tskip]

    # Load all files into the same dataset
    print(f"Loading {len(tskip_files)} SSHA files...")
    ds = xr.open_mfdataset(tskip_files, combine='by_coords')

    return ds

def extract_ssha_timeseries(ds, lat = 30, lon = 135):
    print(f"Extracting SSHA timeseries at {lat = } and {lon = }")

    # Convert the longitude to the range 0-360 if it's in -180 to 180
    if lon < 0:
        lon = 360 + lon

    # Extract SLA values at the given latitude and longitude, and convert it to a numpy array
    # Using method='nearest' to handle case if exact coordinates are not present in the dataset
    measurements = ds[y_key].sel(Latitude=lat, Longitude=lon, method='nearest').values

    # Create xarray dataset
    ds = xr.Dataset(
        {y_key: (x_key, measurements)},
        coords={x_key: ds[x_key].values}
    )

    return ds

def timeseries_to_rgb(timeseries, x_key, y_key, min_period, max_period):
    #print(f"{timeseries.keys() = }")
    timeseries, fits = fancy_detrend(timeseries, x_key, y_key, terms=['constant', 'trend', 'accel'])

    if make_plots: plot_timeseries(timeseries,title="Detrended time series")

    fitted_timeseries = turn_fits_into_timeseries(timeseries, x_key, y_key, fits=fits)

    if make_plots: plot_timeseries(fitted_timeseries,title="Fitted time series")

    # Perform non-uniform FFT to get power spectrum.
    power_spectrum = nfft_power(timeseries)

    power_spectrum = convert_spectrum_from_frequency_to_period(power_spectrum)
    
    #print(f"{min_period = } and {np.max(power_spectrum.period.values) = }")
    if min_period < np.min(power_spectrum.period.values) or min_period >= np.max(power_spectrum.period.values):
        #print(f"!!! WARNING!!! originally {min_period = } but {np.min(power_spectrum.period.values) = } and {np.max(power_spectrum.period.values) = }")
        min_period = np.min(power_spectrum.period.values)
        #print(f"So now {min_period = } which equals {np.min(power_spectrum.period.values) = }")
    #print(f"{max_period = } and {np.max(power_spectrum.period.values) = }")
    if max_period <= np.min(power_spectrum.period.values) or max_period > np.max(power_spectrum.period.values):
        #print(f"!!! WARNING!!! originally {max_period = } but {np.min(power_spectrum.period.values) = } and {np.max(power_spectrum.period.values) = }")
        max_period = np.max(power_spectrum.period.values)
        #print(f"So now {max_period = } which equals {np.max(power_spectrum.period.values) = }")
    signal_period = 365.25
    #power_spectrum = power_spectrum.where((power_spectrum.period > signal_period * 0.2) & (power_spectrum.period < signal_period * 3), drop=True)

    if make_plots: plot_fft_spectrum(power_spectrum,title="FFT Power spectrum")

    cie = load_cie_functions()

    mapped_spectrum = map_power_spectrum(cie, power_spectrum, min_period = min_period, max_period = max_period)
    
    #print(f"{mapped_spectrum = }")
        
    if make_plots: plot_light_spectrum(mapped_spectrum,title="Light spectrum")

    #print(f"{cie = }")
    #spectrum = synthetic_spectrum(cie, 530, 30)
    spectrum = mapped_spectrum
    #print(f"{spectrum = }")
    xyz = spectrum2xyz(spectrum, cie, 1.0)
    #print(f"{xyz = }")
    thepower = 1.0
    #print(f"Before raising y to power {thepower}: {xyz = }")
    raise_y_to_power(xyz, thepower)
    #print(f"After  raising y to power {thepower}: {xyz = }")
    rgb = xyz2rgb(xyz)
    #print(f"Before gamma correction: {rgb = }")
    rgb = gamma_correct_rgb(rgb)
    #print(f"After  gamma correction: {rgb = }")

    if make_plots:
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_color_plot.png")
        plot_color(rgb, filename)
    
    return xr.Dataset(rgb)

# Define your function, that returns a Dataset
def rms_and_mean(x):
    rms_value = np.sqrt(np.mean(x**2))
    mean_value = np.mean(x)
    return xr.Dataset({'rms': rms_value, 'mean': mean_value})

def write_gmt_scripts(plot_options, grid, results):

    # Create GMT script in outputfolder, which should be netcdf_output.
    new_file = os.path.join(plot_options['outputfolder'], 'create_plots.sh')
    try:
        new_fp = open(new_file, 'wb')
    except IOError:
        print("The create_plots.sh GMT script couldn't be created.")

    # Every bash script needs this to be the first line.
    new_fp.write(b"#!/bin/bash\n")
    new_fp.write(b"#set -x #Uncomment to echo these commands.\n")

    # Create flip_backgrounds script in outputfolder, which should be netcdf_output.
    flip_file = os.path.join(plot_options['outputfolder'], 'flip_backgrounds.sh')
    try:
        flip_fp = open(flip_file, 'wb')
    except IOError:
        print("The flip_backgrounds.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    flip_fp.write(b"#!/bin/bash\n")
    flip_fp.write(b"set -x\n")

    # Create trim script in outputfolder, which should be netcdf_output.
    trim_file = os.path.join(plot_options['outputfolder'], 'trim.sh')
    try:
        trim_fp = open(trim_file, 'wb')
    except IOError:
        print("The trim.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    trim_fp.write(b"#!/bin/bash\n")
    trim_fp.write(b"set -x\n")

    # If RGB, only need GMT commands for a single map.
    if 1:#len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
        just_the_filenames = plot_options['just_the_filenames'][:1]

    for i in range(len(just_the_filenames)):
        # If this is the first file...
        if i == 0:
            write_gmt_defs(new_fp)
            new_fp.write(f"color_scheme={plot_options['color_scheme']} #1/2=white/black background\n".encode())
            new_fp.write(f"montage={plot_options['montage']} #1=left-justify titles, add (a),(b), run montage.sh.\n".encode())
            new_fp.write(b"prefixes=('(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' '(u)' '(v)' '(w)' '(x)' '(y)' '(z)')\n")
            new_fp.write(b"index=-1 #Increments on each map, accesses prefixes above for montage.\n")
            new_fp.write(b'png_options=" -P -Tg " #PDF default: -E720, else 300 dpi.\n')
            new_fp.write(b"if [ $montage != 0 ]\nthen\n  png_options=\" -A\"$png_options\nfi\n")
            new_fp.write(b"#Force off-white(dark gray) fore(back)ground color because\n#flip_backgrounds.sh can change the maps' text from\n#black to white, and their backgrounds from white to black.\n")
            new_fp.write(b"gmt set COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=253/253/253\n")
            new_fp.write(f"digits={plot_options['scale_digits']}\n".encode())
            new_fp.write(b"gmt set D_FORMAT=%.${digits}f\n")

        # Define the filenames for PostScript output.
        s = f"{plot_options['output_base']}_{i+1:04d}"
        new_fp.write(b"#######################################################\n")
        if 1:#len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
            new_fp.write(b"data_name=redgreenblue\n")
        else:
            new_fp.write(f'data_name="{just_the_filenames[i]}"\n'.encode())
        new_fp.write(f'plot_base="{s}"\n'.encode())
        new_fp.write(b"let index=$index+1\n")
        new_fp.write(b"#######################################################\n")

        # Record the filename bases in the flip_backgrounds.sh and trim.sh scripts.
        if i == 0:
            if len(just_the_filenames) == 1:
                flip_fp.write(f'all_bases="{s}"\n'.encode())  # First and last base starts and ends the string list.
                trim_fp.write(f'all_bases="{s}"\n'.encode())  # First and last base starts and ends the string list.
            else:
                flip_fp.write(f'all_bases="{s}\n'.encode())  # First base starts the string list.
                trim_fp.write(f'all_bases="{s}\n'.encode())  # First base starts the string list.
        elif i < len(just_the_filenames) - 1:
            flip_fp.write(f'{s}\n'.encode())  # Bases that aren't first or last are unadorned.
            trim_fp.write(f'{s}\n'.encode())  # Bases that aren't first or last are unadorned.
        else:
            flip_fp.write(f'{s}"\n'.encode())  # Last base ends with a " to terminate string list.
            trim_fp.write(f'{s}"\n'.encode())  # Last base ends with a " to terminate string list.

        # Plot data, with title on top and coastlines.
        write_gmt_map_data(results, grid, plot_options, new_fp, results['titles'][i], i)

        new_fp.write(b"convert $png_options $plot_base.ps #Convert PS to PNG format.\n")
        new_fp.write(b"#convert -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n")

        # Delete PS file because it's twice as large (more at -E2000 and 0.25x0.25 global- 2.7MB PDF, 165MB PS!) as the PDF or PNG.
        new_fp.write(b"rm -f $plot_base.ps\n")
        
        # Move cpt file to backup version so it won't be automatically used
        # if this script is executed again... just in case write_gmt_cpt fails:
        # it's good for that failure to be obvious.
        if len(results['rgb']) != 3:
            new_fp.write(b"mv map.cpt Zbackup_cpt_$plot_base.cpt\n")

        # When this plot is finished, store its data_name in previous_data_name
        # so that the next plot (if it's a phase plot with amplitude masking) can
        # access that data.
        new_fp.write(b"previous_data_name=$data_name\n")
    # Trim images.
    new_fp.write(b"#######################################################\n")
    new_fp.write(b". ./trim.sh\n")

    # If requested, change colors using flip_backgrounds.sh.
    new_fp.write(b"#######################################################\n")
    new_fp.write(b"if [ $color_scheme == 2 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  . ./flip_backgrounds.sh\n")
    new_fp.write(b"fi\n")

    # If requested, make a montage using montage.sh.
    new_fp.write(b"#######################################################\n")
    new_fp.write(b"if [ $montage != 0 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  . ./montage.sh\n")
    new_fp.write(b"fi\n")
    new_fp.close()

    # Either way, finish writing flip_backgrounds.sh so it can be used later.
    finish_flip_backgrounds(flip_fp)
    flip_fp.close()

    # Finish writing trim.sh.
    finish_trim(trim_fp)
    trim_fp.close()

    # Create animate script in outputfolder, which should be netcdf_output.
    extra_file = os.path.join(plot_options['outputfolder'], "animate.sh")
    try:
        extra_fp = open(extra_file, 'wb')
    except IOError:
        print("The animate.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    extra_fp.write(b"#!/bin/bash\n")
    extra_fp.write(b"set -x\n")
    extra_fp.write(b"delay=100 #delay in hundredths of a second.\n")
    extra_fp.write(b"#size=\"640x480\"\n")
    extra_fp.write(b"#size=\"800x600\"\n")
    extra_fp.write(b"size=\"1024x768\"\n")
    extra_fp.write(f"output_base=\"{plot_options['output_base']}\"\n".encode())
    extra_fp.write(b"#Imagemagick can also output .mng (animated PNG, not well-supported), but ffmpeg is needed as a delegate for .mp4.\n")
    extra_fp.write(b"#convert -verbose -delay $delay -loop 0 $output_base* -resize $size animation.gif\n")
    extra_fp.write(b"#Or ffmpeg can output .mp4 directly.\n")
    extra_fp.write(b"ffmpeg -f image2 -i $output_base%d.png animation.mp4\n")
    extra_fp.close()

    # Create montage script in outputfolder, which should be netcdf_output.
    extra_file = os.path.join(plot_options['outputfolder'], "montage.sh")
    try:
        extra_fp = open(extra_file, 'wb')
    except IOError:
        print("The montage.sh script couldn't be created.")

    # Every bash script needs this to be the first line.
    extra_fp.write(b"#!/bin/bash\n")
    extra_fp.write(b"set -x\n")
    extra_fp.write(f"output_base=\"{plot_options['output_base']}\"\n".encode())
    extra_fp.write(b"montage $output_base* -geometry +2+2 montage.png\n")
    extra_fp.close()

def write_gmt_coastlines(new_fp):
    """
    Purpose: This function writes the GMT coastlines command(s).
    Input:  
        new_fp - file pointer to current gmt script file
    """
    new_fp.write(b"gmt coast -W$coast_thk/$coast_color $coast_res $range $projection $map_pos $middle >> $plot_base.ps\n")
    new_fp.write(b"if [ $coastlines == 2 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  gmt plot -N $coast_file -: -Sc$coast_thk -W$coast_thk/$coast_color $range $projection $map_pos $middle >> $plot_base.ps\n")
    new_fp.write(b"fi\n")

def is_polar(results: Dict, grid: Dict) -> int:
    """
    Purpose: This function examines results,grid, and returns 0 if global,
            or 1 (2) if it's roughly centered on the north (south) pole.
    Input: results, grid (types are dictionary)
    Output: Return value described above.
    """
    polar = 0  # 0 - global, 1 - north pole, 2 - south pole.
    
    if results['options']['output_choice'] == 1 or results['options']['output_choice'] == 4:
        results['minlat'] = min(grid['lat'])
        results['maxlat'] = max(grid['lat'])
        results['minlon'] = min(grid['lon'])
        results['maxlon'] = max(grid['lon'])
        
    elif results['options']['output_choice'] == 5:
        results['minlat'] = min(results['latlon']['lat'])
        results['maxlat'] = max(results['latlon']['lat'])
        results['minlon'] = min(results['latlon']['lon'])
        results['maxlon'] = max(results['latlon']['lon'])
        
    else:
        print(f"!!!!WARNING!!!!!! results['options']['output_choice'] {results['options']['output_choice']} isn't recognized.")
        
    if results['minlat'] < 0 and results['maxlat'] > 0:
        polar = 0
    elif results['minlat'] > 0 and results['maxlat'] > 0:
        polar = 1
    elif results['minlat'] < 0 and results['maxlat'] < 0:
        polar = 2
    
    return polar

def write_gmt_colorscale(new_fp, results, kml_output):
    """
    Purpose: This function writes the GMT color scale command.
    Input:  
        new_fp - file pointer to current gmt script file
        kml_output - 0 for GMT scale next to plot, 1 for KMZ by itself.
    """
    rgb = {}  # Just to init choice and maxes.

    if len(results['rgb']) != 3:
        new_fp.write(b"cpt_name=\"-Cmap.cpt \"\n")
    else:
        new_fp.write(b"cpt_name=\"-C../../rgb00001.cpt \"\n")
    
    if len(results['rgb']) != 3 or results['rgb_choice'] < 2:
        # If this is intended for KMZ output, don't overlay the scale.
        if kml_output:
            # KMZ version writes to a different file.
            new_fp.write(b"gmt colorbar $cpt_name -L $scale_format $overflow $scale_pos -A $start > scale_$plot_base.ps\n")
            new_fp.write(b"# Print units manually, otherwise they're too close to numbers on scale.\n")
            new_fp.write(b"echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $end >> scale_$plot_base.ps\n")
        else:
            new_fp.write(b"gmt colorbar $cpt_name -L $scale_format $overflow $scale_pos -A $middle >> $plot_base.ps\n")
            new_fp.write(b"# Print units manually, otherwise they're too close to numbers on scale.\n")
            new_fp.write(b"echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n")
    else:
        new_fp.write(f"numwidths={results['max_widths']}\n".encode())
        new_fp.write(b"gmt set TICK_LENGTH 0.3c\n")
        new_fp.write(b"scale_width=$(bc <<< \"scale=5; $scale_width / $numwidths\")\n")
        new_fp.write(b"for (( j=1; j <= $numwidths; j++ ))\n")
        new_fp.write(b"do\n")
        new_fp.write(b"  j_string=$(printf '%%05d' $j)\n")
        new_fp.write(b"  cpt_name=\"-C../rgb$j_string.cpt\"\n")
        new_fp.write(b"  gmt colorbar $cpt_name -L $scale_format $overflow $scale_pos -S -A $middle >> $plot_base.ps\n")
        new_fp.write(b"  # Print units manually, otherwise they're too close to numbers on scale.\n")
        new_fp.write(b"  echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n")
        new_fp.write(b"  # Move next scale to the right and get rid of the tick marks.\n")
        new_fp.write(b"  scale_x=$(bc <<< \"scale=5; $scale_x+$scale_width\")\n")
        new_fp.write(b"  scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n")
        new_fp.write(b"  gmt set TICK_LENGTH 0.0\n")
        new_fp.write(b"done\n")

def write_rgb_colorscale(results, cie, plot_options, verbose):
    # Initialize objects
    i = j = l = 0
    s = [None]*max_length
    new_file = ''
    old = copy = all = results_s()

    if results['rgb_choice'] == 0:
        new_file = plot_options['outputfolder'] + "rgb00001.cpt"
        with open(new_file, 'w') as new_fp:
            new_fp.write(b"# COLOR_MODEL = RGB\n")
            print("Writing RGB colorscale to disk.")
            for i in range(len(results['xy']['x_values'][0][0])):
                copy = results
                create_synthetic_plot(copy,20,0,0,0,results['xy']['x_values'][0][0][i],0.2)
                match_x_axes(copy, cie)
                interpolate(copy, cie)
                spectrum2xyz(copy, cie, verbose)
                xyz2rgb(copy)
                rescale_single_rgb(copy, verbose)
                gamma_correct_single_rgb(copy, verbose)
                for l in range(3): copy.rgb[l] *= 255.0
                if i > 0:
                    new_fp.write(f"{results['xy']['x_values'][0][0][i-1]:12.6e} {int(old.rgb[0]):3d} {int(old.rgb[1]):3d} {int(old.rgb[2]):3d} {results['xy']['x_values'][0][0][i]:12.6e} {int(copy.rgb[0]):3d} {int(copy.rgb[1]):3d} {int(copy.rgb[2]):3d}\n".encode())
                old = copy
    elif results['rgb_choice'] == 1:
        # Create CPT file in outputfolder, which is NOT netcdf_output bc that doesn't exist yet.
        new_file = os.path.join(plot_options['outputfolder'], 'rgb00001.cpt')
        try:
            with open(new_file, 'w') as new_fp:
                new_fp.write(b"# COLOR_MODEL = RGB\n")
                print("Writing RGB colorscale to disk.")
                
                hw = 2 * (results['xy']['x_values'][0][0][1] - results['xy']['x_values'][0][0][0])
                all['options']['output_choice'] = 5
                all['latlon']['lat'] = [0]*len(results['xy']['x_values'][0][0])
                all['latlon']['lon'] = [0]
                j = 0  # Because there's only 1 lon.
                init_latlon(all,3)
                
                # Loop through periods for colorscale divisions.
                for i in range(len(results['xy']['x_values'][0][0])):
                    # Calculate color for each period.
                    copy = deepcopy(results)
                    create_synthetic_plot(copy,20,0,0,0,results['xy']['x_values'][0][0][i],hw) # choice,#pts,x0,delta,param1,2.
                    match_x_axes(copy,cie)
                    interpolate(copy,cie)
                    spectrum2xyz(copy,cie,verbose)
                    xyz2rgb(copy)
                    for l in range(3):
                        all['latlon']['outputs'][l][i][j] = copy.rgb[l]
                
                spectra_end(all,2000,1) # mod_choice,norm_rgb
                
                for i in range(1,len(results['xy']['x_values'][0][0])):
                    new_fp.write(b"%12.6e %3d %3d %3d %12.6e %3d %3d %3d\n" % (
                        results['xy']['x_values'][0][0][i-1], int(all['latlon']['outputs'][0][i-1][j]), int(all['latlon']['outputs'][1][i-1][j]), int(all['latlon']['outputs'][2][i-1][j]),
                        results['xy']['x_values'][0][0][i], int(all['latlon']['outputs'][0][i][j]), int(all['latlon']['outputs'][1][i][j]), int(all['latlon']['outputs'][2][i][j])
                    ))
        except Exception as e:
            print(f"The rgb00001.cpt file couldn't be created due to the following error: {e}")

    else: print(f"!!!WARNING!!! Didn't recognize {results['rgb_choice'] = }")

def write_gmt_map_data(results, grid, plot_options, new_fp, title, i):
    """
    Purpose: This function writes the GMT data plotting command.
            If the title is blank, this is assumed to be a KMZ map so
            the map fills the page in the cylindrical equidistant
            projection that works in Google Earth.
    Input:  new_fp - file handle to current gmt script file
            title - string to be plotting at top of map data. If blank,
                    title is suppressed and Google Earth KMZ is output.
            i - (output) index of parameter being mapped.
    """

    # Is this a polar or global plot?
    polar = is_polar(results, grid)  # 0 - global, 1 - north pole, 2 - south pole.

    if i == 0:
        new_fp.write(b"#Set resolution, coast_file, coast_thickness, and coastlines\n")
        new_fp.write(b"#on first map only because they should be universal.\n")
        coast_file = plot_options['outputfolder'] + "data/ancillary/Rignot/InSAR_GL_Antarctica.txt"
        new_fp.write(f'coast_file="{coast_file}"\n'.encode())

        # Only latlon data determines resolution using delta_lat.
        if results['options']['output_choice'] == 5:
            delta_lat = abs(results['latlon']['lat'][1] - results['latlon']['lat'][0])
            if delta_lat < 0.4:  # Latlon lat spacing controls the resolution.
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Df+ "\n')
            elif delta_lat < 0.9:  # Latlon lat spacing controls the resolution.
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Df+ "\n')
            else:
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Di+ "\n')
        elif results['options']['output_choice'] in [1, 4]:
            new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
            new_fp.write(b'coast_res=" -Di+ "\n')
        else:
            print(f"!!!!WARNING!!!!!! results['options']['output_choice'] {results['options']['output_choice']} isn't recognized.")

        new_fp.write(b'coast_res_orig=$coast_res #Don\'t want USA maps to repeatedly add -N2.\n')
        new_fp.write(b'coast_thk="0.6"\n')
        new_fp.write(b'coast_thk="0.009"\n')
        if polar == 1:
            coastlines = 1  # InSAR is only in Antarctica, so disable for NP plots.
        new_fp.write(f'coastlines={coastlines} #1:coast, 2:coast+InSAR.\n'.encode())

    new_fp.write(b"#coast_color is gray82 for off-white, or gray10 for dark coastlines.\n")

    # Adjust max/min latitudes for mapping points, otherwise points on edge aren't visible.
    if results['options']['output_choice'] in [1, 4]:
        buffer = 5
        if results['maxlat'] <= 90 - buffer:
            results['maxlat'] += buffer
        if results['minlat'] >= -90 + buffer:
            results['minlat'] -= buffer

    # Record text formats, which are the same for all projections and data types.
    new_fp.write(b'title_format="0 0 30 0 0 MC"\n')
    new_fp.write(b'blurb_format="0 0 15 0 1 ML"\n')
    new_fp.write(b'units_format="0 0 13 0 0 MC"\n')
    # Record units for the scale, which are the same for all projections and data types.
    new_fp.write(f"scale_units=\"{results['units'][i]}\"\n".encode())

    new_fp.write(b'misc_range=" -R0/1/0/1 -JX1c "\n')
    new_fp.write(b"#grdcut requires actual limits, but if grdimage uses them: GMT Fatal Error: grdimage could not allocate memory [21.69 Gb, n_items = 5823567396]\n")
    new_fp.write(b'minlon=%.3f\n' % 0.0)  # results['minlon']
    new_fp.write(b'maxlon=%.3f\n' % 360.0)  # results['maxlon']
    new_fp.write(b'minlat=%.3f\n' % -90.0)  # results['minlat']
    new_fp.write(b'maxlat=%.3f\n' % 90.0)  # results['maxlat']

    if title:
        # All projections are always available.
        if i == 0:  # Only print this guide for the first map.
            new_fp.write(b"#Global projections:\n")
            new_fp.write(b"#    1 - Robinson\n")
            new_fp.write(b"#    2 - Winkel Tripel\n")
            new_fp.write(b"#    3 - Mollweide\n")
            new_fp.write(b"#    4 - Miller\n")
            new_fp.write(b"#Polar projections:\n")
            new_fp.write(b"#  101 - N. Azimuthal Equidistant\n")
            new_fp.write(b"#  102 - S. Azimuthal Equidistant\n")
            new_fp.write(b"#Specific regions:\n")
            new_fp.write(b"# 1001 - North America\n")
            new_fp.write(b"# 1002 - South America\n")
            new_fp.write(b"# 1003 - Africa\n")
            new_fp.write(b"# 1004 - Greenland\n")
            new_fp.write(b"# 1005 - South Asia\n")
            new_fp.write(b"# 1006 - Australia\n")
            new_fp.write(b"# 1007 - Europe\n")
            new_fp.write(b"# 1101 - Contiguous United States\n")
            new_fp.write(b"# 1102 - California\n")

        if polar == 0:
            new_fp.write(f"{'#' if i > 0 else ''}projection_choice={plot_options.projection}\n".encode())
        else:
            if polar == 1:  # North pole.
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=101\n".encode())
            else:  # South pole.
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=102\n".encode())

        new_fp.write(b"standard_circle=0 #1=all specific regions use standard circular projection.\n")
        new_fp.write(b"standard_rect=0 #1=all specific regions use standard rectangular projection.\n")
        new_fp.write(b". ./projections.sh\n")
        new_fp.write(b"if [ $projection_choice == 101 ]\n")
        new_fp.write(b"then\n")

        minlat = results['minlat'] if polar == 1 else 0.0
        new_fp.write(f"  minlat={minlat:.3f}\n".encode())
        new_fp.write(b"  actual_range=\" -R0.0/360.0/$minlat/90.0 \"\n")
        polar_radius = 90 - minlat
        new_fp.write(f"  polar_radius={polar_radius}\n".encode())
        new_fp.write(b"  projection=\" -JE0/90.0/${polar_radius}/${map_width}c \" #N. Azimuthal Equidistant\n")
        new_fp.write(b"elif [ $projection_choice == 102 ]\n")
        new_fp.write(b"then\n")

        maxlat = results['maxlat'] if polar == 2 else 0.0
        new_fp.write(f"  maxlat={maxlat:.3f}\n".encode())
        new_fp.write(b"  actual_range=\" -R0.0/360.0/-90.0/$maxlat \"\n")
        polar_radius = 90 + maxlat
        new_fp.write(f"  polar_radius={polar_radius}\n".encode())
        new_fp.write(b"  projection=\" -JE0/-90.0/{polar_radius}/${map_width}c \" #S. Azimuthal Equidistant\n")
        new_fp.write(b"fi\n")

        new_fp.write(b"range=\" -R${minlon}/${maxlon}/${minlat}/${maxlat} \"\n")
        new_fp.write(b"map_pos=\" -Xa${map_x}c -Ya${map_y}c \"\n")

        if len(results['rgb']) == 3 and results['rgb_choice'] >= 2:
            new_fp.write(b"scale_width=1.2 #Override for RGB maps.\n")

        new_fp.write(b"scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n")
        new_fp.write(b"units_x=$(bc <<< \"scale=5; $scale_x+$scale_width/2\")\n");
        new_fp.write(b"units_y=$(bc <<< \"scale=5; $scale_y+$scale_length/2\")\n");
        new_fp.write(b"units_pos=\" -Xa${units_x}c -Ya${units_y}c \"\n")
        new_fp.write(b"blurb_pos=\" -Xa${blurb_x}c -Ya${blurbs_y}c \"\n")
        new_fp.write(b"blurb2_pos=\" -Xa${blurb2_x}c -Ya${blurbs_y}c \"\n")
    else: print("!!!WARNING!!! NO TITLE!")

    # Plot data, with title on top.
    new_fp.write(f"title=\"{title}\"\n".encode())
    if len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
        new_fp.write(b"gmt grdimage red.nc green.nc blue.nc $boundary $resolution $range $projection $map_pos $start > $plot_base.ps\n")
    else:
        new_fp.write(b"gmt grdimage $data_name $boundary $resolution $range $projection $map_pos -Cmap.cpt $start > $plot_base.ps\n")

    write_gmt_coastlines(new_fp)

    # If marker_lats multivectors have the right size, plot markers in the color "sandy brown".
    # Plot markers before mascons because mascons are smaller than markers.
    if len(results['marker_lats']) == len(results['latlon']['outputs']) and len(results['marker_lons']) == len(results['latlon']['outputs']) and results['latlon']['outputs']:
        if results['marker_lats'][i] and results['marker_lons'][i]:
            new_fp.write(b"gmt plot -N $data_name -bcmarker_lons/marker_lats -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n")

    new_fp.write(b"#Uncomment to put a marker at echoed coords, given as lon lat:\n")
    new_fp.write(b"#echo -85.19 -77.36 | gmt plot -N -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n")

    # If requested and mascon_lats/lons vectors aren't empty, plot mascon centers in the color "saddle brown".
    if plot_options['plot_mascons'] != 0 and results['latlon']['mascon_lats']:
        new_fp.write(b"gmt plot $data_name -bcmascon_lons/mascon_lats -Sc0.01c -G139/69/19 $range $projection $map_pos $middle >> $plot_base.ps\n")

    new_fp.write(b"if [ $montage != 0 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  title=${prefixes[$index]}\" \"$title\n")
    new_fp.write(b"  title_format=\"0 0 30 0 0 ML\" #Left-justify so montage titles are uniform.\n")
    new_fp.write(b"  title_x=$(bc <<< \"scale=5; $blurb_x-0.1\")\n")
    new_fp.write(b"else\n")
    new_fp.write(b"  title_x=$(bc <<< \"scale=5; $map_x+$map_width/2\")\n")
    new_fp.write(b"fi\n")
    new_fp.write(b"title_pos=\" -Xa${title_x}c -Ya${title_y}c \"\n")
    new_fp.write(b"echo $title_format $title | gmt pstext -N $title_pos $misc_range $middle >> $plot_base.ps\n")

    # Draw color scale with units printed above.
    write_gmt_colorscale(new_fp, results, 0)  # 0 = next to map for GMT PDF output.

    # Print blurb about data range, or masked amplitudes for phase plots.
    # Don't print blurb if disabled or if this is an RGB map.
    if len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
        plot_options['blurb_disabled'] = 1
    if plot_options['blurb_disabled']:
        new_fp.write(b"blurb_contents=\"\"\n")

    new_fp.write(b"echo $blurb_format $blurb_contents | gmt pstext -N $blurb_pos $misc_range $middle >> $plot_base.ps\n")

    # Only print error bars if they're the right length, and this one is > 0.0.
    # Need to have separate copies for unstructured and latlon maps because
    # in each case the "right length" is defined by a different multivector.
    blurb2_written = 0
    if len(results['latlon']['outputs']) == len(results['error_bars']) and results['latlon']['outputs']:
        if results['error_bars'][i] > 0.0:
            blurb2_written = 1
            new_fp.write(f"blurb2_contents=\"Error bar: {results['error_bars'][i]:.1f} $scale_units\"\n".encode())

    if not blurb2_written:
        new_fp.write(b"blurb2_contents=\"\" #Error bar: N/A $scale_units\n")

    new_fp.write(b"echo $blurb_format $blurb2_contents | gmt pstext -N $blurb2_pos $misc_range $end >> $plot_base.ps\n")

def write_gmt_defs(new_fp):
    """
    Purpose: This function writes some clarifying definitions that help me to consistently write working GMT data plotting commands.
    Input:  new_fp - file handle to current gmt script file
    """
    new_fp.write(b"#######################################################\n")
    new_fp.write(b"#Clarifying definitions. Do not change!################\n")
    new_fp.write(b"start=\" -K \" #Should always redirect using > to write new PS.\n")
    new_fp.write(b"middle=\" -O -K \" #Should always redirect using >> to append to PS.\n")
    new_fp.write(b"end=\" -O \" #Should always redirect using >> to append to PS.\n")
    new_fp.write(b"#######################################################\n")

def finish_flip_backgrounds(flip_fp):
    """
    Purpose: This function finishes writing the flip_backgrounds.sh script.

    Input:  
        flip_fp - file pointer to flip_backgrounds.sh script file
    """
    flip_fp.write(b"for current_base in $all_bases\n")
    flip_fp.write(b"do\n")
    flip_fp.write(b"  #Change black to cyan temporarily.\n")
    flip_fp.write(b"  convert $current_base.png -fill cyan -opaque black $current_base.png\n")
    flip_fp.write(b"  #Change white to black.\n")
    flip_fp.write(b"  convert $current_base.png -fill black -opaque white $current_base.png\n")
    flip_fp.write(b"  #Change temporary cyan to white.\n")
    flip_fp.write(b"  convert $current_base.png -fill white -opaque cyan $current_base.png\n")
    flip_fp.write(b"done\n")

def finish_trim(trim_fp):
    """
    Purpose: This function finishes writing the trim.sh script.

    Input:  
        trim_fp - file pointer to trim.sh script file
    """
    trim_fp.write(b"for current_base in $all_bases\n")
    trim_fp.write(b"do\n")
    trim_fp.write(b"  convert $current_base.png -trim $current_base.png\n")
    trim_fp.write(b"done\n")

if __name__ == "__main__":
    plt.style.use('dark_background')
    plt.rcParams['font.size'] = 14  # Change the global font size
    plt.rcParams['axes.linewidth'] = 2  # Change the global linewidth
    myfigsize=(10,5)

    make_plots = 0

    x_key = 'Time'
    y_key = 'SLA'
    min_period = 220
    max_period = 2000

    # Calculate the wavelength ratio
    #wavelength_ratio = cie['wavelength'].max() / cie['wavelength'].min()
    # Calculate the max_period
    #max_period = min_period * wavelength_ratio

    ds = load_ssha_files(tskip=1)
    xskip = 6#96#192
    print(f"Grabbing one lat/lon point in every {xskip**2} points...",end="")
    ds = ds.isel(Latitude=slice(None, None, xskip), Longitude=slice(None, None, xskip))
    print(" done.")
    #breakpoint()
    
    # Check if the size of x_key dimension is odd
    if ds.dims[x_key] % 2 == 1:
        print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {ds.dims[x_key] = }")
        print(f"!!! DELETING LAST DATA POINT!")
        # If it is, select all elements up to the second last one
        ds = ds.isel({x_key: slice(None, -1)})

    # Time the code
    start_time = timeit.default_timer()
    print("Starting map operation WITHOUT dask...")

    # Stack 'Latitude' and 'Longitude' into a new single dimension 'position'
    stacked = ds.stack(position=['Latitude', 'Longitude'])

    # Create a new function with min_period and max_period filled
    timeseries_to_rgb_partial = functools.partial(timeseries_to_rgb, x_key=x_key, y_key=y_key, min_period=min_period, max_period=max_period)

    # Apply the function to the 'SLA' variable of the stacked dataset
    result = stacked.groupby('position').map(timeseries_to_rgb_partial)

    # Now `result` is a Dataset of Datasets, create a dictionary of unstacked datasets
    datasets = {key: result[key].unstack('position') for key in result.data_vars}

    # Stop the timer
    end_time = timeit.default_timer()

    # Calculate the time taken
    time_taken = end_time - start_time

    print(f"Time taken WITHOUT DASK: {time_taken:.2f} seconds")

    if 0:
        # Assuming `ds` is your xarray Dataset
        sla_rms = ds['SLA_RMS']
        title = 'SLA RMS'
        # Create the figure and axes objects
        fig = plt.figure(figsize=myfigsize)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        # Make the plot
        sla_rms.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis', cbar_kwargs={'label': 'RMS', 'orientation': 'horizontal', 'pad': 0.1})
        # Add gridlines and labels
        ax.coastlines()
        ax.gridlines(draw_labels=True)
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder, f"{date_str}_{title.replace(' ','_')}.png")
        # Save the figure with the desired options
        plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    for thekey in datasets.keys():
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder, f"{thekey.replace(' ', '_')}.nc")
        
        # Save the dataset to a NetCDF file
        print(f"Saving {filename}...")
        datasets[thekey].to_netcdf(filename)
    print("Finished saving.")
    
    write_gmt_scripts(plot_options, grid, results)
    
    #Internal Docker folder:
    docker_internal_folder = '/home/jovyan'
    # Define Docker internal script
    docker_internal_filename = 'docker_internal.sh'
    docker_internal_script = "#!/bin/bash\n" \
                             "cp create_plots.sh Zbackup_create_plots.sh\n" \
                             "cp projections.sh Zbackup_projections.sh\n" \
                             "./create_plots.sh\n"
    docker_internal_script = "#!/bin/bash\n" \
                             "echo '----- Environment Variables -----'\n" \
                             "echo 'PATH: '\n" \
                             "echo $PATH\n" \
                             "echo 'LD_LIBRARY_PATH: '\n" \
                             "echo $LD_LIBRARY_PATH\n" \
                             "echo '----- Conda Environments -----'\n" \
                             "/opt/conda/bin/conda env list\n" \
                             "echo '----- Working Directory -----'\n" \
                             "pwd\n" \
                             "cp create_plots.sh Zbackup_create_plots.sh\n" \
                             "cp projections.sh Zbackup_projections.sh\n" \
                             "./create_plots.sh\n"
    
    # Define Docker external command
    docker_command = f'docker run -e TZ=America/Los_Angeles -v {outputfolder}:{docker_internal_folder} -w {docker_internal_folder} grace/testing-bpr-grace2 {docker_internal_folder}/{docker_internal_filename}'
    docker_command = f'docker run -it -e TZ=America/Los_Angeles -v .:{docker_internal_folder} -w {docker_internal_folder} grace/testing-bpr-grace2 /bin/bash'

    docker_external_filename = 'docker_external.bat'
    #Used to have this line in the docker_external_script before the docker command: @echo off
    docker_external_script = f"""
    {docker_command}
    """

    # Write Docker internal script
    with open(os.path.join(outputfolder, docker_internal_filename), 'wb') as file:
        file.write(docker_internal_script.encode())

    # Write Docker external script
    with open(os.path.join(outputfolder,docker_external_filename), 'w') as file:
        file.write(docker_external_script)

    # Change permissions of the Docker internal script
    os.chmod(os.path.join(outputfolder, docker_internal_filename), 0o755)

    # Run Docker external script
    print(f"Opening an interactive shell. Run {docker_external_filename}, then run {docker_internal_filename}")
    # Go to the outputfolder directory
    os.chdir(outputfolder)
    # Start an interactive shell
    os.system('cmd')

    print("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
