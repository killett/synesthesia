import csv
import os
import xarray as xr
import numpy as np
import pandas as pd
import nfft # pip install nfft From: https://github.com/jakevdp/nfft
import matplotlib.pyplot as plt
import datetime
import time

outputfolder = os.path.join('.', 'output')

dpi_choice = 300

def load_cie_functions():
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
    rgb = {'r': rgb[0], 'g': rgb[1], 'b': rgb[2]}
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
        {'measurements': ('time', measurements)},
        coords={'time': dates}
    )

    return ds

def nfft_power(ds):
    # Convert datetime index to numeric (we use 'day' as the unit)
    print("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
    x = (ds.time - ds.time[0]).values.astype(float) / (24*3600*1e9)

    #If timeseries has an odd number of points,
    #remove the last data point, then calculate min, range.
    N = 1
    while N % 2:
        N = len(x)
        x_min = np.min(x)
        x_range = np.max(x) - np.min(x)
        x_norm = (x - x_min) / x_range - 0.5
        print(f"{N = }, {x_min = }, {x_range = }")
        if N % 2:
            print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(x) = }")
            print(f"!!! DELETING LAST DATA POINT!")
            x = np.delete(x, -1)

    # Define Fourier modes
    k = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range
    
    #Perform NFFT.
    f_k = nfft.nfft(x_norm,ds.measurements)

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

def old_map_power_spectrum(cie, power_spectrum, min_period):
    # Calculate the wavelength ratio
    wavelength_ratio = cie['wavelength'].max() / cie['wavelength'].min()
    
    # Calculate the max_period
    max_period = min_period * wavelength_ratio
    
    # Normalize the 'period' coordinate in power_spectrum to range from min_period to max_period
    power_spectrum['period'] = (power_spectrum['period'] - power_spectrum['period'].min()) / (power_spectrum['period'].max() - power_spectrum['period'].min()) * (max_period - min_period) + min_period

    # Map 'power' from power_spectrum onto a new wavelength scale before interpolation
    mapped_power_values = np.interp(np.linspace(min_period, max_period, len(cie['wavelength'])), power_spectrum['period'], power_spectrum['power'])

    # Interpolate the power values onto the new wavelength scale
    interpolated_power_values = np.interp(cie['wavelength'].values, np.linspace(min_period, max_period, len(power_spectrum['power'])), power_spectrum['power'])

    # Create a new Dataset that shares the 'wavelength' coordinate with cie, with 'power' as its data variable
    new_power_spectrum = xr.Dataset(
        {'power': (('wavelength',), interpolated_power_values)},
        coords={'wavelength': cie['wavelength']}
    )

    # Plotting the 'power' values before and after interpolation
    plt.figure(figsize=myfigsize)
    # Plot the mapped 'power' values before interpolation
    plt.plot(cie['wavelength'], mapped_power_values, color='lime', linewidth=2.0, label='Before Interpolation')
    plt.scatter(cie['wavelength'], mapped_power_values, marker='s', color='cyan', s=10)
    # Plot the 'power' values after interpolation
    plt.plot(cie['wavelength'], new_power_spectrum['power'], color='red', linewidth=2.0, label='After Interpolation')
    plt.scatter(cie['wavelength'], new_power_spectrum['power'], marker='o', color='orange', s=10)

    plt.title('Power spectrum')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Power')
    plt.grid(True, color='gray')  # set grid color to gray for visibility
    plt.legend()

    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_interpolating_spectrum.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    return new_power_spectrum

def map_power_spectrum(cie, power_spectrum, min_period):
    # calculate the wavelength ratio
    wavelength_ratio = cie['wavelength'].max() / cie['wavelength'].min()
    
    # calculate the max_period
    max_period = min_period * wavelength_ratio

    # Normalize the period in power_spectrum to be within the range min_period to max_period
    period_norm = (power_spectrum['period'] - power_spectrum['period'].min()) / (power_spectrum['period'].max() - power_spectrum['period'].min())
    period_scaled = period_norm * (max_period - min_period) + min_period
    
    # Now scale this normalized period to match the wavelength range in cie
    period_to_wavelength = period_scaled * (cie['wavelength'].max() - cie['wavelength'].min()) + cie['wavelength'].min()
    
    # Plot the power_spectrum['power'] values mapped onto cie wavelengths BEFORE interpolation
    plt.figure(figsize=myfigsize)
    plt.plot(period_to_wavelength, power_spectrum['power'], color='lime', linewidth=2.0, label='Before Interpolation')
    plt.scatter(period_to_wavelength, power_spectrum['power'], marker='s', color='cyan', s=10)
    
    # Interpolate the power_spectrum['power'] onto the cie wavelengths
    power_interpolated = xr.DataArray(
        np.interp(cie['wavelength'], period_to_wavelength, power_spectrum['power']),
        dims=['wavelength'],
        coords={'wavelength': cie['wavelength']}
    )
    new_power_spectrum = xr.Dataset({'power': power_interpolated})

    # Plot the 'power' values after they have been interpolated to the cie wavelengths
    plt.plot(new_power_spectrum['wavelength'], new_power_spectrum['power'], color='red', linewidth=2.0, label='After Interpolation')
    plt.scatter(new_power_spectrum['wavelength'], new_power_spectrum['power'], marker='o', color='orange', s=10)
    
    plt.title('Power spectrum')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Power')
    plt.grid(True, color='gray')  # set grid color to gray for visibility
    plt.legend()

    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_light_spectrum.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    return new_power_spectrum

def plot_timeseries(ds):
    plt.figure(figsize=myfigsize)
    plt.plot(ds['time'], ds['measurements'], color='lime') # using a bright color for visibility
    plt.scatter(ds['time'], ds['measurements'], marker='s', color='cyan', s=10)
    plt.title("Time series", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_timeseries.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

def plot_power_spectrum(power_spectrum):
    plt.figure(figsize=myfigsize)
    plt.plot(power_spectrum.period, power_spectrum['power'], color = 'lime', linewidth=2.0)
    plt.scatter(power_spectrum.period, power_spectrum['power'], marker='s', color='cyan', s=10)
    plt.title('Power spectrum')
    plt.xlabel('Period')
    plt.ylabel('Power')
    plt.grid(True, color='gray')  # set grid color to gray for visibility
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_fft_test_power.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

def plot_light_spectrum(power_spectrum):
    plt.figure(figsize=myfigsize)
    plt.plot(power_spectrum.wavelength, power_spectrum['power'], color = 'lime', linewidth=2.0)
    plt.scatter(power_spectrum.wavelength, power_spectrum['power'], marker='s', color='cyan', s=10)
    plt.title('Power spectrum')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Power')
    plt.grid(True, color='gray')  # set grid color to gray for visibility
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_light_spectrum.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

def plot_color(rgb, filename):
    fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=dpi_choice)

    # Set the facecolor using the normalized RGB values
    ax.set_facecolor(tuple(rgb[key] for key in ['r', 'g', 'b']))

    # Remove all axes and labels
    ax.axis('off')

    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

if __name__ == "__main__":
    plt.style.use('dark_background')
    plt.rcParams['font.size'] = 14  # Change the global font size
    plt.rcParams['axes.linewidth'] = 2  # Change the global linewidth
    myfigsize=(10,5)
    
    min_period = 200

    timeseries = synthetic_timeseries(signal='annual', signal_amplitude=1.0, noise='white', noise_level=0.5, 
                         temporal_resolution='monthly', time_start=datetime.datetime(2001, 1, 1), 
                         time_stop=datetime.datetime(2020, 1, 1))

    N = len(timeseries['time'])
    if N % 2: print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(timeseries['time']) = }")

    plot_timeseries(timeseries)
    
    print("!!!! DETREND BY REMOVING CONSTANT, TREND, ACCEL, AND POSSIBLY ANNUAL?")

    # Perform non-uniform FFT to get power spectrum.
    power_spectrum = nfft_power(timeseries)

    power_spectrum = convert_spectrum_from_frequency_to_period(power_spectrum)
    
    signal_period = 365.25
    power_spectrum = power_spectrum.where((power_spectrum.period > signal_period * 0.5) & (power_spectrum.period < signal_period * 2), drop=True)

    plot_power_spectrum(power_spectrum)

    cie = load_cie_functions()

    mapped_spectrum = map_power_spectrum(cie, power_spectrum, min_period)
    
    print(f"{mapped_spectrum = }")
    
    plot_light_spectrum(mapped_spectrum)

    if 0:
        print(f"{cie = }")    
        spectrum = synthetic_spectrum(cie, 530, 30)
        print(f"{spectrum = }")
        xyz = spectrum2xyz(spectrum, cie, 1.0)
        print(f"{xyz = }")
        thepower = 1.0
        print(f"Before raising y to power {thepower}: {xyz = }")
        raise_y_to_power(xyz, thepower)
        print(f"After  raising y to power {thepower}: {xyz = }")
        rgb = xyz2rgb(xyz)
        print(f"Before gamma correction: {rgb = }")
        rgb = gamma_correct_rgb(rgb)
        print(f"After  gamma correction: {rgb = }")

        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_color_plot.png")
        #plot_color(rgb, filename)
