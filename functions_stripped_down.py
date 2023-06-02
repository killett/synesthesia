import csv
import os
import xarray as xr
import numpy as np
import pandas as pd
import nfft # pip install nfft From: https://github.com/jakevdp/nfft
import matplotlib.pyplot as plt
import datetime

outputfolder = os.path.join('.', 'output')

dpi_choice = 300

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
    x = (ds.time - ds.time[0]).values.astype(float) / (24*3600*1e9)

    # number of sample points
    N = len(x)
    if N % 2: print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(x) = }")

    x_min = np.min(x)
    x_range = np.max(x) - np.min(x)
    x_norm = (x - x_min) / x_range - 0.5

    # Define Fourier modes
    k = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range
    
    #Perform NFFT.
    f_k = nfft.nfft(x_norm,ds.measurements)

    #Compute power spectrum, which is the square of the absolute value of the Fourier Transform
    power_spectrum = np.abs(f_k)**2

    #Only take the positive frequencies. Since the output is symmetric, this will not lose any information.
    power_spectrum = power_spectrum[N//2:]
    xf_half = xf[N//2:]

    # Create xarray DataArray with coordinates
    power_spectrum_da = xr.DataArray(power_spectrum, coords=[('frequency', xf_half)], name='power_spectrum')

    # Convert this DataArray to a Dataset
    ds = power_spectrum_da.to_dataset()
    ds['frequency'].attrs['units'] = '1/days'
    
    return ds

def convert_spectrum_from_frequency_to_period(ds):
    freq_units = ds['frequency'].attrs.get('units', None)

    # Map of frequency units to period units
    units_map = {'1/days': 'days', 'Hz': 'seconds', '1/years': 'years'}

    # Calculate period as reciprocal of frequency.
    epsilon = 1e-9
    period = 1.0 / (ds['frequency'] + epsilon)
    
    # Replace the 'frequency' coordinate with 'period'
    ds = ds.rename({'frequency': 'period'})
    ds['period'] = period

    # Handle units attribute
    if freq_units in units_map:
        ds['period'].attrs['units'] = units_map[freq_units]
    else:
        print(f"!!!WARNING!!! DID NOT RECOGNIZE {freq_units = }")

    # Reorder the dataset so that period is increasing
    ds = ds.sortby('period')

    return ds

if __name__ == "__main__":

    timeseries = synthetic_timeseries(signal='annual', signal_amplitude=10, noise='white', noise_level=0.0, 
                         temporal_resolution='monthly', time_start=datetime.datetime(2001, 1, 1), 
                         time_stop=datetime.datetime(2055, 1, 1))

    N = len(timeseries['time'])
    if N % 2: print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(timeseries['time']) = }")

    # Perform non-uniform FFT to get power spectrum.
    power_spectrum = nfft_power(timeseries)

    power_spectrum = convert_spectrum_from_frequency_to_period(power_spectrum)

    plt.figure(figsize=(10,5))
    power_spectrum['power_spectrum'].plot.line('o-') # 'o-' will create a line plot with markers at data points
    plt.title('Power spectrum')
    plt.xlabel(power_spectrum['period'].attrs['units'])
    plt.ylabel('Power')
    plt.grid(True)
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_fft_test_power_vs_period.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')
