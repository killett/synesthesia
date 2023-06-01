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

def synthetic_plot(cie,mu,sig):

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
                         time_stop=datetime.datetime(2020, 1, 1)):
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
        coords={
            'time': dates
        }
    )

    return ds

def nfft_power(ds):
    # Convert datetime index to numeric (we use 'day' as the unit)
    x = (ds.time - ds.time[0]).values.astype(float) / (24*3600*1e9)
    print(f"{x = }")

    # number of sample points
    N = len(x)

    x_min = np.min(x)
    x_range = np.max(x) - np.min(x)
    x_norm = (x - x_min) / x_range - 0.5

    # Define Fourier modes
    k = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range
    
    f_k = nfft.nfft(x,ds.measurements)

    #Compute power spectrum, which is the square of the absolute value of the Fourier Transform
    power_spectrum = np.abs(f_k)**2

    #Only take the positive frequencies. Since the output is symmetric, this will not lose any information.
    power_spectrum = power_spectrum[N//2:]
    xf_half = xf[N//2:]

    # Create xarray DataArray with coordinates
    power_spectrum_da = xr.DataArray(power_spectrum, coords=[('frequency', xf_half)], name='power_spectrum')

    # Convert this DataArray to a Dataset
    ds = power_spectrum_da.to_dataset()

    return ds

def plot_timeseries(ds,filename):
    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot measurements
    ds.measurements.plot(ax=ax)

    plt.title('Time series')
    plt.xticks(rotation=30)

    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

def plot_color(rgb, filename):
    fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=dpi_choice)

    # Set the facecolor using the normalized RGB values
    ax.set_facecolor(tuple(rgb[key] for key in ['r', 'g', 'b']))

    # Remove all axes and labels
    ax.axis('off')

    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

if __name__ == "__main__":

    if 1:
        # number of sample points
        N = 1000
        timeseries_length = 10.0

        # Generate N random x values between 0 and timeseries_length
        #x = np.sort(np.random.uniform(0, timeseries_length, N))
        #COMPLETELY UNIFORM POINTS!
        x =  np.linspace(0, timeseries_length, N)
        #NEARLY uniform points!
        x =  np.sort(np.linspace(0, timeseries_length, N) + 2*np.random.random(N))

        x_min = np.min(x)
        x_range = np.max(x) - np.min(x)
        x_norm = (x - x_min) / x_range - 0.5

        # Define Fourier modes
        k = -(N // 2) + np.arange(N)
        # Convert Fourier modes to frequencies
        xf = k / x_range
        
        #print(f"{xf = }")

        #Define based on original x:
        y = 10*np.sin(1.0 * 2.0 * np.pi * x)

        plt.figure(figsize=(10,5))
        plt.plot(x, y, color='red')
        plt.title("Time series")
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_timeseries.png")
        # Save the figure with the desired options
        plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

        yf = np.abs(nfft.nfft(x_norm, y))
        f_k = nfft.nfft(x_norm, y)

        plt.figure(figsize=(10,5))
        plt.plot(xf, f_k.real, label='real', color = 'green')
        plt.plot(xf, f_k.imag, label='imag', color = 'orange')
        plt.legend()
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_fft_test_complex.png")
        # Save the figure with the desired options
        plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

        #Compute power spectrum, which is the square of the absolute value of the Fourier Transform
        power_spectrum = np.abs(f_k)**2

        #Only take the positive frequencies. Since the output is symmetric, this will not lose any information.
        power_spectrum = power_spectrum[N//2:]
        xf_half = xf[N//2:]

        #Plotting the Power Spectrum
        plt.figure(figsize=(10,5))
        plt.plot(xf_half, power_spectrum)
        plt.title('Power spectrum')
        plt.xlabel('Frequency')
        plt.ylabel('Power')
        plt.grid(True)
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_fft_test_power.png")
        # Save the figure with the desired options
        plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

    if 0:
        timeseries = synthetic_timeseries(signal='annual', signal_amplitude=1, noise='white', noise_level=0.1, 
                             temporal_resolution='monthly', time_start=datetime.datetime(2001, 1, 1), 
                             time_stop=datetime.datetime(2020, 1, 1))

        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_timeseries.png")
        plot_timeseries(timeseries,filename)
        # Perform non-uniform FFT to get power spectrum.
        power_spectrum = nfft_power(timeseries)

        plt.figure(figsize=(10,5))
        power_spectrum['power_spectrum'].plot.line('o-') # 'o-' will create a line plot with markers at data points
        plt.title('Power spectrum')
        plt.xlabel('Frequency')
        plt.ylabel('Power')
        plt.grid(True)
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_fft_test_power.png")
        # Save the figure with the desired options
        plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

        cie = load_cie_functions()
        print(f"{cie = }")    
        spectrum = synthetic_plot(cie, 530, 30)
        print(f"{spectrum = }")
        xyz = spectrum2xyz(spectrum,cie, 1.0)
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
