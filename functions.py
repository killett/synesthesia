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

def perform_nfft(ds):
    # Convert datetime index to numeric (we use 'day' as the unit)
    time_numeric = (ds.time - ds.time[0]).values.astype(float) / (24*3600*1e9)

    fft = nfft.nfft(time_numeric,ds.measurements)
    
    return np.abs(fft)

def plot_timeseries(ds,filename):
    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot measurements
    ds.measurements.plot(ax=ax)

    # Rotate x-axis labels
    plt.xticks(rotation=30)

    # Show the plot
    #plt.show()
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

    if 0:
        N = len(x)
        x_min = np.min(x)
        x_range = np.max(x) - np.min(x)
        x_norm = (x - x_min) / x_range - 0.5

        xf_norm = np.linspace(-0.5, 0.5, M)

        # Define Fourier modes

        k = np.arange(-N/2, N/2)
        # Convert Fourier modes to frequencies
        # Here, T is the total time span of your data
        T = np.max(x) - np.min(x)
        frequencies = k / T

    if 1:
        # number of sample points
        N = 400

        # Simulated non-uniform data
        x = np.linspace(0.0,0.5-0.02, N) + np.random.random((N)) * 0.001
        #print(x)

        y = np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
        yf = np.abs(nfft.nfft(x, y))

        fig, axs = plt.subplots(1)
        fig_f, axs_f = plt.subplots(1)

        axs.plot(x, y, '.', color='red')
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_fft_test_x.png")
        # Save the figure with the desired options
        fig.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

        xf = np.fft.fftfreq(N,1./N)

        axs_f.plot(xf[:int(N/2)], yf[:int(N/2)], color='blue')
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_fft_test_f.png")
        # Save the figure with the desired options
        fig.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

    if 0:
        # number of sample points
        N = 400

        # Simulated non-uniform data
        x = np.linspace(0.0,0.5-0.02, N) + np.random.random((N)) * 0.001
        #print(x)

        #print( 'random' )
        #print( np.random.random((N)) * 0.001 )

        y = np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
        yf = np.abs(nfft.nfft(x, y))

        fig, axs = plt.subplots(1)
        fig_f, axs_f = plt.subplots(1)

        axs.plot(x, y, '.', color='red')

        xf = np.fft.fftfreq(N,1./N)

        axs_f.plot(xf[:int(N/2)], yf[:int(N/2)], color='red')
        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_fft_test_f.png")
        # Save the figure with the desired options
        plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

    if 0:
        x = -0.5 + np.random.rand(1000)
        f = np.sin(10 * 2 * np.pi * x)

        #k = -20 + np.arange(40)
        f_k = nfft.nfft(x, f)

        plt.plot(f_k.real, label='real')
        plt.plot(f_k.imag, label='imag')
        plt.legend()

        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_fft_test1.png")
        # Save the figure with the desired options
        plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight')

    if 0:
        timeseries = synthetic_timeseries(signal='annual', signal_amplitude=1, noise='white', noise_level=0.1, 
                             temporal_resolution='monthly', time_start=datetime.datetime(2001, 1, 1), 
                             time_stop=datetime.datetime(2020, 1, 1))

        # Create filename with current date
        date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        filename = os.path.join(outputfolder,f"{date_str}_timeseries.png")
        #plot_timeseries(timeseries,filename)
        # Perform non-uniform FFT
        fft_results = perform_nfft(timeseries)

        # Print results
        print(fft_results)

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
