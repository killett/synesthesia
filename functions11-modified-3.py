import os
import xarray as xr
import numpy as np
import nfft # pip install nfft From: https://github.com/jakevdp/nfft
import matplotlib.pyplot as plt
import datetime

outputfolder = os.path.join('.', 'output')

dpi_choice = 300

if __name__ == "__main__":
    # number of sample points
    N = 100
    timeseries_length = 10.0
    frequency = 2.0

    #COMPLETELY UNIFORM POINTS!
    x =  np.linspace(0, timeseries_length, N, endpoint=False)

    x_min = np.min(x)
    x_range = np.max(x) - np.min(x)
    x_norm = (x - x_min) / x_range - 0.5
    N = len(x)
    if N % 2: print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(x) = }")

    # Define Fourier modes
    k = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range
    
    #Define based on original x:
    y = np.sin(frequency * 2.0 * np.pi * x)

    plt.style.use('dark_background')
    plt.rcParams['font.size'] = 14  # Change the global font size
    plt.rcParams['axes.linewidth'] = 2  # Change the global linewidth

    plt.figure(figsize=(10,5))
    plt.plot(x, y, color='lime') # using a bright color for visibility
    plt.scatter(x, y, marker='s', color='cyan', s=10) # using a bright color for visibility
    plt.title(f"Time series, {frequency = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_timeseries.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    #Take NFFT
    f_k = nfft.nfft(x_norm, y)

    plt.style.use('dark_background')
    plt.figure(figsize=(10,5))
    plt.plot(xf, f_k.real, label='real', color = 'lime', linewidth=2.0) # set linewidth to increase thickness
    plt.plot(xf, f_k.imag, label='imag', color = 'cyan', linewidth=2.0) # set linewidth to increase thickness
    plt.legend()
    plt.title(f"Complex NFFT, {frequency = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_nfft_test_complex.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    # Convert Fourier modes to frequencies
    xf = np.fft.fftfreq(N, d=x[1]-x[0])
    # Take FFT
    f_k2 = np.fft.fft(y)

    plt.style.use('dark_background')
    plt.figure(figsize=(10,5))
    plt.plot(xf, f_k2.real, label='real', color = 'lime', linewidth=2.0) # set linewidth to increase thickness
    plt.plot(xf, f_k2.imag, label='imag', color = 'cyan', linewidth=2.0) # set linewidth to increase thickness
    plt.legend()
    plt.title(f"Complex FFT, {frequency = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_fft_test_complex.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')
