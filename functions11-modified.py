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
    N = 6940#1000
    timeseries_length = float(N)#10.0
    signal_period = 365.25

    # Generate N random x values between 0 and timeseries_length
    #x = np.sort(np.random.uniform(0, timeseries_length, N))
    #COMPLETELY UNIFORM POINTS!
    x =  np.linspace(0, timeseries_length, N)
    #NEARLY uniform points!
    #x =  np.sort(np.linspace(0, timeseries_length, N) + 0.5*np.random.random(N))

    x_min = np.min(x)
    x_range = np.max(x) - np.min(x)
    x_norm = (x - x_min) / x_range - 0.5
    N = len(x)
    if N % 2: print(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(x) = }")

    # Define Fourier modes
    k = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range
    
    #print(f"{xf = }")

    #Define based on original x:
    y = 100*np.sin(1.0/signal_period * 2.0 * np.pi * x)

    plt.style.use('dark_background')
    plt.rcParams['font.size'] = 14  # Change the global font size
    plt.rcParams['axes.linewidth'] = 2  # Change the global linewidth

    plt.figure(figsize=(10,5))
    plt.plot(x, y, color='lime') # using a bright color for visibility
    plt.scatter(x, y, marker='s', color='cyan', s=10) # using a bright color for visibility
    plt.title("Time series", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_timeseries.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    yf = np.abs(nfft.nfft(x_norm, y))
    f_k = nfft.nfft(x_norm, y)

    plt.style.use('dark_background')
    plt.figure(figsize=(10,5))
    plt.plot(xf, f_k.real, label='real', color = 'lime', linewidth=2.0)  # set linewidth to increase thickness
    plt.plot(xf, f_k.imag, label='imag', color = 'cyan', linewidth=2.0)  # set linewidth to increase thickness
    plt.legend()
    plt.title("Complex NFFT", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_fft_test_complex.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    #Compute power spectrum, which is the square of the absolute value of the Fourier Transform
    power_spectrum = np.abs(f_k)**2

    #Only take the positive frequencies. Since the output is symmetric, this will not lose any information.
    power_spectrum = power_spectrum[N//2:]
    xf_half = xf[N//2:]

    plt.figure(figsize=(10,5))
    plt.plot(xf_half, power_spectrum, color = 'lime', linewidth=2.0)  # set linewidth to increase thickness
    plt.title('Power spectrum')
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.grid(True, color='gray')  # set grid color to gray for visibility
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_fft_test_power.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')
