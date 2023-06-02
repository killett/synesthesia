import os
import xarray as xr
import numpy as np
import nfft # pip install nfft From: https://github.com/jakevdp/nfft
import matplotlib.pyplot as plt
import datetime


import timeit

outputfolder = os.path.join('.', 'output')

dpi_choice = 300

if __name__ == "__main__":
    # number of sample points
    N = 10000
    timeseries_length = 100.0
    signal_period = 2

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
    y = np.sin(1/signal_period * 2.0 * np.pi * x)

    plt.style.use('dark_background')
    plt.rcParams['font.size'] = 14  # Change the global font size
    plt.rcParams['axes.linewidth'] = 2  # Change the global linewidth

    plt.figure(figsize=(10,5))
    plt.plot(x, y, color='lime') # using a bright color for visibility
    plt.scatter(x, y, marker='s', color='cyan', s=10) # using a bright color for visibility
    plt.title(f"Time series, {signal_period = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_timeseries.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    # Time the code
    start_time = timeit.default_timer()

    #Take NFFT
    f_k = nfft.nfft(x_norm, y)
    
    # Stop the timer
    end_time = timeit.default_timer()

    # Calculate the time taken
    time_taken = end_time - start_time

    print(f"Time taken: {time_taken} seconds")

    # Time the code
    start_time = timeit.default_timer()

    #Take NFFT
    f_k = nfft.nfft(x_norm, y, tol = 1e-15)
    
    # Stop the timer
    end_time = timeit.default_timer()

    # Calculate the time taken
    time_taken2 = end_time - start_time

    print(f"Time taken: {time_taken2} seconds")
    print(f"Slowdown factor: {time_taken2/time_taken}")
    
    # Ignore zero-frequency term
    xf = np.delete(xf, N // 2)
    f_k = np.delete(f_k, N // 2)
    print(f"{xf = }")
    print(f"{xf[N//2] = }")
    nfft_periods = 1.0 / xf

    if 0:
        print(f"BEFORE SORTING: {nfft_periods = }")
        # Create an array of tuples and sort by period
        sorted_tuples = sorted(zip(nfft_periods, f_k))
        # Separate the sorted tuples back into two arrays
        nfft_periods, f_k = zip(*sorted_tuples)
        # Convert the results back to numpy arrays
        nfft_periods = np.array(nfft_periods)
        f_k = np.array(f_k)
        print(f"AFTER  SORTING: {nfft_periods = }")

    plt.style.use('dark_background')
    plt.figure(figsize=(10,5))
    plt.plot(nfft_periods, f_k.real, label='real', color = 'lime', linewidth=2.0) # set linewidth to increase thickness
    plt.plot(nfft_periods, f_k.imag, label='imag', color = 'cyan', linewidth=2.0) # set linewidth to increase thickness
    plt.legend()
    plt.title(f"Complex NFFT, {signal_period = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_nfft_complex_vs_period.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    # Convert Fourier modes to frequencies
    xf = np.fft.fftfreq(N, d=x[1]-x[0])
    # Take FFT
    f_k2 = np.fft.fft(y)
    # Ignore zero-frequency term
    xf = xf[1:]
    f_k2 = f_k2[1:]
    fft_periods = 1.0 / xf

    if 0:
        print(f"BEFORE SORTING: {fft_periods = }")
        # Create an array of tuples and sort by period
        sorted_tuples = sorted(zip(fft_periods, f_k2))
        # Separate the sorted tuples back into two arrays
        fft_periods, f_k2 = zip(*sorted_tuples)
        # Convert the results back to numpy arrays
        fft_periods = np.array(fft_periods)
        f_k2 = np.array(f_k2)
        print(f"AFTER  SORTING: {fft_periods = }")

    plt.style.use('dark_background')
    plt.figure(figsize=(10,5))
    plt.plot(fft_periods, f_k2.real, label='real', color = 'lime', linewidth=2.0) # set linewidth to increase thickness
    plt.plot(fft_periods, f_k2.imag, label='imag', color = 'cyan', linewidth=2.0) # set linewidth to increase thickness
    plt.legend()
    plt.title(f"Complex FFT, {signal_period = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_fft_complex_vs_period.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')
