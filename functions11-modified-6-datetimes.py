import os
import xarray as xr
import numpy as np
import pandas as pd
import nfft # pip install nfft From: https://github.com/jakevdp/nfft
import matplotlib.pyplot as plt
import datetime

import timeit

outputfolder = os.path.join('.', 'output')

dpi_choice = 300

if __name__ == "__main__":

    plt.style.use('dark_background')
    plt.rcParams['font.size'] = 14  # Change the global font size
    plt.rcParams['axes.linewidth'] = 2  # Change the global linewidth

    # number of sample points
    N = 1000
    timeseries_length = 100.0
    signal_period = 365.25

    #COMPLETELY UNIFORM POINTS!
    #x =  np.linspace(0, timeseries_length, N, endpoint=False)

    time_start=datetime.datetime(2001, 1, 1)
    time_stop=datetime.datetime(2022, 1, 1)
    dates = pd.date_range(start=time_start, end=time_stop, freq='D')

    t = (dates - time_start).days / signal_period
    signal_values = 1.0 * np.sin(2 * np.pi * t)    

    plt.figure(figsize=(10,5))
    plt.plot(range(len(t)),t, color='lime') # using a bright color for visibility
    plt.scatter(range(len(t)),t, marker='s', color='cyan', s=10) # using a bright color for visibility
    plt.title(f"'t', {signal_period = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_t_in_timeseries.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    plt.figure(figsize=(10,5))
    plt.plot(t,signal_values, color='lime') # using a bright color for visibility
    plt.scatter(t,signal_values, marker='s', color='cyan', s=10) # using a bright color for visibility
    plt.title(f"signal_values vs t, {signal_period = }", color='white')
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_signal_values_vs_t_in_timeseries.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')

    print("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
    x = (dates - dates[0]).values.astype(float) / (24*3600*1e9)

    print("!!!! DETREND BY REMOVING CONSTANT, TREND, ACCEL, AND POSSIBLY ANNUAL?")

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
    #f_k = nfft.nfft(x_norm, y, tol = 1e-15) #SLOWER BUT NO DIFFERENCE
    
    # Stop the timer
    end_time = timeit.default_timer()

    # Calculate the time taken
    time_taken = end_time - start_time
    print(f"nfft() took {time_taken} seconds")
        
    # Ignore zero-frequency term
    xf = np.delete(xf, N // 2)
    f_k = np.delete(f_k, N // 2)
    print(f"{xf = }")
    print(f"{xf[N//2] = }")
    nfft_periods = 1.0 / xf

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

    nfft_periods = nfft_periods[N//2:]
    f_k = f_k[N//2:]
    power_spectrum = np.abs(f_k)**2
    
    window = (nfft_periods > signal_period*0.5) & (nfft_periods < signal_period*2)
    #window = nfft_periods < 1e99

    plt.figure(figsize=(10,5))
    plt.plot(nfft_periods[window], power_spectrum[window], color = 'lime', linewidth=2.0)  # set linewidth to increase thickness
    plt.title('Power spectrum')
    plt.xlabel('Period')
    plt.ylabel('Power')
    plt.grid(True, color='gray')  # set grid color to gray for visibility
    # Create filename with current date
    date_str = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder, f"{date_str}_fft_test_power.png")
    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=False, bbox_inches='tight', facecolor='black')
