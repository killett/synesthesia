import csv
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

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

import numpy as np

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

def plot_color(rgb, filename):
    fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=dpi_choice)

    # Set the facecolor using the normalized RGB values
    ax.set_facecolor(tuple(rgb[key] for key in ['r', 'g', 'b']))

    # Remove all axes and labels
    ax.axis('off')

    # Save the figure with the desired options
    plt.savefig(filename, dpi=dpi_choice, format='png', transparent=True, bbox_inches='tight')

if __name__ == "__main__":
    cie = load_cie_functions()
    print(f"{cie = }")    
    spectrum = synthetic_plot(cie, 530, 30)
    print(f"{spectrum = }")
    xyz = spectrum2xyz(spectrum,cie, 1.0)
    print(f"{xyz = }")
    rgb = xyz2rgb(xyz)
    print(f"{rgb = }")

    # Create filename with current date
    date_str = datetime.now().strftime('%Y%m%d-%H%M%S')
    filename = os.path.join(outputfolder,f"{date_str}_color_plot.png")
    
    plot_color(rgb, filename)
