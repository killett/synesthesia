import numpy as np
from colour import SpectralDistribution
from colour.colorimetry import get_cmfs, sd_to_XYZ

# Define your spectral power distribution data.
wavelengths = np.arange(360, 781)  # For example, from 360 to 780 nm
values = np.random.random_sample(len(wavelengths))  # Your actual SPD data here

# Create a SpectralDistribution object
spd = SpectralDistribution(values, wavelengths)

# Define the standard observer color matching functions.
cmfs = get_cmfs('CIE 1931 2 Degree Standard Observer')

# Compute the XYZ tristimulus values.
XYZ = sd_to_XYZ(spd, cmfs)

print(XYZ)
