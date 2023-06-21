import numpy as np
from colour.colorimetry import MSDS_CMFS_STANDARD_OBSERVER
from colour import SpectralDistribution, sd_to_XYZ

# Define your spectral power distribution data.
wavelengths = np.arange(360, 781)  # For example, from 360 to 780 nm
values = np.random.random_sample(len(wavelengths))  # Your actual SPD data here

# Create a SpectralDistribution object
spd = SpectralDistribution(values, wavelengths)

# Define the standard observer color matching functions.
cmfs = MSDS_CMFS_STANDARD_OBSERVER['CIE 1931 2 Degree Standard Observer']

# Compute the XYZ tristimulus values.
XYZ = sd_to_XYZ(spd, cmfs)

print(XYZ)
