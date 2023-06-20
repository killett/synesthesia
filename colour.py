import numpy as np
import colour

# Define your spectral power distribution data.
wavelengths = np.arange(360, 781)  # For example, from 360 to 780 nm
values = np.random.random_sample(len(wavelengths))  # Your actual SPD data here

# Create a dictionary that maps wavelengths to values
spd_data = dict(zip(wavelengths, values))

# Create a SpectralDistribution object
spd = colour.SpectralDistribution(spd_data)

# Define the standard observer color matching functions.
cmfs = colour.STANDARD_OBSERVERS_CMFS['CIE 1931 2 Degree Standard Observer']

# Compute the XYZ tristimulus values.
XYZ = colour.spectral_to_XYZ(spd, cmfs)

print(XYZ)
