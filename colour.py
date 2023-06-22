from PIL import Image
import numpy as np
import colour

# Open image file
im = Image.open('image.png')

# Convert image data to an array
data = np.array(im)

# Normalise the RGB values
rgb = data / 255.0

# Convert RGB to XYZ
xyz = colour.sRGB_to_XYZ(rgb)

# Convert XYZ to Spectral Distribution
sd = colour.recovery.XYZ_to_sd_Rayleigh_Scattering(xyz)
