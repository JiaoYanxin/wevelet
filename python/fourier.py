from scipy.fft import fft2, fftshift
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Load the image from the file system
image_path = '/home/jiaoyx/图片/截图/截图 2024-04-20 17-12-33.png'
image = Image.open(image_path).convert('L')  # convert to grayscale

# Perform a 2D Fourier transform
f_transform = fft2(np.array(image))
f_shifted = fftshift(f_transform)  # Shift the zero frequency component to the center
magnitude_spectrum = 20*np.log(np.abs(f_shifted))

# Plot the magnitude spectrum
plt.figure(figsize=(6, 6))
plt.imshow(magnitude_spectrum, cmap='gray')
plt.title('Magnitude Spectrum')
plt.axis('off')  # Turn off the axis
plt.show()