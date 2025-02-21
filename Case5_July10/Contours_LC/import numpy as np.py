import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.ndimage import gaussian_filter

# Generate synthetic data (or load an image as a 2D NumPy array)
np.random.seed(0)
image = np.random.rand(100, 100)  # Replace with your image array

# Smooth the image to create more meaningful contours
image_smooth = gaussian_filter(image, sigma=2)

# Define contour levels
levels = np.linspace(image_smooth.min(), image_smooth.max(), 5)

# Create the contour plot
fig, ax = plt.subplots()
contour = ax.contour(image_smooth, levels=levels, colors='r')

# Get contour paths
contour_paths = contour.collections

# Create a grid of pixel coordinates
y, x = np.indices(image.shape)
points = np.vstack((x.ravel(), y.ravel())).T  # (N, 2) shape

# Count pixels inside each contour
counts = []
for collection in contour_paths:
    for path in collection.get_paths():
        mask = np.array([path.contains_point(pt) for pt in points])
        count = np.sum(mask)
        counts.append(count)

# Show the image and contours
ax.imshow(image_smooth, cmap="gray")
plt.colorbar()
plt.show()

print(f"Pixel counts inside contours: {counts}")
