import imageio.v2 as imageio  # fixes DeprecationWarning
import os
import numpy as np

# Settings
folder = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/localize_brightenings/Enhancment/NB04/"
output_file = "c7_nb04_video.mp4"
fps = 10

def make_even(image):
    """Crop 1px from right/bottom if dimensions are odd (libx264 requires even dims)."""
    h, w = image.shape[:2]
    return image[:h - (h % 2), :w - (w % 2)]

# Get all PNG files sorted
images = sorted([img for img in os.listdir(folder) if img.endswith(".png")])

if not images:
    raise FileNotFoundError(f"No PNG files found in: {folder}")

print(f"Found {len(images)} frames. Writing video...")

with imageio.get_writer(
    output_file,
    fps=fps,
    codec="libx264",
    pixelformat="yuv420p",
) as writer:
    for filename in images:
        image = imageio.imread(os.path.join(folder, filename))
        image = make_even(image)
        writer.append_data(image)

print(f"Saved video to {output_file}")