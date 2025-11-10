import imageio
import os

# Settings
folder = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/localize_brightenings/Enhancment/NB04/movie/"

output_file = "c8_nb04_video.mp4"
fps = 5

# Get all png files sorted
images = sorted([img for img in os.listdir(folder) if img.endswith(".png")])

# Create video
with imageio.get_writer(output_file, fps=fps) as writer:
    for filename in images:
        image = imageio.imread(os.path.join(folder, filename))
        writer.append_data(image)

print(f"Saved video to {output_file}")

