import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
from glob import glob

lc_array=[]
tm_array=[]

def load_aia131_images(file_pattern):
    """Load AIA 131 images from FITS files."""
    print('Listing images')
    file_list = sorted(glob(file_pattern))
    print(len(file_list))
    maps = [sunpy.map.Map(f) for f in file_list]
    return maps

def running_difference(maps):
    """Compute running difference images."""
    print('Running diff ..')
    diff_maps = []
    for i in range(2, len(maps)):
        diff_data = maps[i].data - maps[i-1].data
        diff_data[diff_data<0]=0
        lc_array.append(np.sum(diff_data))
        tm_array.append(maps[i].date)
        diff_maps.append(sunpy.map.Map(diff_data, maps[i].meta))
    print('done')
    return diff_maps

def plot_running_difference(diff_maps, output_folder="running_diff_output"):
    """Plot and save running difference images."""
    import os
    os.makedirs(output_folder, exist_ok=True)
    
    for i, diff_map in enumerate(diff_maps):
        plt.figure(figsize=(8, 8))
        plt.imshow(diff_map.data, origin='lower', vmin=0, vmax=100)
        plt.colorbar(label="Intensity Difference")
        plt.title(f"AIA 131 ({diff_map.date})")
        plt.savefig(f"{output_folder}/f{diff_map.date}.png")
        plt.close()

# Example usage
file_pattern = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/AIA_Images/aia.*131.image_lev1.fits"  # Adjust based on your file naming pattern
maps = load_aia131_images(file_pattern)
diff_maps = running_difference(maps)
plot_running_difference(diff_maps)
print(lc_array)
np.savetxt('aia131_running_diff.csv',np.c_[tm_array,lc_array],delimiter=',',fmt='%s')
