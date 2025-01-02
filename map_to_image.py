import sunpy.map
import matplotlib.pyplot as plt

def sunpy_map_to_jpg(sunpy_map, output_file):
    """
    Convert a SunPy map to a JPEG image.

    Parameters:
    sunpy_map (sunpy.map.Map): The SunPy map object to convert.
    output_file (str): The path to the output JPEG file.
    """
    # Check if the input is a SunPy map
    if not isinstance(sunpy_map, sunpy.map.Map):
        raise TypeError("Input must be a SunPy map object")
    
    # Create a figure
    fig = plt.figure()
    
    # Plot the SunPy map
    ax = fig.add_subplot(111, projection=sunpy_map)
    sunpy_map.plot(ax)
    
    # Remove the colorbar and axes for a cleaner image
    plt.colorbar().remove()
    ax.set_axis_off()
    
    # Save the figure to a JPEG file
    plt.savefig(output_file, format='jpeg', bbox_inches='tight', pad_inches=0)
    
    # Close the figure to free memory
    plt.close(fig)

# Example usage:
# map = sunpy.map.Map('path_to_fits_file.fits')
# sunpy_map_to_jpg(map, 'output_image.jpg')
