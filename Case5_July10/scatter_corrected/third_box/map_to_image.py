import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u

def sunpy_map_to_jpg(sunpy_map, output_file):
    """
    Convert a SunPy map to a JPEG image.

    Parameters:
    sunpy_map (sunpy.map.Map): The SunPy map object to convert.
    output_file (str): The path to the output JPEG file.
    """
    # Check if the input is a SunPy map
  
    # Create a figure
    fig = plt.figure()
    
    # Plot the SunPy map
    ax = fig.add_subplot(projection=sunpy_map)
    sunpy_map.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    
    # Remove the colorbar and axes for a cleaner image
    plt.colorbar().remove()
    ax.set_axis_off()
    
    # Save the figure to a JPEG file
    plt.savefig(output_file, format='jpeg', bbox_inches='tight', pad_inches=0,dpi=300)
    
    # Close the figure to free memory
    plt.close(fig)

# Example usage:
# map = sunpy.map.Map('path_to_fits_file.fits')
# sunpy_map_to_jpg(map, 'output_image.jpg')
