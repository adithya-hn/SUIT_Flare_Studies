import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u

def save_jpg(sunpy_map, output_file):
    """
    Convert a SunPy map to a JPEG image.

    Parameters:
    sunpy_map (sunpy.map.Map): The SunPy map object to convert.
    output_file (str): The path to the output JPEG file.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection=sunpy_map)
    sunpy_map.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    plt.colorbar().remove()
    ax.set_axis_off()
    plt.savefig(output_file, format='jpeg', bbox_inches='tight', pad_inches=0,dpi=300)
    
    # Close the figure to free memory
    plt.close(fig)

def View_map(sunpy_map,close_map=False):
    fig = plt.figure()
    ax = fig.add_subplot(projection=sunpy_map)
    sunpy_map.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    output_file=sunpy_map.meta.get('F_NAME')[:-4]+'jpg'
    plt.savefig(output_file, format='jpeg',dpi=300)
    if close_map:
        plt.close(fig)
    if close_map is False:
        plt.show()
    

if __name__=="__main__":
    map = sunpy.map.Map('path_to_fits_file.fits')
    save_jpg(map, 'output_image.jpg')
