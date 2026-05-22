
import numpy as np
from astropy.io import fits
from scipy.ndimage import zoom

iris=fits.open('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/IRIS_flux_calib/sunpy_fits/iris_real_units_20241007.fits')

image=iris[0].data
iris_header=iris[0].header

wvl = iris[1].data.astype(np.float64)
mask = (wvl >= lower_wvl) & (wvl <= upper_wvl)
iris_data_masked = iris_data[mask,:,:]
iris_data_sum = np.sum(iris_data_masked,axis=0)

def rebin_by_pixel_scale_flux_conserve(image,  old_scale_x, old_scale_y, new_scale_x, new_scale_y):

    # Compute zoom factors
    zoom_x = old_scale_x / new_scale_x
    zoom_y = old_scale_y / new_scale_y

    # Interpolate
    rebinned = zoom(image, (zoom_y, zoom_x), order=1)

    # Flux conservation correction
    old_pixel_area = old_scale_x * old_scale_y
    new_pixel_area = new_scale_x * new_scale_y

    area_ratio = new_pixel_area / old_pixel_area

    rebinned_flux_conserved = rebinned * area_ratio

    return rebinned_flux_conserved

print(iris_header['CDELT1'],iris_header['CDELT2'],)
# Example
rebinned = rebin_by_pixel_scale_flux_conserve(
    image,
    old_scale_x=iris_header['CDELT1'],
    old_scale_y=iris_header['CDELT2'],
    new_scale_x=0.698,
    new_scale_y=0.698)



# choose a center pixel in original image
cx, cy = 200, 200   # example center

# original 10x10 box
half_size = 5
orig_cut = image[cy-half_size:cy+half_size,
                 cx-half_size:cx+half_size]

# compute scaling
scale_x = abs(iris_header['CDELT1'] / 0.698)
scale_y = abs(iris_header['CDELT2'] / 0.698)

# corresponding center in rebinned image
cx_new = int(cx * scale_x)
cy_new = int(cy * scale_y)

# new half size (scaled)
half_size_x = int(half_size * scale_x)
half_size_y = int(half_size * scale_y)


# extract same physical area
rebinned_cut = rebinned[cy_new-half_size_y:cy_new+half_size_y,cx_new-half_size_x:cx_new+half_size_x]
