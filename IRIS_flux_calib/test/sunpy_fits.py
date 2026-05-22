from astropy.io import fits
from astropy.time import Time
import numpy as np
from scipy.io import readsav
import os
from glob import glob

def write_iris_fits_cube(s, filename):
    # -------------------------------------------------
    # Original data
    # INT shape: (ny, nx, nlam)
    # -------------------------------------------------
    ny, nx, nlam = s['INT'].shape

    # Reorder to FITS standard: (nlam, ny, nx)
    data = np.transpose(s['INT'], (2, 0, 1))

    solar_x = s['SOLAR_X']          # (nx,)
    solar_y = s['SOLAR_Y']          # (ny,)
    scale_x, scale_y = s['SCALE']

    hdr = fits.Header()

    # -------------------------------------------------
    # Axis definitions
    # -------------------------------------------------
    hdr['NAXIS']  = 3
    hdr['NAXIS1'] = nx     # X
    hdr['NAXIS2'] = ny     # Y
    hdr['NAXIS3'] = nlam   # Wavelength

    # -------------------------------------------------
    # Spatial WCS
    # -------------------------------------------------
    hdr['CTYPE1'] = 'HPLN-TAN'
    hdr['CTYPE2'] = 'HPLT-TAN'
    hdr['CUNIT1'] = 'arcsec'
    hdr['CUNIT2'] = 'arcsec'
    hdr['CDELT1'] = scale_x
    hdr['CDELT2'] = scale_y
    hdr['CRPIX1'] = (nx + 1) / 2
    hdr['CRPIX2'] = (ny + 1) / 2
    hdr['CRVAL1'] = solar_x[nx // 2]
    hdr['CRVAL2'] = solar_y[ny // 2]

    hdr['HGLN_OBS'] = 0.0
    hdr['HGLT_OBS'] = 0.0
    hdr['DSUN_OBS'] = 1.496e11   # meters
    hdr['RSUN_REF'] = 6.96e8     # meters
    hdr['RSUN_OBS'] = 959.63     # arcsec
    hdr['TELESCOP'] = 'IRIS'
    hdr['INSTRUME'] = 'SJI/Raster'

    hdr['WAVELNTH'] = float(np.mean(s['WVL'])) #for sunpy readability
    hdr['WAVEUNIT'] = 'Angstrom' #for sunpy
    # -------------------------------------------------
    # Spectral axis
    # -------------------------------------------------
    hdr['CTYPE3'] = 'WAVE'
    hdr['CUNIT3'] = 'Angstrom'
    hdr['CRPIX3'] = 1.0
    hdr['CRVAL3'] = s['WVL'][0]
    hdr['CDELT3'] = s['WVL'][1] - s['WVL'][0]
    hdr['WCSAXES'] = 3
    hdr['PC1_1'] = 1.0
    hdr['PC2_2'] = 1.0
    hdr['PC3_3'] = 1.0

    # -------------------------------------------------
    # Metadata
    # -------------------------------------------------
    hdr['BUNIT'] = s['UNITS'].decode('utf-8')
    hdr['DATE_OBS'] =s['TIME_CCSDS'][0].decode('utf-8')
    # hdr['DATE-OBS'] = s['DATE_OBS'].decode('utf-8')
    hdr['EXPTIME'] = float(np.mean(s['EXPOSURE_TIME']))

    # -------------------------------------------------
    # Primary HDU with data cube
    # -------------------------------------------------
    hdu_primary = fits.PrimaryHDU(data=data, header=hdr)

    # -------------------------------------------------
    # Wavelength table HDU
    # -------------------------------------------------
    wavelengths = s['WVL']
    col = fits.Column(name='WAVELENGTH', format='D', array=wavelengths)
    hdu_table =fits.BinTableHDU.from_columns([col], name='WAVELENGTH_TABLE')

    # -------------------------------------------------
    # Combine and write to FITS
    # -------------------------------------------------
    hdul = fits.HDUList([hdu_primary, hdu_table])
    hdul.writeto(filename, overwrite=True)


# -------------------------
# Input and output folders
# -------------------------
input_folder = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/IRIS_flux_calib/iris_images"
output_folder = "./sunpy_fits/"
os.makedirs(output_folder, exist_ok=True)

# -------------------------
# Process all .sav files
# -------------------------
sav_files = glob(os.path.join(input_folder, "*.sav"))

for sav_file in sav_files:
    data = readsav(sav_file)
    iris = data['iris_real_units_per_ang']
    s = iris[0]

    base_name = os.path.basename(sav_file).replace(".sav", ".fits")
    out_file = os.path.join(output_folder, base_name)

    write_iris_fits_cube(s, out_file)
    print(f"Saved {out_file}")

print("All SAV files have been converted to FITS with wavelength tables.")
# import sunpy.map
# m = sunpy.map.Map(out_file)
# m.peek()
