from astropy.io import fits
import pandas as pd
import os

# Load the query results from the CSV file
csv_file = 'drms_query_results_hmi.M_45s[2024.11.12_22:15:00_TAI-2024.11.13_00:40:00_TAI].csv'
query_results = pd.read_csv(csv_file)

# Path to the directory containing the downloaded FITS files
fits_dir = '/Analysis/Projects_Data/Flare_Data/Nov13_Flare_Data1/HMI/'

# Loop through the query results and update FITS headers
for index, row in query_results.iterrows():
    # Construct the file path (assuming filenames are based on T_REC or other metadata)
    fits_file = os.path.join(fits_dir, f"hmi.m_45s.{row['T_REC'].replace(':', '').replace('-', '').replace('.','')}.2.magnetogram.fits")
    #print(fits_file)
    

    # Check if the FITS file exists
    if not os.path.exists(fits_file):
        print(f"FITS file not found: {fits_file}")
        continue

    # Open and update the FITS file
    with fits.open(fits_file, mode='update') as hdul:
        header = hdul[1].header

        # Update header with values from the query results
        header['CRPIX1'] = row.get('CRPIX1', 2048)  # Default to 2048 if missing
        header['CRPIX2'] = row.get('CRPIX2', 2048)
        header['CRVAL1'] = row.get('CRVAL1', 0.0)
        header['CRVAL2'] = row.get('CRVAL2', 0.0)
        header['CDELT1'] = row.get('CDELT1', 0.6)
        header['CDELT2'] = row.get('CDELT2', 0.6)
        header['CUNIT1'] = row.get('CUNIT1', 'arcsec')
        header['CUNIT2'] = row.get('CUNIT2', 'arcsec')
        header['CTYPE1'] = row.get('CTYPE1', 'HPLN-TAN')
        header['CTYPE2'] = row.get('CTYPE2', 'HPLT-TAN')
        header['DATE-OBS'] = row.get('DATE-OBS', '2023-01-01T12:00:00')
        header['T_REC'] = row.get('T_REC', 'UTC')
        header['TELESCOP'] = row.get('TELESCOP', 'SDO/HMI')
        header['INSTRUME'] = row.get('INSTRUME', 'HMI')
        header['BUNIT'] = row.get('BUNIT', 'Gauss')
        header['RSUN_OBS'] = row.get('RSUN_OBS', 975.0)
        header['DSUN_OBS'] = row.get('DSUN_OBS', 1.496e11)
        header['CRLN_OBS'] = row.get('CRLN_OBS', 0.0)
        header['CRLT_OBS'] = row.get('CRLT_OBS', 0.0)
        header['CROTA2'] = row.get('CROTA2', 0.0)
        header['WCSNAME']= 'Helioprojective-cartesian'

        # Save changes
        hdul.flush()

    print(f"Updated FITS file: {fits_file}")
