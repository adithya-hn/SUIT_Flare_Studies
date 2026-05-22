#####################################################
# @Author: Abhilash Sarwade
# @Date:   2026-01-01 07:22:26 pm
# @email: sarwade@ursc.gov.in
# @File Name: deadtime.py
# @Project: solexs_tools
#
# @Last Modified time: 2026-02-04 09:38:53 pm
#####################################################

import os, argparse
import numpy as np
from astropy.io import fits
from .caldb_utils import get_caldb_file
from scipy.special import lambertw 

def get_deadtime_params(filter_sdd,obs_date=None):

    dt_file = get_caldb_file('deadtime',filter_sdd,obs_date)
    
    with fits.open(dt_file) as hdul:
        
        # eff_factor = hdul[1].data['EFFICIENCY-FACTOR'][0]
        tau_temporal = hdul[1].data['DEAD-TIME-TEMPORAL'][0]
        offset_cr1 = hdul[1].data['OFFSET-COUNTRATE-1'][0]
        offset_cr2 = hdul[1].data['OFFSET-COUNTRATE-2'][0]
        
    return tau_temporal, offset_cr1, offset_cr2, dt_file

def apply_deadtime_correction(pi_file, hk_file, output_file=None,clobber=True):
    """
    Applies deadtime correction to a Type II PI file.
    Updates the EXPOSURE column based on HK data.
    """
    hdu1 = fits.open(pi_file)

    if hdu1[0].header['CONTENT'] != 'Type II PHA file':
        raise TypeError('Input File is not Type II PHA file.')


    filter_sdd = hdu1[1].header['FILTER']
    obs_date = hdu1[0].header['OBS_DATE']
    tau_temporal, offset_cr1, offset_cr2, dt_file = get_deadtime_params(filter_sdd,obs_date)

    hk_hdul = fits.open(hk_file)
    hk_data = hk_hdul[1].data

    slow_cr = hk_data['SLOW_COUNTS']
    fast_cr = hk_data['FAST_COUNTS']
    
    if np.nanmin(fast_cr) < 500:
        offset_cr = offset_cr2
    else:
        offset_cr = offset_cr1

    fast_cr_dt_corr = -lambertw(-fast_cr * tau_temporal).real / tau_temporal

    corr_factor = (fast_cr_dt_corr - offset_cr) / slow_cr
    new_exposures = 1/corr_factor
    new_exposures = np.clip(new_exposures,0,1)
    

    data = hdu1[1].data
    data['EXPOSURE'] = new_exposures

    header = hdu1[1].header

    header['DTCORR'] = (True, 'Deadtime correction applied')

    header['HISTORY'] = f"Deadtime corrected using {os.path.basename(dt_file)}"
    header['HISTORY'] = f"Deadtime correction Offset Count Rate ={offset_cr}"

    if output_file is None:
        pi_file_basename = os.path.basename(pi_file)
        pi_file_basename = pi_file_basename.split('.')[0]
        output_file = f"{pi_file_basename}_dt_corr.pi.gz"

    hdu1.writeto(output_file,overwrite=clobber,checksum=True)

    return output_file


def solexs_deadtime_correction_cli():
    parser = argparse.ArgumentParser(description="Apply Deadtime Correction to Level 1 PI spectrogram file (Type II)")

    parser.add_argument("-i", "--infile", required=True, help="Path to the Level 1 PI spectrogram file (Type II)")
    parser.add_argument("-hk", "--hkfile", required=True, help="Path to the Level 1 Housekeeping file")
    parser.add_argument("-o", "--outfile", help="Output filename (default: <input>_dt_corr.pi)")
    parser.add_argument('-c','--clobber', type=bool, default= False, help='Overwrite existing file if it exists')

    args = parser.parse_args()

    try:
        output_file = apply_deadtime_correction(args.infile, args.hkfile, args.outfile, args.clobber)
        print(f"Output written to {output_file}.")
    except Exception as e:
        print(f"Error: {e}")
        exit(1)
