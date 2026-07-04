

"""
File Name   : deadtimme_correction_user.py
Project     : Pre-flare study
Author      : Adithya H N
Created On  : 2026-02-12
Last Updated: 2026-02-12
Version     : 1.0

Description
-----------
prep for solex fitting.
user input file for deadtime correction
"""

import deadtime
import numpy as np
from astropy.io import fits
import subprocess   

pi_file='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/solexs/AL1_SLX_L1_20241009_v1.0/SDD2/AL1_SOLEXS_20241009_SDD2_L1.pi.gz'
hk_file='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/data/solexs/AL1_SLX_L1_20241009_v1.0/SDD2/AL1_SOLEXS_20241009_SDD2_L1.hk'
pi_data,pi_header=fits.open(pi_file)

dead_time_corrected=deadtime.apply_deadtime_correction(pi_file, hk_file, output_file='dt_corr.pi',clobber=True)

subprocess.run(['gzip' , 'dt_corr.pi'])


# fits.writeto('dt_corr.pi',dead_time_corrected, pi_header,overwrite=True)