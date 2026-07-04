

import numpy as np
import subprocess
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import deadtime

command=[]


date_str = "20241101"
dir_path = f'AL1_SLX_L1_{date_str}_v2.0/SDD2/'

pi_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.pi.gz")
lc_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.lc.gz")
hk_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.hk.gz")
gti_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.gti.gz")


# dead_time_corrected=deadtime.apply_deadtime_correction(pi_file, hk_file, output_file='dt_corr.pi',clobber=True)
# solexs-dtcorr -i <l1_pi_file> -hk <l1_hk_file> [-o <outfile>] [--clobber <True/False>]

subprocess.run(['solexs-dtcorr' , 'dt_corr.pi'])

subprocess.run(['gzip' , 'dt_corr.pi'])

solexs-utc2time 2024-02-12T11:00:00

# solexs-genspec -i <l1_pi_file> -tstart <tstart> -tstop <tstop> -gti <l1_gti_file> [-o <outfile>] [--clobber <True/False>]