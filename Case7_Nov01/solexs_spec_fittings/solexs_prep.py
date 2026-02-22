"""
@Author      : Adithya H N
@Created On  : 2026-02-14
@Last Updated: 2026-02-14
@Project     : Project Name
@Version     : 1.0

@Description
-----------
Brief description
"""


import numpy as np
import subprocess
from astropy.io import fits
import matplotlib.pyplot as plt
import os

date_str = "20241101"
dir_path = f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/solexs/AL1_SLX_L1_{date_str}_v1.0/SDD2/'

pi_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.pi.gz")
lc_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.lc.gz")
hk_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.hk.gz")
gti_file = os.path.join(dir_path,f"AL1_SOLEXS_{date_str}_SDD2_L1.gti.gz")

dt_start='2024-11-01T01:33:29' 
dt_end='2024-11-01T01:34:29'


# subprocess.run(['solexs-dtcorr' ,'-i',pi_file,'-hk',hk_file,'-o','dt_corr.pi','--clobber','True'])
# subprocess.run(['gzip' , 'dt_corr.pi'])
t_st=subprocess.run(['solexs-utc2time' , dt_start],capture_output=True, text=True)
t_end=subprocess.run(['solexs-utc2time' , dt_end],capture_output=True,text=True)

# output = result.stdout.strip()

# Extract timestamp
timestamp1 = (t_st.stdout.strip()).split(":")[1].strip()
timestamp2 = (t_end.stdout.strip()).split(":")[1].strip()
# Convert to integer if needed
start_t = timestamp1
end_t=timestamp2


print(start_t,end_t)


dt_corr_fl='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/solexs_spec_fittings/dt_corr.pi.gz'
spec_file='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/solexs_spec_fittings/solexs_spec_2.pi'
subprocess.run(['solexs-genspec', '-i',dt_corr_fl,  '-tstart',start_t, '-tstop', end_t, '-gti', gti_file,'-o','spec_file_1', '--clobber', 'True'])
