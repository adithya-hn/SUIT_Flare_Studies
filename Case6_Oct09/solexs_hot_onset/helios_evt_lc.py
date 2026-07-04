import subprocess
import os
from astropy.io import fits
import glob
from datetime import datetime,timedelta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from scipy.signal import find_peaks



aa = fits.open("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/suit_helios_pks/evt.fits");
c1 = aa[1].data;        # Cdte1 event data
c2 = aa[2].data;    

t_start = pd.to_datetime(c1['utc-isot'][0]); t_end =  pd.to_datetime(c1['utc-isot'][-1]);
print("Start time:", t_start)
print("End time:", t_end)

command = [
    "python3",
    "LcFromEvtFile_V06.py",
    "--fits_file", "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/suit_helios_pks/evt.fits",  
    "--mode", "CDTE_TOTAL",
    "--start_time", "2024-10-09T01:10:00.000",
    "--end_time",   "2024-10-09T02:11:00.000",
    "--time_bin_size", "30",
    "--output", "LightCurve.csv"
]


# Inputs that the script expects after running
user_inputs = "1\n22 30\n n\n"

result = subprocess.run(
    command,
    input=user_inputs,
    text=True,
    capture_output=True
)
