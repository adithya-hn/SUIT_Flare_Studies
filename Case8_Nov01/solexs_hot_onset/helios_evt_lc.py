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

command = [
    "python3",
    "LcFromEvtFile_V06.py",
    "--fits_file", "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/helios/2024/11/01/HLS_20241101_120000_43198sec_lev1_V111/events/evt.fits",
    "--mode", "CDTE_TOTAL",
    "--start_time", "2024-11-01T14:10:00.000",
    "--end_time", "2024-11-01T14:40:00.000",
    "--time_bin_size", "30",
    "--output", "LightCurve.csv"
]

# Inputs that the script expects after running
user_inputs = "1\n22 30\nn \n"

result = subprocess.run(
    command,
    input=user_inputs,
    text=True,
    capture_output=True
)
