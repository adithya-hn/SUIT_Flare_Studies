"""
@Author      : Adithya H N
@Created On  : 2026-06-30
@Last Updated: 2026-06-30
@Project     : Pre-flare study using Aditya-L1
@Version     : 1.0

@Description
-----------
Brief description: Plot SoLEXS light curves for finding quiet background for hot-onset analysis.

"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime,timedelta
# from astropy.time import Time, TimeDelta
import numpy as np
import astropy.units as u
import astropy.io.fits as fits
import scienceplots
plt.style.use('science')


lc1_file ="AL1_SOLEXS_20240710_SDD2_L1_2.03_11.99keV_30sec.lc" #solexs 2-12 kev
lc2_file= "AL1_SOLEXS_20240710_SDD2_L1_12.08_21.99keV_30sec.lc" #solexs 12-22 kev
lc3_file= 'LightCurve.csv' #helios above 22 kev
lc1_data=fits.open(lc1_file)
lc2_data=fits.open(lc2_file)
lc3_data=pd.read_csv(lc3_file)

lc1_dt=np.array(lc1_data[1].data['TIME'],dtype='datetime64[s]')
lc1_rate=np.array(lc1_data[1].data['COUNTS'],dtype='float')

lc2_dt=np.array(lc2_data[1].data['TIME'],dtype='datetime64[s]')
lc2_rate=np.array(lc2_data[1].data['COUNTS'],dtype='float')

lc3_dt=np.array(lc3_data['TIME'],dtype='datetime64[s]')
lc3_rate=np.array(lc3_data['COUNTS'],dtype='float')

flare_start = datetime.fromisoformat("2024-07-10T15:25:00")
flare_peak  = datetime.fromisoformat("2024-07-10T15:37:00")
background_start = datetime.fromisoformat("2024-07-10T15:20:30")
background_end   = datetime.fromisoformat("2024-07-10T15:23:30")

plt.figure(figsize=(10, 6))
plt.plot(lc1_dt, lc1_rate, marker='o', linestyle='-', markersize=2, label='SoLEXS 2-12 keV',color='blue')
# plt.plot(lc2_dt, lc2_rate, marker='o', linestyle='-', markersize=2)
# plt.plot(lc3_dt, lc3_rate, marker='o', linestyle='-', markersize=2,label=r'HEL1OS $>$ 22 keV',color='orange')

plt.axvline(flare_start, color='red', linestyle='--', label='Flare Start (15:25:00 UT)')
plt.axvline(flare_peak, color='green', linestyle='--', label='Flare Peak (15:37:00 UT)')
# plt.axvline(background_start, color='blue', linestyle='--', label='Background Start')
# plt.axvline(background_end, color='purple', linestyle='--', label='Background End')
plt.axvspan(background_start, background_end, color='gray', alpha=0.3, label='Background Interval (15:20:30-15:23:30 UT)')

plt.legend()
plt.yscale('log')
plt.xlabel('Time')
plt.ylabel('Count rate')
plt.title('SoLEXS-HEL1OS flare onset light curves')
plt.xlim(background_start-timedelta(minutes=60), flare_peak+timedelta(minutes=10))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
plt.gcf().autofmt_xdate()
# plt.grid()
plt.savefig('solexs_lc.png', dpi=300, bbox_inches='tight')
plt.show()