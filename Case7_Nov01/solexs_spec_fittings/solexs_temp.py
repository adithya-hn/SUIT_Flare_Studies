from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ImageNormalize, AsinhStretch
import datetime
import astropy.units as u

data=fits.open('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/light_curve/csv_files/AL1_SOLEXS_20241009_SDD2_L1_puc_tb_fit_results_T_EM.fits')


data.info()

print(data[1].columns)

temp = data[1].data['temperature']
temp_er= data[1].data['TEMPERATURE_ERR']


time = data[1].data['TIME']
emission_measure = data[1].data['em']


Time=[]
for t in time:
    Time.append(datetime.datetime.fromtimestamp(t)-datetime.timedelta(hours=5,minutes=30))

plt.errorbar(Time, temp, yerr=temp_er, fmt='o',markersize=1, ecolor='r', capthick=2)
plt.show()

# print(time)


