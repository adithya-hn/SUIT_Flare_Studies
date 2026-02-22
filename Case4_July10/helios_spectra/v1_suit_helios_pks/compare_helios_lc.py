import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


spec_lc=np.genfromtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/helios_spectra/binned_counts_10_30kev_lc.csv',delimiter=',',names=True, dtype=None, encoding='utf-8')
evt_lc =np.genfromtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/helios_spectra/helios_CdTe_c4.csv',delimiter=',',names=True, dtype=None, encoding='utf-8')

spec_dt=np.array(spec_lc['Date'],dtype='datetime64')
spec_count=np.array(spec_lc['counts'],dtype=float)
spec_err=np.sqrt(spec_count)
evt_dt=np.array(evt_lc['Time'],dtype='datetime64')
evt_ct=np.array(evt_lc['Total'],dtype=float)
evt_er=np.array(evt_lc['CdTeEr'],dtype=float)
print(spec_lc)
fig,ax=plt.subplots(figsize=(14,8))
ax1=ax.twinx()
ax.errorbar(evt_dt,evt_ct,yerr=evt_er,alpha=0.6)
#ax1.errorbar(spec_dt,spec_count,yerr=spec_err,color='g',alpha=0.6)
ax.set_yscale('log')
ax1.set_yscale('log')
plt.show()


