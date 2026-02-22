import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

spec_lc=np.genfromtxt('binned_counts_10_30kev_lc.csv',delimiter=',', names=True,dtype=None)
srk_lc=np.genfromtxt('helios_CdTe_c5.csv',delimiter=',', names=True,dtype=None)


spec_dt=np.array(spec_lc['Date'],dtype='datetime64')
srk_dt =np.array(srk_lc['Time'], dtype='datetime64')

spec_count=np.array(spec_lc['counts'],dtype=float)
srk_count =np.array(srk_lc['Total'],dtype=float)
srk_er    =np.array(srk_lc['CdTeEr'],dtype=float)
spec_er   =np.sqrt(spec_count)

fig,ax=plt.subplots(figsize=(14,8))
plt.title('Helios light curve comparison')
ax.errorbar(spec_dt,spec_count,yerr=spec_er,label='spectra file light curve')
ax.errorbar(srk_dt,srk_count,yerr=srk_er,label='events light curve')
ax.set_yscale('log')
ax.set_xlabel('Time')
ax.set_ylabel('Counts')
ax.legend()
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('light curve.png',dpi=300)
plt.show()