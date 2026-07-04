import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates


header_data=pd.read_csv("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/SHARP_DATA/sharp_headers.csv")

dt=pd.to_datetime(header_data['T_REC'], format="%Y.%m.%d_%H:%M:%S")
us_flux=np.array(header_data['USFLUX'])
mn_shear=np.array(header_data['MEANSHR'])
mn_pot=np.array(header_data['MEANPOT'])
tot_pot=np.array(header_data['TOTPOT'])
shear=np.array(header_data['SHRGT45'])
mn_jz=np.array(header_data['MEANJZD'])
tot_jz=np.array(header_data['TOTUSJZ'])
mn_bh=np.array(header_data['MEANGBH'])

#print(dt)
fig,ax=plt.subplots(figsize=(14,8))
ax2=ax.twinx()
ax.plot(dt,us_flux/np.max(us_flux),label='USFLUX')
ax.plot(dt,mn_shear/np.max(mn_shear),label='MEANSHR')
ax.plot(dt,mn_pot/np.max(mn_pot),label='MEANPOT')
ax.plot(dt,tot_pot/np.max(tot_pot),label='TOTPOT')
ax.plot(dt,shear/np.max(shear),label='SHRGT45')
ax2.plot(dt,mn_jz/np.max(mn_jz),label='MEANJZD',color='k')
ax.plot(dt,tot_jz/np.max(tot_jz),label='TOTUSJZ')
ax.plot(dt,mn_bh/np.max(mn_bh),label='MEANGBH')
ax.legend()
ax2.legend()
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.title('Normalized SHARP params')
plt.savefig('SHARP_params.png',dpi=300)
plt.show()