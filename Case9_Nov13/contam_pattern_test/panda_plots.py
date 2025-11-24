import pandas as pd
import matplotlib.pyplot as plt
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

filter='NB04'
df = pd.read_csv(f'{filter}_ar_qs_th_count.csv', parse_dates=['Time'])
drift=pd.read_csv('case9_crpix_values_2k.csv', parse_dates=['time'])

# Sort and calculate time difference
df = df.sort_values('Time')
df['dt'] = df['Time'].diff().dt.total_seconds()

# drift=drift.sort_values('Time')
# drift['dt'] = drift['Time'].diff().dt.total_seconds()

# Split where the gap is large (>120 sec)
gap_indices = df[df['dt'] > 120].index
# drift_gap_idx=drift[drift['dt']>120].index
# Plot each continuous segment
params=['QS1_tot','QS2_tot','QS3_tot','AR_tot','Meas_Exp','Threshold_count']

for par in params:
    start = 0
    fig,ax=plt.subplots(figsize=(12,5))
    #ax1=ax.twinx()
    ax2=ax.twinx()
    for idx in gap_indices:
        ax.plot(df['Time'][start:idx], df[par][start:idx], color='tab:blue')
        #ax1.plot(df['Time'][start:idx], df['QS2mean'][start:idx], color='tab:red')
        start = idx
    ax.plot(df['Time'][start:], df[par][start:], color='tab:blue',label=f'{par}')
    #ax1.plot(df['Time'][start:], df['QS2mean'][start:], color='tab:red')
    ax2.plot(drift['time'],drift['crpix1'],color='tab:green',label='2K Drift pattern')
    ax2.set_ylim(1270,1290)
    plt.xlabel('Time')
    plt.ylabel('Counts')
    plt.title(f'{filter} Light Curve')
    ax.legend(loc='upper right')
    ax2.legend(loc='upper left')
    #plt.grid(True)
    plt.savefig(f'{filter}_{par}.png',dpi=300)
    plt.close()
