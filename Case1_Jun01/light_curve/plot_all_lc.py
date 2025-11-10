import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

files = {
    "NB01": "csv_files/c1_NB01_lc_data.csv",
    "NB02": "csv_files/c1_NB02_lc_data.csv",
    "NB03": "csv_files/c1_NB03_lc_data.csv",
    "NB04": "csv_files/c1_NB04_lc_data.csv",
    "NB05": "csv_files/c1_NB05_lc_data.csv",
    "NB06": "csv_files/c1_NB06_lc_data.csv",
    "NB08": "csv_files/c1_NB08_lc_data.csv",
    "HEL1OS":"csv_files/helios_CdTe_c1.csv"
}
#"NB07": "csv_files/c1_NB07_lc_data.csv",
fig,ax1=plt.subplots(figsize=(12,5))
ax2=ax1.twinx()
for filt, fname in files.items():
    # fig,ax1=plt.subplots(figsize=(12,5))
    # ax2=ax1.twinx()
    df = pd.read_csv(fname)
    df['Time'] = pd.to_datetime(df['Time'])


    # Normalize intensity between 0 and 1
    checked=(df['Total']>0.2) #for helios
    
    #df = df[checked] #df['Total']
    y=df['Total']
    y_norm = (y - y.min()) / (y.max() - y.min())
    if filt=='HEL1OS':
        ax2.errorbar(df['Time'], df['Total'],yerr=df["CdTe1+2Er"] ,label=filt,fmt='ko-',capsize=2,markersize=2, linewidth=1.5,alpha=0.4)
    else:
        ax1.plot(df['Time'], y_norm, label=filt, linewidth=1.5)
    

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xlabel('Time (UT)')
ax1.set_ylabel('Intensity (Counts/s)')
ax2.set_ylabel('HEL1OS count rate (c/min)')
ax2.set_yscale('log')
ax2.legend(loc='upper right')
ax1.legend(loc='upper left')

plt.title('Multi-filter Light Curves')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('All_filter_lc.png',dpi=300)
plt.show()
