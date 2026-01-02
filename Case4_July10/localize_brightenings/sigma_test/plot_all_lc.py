import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
files = {
   
    "3 Sigma": "3sig_Diff_img_data_NB04.csv",
    "4 Sigma": "4sig_Diff_img_data_NB04.csv",
    "5 Sigma": "5sig_Diff_img_data_NB04.csv",
    
}#"HEL1OS":  "helios_CdTe_c4.csv"
#"NB07": "csv_files/c1_NB07_lc_data.csv", "NB01": "csv_files/c5_NB01_lc_data.csv",
fig,ax1=plt.subplots(figsize=(12,5))
ax2=ax1.twinx()
lc=[]
a=[]
for filt, fname in files.items():
    df = pd.read_csv(fname)
    df['Date'] = pd.to_datetime(df['Date'])

    y=df['Img_count']
    area=df['area']
    lc.append(y)
    a.append(area)
    y_norm = y#(y - y.min()) / (y.max() - y.min())
    #ax1.plot(df['Date'], y_norm, label=filt, linewidth=1.5)
    ax2.plot(df['Date'],y/area,label=filt)
    
    

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xlabel('Time (UT)')
#ax1.set_ylabel('Intensity (Counts/s)')
ax2.set_ylabel('Area normalized intensity')
ax2.set_yscale('log')
ax2.legend(loc='upper right')
ax1.legend(loc='upper left')

plt.title('Multi-threshold Light Curves')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('All_filter_lc.png',dpi=300)
plt.show()

lc=np.array(lc)
print(((lc[1]/a[1])-(lc[2]/a[2]))/(lc[1]/a[1]))
