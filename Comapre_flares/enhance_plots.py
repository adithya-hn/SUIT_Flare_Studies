import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob("brightness_csv/*.csv"))
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

# plt.figure(figsize=(10,6))

for f in files:
    data = np.loadtxt(f, skiprows=1,delimiter=',',dtype='str')
    plt.figure(figsize=(10,6))
    t = np.array(data[:,0],dtype='datetime64[s]')     # can be different for each file
    t = (t - t[0]) / np.timedelta64(1,'s')  # convert to seconds
    y = np.array(data[:,3],dtype=float)
    #y = y / np.max(y)
    t_peak = t[np.argmax(y)]
    t_shifted = t - t_peak
    # --- Min–Max normalization ---
    y_min = np.min(y)
    y_max = np.max(y)
    y_norm = (y - y_min) / (y_max - y_min)

    plt.plot(t_shifted, y_norm, label=f[15:])

    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Value")
    #plt.yscale('log')
    plt.tight_layout()
    plt.savefig(f'suit_flares_norm_{f[15:18]}.png',dpi=300)
    plt.close()
