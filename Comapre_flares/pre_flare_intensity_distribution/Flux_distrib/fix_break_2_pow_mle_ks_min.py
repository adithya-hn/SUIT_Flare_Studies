import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma
import seaborn as sns
import glob
import os
from scipy.optimize import curve_fit
import scienceplots
plt.style.use('science')
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, curve_fit
from scipy.stats import kstest
from scipy.optimize import minimize_scalar
from itertools import product
from scipy.optimize import minimize

scol =sns.color_palette("colorblind")

files = sorted(glob.glob("csv_files/*NB04.csv"))
time_fls=sorted(glob.glob("x-ray_peaks/*.csv"))
norm =np.genfromtxt('norm_vals.csv',delimiter=',')


all_times = []
all_ints = []
t_peaks=[]
for t in time_fls:
    tp=np.genfromtxt(t, delimiter=',', names=True, dtype=None, encoding='utf-8')
    t_peaks.append(tp['date_time'])

t_peaks=np.concatenate(t_peaks)
# print(t_peaks)
l=0
for f in files:

    data = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding='utf-8')
    hpks=(np.genfromtxt(time_fls[l], delimiter=',', names=True, dtype=None, encoding='utf-8'))['date_time']
    h_counts=(np.genfromtxt(time_fls[l], delimiter=',', names=True, dtype=None, encoding='utf-8'))['helios_counts']
    l+=1
    h_pks=np.array(hpks,dtype='datetime64[ms]')
    # Convert time to datetime64
    t = np.array(data['Date'], dtype='datetime64[ms]')
    area=np.array(data['area'],dtype=float)
    # x = np.array(data['diff_count'], dtype=float)#/area
    x = np.array(data['Img_count'], dtype=float)#*(area*11.444e-6)#(np.pi/3600/180)**2*0.698*0.698)
    # print(x)
    case_id=os.path.basename(f)[0:3]
    flt=os.path.basename(f)[-8:-4]
    #print(flt)
    
    if case_id=='c04':
        t_start =  np.datetime64("2024-07-10T03:59:00")
        t_end   =  np.datetime64("2024-07-10T05:44:00")
        m=1.87e3
        c=3.39e6
        

    elif case_id=='c05':
        t_start =  np.datetime64("2024-07-10T13:37:00")
        t_end   =  np.datetime64("2024-07-10T15:25:00")
        m=1.87e3
        c=3.39e6

    elif case_id=='c06':
        t_start =  np.datetime64("2024-10-08T23:56:00")
        t_end   =  np.datetime64("2024-10-09T01:25:00")
        m=2.32e3
        c=-2.77e6
        

    elif case_id=='c07':
        t_start =  np.datetime64("2024-11-01T00:16:00")
        t_end   =  np.datetime64("2024-11-01T02:05:00")
        m=2.84e3
        c=-4.78e6
        

    elif case_id=='c08':
        t_start =  np.datetime64("2024-11-01T12:31:00")
        t_end   =  np.datetime64("2024-11-01T14:18:00")
        m=2.84e3
        c=-4.78e6

    elif case_id=='c09':
        t_start =  np.datetime64("2024-11-12T22:22:00")
        t_end   =  np.datetime64("2024-11-13T00:10:00")
        m=3.28e3
        c=-4.43e6

    elif case_id=='c10':
        t_start =  np.datetime64("2024-11-13T15:08:00")
        t_end   =  np.datetime64("2024-11-13T16:57:00")
        m=3.28e3
        c=-4.43e6
    
    mask = (t >= t_start) & (t <= t_end)
    t_cut = t[mask]
    img_count_cut  = (x[mask]*m+c*area[mask])*((np.pi/3600/180)**2*0.698*0.698)
    # print(area[mask])
    mid_cut=np.median(img_count_cut)

    fig,ax=plt.subplots(figsize=(14,6))

    threshold=0#2.5e6 #
    peaks_=find_peaks(img_count_cut, height=threshold)
    dt_=np.array(t_cut[peaks_[0]])
    pk_=np.array(img_count_cut[peaks_[0]])
    for i in range(len(peaks_[0])):
        plt.axvline(t_cut[peaks_[0][i]], color='g',alpha=0.2)

    peaks,props = find_peaks(img_count_cut, prominence=0)#,width=3)
    dt=dt_#np.array(t_cut[peaks])
    pk=pk_#np.array(img_count_cut[peaks])

    all_times.append(dt)
    all_ints.append(pk)

    df = pd.DataFrame({"Peak time": dt,  "Peak total count": pk})
    df.to_csv(f"pks/{case_id}_{flt}_peaks.csv", index=False)
    ax.plot(t_cut,img_count_cut)

    for k in range(len(h_pks)):
        plt.axvline(h_pks[k], color='r',alpha=0.2)
        # print(h_pks[k],h_counts[k])
    plt.scatter(h_pks,h_counts,marker='+',c='magenta')
    #plt.yscale('log')
    plt.title(f'{case_id}')
    plt.axhline(mid_cut,color='r',linestyle='-',label='Median')
    plt.axhline(threshold,color='g', linestyle=':',label='Threshold')
    plt.legend()
    plt.yscale('log')
    plt.savefig(f'pks/{case_id}_{flt}_local_peaks.png',dpi=300)
    plt.close()


# Convert to datetime64 if not already
t_peaks = np.array(t_peaks, dtype='datetime64[ms]')
tol = np.timedelta64(30, 's')
matched_ints = []
all_c = np.concatenate(all_ints)
all_t = np.concatenate(all_times)

print('number peaks: ',len(all_c),len(t_peaks))

for dt_arr, pk_arr in zip(all_times, all_ints):
    dt_arr = np.array(dt_arr, dtype='datetime64[ms]')
    
    for t_x in t_peaks:
        mask = np.abs(dt_arr - t_x) <= tol
        matched_ints.extend(pk_arr[mask])

matched_ints = np.array(matched_ints)
matched_ints = matched_ints#[matched_ints >= Emin_opt]
matched_ints = np.sort(matched_ints)

E =np.sort(np.array(all_c))


eb1 = 0.36
eb2 = 0.99

seg1 = E[(E >= eb1) & (E < eb2)]
seg2 = E[E >= eb2]

def mle_alpha(E, Emin):
    return 1.0 + len(E) / np.sum(np.log(E / Emin))

def neg_logL_truncated(alpha, data, xmin, xmax):
    if alpha <= 1:
        return np.inf

    N = len(data)

    norm = (xmax**(1-alpha) - xmin**(1-alpha)) / (1-alpha)

    return alpha * np.sum(np.log(data)) + N * np.log(norm)


def fit_alpha_truncated(data, xmin, xmax):
    res = minimize_scalar(
        neg_logL_truncated,
        bounds=(1.01, 10),
        args=(data, xmin, xmax),
        method='bounded'
    )
    return res.x


def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

def model_ccdf(E, Emin, alpha):
    return (E / Emin) ** (-(alpha - 1))

def get_ccdf_y(E_sorted, ccdf_full, values):
    y_vals = []
    
    for v in values:
        idx = np.searchsorted(E_sorted, v, side='left')
        if idx < len(ccdf_full):
            y_vals.append(ccdf_full[idx])
        else:
            y_vals.append(np.nan)
    
    return np.array(y_vals)

seg1 = E[(E >= eb1) & (E < eb2)]
seg2 = E[(E >= eb2)&(E<1.7)]

alpha1 = fit_alpha_truncated(seg1, eb1, eb2)   # truncated since upper is cutoff
alpha2 = mle_alpha(seg2, eb2)                  # third segment

# alpha1 = fit_alpha_truncated(seg1, eb1, eb2)
# alpha2 = fit_alpha_truncated(seg2, eb2, max(E))

print(alpha1, alpha2)
alpha = fit_alpha_truncated(E, eb1, max(E))
print(alpha)

E1fit=seg1
E2fit=seg2
# ccdf_emp = empirical_ccdf(E1fit)
# ccdf_model = model_ccdf(E1fit, eb1, alpha1)

E_all = np.sort(E)                     # all base-subtracted peak counts
N_all = np.arange(len(E_all), 0, -1)
N0 = len(E1fit)  # number of events above Emin
model_ccdf_seg1 = 0.64 * (E1fit / eb1) ** (-(alpha1 - 1))
model_ccdf_seg2 = 0.3 * (E2fit / eb2) ** (-(alpha2 - 1))

ccdf = empirical_ccdf(E)
y_matched = get_ccdf_y(E, ccdf, matched_ints)


alpha_sig1 = (alpha1 - 1) / np.sqrt(len(E1fit))
alpha_sig2 = (alpha2 - 1) / np.sqrt(len(E2fit))


plt.figure(figsize=(12,8))
plt.title('CCDF of all pre-flare transients',fontsize=24)

# Power-law fit (only above Emin)
plt.loglog(E, ccdf, '+',color='k', label='SUIT transeint events')
plt.loglog(E1fit, model_ccdf_seg1, '-',linewidth=2,label=fr'Power-law fit (alpha = {alpha1:.2f}$\pm$ {alpha1:.2f} )')
plt.loglog(E2fit, model_ccdf_seg2, '-',linewidth=2,label=fr'Power-law fit (alpha = {alpha2:.2f}$\pm$ {alpha2:.2f} )')

plt.scatter(matched_ints, y_matched,  facecolors='none', edgecolors='red',s=70, linewidth=1.5,  label='Co-temporal HEL1OS')

# Emin marker
plt.axvline(eb1, color='b', linestyle='--',label=r'$E_{\min}$')

# y_matched = get_ccdf_y(E, ccdf, matched_ints)

plt.title('CCDF of pre-flare transients')
plt.xlabel(r'Energy ($erg^-1 \times cm^-2$)', fontsize=18)
plt.ylabel('CCDF', fontsize=18)
plt.legend(fontsize=14)
plt.grid(True, which='both', alpha=0.3)
plt.savefig('ccdf_mle_ks_min2.pdf',dpi=300)
plt.show()
