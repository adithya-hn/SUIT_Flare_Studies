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
from scipy.optimize import minimize
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

scol =sns.color_palette("colorblind")

files = sorted(glob.glob("csv_comp_detections/*.csv"))
norm =np.loadtxt('norm_vals.csv',delimiter=',')


all_times = []
all_ints = []

for f in files:
    data = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding='utf-8')
    # Convert time to datetime64
    t = np.array(data['date_time'], dtype='datetime64[ms]')
    area=np.array(data['area'],dtype=float)
    #x = np.array(data['diff_count'], dtype=float)#/area
    x = np.array(data['diff_count'], dtype=float)#/area
    
    case_id=os.path.basename(f)[0:3]
    flt=os.path.basename(f)[-8:-4]
    #print(flt)
    if case_id=='c02':
        t_start =  np.datetime64("2024-06-02T02:50:00")
        t_end   =  np.datetime64("2024-06-02T04:41:00")
        norm_v=norm[0]
        
    
    if case_id=='c03':
        t_start =  np.datetime64("2024-06-02T06:50:00")
        t_end   =  np.datetime64("2024-06-02T08:40:00")
        norm_v=norm[0]
    
    if case_id=='c04':
        t_start =  np.datetime64("2024-07-10T03:59:00")
        t_end   =  np.datetime64("2024-07-10T05:44:00")
        norm_v=norm[0]
        

    elif case_id=='c05':
        t_start =  np.datetime64("2024-07-10T13:37:00")
        t_end   =  np.datetime64("2024-07-10T15:25:00")
        norm_v=norm[1]

    elif case_id=='c06':
        t_start =  np.datetime64("2024-10-08T23:56:00")
        t_end   =  np.datetime64("2024-10-09T01:25:00")
        norm_v=norm[2]
        

    elif case_id=='c07':
        t_start =  np.datetime64("2024-11-01T00:16:00")
        t_end   =  np.datetime64("2024-11-01T02:05:00")
        norm_v=norm[3]
        print(x)
        

    elif case_id=='c08':
        t_start =  np.datetime64("2024-11-01T12:31:00")
        t_end   =  np.datetime64("2024-11-01T14:18:00")
        norm_v=norm[4]

    elif case_id=='c09':
        t_start =  np.datetime64("2024-11-12T22:22:00")
        t_end   =  np.datetime64("2024-11-13T00:10:00")
        norm_v=norm[5]

    elif case_id=='c10':
        t_start =  np.datetime64("2024-11-13T15:08:00")
        t_end   =  np.datetime64("2024-11-13T16:57:00")
        norm_v=norm[6]
    
    mask = (t >= t_start) & (t <= t_end)
    t_cut = t[mask]
    img_count_cut  = x[mask]/norm_v
    mid_cut=np.median(img_count_cut)

    fig,ax=plt.subplots(figsize=(14,6))

    threshold=2.0e3 #5661*4 #2.5e6 #
    peaks=find_peaks(img_count_cut, height=threshold)
    
    dt=np.array(t_cut[peaks[0]])
    pk=np.array(img_count_cut[peaks[0]])

    all_times.append(t)
    all_ints.append(x)
    

    df = pd.DataFrame({"Peak time": dt,  "Peak total count": pk})
    df.to_csv(f"pks/{case_id}_{flt}_peaks.csv", index=False)
    ax.plot(t_cut,img_count_cut)

    for i in range(len(peaks[0])):
        plt.axvline(t_cut[peaks[0][i]], color='b',alpha=0.2)
    #plt.yscale('log')
    plt.title(f'{case_id}')
    plt.axhline(mid_cut,color='r',linestyle='-',label='Median')
    plt.axhline(threshold,color='g', linestyle=':',label='Threshold')
    plt.legend()
    plt.savefig(f'pks/{case_id}_{flt}_local_peaks.png',dpi=300)
    plt.close()


all_c = np.concatenate(all_ints)
all_t = np.concatenate(all_times)
plt.scatter(all_t, all_c)
plt.ylim(1e3,1e5)
plt.show()
# print(all_c)
E = np.array(all_c)
E = E[E > 0]
E = np.sort(E)
print('Total events detected:', len(E))
def mle_alpha(E, Emin):
    return 1.0 + len(E) / np.sum(np.log(E / Emin))

def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

def model_ccdf(E, Emin, alpha):
    return (E / Emin) ** (-(alpha - 1))

def ks_distance(E, Emin):
    Efit = E[E >= Emin]
    n = len(Efit)

    # Require enough statistics
    if n < 30:
        return np.nan, np.nan

    alpha = mle_alpha(Efit, Emin)

    ccdf_emp = empirical_ccdf(Efit)
    ccdf_mod = model_ccdf(Efit, Emin, alpha)

    D = np.max(np.abs(ccdf_emp - ccdf_mod))
    return D, alpha
Emins = np.unique(E)
ks_vals = []
alphas = []

for Emin in Emins:
    D, a = ks_distance(E, Emin)
    ks_vals.append(D)
    alphas.append(a)

ks_vals = np.array(ks_vals)
alphas = np.array(alphas)
valid = ~np.isnan(ks_vals)

Emins_valid = Emins[valid]
ks_valid = ks_vals[valid]
alphas_valid = alphas[valid]

idx = np.argmin(ks_valid)

Emin_opt = Emins_valid[idx]

alpha_opt = alphas_valid[idx]

print("Optimal Emin =", Emin_opt)
print("MLE alpha =", alpha_opt)
# plt.figure(figsize=(6,5))
# plt.loglog(Emins_valid, ks_valid, 'o-')
# plt.axvline(Emin_opt, color='r', linestyle='--', label=r'Optimal $E_{\min}$')
# plt.xlabel(r'$E_{\min}$')
# plt.ylabel('KS distance')
# plt.legend()
# plt.grid(True, which='both', alpha=0.3)
# plt.show()

Efit = E[E >= Emin_opt]
alpha_final = mle_alpha(Efit, Emin_opt)

print("Final alpha =", alpha_final)
print("Number of events =", len(Efit))

sigma_alpha = (alpha_final - 1) / np.sqrt(len(Efit))
print("Alpha uncertainty ~", sigma_alpha)

Efit = np.sort(Efit)
ccdf_emp = empirical_ccdf(Efit)

ccdf_model = model_ccdf(Efit, Emin_opt, alpha_final)
E_all = np.sort(E)                     # all base-subtracted peak counts
N_all = np.arange(len(E_all), 0, -1)
N0 = len(Efit)  # number of events above Emin
model_ccdf = N0 * (Efit / Emin_opt) ** (-(alpha_final - 1))

plt.figure(figsize=(12,8))
plt.title('CCDF of pre-flare transients',fontsize=24)
# Full CCDF (all events)
plt.loglog(E_all, N_all, 'o', markersize=4,
           label='All detected events')

# Power-law fit (only above Emin)
plt.loglog(Efit, model_ccdf, '-',
           linewidth=2,
           label=f'Power-law fit (alpha = {alpha_final:.2f})')

# Emin marker
plt.axvline(Emin_opt, color='k', linestyle='--',label=r'$E_{\min}$')

plt.xlabel('Peak excess Mg II h counts',fontsize=18)
plt.ylabel(r'N($\geq$Counts)',fontsize=18)
plt.legend(fontsize=18)
#plt.grid(True, which='both', alpha=0.3)
#plt.tight_layout()
plt.savefig('MLE_Power_law_fit.pdf',dpi=300)
plt.show()

E = np.array(all_c)
E = E[E > 0]
E = np.sort(E)

def neg_log_likelihood(params, E, Emin):
    alpha1, alpha2, Eb = params

    if alpha1 <= 1 or alpha2 <= 1 or Eb <= Emin:
        return np.inf

    # Split data
    E1 = E[(E >= Emin) & (E < Eb)]
    E2 = E[E >= Eb]

    # Normalization
    Cinv = (
        (Eb**(1-alpha1) - Emin**(1-alpha1)) / (1-alpha1)
        + Eb**(alpha2-alpha1) * (Eb**(1-alpha2)) / (alpha2-1)
    )
    C = 1 / Cinv

    logL = 0.0
    if len(E1) > 0:
        logL += np.sum(np.log(C) - alpha1*np.log(E1))
    if len(E2) > 0:
        logL += np.sum(np.log(C) + (alpha2-alpha1)*np.log(Eb)
                       - alpha2*np.log(E2))

    return -logL

#Emin_opt=E[6]
Efit = E[E >= Emin_opt]

init = [2.0, 2.5, np.median(Efit)]
bounds = [(1.01, 5), (1.01, 5), (Emin_opt*1.1, Efit.max())]

res = minimize(neg_log_likelihood,init,args=(Efit, Emin_opt),bounds=bounds)

alpha1_mle, alpha2_mle, Eb_mle = res.x

print(alpha1_mle, alpha2_mle, Eb_mle)

E_sorted = np.sort(E)
N_emp = np.arange(len(E_sorted), 0, -1)

def broken_ccdf(E, Emin, alpha1, alpha2, Eb, N0):
    ccdf = np.zeros_like(E, dtype=float)

    mask1 = (E >= Emin) & (E < Eb)
    mask2 = E >= Eb

    ccdf[mask1] = N0 * (E[mask1]/Emin)**(-(alpha1-1))
    ccdf[mask2] = (
        N0 * (Eb/Emin)**(-(alpha1-1))
        * (E[mask2]/Eb)**(-(alpha2-1))
    )

    return ccdf

#Emin_opt=Emin
# N0 = np.sum(E >= Emin_opt)
# E_plot = E_sorted[E_sorted >= Emin_opt]

# ccdf_model = broken_ccdf(E_plot, Emin_opt,alpha1_mle, alpha2_mle, Eb_mle, N0)

# plt.figure(figsize=(12,8))
# plt.loglog(E_sorted, N_emp, 'o', ms=4, label='All events')
# plt.loglog(E_plot, ccdf_model, 'r-', lw=2,label='Broken power-law MLE')

# plt.axvline(Emin_opt, ls='--', c='k', label=r'$E_{\min}$')
# plt.axvline(Eb_mle, ls=':', c='r', label=r'$E_b$')

# plt.xlabel('Peak excess Mg II h counts')
# plt.ylabel(r'N($\geq$E)')
# plt.legend()
# plt.grid(True, which='both', alpha=0.3)
# plt.xlim(1e4, 1e8)
# plt.savefig('Broken_MLE_Power_law_fit.png',dpi=300)
# plt.show()

# Logarithmic bins
bins = np.logspace(np.log10(Emin_opt),
                   np.log10(E.max()), 12)

hist, edges = np.histogram(E[E >= Emin_opt], bins=bins)
bin_centers = np.sqrt(edges[:-1] * edges[1:])
dNdF = hist / np.diff(edges)

plt.loglog(bin_centers, dNdF, 'o')
plt.xlabel('Peak excess Mg II h counts')
plt.ylabel('dN/dF')
plt.title('Differential size distribution')
#plt.grid(True, which='both', alpha=0.3)
plt.savefig('dNdF_plot.png',dpi=300)
plt.show()
