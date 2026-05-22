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
from scipy.optimize import minimize_scalar

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
    x = np.array(data['Img_count'], dtype=float)#/area
    
    case_id=os.path.basename(f)[0:3]
    flt=os.path.basename(f)[-8:-4]
    #print(flt)
    
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
    # df.to_csv(f"pks/{case_id}_{flt}_peaks.csv", index=False)
    ax.plot(t_cut,img_count_cut)

    # for i in range(len(peaks)):
    #     plt.axvline(t_cut[peaks[i]], color='b',alpha=0.2)
    for k in range(len(h_pks)):
        plt.axvline(h_pks[k], color='r',alpha=0.2)
        # print(h_pks[k],h_counts[k])
    plt.scatter(h_pks,h_counts*1e7,marker='+',c='magenta')
    #plt.yscale('log')
    plt.title(f'{case_id}')
    plt.axhline(mid_cut,color='r',linestyle='-',label='Median')
    plt.axhline(threshold,color='g', linestyle=':',label='Threshold')
    plt.legend()
    # plt.savefig(f'pks/{case_id}_{flt}_local_peaks.png',dpi=300)
    plt.close()
    # print('sigma_thresh : ', np.std(img_count_cut)*0.5)


# Convert to datetime64 if not already
t_peaks = np.array(t_peaks, dtype='datetime64[ms]')

# tolerance
tol = np.timedelta64(30, 's')

matched_ints = []
all_c = np.concatenate(all_ints)
all_t = np.concatenate(all_times)

print('number peaks: ',len(all_c),len(t_peaks))

# for dt_arr, pk_arr in zip(all_t, all_c):
#     dt_arr = np.array(dt_arr, dtype='datetime64[ms]')
    
#     for t_x in t_peaks:
#         mask = np.abs(dt_arr - t_x) <= tol
#         matched_ints.extend(pk_arr[mask])

# matched_ints = np.array(matched_ints)

for dt_arr, pk_arr in zip(all_times, all_ints):
    dt_arr = np.array(dt_arr, dtype='datetime64[ms]')
    
    for t_x in t_peaks:
        mask = np.abs(dt_arr - t_x) <= tol
        matched_ints.extend(pk_arr[mask])

matched_ints = np.array(matched_ints)


E = np.array(all_c)


E = E[E > 0]
E = np.sort(E)
print('---------->',len(E))
# E = E[60:]
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
# '''
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
plt.figure(figsize=(6,5))
plt.loglog(Emins_valid, ks_valid, 'o-')
plt.axvline(Emin_opt, color='r', linestyle='--', label=r'Optimal $E_{\min}$')
plt.xlabel(r'$E_{\min}$')
plt.ylabel('KS distance')
plt.legend()
plt.grid(True, which='both', alpha=0.3)
plt.close()

Efit = E[E >= Emin_opt]
ccdf_full = empirical_ccdf(E)
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
model_CCDF = N0 * (Efit / Emin_opt) ** (-(alpha_final - 1))

matched_ints = matched_ints#[matched_ints >= Emin_opt]
matched_ints = np.sort(matched_ints)
def mle_truncated(E, Emin, Emax):
    Eseg = E[(E >= Emin) & (E < Emax)]
    N = len(Eseg)
    
    if N < 20:
        return np.nan, -np.inf
    
    sum_logE = np.sum(np.log(Eseg))
    
    def neg_logL(alpha):
        if alpha <= 1:
            return np.inf
        
        C = (alpha - 1) / (Emin**(1-alpha) - Emax**(1-alpha))
        logL = N*np.log(C) - alpha*sum_logE
        return -logL
    
    res = minimize_scalar(neg_logL, bounds=(1.01, 5), method='bounded')
    
    return res.x, -res.fun  # alpha, logL

def get_ccdf_y(E_sorted, ccdf_full, values):
    y_vals = []
    
    for v in values:
        idx = np.searchsorted(E_sorted, v, side='left')
        if idx < len(ccdf_full):
            y_vals.append(ccdf_full[idx])
        else:
            y_vals.append(np.nan)
    
    return np.array(y_vals)

y_matched = get_ccdf_y(E, ccdf_full, matched_ints)

plt.figure(figsize=(6,5))

# FULL CCDF (all data)
plt.loglog(E, ccdf_full, '.', label='All peaks')

# Power-law fit ONLY above Emin
Efit = E[E >= Emin_opt]
plt.loglog(Efit,
           model_ccdf(Efit, Emin_opt, alpha_final),
           '-', label=f'Fit (E ≥ Emin), α={alpha_final:.2f}')

# Overlay ALL co-temporal points (no filtering)
plt.scatter(matched_ints, y_matched,
            facecolors='none', edgecolors='red',
            s=70, linewidth=1.5,
            label='X-ray co-temporal')

plt.axvline(Emin_opt, color='k', linestyle='--', label='Emin')

plt.xlabel('E')
plt.ylabel('CCDF')
plt.legend()
plt.grid(True, which='both', alpha=0.3)
plt.show()


plt.figure(figsize=(12,8))
# Full CCDF (all events)
plt.loglog(E_all, N_all, 'o', markersize=4,label='All detected events')

# Power-law fit (only above Emin)
plt.loglog(Efit, model_CCDF, '-',linewidth=2, label=f'Power-law fit (α = {alpha_final:.2f})')

# Emin marker
plt.title('CCDF of all pre-flare transients',fontsize=24)
plt.axvline(Emin_opt, color='k', linestyle='--',label=r'$E_{\min}$')

plt.xlabel('Peak Mg II h counts',fontsize=18)
plt.ylabel('N(≥E)',fontsize=18)
plt.legend(fontsize=18)
plt.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.savefig('MLE_Power_law_fit.png',dpi=300)
plt.show()

E_sorted = np.sort(E)

def mle_alpha_segment(E_segment, Emin):
    return 1 + len(E_segment) / np.sum(np.log(E_segment / Emin))

Eb_candidates = np.logspace(np.log10(min(E)), np.log10(max(E)), 30)
best_score = np.inf
best_params = None

for i in range(5, len(Eb_candidates)-10):
    for j in range(i+5, len(Eb_candidates)-5):
        
        Eb1 = Eb_candidates[i]
        Eb2 = Eb_candidates[j]
        
        seg1 = E[(E >= E[0]) & (E < Eb1)]
        seg2 = E[(E >= Eb1) & (E < Eb2)]
        seg3 = E[(E >= Eb2)]
        
        if len(seg1) < 30 or len(seg2) < 30 or len(seg3) < 30:
            continue
        
        a1 = mle_alpha_segment(seg1, seg1.min())
        a2 = mle_alpha_segment(seg2, Eb1)
        a3 = mle_alpha_segment(seg3, Eb2)
        
        # simple error metric: continuity mismatch
        err = abs(a1 - a2) + abs(a2 - a3)
        
        if err < best_score:
            best_score = err
            best_params = (Eb1, Eb2, a1, a2, a3)


ccdf = empirical_ccdf(E_sorted)

plt.figure(figsize=(6,5))
plt.loglog(E_sorted, ccdf, '.', label='Data')
# Segment masks
m1 = (E_sorted < Eb1)
m2 = (E_sorted >= Eb1) & (E_sorted < Eb2)
m3 = (E_sorted >= Eb2)

plt.loglog(E_sorted[m1],(E_sorted[m1]/E_sorted[m1][0])**(-(a1-1)),'-', label=f'α1={a1:.2f}')

plt.loglog(E_sorted[m2], (E_sorted[m2]/Eb1)**(-(a2-1)),'-', label=f'α2={a2:.2f}')

plt.loglog(E_sorted[m3],(E_sorted[m3]/Eb2)**(-(a3-1)),'-', label=f'α3={a3:.2f}')

plt.axvline(Eb1, color='k', linestyle='--')
plt.axvline(Eb2, color='k', linestyle='--')

plt.xlabel('E')
plt.ylabel('CCDF')
plt.legend()
plt.grid(True, which='both', alpha=0.3)
plt.show()

def fit_single(E):
    Emin = E.min()
    Emax = E.max()
    
    alpha, logL = mle_truncated(E, Emin, Emax)
    k = 1
    
    return logL, k, (alpha,)

def fit_two(E):
    E = np.sort(E)
    Eb_candidates = np.logspace(np.log10(E.min()), np.log10(E.max()), 30)
    
    best_L = -np.inf
    best_params = None
    
    for Eb in Eb_candidates[5:-5]:
        
        a1, L1 = mle_truncated(E, E.min(), Eb)
        a2, L2 = mle_truncated(E, Eb, E.max())
        
        if np.isnan(a1) or np.isnan(a2):
            continue
        
        L = L1 + L2
        
        if L > best_L:
            best_L = L
            best_params = (a1, a2, Eb)
    
    k = 3  # α1, α2, Eb
    return best_L, k, best_params

def fit_three(E):
    E = np.sort(E)
    Eb_candidates = np.logspace(np.log10(E.min()), np.log10(E.max()), 30)
    
    best_L = -np.inf
    best_params = None
    
    for i in range(5, len(Eb_candidates)-10):
        for j in range(i+5, len(Eb_candidates)-5):
            
            Eb1 = Eb_candidates[i]
            Eb2 = Eb_candidates[j]
            
            a1, L1 = mle_truncated(E, E.min(), Eb1)
            a2, L2 = mle_truncated(E, Eb1, Eb2)
            a3, L3 = mle_truncated(E, Eb2, E.max())
            
            if np.isnan(a1) or np.isnan(a2) or np.isnan(a3):
                continue
            
            L = L1 + L2 + L3
            
            if L > best_L:
                best_L = L
                best_params = (a1, a2, a3, Eb1, Eb2)
    
    k = 5  # α1, α2, α3, Eb1, Eb2
    return best_L, k, best_params

def compute_scores(logL, k, N):
    AIC = -2*logL + 2*k
    BIC = -2*logL + k*np.log(N)
    return AIC, BIC

N = len(E)

L1, k1, p1 = fit_single(E)
L2, k2, p2 = fit_two(E)
L3, k3, p3 = fit_three(E)

AIC1, BIC1 = compute_scores(L1, k1, N)
AIC2, BIC2 = compute_scores(L2, k2, N)
AIC3, BIC3 = compute_scores(L3, k3, N)

print("Model   AIC      BIC")
print("1-seg:", AIC1, BIC1)
print("2-seg:", AIC2, BIC2)
print("3-seg:", AIC3, BIC3)